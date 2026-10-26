"""Tests of the test harness itself: random numbers, genomes, the read
simulator's truth, tag recomputation, the parsers, the process runner and
the out-of-tree build."""

import gzip
import os
import shutil
import sys

import pytest

import conftest
import minimap2_build
import minimap2_helpers as H


# --- random numbers and genomes -------------------------------------------

def test_rng_matches_reference_splitmix64():
	# Reference outputs of splitmix64 seeded with 0 (Vigna's splitmix64.c).
	r = H.Rng(0)
	assert [r.u64() for _ in range(3)] == [0xE220A8397B1DCDAF, 0x6E789E6AA1B965F4, 0x06C45D188009454F]


def test_rng_helpers_stay_in_range():
	r = H.Rng(5)
	assert all(0 <= r.below(7) < 7 for _ in range(1000))
	assert all(3 <= r.randint(3, 5) <= 5 for _ in range(1000))
	assert all(0.0 <= r.random() < 1.0 for _ in range(1000))
	x = [r.normal_int(300, 30) for _ in range(4000)]
	mean = sum(x) / float(len(x))
	sd = (sum((v - mean) ** 2 for v in x) / len(x)) ** 0.5
	assert abs(mean - 300) < 3 and abs(sd - 30) < 3


def test_rng_dna_alphabets():
	# four letters: 2 bits per letter, low bits first (the stream every pin depends on)
	x = H.Rng(9).u64()
	assert H.Rng(9).dna(32) == "".join("ACGT"[(x >> (2 * i)) & 3] for i in range(32))
	assert H.Rng(9).dna(40, "TGCA") == H.Rng(9).dna(40).translate(str.maketrans("ACGT", "TGCA"))
	# other sizes: every letter of the alphabet, nothing else
	for alphabet in ("N", "AC", "ACG", "ACGTN", "ACGTURYKMSWBDHVN"):
		s = H.Rng(3).dna(3000, alphabet)
		assert len(s) == 3000 and set(s) == set(alphabet)
	with pytest.raises(ValueError):
		H.Rng(1).dna(5, "")


def test_standard_genomes_are_stable():
	"""The standard genomes must not change across Python versions or platforms: pins depend on them."""
	g = conftest.make_synth()
	assert g.names == ["chr1", "chr2", "short1", "short2"]
	assert [g.length(n) for n in g.names] == [80000, 50000, 125, 96]
	assert g.digest() == "fa3f30588eaf482002964f679d65406bccb929ccd99d27a1d6f07cd87e58de91"
	t = conftest.make_synth_tx()
	assert t.names == ["chr1", "chr2"] and len(t.genes) == 8
	assert t.digest() == "66f11a56ccc351c5e62c1d0b0370980fc9a88fa4a56b97d1b1162a45795e2e39"


def test_genome_features_are_what_they_claim():
	g = conftest.make_synth()
	kinds = {}
	for f in g.features:
		kinds[f.kind] = kinds.get(f.kind, 0) + 1
		s = g[f.contig][f.start:f.end]
		if f.kind == "N":
			assert set(s) == {"N"}
		elif f.kind == "lowcomp":
			assert s == (f.info * len(s))[:len(s)]
	assert kinds == {"repeat": 6, "inversion": 2, "lowcomp": 3, "N": 2}
	fams = {}
	for f in g.features:
		if f.kind == "repeat":
			s = g[f.contig][f.start:f.end]
			fams.setdefault(f.info[0], []).append(s if f.info[2] == "+" else H.revcomp(s))
	for copies in fams.values():  # divergence 0: identical copies
		assert len(set(copies)) == 1
	arms = [f for f in g.features if f.kind == "inversion"]
	assert arms[0].contig == arms[1].contig and arms[0].end < arms[1].start
	assert g[arms[0].contig][arms[0].start:arms[0].end] == H.revcomp(g[arms[1].contig][arms[1].start:arms[1].end])
	assert not g.is_unique(g.features[0].contig, g.features[0].start, g.features[0].start + 1)


def test_genes_have_gt_ag_introns():
	g = conftest.make_synth_tx()
	for gene in g.genes:
		assert 2 <= len(gene.exons) <= 6
		assert all(a < b for a, b in gene.exons)
		assert all(gene.exons[i][1] < gene.exons[i + 1][0] for i in range(len(gene.exons) - 1))
		tx = gene.transcript(g)
		assert len(tx) == gene.length
		for s, e in gene.introns:
			intron = g.fetch(gene.contig, s, e, gene.strand)
			assert intron[:2] == "GT" and intron[-2:] == "AG"
			assert 80 <= e - s <= 2000
		# no gene overlaps a repeat
		assert g.is_unique(gene.contig, gene.start, gene.end)
	assert set(gn.strand for gn in g.genes) == {"+", "-"}


def test_annotation_writers_hand_worked():
	gene = H.Gene("g1", "c1", "-", [(100, 150), (300, 420), (500, 510)])
	assert H.genes_bed12([gene]) == "c1\t100\t510\tg1\t0\t-\t100\t510\t0\t3\t50,120,10,\t0,200,400,\n"
	gtf = H.genes_gtf([gene]).splitlines()
	assert gtf[0] == 'c1\tsim\ttranscript\t101\t510\t.\t-\t.\tgene_id "g1"; transcript_id "g1";'
	assert gtf[2].split("\t")[2:5] == ["exon", "301", "420"]
	assert H.junctions_bed([gene]) == "c1\t150\t300\tg1\t0\t-\nc1\t420\t500\tg1\t0\t-\n"
	assert gene.introns == [(150, 300), (420, 500)] and gene.length == 180


# --- read simulator ------------------------------------------------------

def _apply_truth(ref, seg):
	"""Rebuild the read bases of a segment (reference orientation) from the reference, CIGAR and edits."""
	subs = dict((e[0], e[3]) for e in seg.edits if e[1] == "X")
	ins = dict((e[0], e[3]) for e in seg.edits if e[1] == "I")
	out, r = [], seg.pos
	for n, op in H.parse_cigar(seg.cigar):
		if op == "M":
			for _ in range(n):
				out.append(subs.get(r, ref[r].upper() if ref[r].upper() in "ACGT" else "N"))
				r += 1
		elif op == "I":
			out.append(ins[r])
		elif op in "DN":
			r += n
	assert r == seg.end
	return "".join(out)


def _segment_seq(read, seg):
	part = read.seq[seg.qstart:seg.qend]
	return part if seg.strand == "+" else H.revcomp(part)


def test_simulated_reads_match_their_truth():
	g = conftest.make_synth()
	reads = H.simulate_reads(g, 300, seed=2, sub_rate=0.02, indel_rate=0.002)
	assert len(set(r.name for r in reads)) == 300
	n_ind = 0
	for r in reads:
		assert len(r.seq) == len(r.qual) == 150
		assert H.cigar_query_span(r.cigar) == 150
		assert H.cigar_ref_span(r.cigar) == r.end - r.pos
		assert _apply_truth(g[r.contig], r.segments[0]) == r.fwd_seq
		assert r.seq == (r.fwd_seq if r.strand == "+" else H.revcomp(r.fwd_seq))
		# bases copied from an N run are not edits, but the SAM spec counts them in NM
		assert H.compute_nm(g[r.contig], r.pos, r.fwd_seq, r.cigar) <= r.nm + r.fwd_seq.count("N")
		n_ind += "I" in r.cigar or "D" in r.cigar
	assert n_ind > 10


def test_unique_reads_avoid_features():
	g = conftest.make_synth()
	for r in H.simulate_reads(g, 200, seed=3, unique=True):
		assert g.is_unique(r.contig, r.pos, r.end)


@pytest.mark.parametrize("profile,lo,hi", [("ont", 0.05, 0.10), ("hifi", 0.001, 0.01), ("clr", 0.05, 0.14)])
def test_long_read_profiles_have_their_error_rates(profile, lo, hi):
	g = conftest.make_synth()
	reads = H.simulate_long_reads(g, 12, profile=profile, seed=4)
	n_err = sum(r.nm for r in reads)
	n_base = sum(len(r) for r in reads)
	assert lo < n_err / float(n_base) < hi
	for r in reads:
		assert _apply_truth(g[r.contig], r.segments[0]) == r.fwd_seq
		lr, hr = H.PROFILES[profile]["length"]
		assert lr <= len(r) <= hr


def test_simulated_pairs_are_fr_with_the_drawn_insert():
	g = conftest.make_synth()
	pairs = H.simulate_pairs(g, 200, seed=4, insert_mean=350, insert_sd=40)
	for p in pairs:
		assert p.r1.name == p.r2.name == p.name
		assert p.r1.contig == p.r2.contig == p.contig
		left, right = (p.r1, p.r2) if p.strand == "+" else (p.r2, p.r1)
		assert left.strand == "+" and right.strand == "-"
		assert left.pos == p.frag_start and right.end == p.frag_end
		for r in p:
			assert _apply_truth(g[r.contig], r.segments[0]) == r.fwd_seq
	ins = [p.insert for p in pairs]
	mean = sum(ins) / float(len(ins))
	assert 330 < mean < 370


def test_spliced_reads_follow_their_gene():
	g = conftest.make_synth_tx()
	reads = H.simulate_spliced_reads(g, 120, seed=5, full_length_frac=0.3, min_len=200)
	n_full = n_anti = 0
	for r in reads:
		gene = r.gene
		seg = r.segments[0]
		assert _apply_truth(g[r.contig], seg) == r.fwd_seq
		assert r.introns == [x for x in gene.introns if r.exons[0][0] < x[0] and x[1] < r.exons[-1][1]]
		assert all(any(s <= a and b <= e for s, e in gene.exons) for a, b in r.exons)
		assert r.exons[0][1] - r.exons[0][0] >= 30 or len(r.exons) == 1
		assert r.exons[-1][1] - r.exons[-1][0] >= 30 or len(r.exons) == 1
		# the read strand on the genome is the gene strand, flipped for antisense reads
		assert r.strand == (gene.strand if r.ts_strand == "+" else "+-"[gene.strand == "+"])
		n_full += r.exons[0][0] == gene.start and r.exons[-1][1] == gene.end
		n_anti += r.ts_strand == "-"
		if r.nm == 0:  # an error-free sense read is a substring of the transcript
			tx = gene.transcript(g)
			assert (r.seq if r.ts_strand == "+" else H.revcomp(r.seq)) in tx
	assert 20 < n_full and 30 < n_anti < 90


def test_transcript_blocks_hand_worked():
	plus = H.Gene("g", "c", "+", [(100, 150), (300, 420), (500, 510)])  # transcript 50 + 120 + 10 = 180 bp
	assert H.transcript_blocks(plus, 0, 180) == [(100, 150), (300, 420), (500, 510)]
	assert H.transcript_blocks(plus, 40, 60) == [(140, 150), (300, 310)]
	assert H.transcript_blocks(plus, 40, 60, min_end_exon=11) == [(300, 310)]
	minus = H.Gene("g", "c", "-", [(100, 150), (300, 420), (500, 510)])  # 5' end at 510
	assert H.transcript_blocks(minus, 0, 15) == [(415, 420), (500, 510)]
	assert H.transcript_blocks(minus, 0, 15, min_end_exon=6) == [(500, 510)]


def test_fixed_length_spliced_reads_and_spliced_pairs():
	g = conftest.make_synth_tx()
	reads = H.simulate_spliced_reads(g, 40, seed=8, profile="short", length=100, min_end_exon=0)
	assert all(len(r) - sum(len(e[3]) - len(e[2]) for e in r.edits if e[1] in "ID") == 100 for r in reads)
	pairs = H.simulate_spliced_pairs(g, 60, seed=9, read_length=100, insert_mean=250, stranded=True)
	for p in pairs:
		assert p.r1.gene is p.r2.gene and p.contig == p.r1.contig == p.r2.contig
		assert p.r1.ts_strand == "+" and p.r2.ts_strand == "-" and p.r1.strand != p.r2.strand
		for r in p:
			assert _apply_truth(g[r.contig], r.segments[0]) == r.fwd_seq
			assert all(any(s <= a and b <= e for s, e in r.gene.exons) for a, b in r.exons)
	assert sum(1 for p in pairs if p.r1.introns or p.r2.introns) > 10


def test_chimeric_reads_join_distant_pieces():
	g = conftest.make_synth()
	reads = H.simulate_chimeric_reads(g, 20, seed=6, parts=2, min_distance=10000)
	for r in reads:
		s1, s2 = r.segments
		assert s1.qstart == 0 and s1.qend == s2.qstart and s2.qend == len(r.seq)
		assert s1.contig != s2.contig or abs(s1.pos - s2.pos) >= 10000
		for s in (s1, s2):
			assert _apply_truth(g[s.contig], s) == _segment_seq(r, s)


def test_contigs_with_inversions():
	g = conftest.make_synth()
	ctgs = H.simulate_contigs(g, 4, seed=7, inversions=1, length=(15000, 20000))
	for c in ctgs:
		assert len(c.segments) == 3 and c.qual is None
		segs = sorted(c.segments, key=lambda s: s.pos)
		assert segs[0].end == segs[1].pos and segs[1].end == segs[2].pos  # adjacent on the reference
		assert segs[0].strand == segs[2].strand != segs[1].strand
		assert [s.qstart for s in c.segments] == [0, c.segments[0].qend, c.segments[1].qend]
		for s in c.segments:
			assert _apply_truth(g[s.contig], s) == _segment_seq(c, s)


# --- FASTA/FASTQ writers and readers --------------------------------------

def test_fastx_round_trip_plain_and_gzip(tmp_path):
	g = conftest.make_synth()
	reads = H.simulate_reads(g, 10, seed=7)
	for name in ("r.fq", "r.fq.gz", "r.fa", "r.fa.gz"):
		p = H.write_fastx(str(tmp_path / name), reads)
		back = H.read_fastx(p)
		assert [x[0] for x in back] == [r.name for r in reads]
		assert [x[1] for x in back] == [r.seq for r in reads]
		if ".fq" in name:
			assert [x[2] for x in back] == [r.qual for r in reads]
	with open(str(tmp_path / "r.fq.gz"), "rb") as f:
		assert f.read(2) == b"\x1f\x8b"
	pairs = H.simulate_pairs(g, 5, seed=8)
	H.write_interleaved(str(tmp_path / "i.fq"), pairs, suffix=True)
	names = [x[0] for x in H.read_fastx(str(tmp_path / "i.fq"))]
	assert names == [n for p in pairs for n in (p.name + "/1", p.name + "/2")]
	H.write_fasta(str(tmp_path / "c.fa"), [("a", "ACGTACGT", None, "some comment")], line_width=3)
	assert H.read_text(str(tmp_path / "c.fa")) == ">a some comment\nACG\nTAC\nGT\n"


def test_gzip_output_is_byte_stable(tmp_path):
	recs = [("a", "ACGT", "IIII")]
	H.write_fastq(str(tmp_path / "x.fq.gz"), recs)
	H.write_fastq(str(tmp_path / "y.fq.gz"), recs)
	with open(str(tmp_path / "x.fq.gz"), "rb") as f, open(str(tmp_path / "y.fq.gz"), "rb") as g:
		assert f.read() == g.read()


# --- tag recomputation, worked by hand ------------------------------------

REF = "ACGTACGTAC"


@pytest.mark.parametrize("seq,cigar,nm,md,cs,cs_long", [
	("ACGTACGTAC", "10M", 0, "10", ":10", "=ACGTACGTAC"),
	("ACGAACGTAC", "10M", 1, "3T6", ":3*ta:6", "=ACG*ta=ACGTAC"),
	("TCGTACGTAC", "10M", 1, "0A9", "*at:9", "*at=CGTACGTAC"),
	("ACGTACGTAG", "10M", 1, "9C0", ":9*cg", "=ACGTACGTA*cg"),
	("ACGCGTAC", "3M2D5M", 2, "3^TA5", ":3-ta:5", "=ACG-ta=CGTAC"),
	("ACGCGTAA", "3M2D5M", 3, "3^TA4C0", ":3-ta:4*ca", "=ACG-ta=CGTA*ca"),
	("ACGTTTACGTAC", "4M2I6M", 2, "10", ":4+tt:6", "=ACGT+tt=ACGTAC"),
	("GGACGTACGTAC", "2S10M", 0, "10", ":10", "=ACGTACGTAC"),
	("ACGTACGTNC", "10M", 1, "8A1", ":8*an:1", "=ACGTACGT*an=C"),
	("ACGAGTAC", "3M2D5M", 3, "3^TA0C4", ":3-ta*ca:4", "=ACG-ta*ca=GTAC"),
	("ACGTTCGTAC", "4=1X5=", 1, "4A5", ":4*at:5", "=ACGT*at=CGTAC"),
	("acgtacgtac", "10M", 0, "10", ":10", "=ACGTACGTAC"),
])
def test_tag_recomputation_hand_worked(seq, cigar, nm, md, cs, cs_long):
	assert H.compute_nm(REF, 0, seq, cigar) == nm
	assert H.compute_md(REF, 0, seq, cigar) == md
	assert H.compute_cs(REF, 0, seq, cigar) == cs
	assert H.compute_cs(REF, 0, seq, cigar, long_form=True) == cs_long


def test_spliced_tag_recomputation():
	ref = "AAGGTCCCCCCCAGTTT"  # exon AAG, intron GTCCCCCCCAG (11 bp), exon TTT
	assert H.compute_cs(ref, 0, "AAGTTT", "3M11N3M") == ":3~gt11ag:3"
	assert H.compute_cs(ref, 0, "AAGTTT", "3M11N3M", long_form=True) == "=AAG~gt11ag=TTT"
	assert H.compute_md(ref, 0, "AAGTAT", "3M11N3M") == "4T1"
	assert H.compute_nm(ref, 0, "AAGTAT", "3M11N3M") == 1
	assert H.compute_mlen_blen(ref, 0, "AAGTAT", "3M11N3M") == (5, 6, 0)


def test_md_conventions():
	assert H.md_as_minimap2("149C0") == "149C" and H.md_as_minimap2("3^TA0") == "3^TA"
	assert H.md_as_minimap2("10") == "10" and H.md_as_minimap2("5C0A0") == "5C0A" and H.md_as_minimap2("0") == "0"
	# N in the read against N in the reference: minimap2's MD treats it as a match, calmd does not
	assert H.compute_md("ACNT", 0, "ACNT", "4M") == "4"
	assert H.compute_md("ACNT", 0, "ACNT", "4M", ambiguous_equal=False) == "2N1"
	# IUPAC codes are stored as N in the index
	assert H.compute_md("ACRT", 0, "ACGT", "4M") == "2N1"
	assert H.compute_cs("ACRT", 0, "ACGT", "4M") == ":2*ng:1"


def test_nm_counts_every_ambiguous_column():
	assert H.compute_nm("ACNT", 0, "ACNT", "4M") == 1
	assert H.compute_nm("ACGT", 0, "ACNT", "4M") == 1
	assert H.compute_cs("ACNT", 0, "ACNT", "4M") == ":4"
	assert H.compute_cs("ACNT", 0, "ACNT", "4M", long_form=True) == "=ACNT"


def test_mlen_blen_de_and_score_hand_worked():
	# 8 columns (one mismatch, one N), a 2-bp deletion and a 1-bp insertion
	ref = "ACGTACGTACGT"
	seq = "ACGANCATCG"  # A C G A(mismatch) N C [del GT] A [ins T] C G
	cigar = "6M2D1M1I2M"
	assert H.cigar_query_span(cigar) == 10 and H.cigar_ref_span(cigar) == 11
	mlen, blen, ambi = H.compute_mlen_blen(ref, 0, seq, cigar)
	assert (mlen, blen, ambi) == (7, 11, 1)  # 9 M columns minus 1 N = 8 (7 match, 1 mismatch) + 2 D + 1 I
	assert H.compute_nm(ref, 0, seq, cigar) == 5  # mismatch, N, two deleted bases, one inserted
	assert H.compute_de(ref, 0, seq, cigar) == pytest.approx(1.0 - 7.0 / (9 + 2))
	# +2 per match, -4 per mismatch, -1 for N, gaps min(4+2k, 24+k)
	assert H.cigar_score(ref, 0, seq, cigar) == 7 * 2 - 4 - 1 - (4 + 2 * 2) - (4 + 2 * 1)
	assert H.cigar_score("A" * 40, 0, "A" * 10, "5M30D5M") == 20 - (24 + 30)
	# -b: A<->G and C<->T mismatches cost the transition penalty, transversions still -B
	assert H.cigar_score("ACGTAC", 0, "GTCAAC", "6M", transition=1) == 2 * 2 - 1 - 1 - 4 - 4
	assert H.cigar_score("ACGTAC", 0, "GTCAAC", "6M", transition=0) == 2 * 2 - 4 * 4


def test_reference_minimizers_brute_force():
	# kmer_hash is a bijection on 2k-bit keys
	assert sorted(H.kmer_hash(x, 3) for x in range(64)) == list(range(64))
	# minimizers(): every k-mer whose hash is the smallest of some window of w consecutive k-mers
	def code(t):
		return sum(H.nt4(c) << 2 * (len(t) - 1 - j) for j, c in enumerate(t))
	rng = H.Rng(21)
	for trial in range(40):
		seq = rng.dna(rng.randint(20, 200))
		k, w = rng.choice([5, 7, 11]), rng.randint(1, 10)
		km = []
		for i in range(k - 1, len(seq)):
			s = seq[i - k + 1:i + 1]
			fwd, rev = code(s), code(H.revcomp(s))
			km.append((i, H.kmer_hash(min(fwd, rev), k), 0 if fwd < rev else 1))
		want = set()
		for st in range(max(1, len(km) - w + 1)):
			win = km[st:st + w]
			m = min(h for _, h, _ in win)
			want.update(x for x in win if x[1] == m)
		assert H.minimizers(seq, w, k) == sorted(want), (trial, seq, w, k)
	# a window with its minimum twice: kept by default, refused with ties=False
	with pytest.raises(ValueError):
		H.minimizers("ACGTTACGTTACGTTACGTT", 6, 5, ties=False)
	# an ambiguous base splits the k-mers into runs
	assert [len(r) for r in H.kmer_runs("ACGTACGNACGTACGTA", 5)] == [3, 5]
	# -H: k-mers of the homopolymer-compressed sequence, spans in the original
	q = H.kmer_runs("AAACCCGGGTTTA", 3, hpc=True)[0]
	assert [(x.x & 0xff, x.pos) for x in q] == [(9, 8), (9, 11), (7, 12)]


def test_mmi_reader_round_trip(index_of, tmp_path):
	seqs = [("c1", H.Rng(31).dna(3000)), ("c2", H.Rng(32).dna(1500))]
	fa = H.write_fasta(str(tmp_path / "g.fa"), seqs)
	parts = H.read_mmi(index_of(fa))
	assert len(parts) == 1
	p = parts[0]
	assert (p.w, p.k, p.flag, p.names, p.lengths) == (10, 15, 0, ["c1", "c2"], [3000, 1500])
	assert p.minimizers == H.expected_index([s for _, s in seqs], 10, 15)
	assert [bytes(H.nt4(c) for c in s) for _, s in seqs] == p.seq4
	h = H.read_mmi(index_of(fa), header_only=True)
	assert (h.w, h.k, h.names, h.lengths) == (10, 15, ["c1", "c2"], [3000, 1500])


def test_cigar_helpers():
	assert H.cigar_ref_span("5S10M2D3M1I4M3H") == 19
	assert H.cigar_query_span("5S10M2D3M1I4M3H") == 23
	assert H.cigar_query_span("5S10M2D3M1I4M3H", include_hard=True) == 26
	assert H.cigar_clips("5S10M2D3M1I4M3H") == (5, 3)
	assert H.strip_clips("5S10M2D3M1I4M3H") == "10M2D3M1I4M"
	assert H.parse_cigar("*") == []
	with pytest.raises(H.SamFormatError):
		H.parse_cigar("10M2")


# --- parsers -------------------------------------------------------------

SAM = (
	"@HD\tVN:1.6\tSO:unsorted\tGO:query\n"
	"@SQ\tSN:c1\tLN:100\n"
	"@RG\tID:g1\tSM:s\n"
	"@PG\tID:minimap2\tPN:minimap2\tVN:x\tCL:minimap2 -a a b\n"
	"q1\t99\tc1\t5\t60\t2S8M\t=\t40\t50\tAACGTACGTA\tIIIIIIIIII\tNM:i:0\tms:i:16\tAS:i:16\tnn:i:0\ttp:A:P\tcm:i:2\ts1:i:8\ts2:i:0\tde:f:0\tSA:Z:c1,50,+,4S6M,60,0;\tcs:Z::8\n"
	"q1\t147\tc1\t40\t60\t10M\t=\t5\t-50\tACGTACGTAC\t*\tCG:B:I,160,64\n"
	"q2\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n"
)


def test_sam_parser_types_and_flags():
	s = H.parse_sam(SAM)
	assert s.hd == {"VN": "1.6", "SO": "unsorted", "GO": "query"}
	assert s.sq == [("c1", 100)]
	assert s.rg == [{"ID": "g1", "SM": "s"}]
	assert s.pg[0]["CL"] == "minimap2 -a a b"
	r = s.records[0]
	assert (r.qname, r.flag, r.rname, r.pos, r.mapq, r.cigar, r.rnext, r.pnext, r.tlen) == ("q1", 99, "c1", 5, 60, "2S8M", "=", 40, 50)
	assert r.tags["tp"] == "P" and r.tags["de"] == 0.0 and r.tags["cs"] == ":8" and r.tags["SA"] == "c1,50,+,4S6M,60,0;"
	assert r.tag_types["de"] == "f" and r.tag_order[0] == "NM"
	assert r.is_paired and r.is_proper_pair and r.is_read1 and r.mate_reverse and not r.is_reverse and r.is_primary
	assert r.pos0 == 4 and r.ref_end0 == 12 and r.left_clip == 2 and r.unclipped_start0 == 2 and r.query_length == 10
	r2 = s.records[1]
	assert r2.is_reverse and r2.is_read2 and r2.segment == 2 and r2.strand == "-"
	assert r2.tags["CG"] == ("I", [160, 64])
	assert s.records[2].is_unmapped and s.records[2].segment == 0
	assert set(s.primaries()) == {("q1", 1), ("q1", 2), ("q2", 0)}
	assert H.strip_pg(SAM).count("@PG") == 0 and H.sam_body(SAM).count("\n") == 3


def test_sam_introns():
	r = H.parse_sam("q\t0\tc\t11\t60\t2S5M100N3M2D4M50N6M\t*\t0\t0\t" + "A" * 20 + "\t*\n").records[0]
	assert r.introns == [(15, 115), (124, 174)]


@pytest.mark.parametrize("bad", [
	"q\t0\tc1\t1\t60\t4M\t*\t0\t0\tACGT\t*\tNM:i:0\tNM:i:1\n",  # duplicate tag
	"q\t0\tc1\t1\t60\t5M\t*\t0\t0\tACGT\t*\n",  # CIGAR and SEQ lengths differ
	"q\t0\tc1\t1\t60\t4M\t*\t0\t0\tACGT\tIII\n",  # QUAL length
	"q\t0\tc1\t1\t300\t4M\t*\t0\t0\tACGT\t*\n",  # MAPQ range
	"q\t0\tc1\t1\t60\t4M\t*\t0\t0\tACGT\t*\tNM:i:x\n",  # bad integer
	"q\t0\tc1\t1\t60\t4M\t*\t0\t0\tACGT\n",  # too few fields
	"q\t0\tc1\t1\t60\t4M\t*\t0\t0\tACGT\t*",  # no final newline
	"q\t0\tc1\t1\t60\t4Q\t*\t0\t0\tACGT\t*\n",  # bad CIGAR op
	"q\t0\tc1\t1\t60\t4M\t*\t0\t0\tACGT\t*\tnot a tag\n",  # comment copied by -y that is not a tag
])
def test_sam_parser_rejects_malformed(bad):
	with pytest.raises(H.SamFormatError):
		H.parse_sam(bad)


def test_non_strict_parsing_keeps_untyped_fields():
	r = H.parse_sam("q\t0\tc1\t1\t60\t4M\t*\t0\t0\tACGT\t*\tNM:i:0\tnot a tag\n", strict=False).records[0]
	assert r.tags == {"NM": 0} and r.extra == ["not a tag"]
	p = H.parse_paf("q\t4\t0\t4\t+\tc\t10\t0\t4\t4\t4\t60\ttp:A:P\tBC:Z:x y\n")[0]
	assert p.tags["BC"] == "x y"


def test_paf_parser():
	text = ("r1\t150\t0\t150\t+\tc1\t1000\t10\t160\t148\t150\t60\ttp:A:P\tcm:i:5\ts1:i:100\tdv:f:0.0012\tcg:Z:150M\n"
		"r2\t100\t0\t0\t*\t*\t0\t0\t0\t0\t0\t0\trl:i:0\n")
	p = H.parse_paf(text)
	assert (p[0].qname, p[0].qlen, p[0].qstart, p[0].qend, p[0].strand, p[0].tname, p[0].tlen, p[0].tstart, p[0].tend,
		p[0].nmatch, p[0].alen, p[0].mapq) == ("r1", 150, 0, 150, "+", "c1", 1000, 10, 160, 148, 150, 60)
	assert p[0].is_primary and p[0].cigar == "150M" and p[0].tags["dv"] == pytest.approx(0.0012)
	assert p[1].is_unmapped and p[1].qlen == 100 and p[1].strand is None and p[1].tags == {"rl": 0}
	with pytest.raises(H.SamFormatError):
		H.parse_paf("r\t1\t0\t1\t?\tc\t1\t0\t1\t1\t1\t0\n")


def test_invariant_checkers_catch_errors():
	genome = {"c": "ACGTACGTACGTTTGACCAG"}
	good = "@SQ\tSN:c\tLN:20\nq\t0\tc\t3\t60\t2S6M\t*\t0\t0\tTTGTACGT\t*\tNM:i:0\tMD:Z:6\tcs:Z::6\n"
	assert H.sam_invariant_errors(H.parse_sam(good), genome, {"q": ["TTGTACGT"]}) == []
	for old, new in (("NM:i:0", "NM:i:1"), ("MD:Z:6", "MD:Z:5A0"), ("cs:Z::6", "cs:Z::5*ac"), ("\t3\t", "\t4\t")):
		assert H.sam_invariant_errors(H.parse_sam(good.replace(old, new)), genome, {"q": ["TTGTACGT"]}) != []
	paf = H.parse_paf("q\t8\t2\t8\t+\tc\t20\t2\t8\t6\t6\t60\ttp:A:P\tNM:i:0\tde:f:0\tcg:Z:6M\tcs:Z::6\n")
	assert H.paf_invariant_errors(paf, genome, {"q": "TTGTACGT"}) == []
	paf[0].nmatch = 5
	assert H.paf_invariant_errors(paf, genome, {"q": "TTGTACGT"}) != []
	# a long-form cs:Z may start with a substitution: recognised by its '=' runs, not its first character
	assert H.cs_is_long("*ag=TACGT") and H.cs_is_long("=GTACG*ta") and not H.cs_is_long("*ag:5")
	mis = "@SQ\tSN:c\tLN:20\nq\t0\tc\t3\t60\t2S6M\t*\t0\t0\tTTATACGT\t*\tcs:Z:*ga=TACGT\n"
	assert H.sam_invariant_errors(H.parse_sam(mis), genome, {"q": ["TTATACGT"]}) == []
	assert H.sam_invariant_errors(H.parse_sam(mis.replace("*ga=TACGT", "*ga=TACGA")), genome, {"q": ["TTATACGT"]}) != []
	paf = H.parse_paf("q\t8\t2\t8\t+\tc\t20\t2\t8\t5\t6\t60\ttp:A:P\tcg:Z:6M\tcs:Z:*ga=TACGT\n")
	assert H.paf_invariant_errors(paf, genome, {"q": "TTATACGT"}) == []


def test_pair_and_tlen_helpers():
	s = H.parse_sam(
		"p\t99\tc\t11\t60\t10M\t=\t41\t40\tACGTACGTAC\t*\n"
		"p\t147\tc\t41\t60\t10M\t=\t11\t-40\tACGTACGTAC\t*\n")
	assert H.pair_invariant_errors(s) == []
	a, b = s.records
	assert H.tlen_5prime(a, b) == 40 and H.tlen_5prime(b, a) == -40 and H.tlen_outer(a, b) == 40
	bad = H.parse_sam(
		"p\t99\tc\t11\t60\t10M\t=\t41\t40\tACGTACGTAC\t*\n"
		"p\t131\tc\t41\t60\t10M\t=\t11\t-40\tACGTACGTAC\t*\n")
	assert H.pair_invariant_errors(bad) != []


def test_pin_helper(tmp_path, monkeypatch):
	monkeypatch.setattr(H, "PINS_DIR", str(tmp_path))
	monkeypatch.setenv("MINIMAP2_TEST_UPDATE_PINS", "1")
	H.check_pinned("a.txt", "x\ny\n")
	H.check_pinned("b.txt.gz", "x\ny\n")
	monkeypatch.delenv("MINIMAP2_TEST_UPDATE_PINS")
	H.check_pinned("a.txt", "x\ny\n")
	H.check_pinned("b.txt.gz", "x\ny\n")
	with gzip.open(str(tmp_path / "b.txt.gz"), "rt") as f:
		assert f.read() == "x\ny\n"
	with pytest.raises(AssertionError, match="differs from pin"):
		H.check_pinned("a.txt", "x\nz\n")
	with pytest.raises(AssertionError, match="missing pin"):
		H.check_pinned("c.txt", "x\n")


# --- runner and build ----------------------------------------------------

def test_run_mm2_reports_exit_status(mm2_bin):
	r = H.run_mm2(["--no-such-option"])
	assert r.returncode != 0 and r.messages
	with pytest.raises(AssertionError, match="exited with"):
		r.check()
	assert H._has_thread_opt(["-t", "4"]) and H._has_thread_opt(["-t4"]) and not H._has_thread_opt(["-x", "sr"])


def test_processes_run_without_core_dumps_in_a_scratch_directory(tmp_path):
	code = "import os, resource; print(os.getcwd()); print(resource.getrlimit(resource.RLIMIT_CORE)[0]); print(os.environ['LC_ALL'])"
	r = H.run_process([sys.executable, "-c", code], check=True)
	cwd, core, lc = r.stdout.split()
	assert core in ("0", "1") and lc == "C"
	assert not os.path.exists(cwd)  # the temporary working directory is removed
	assert not os.path.realpath(cwd).startswith(os.path.realpath(minimap2_build.src_dir()) + os.sep)
	r = H.run_process(["sh", "-c", "kill -ABRT $$"])
	assert r.returncode == -6


def test_map_reads_adds_the_thread_count(synth, tmp_path, monkeypatch):
	reads = H.write_fasta(str(tmp_path / "r.fa"), H.simulate_reads(synth.genome, 3, seed=9, qual=None))
	monkeypatch.setenv("MINIMAP2_TEST_THREADS", "3")
	r = H.map_reads(synth.fasta, reads, ["-a"])
	assert r.sam().pg[0]["CL"].split()[-4:-2] == ["-t", "3"]
	r = H.map_reads(synth.fasta, reads, ["-a"], threads=2, preset="sr")
	assert r.sam().pg[0]["CL"].split()[1:3] == ["-x", "sr"] and " -t 2 " in r.sam().pg[0]["CL"]
	with pytest.raises(ValueError):
		H.map_reads(synth.fasta, reads, ["-t", "2"], threads=2)


def test_index_genome_is_cached(index_of, synth, index_dir):
	a = index_of(synth.genome, "sr")
	b = index_of(synth.fasta, "sr")
	assert a == b and os.path.dirname(a) == index_dir
	assert os.path.getsize(a) > 0 and index_of(synth.genome, "map-ont") != a


def test_source_hash_tracks_sources_only(tmp_path):
	src = minimap2_build.src_dir()
	dst = str(tmp_path / "src")
	os.makedirs(dst)
	for rel in minimap2_build.source_files(src):
		os.makedirs(os.path.dirname(os.path.join(dst, rel)), exist_ok=True)
		shutil.copy2(os.path.join(src, rel), os.path.join(dst, rel))
	h0 = minimap2_build.source_hash(dst)
	assert h0 == minimap2_build.source_hash(src)
	os.makedirs(os.path.join(dst, "tests"))
	with open(os.path.join(dst, "tests", "x.py"), "w") as f:
		f.write("# not a source file\n")
	assert minimap2_build.source_hash(dst) == h0
	with open(os.path.join(dst, "seed.c"), "a") as f:
		f.write("\n")
	assert minimap2_build.source_hash(dst) != h0


def test_make_args_change_the_build_key(monkeypatch):
	"""MINIMAP2_TEST_MAKE_ARGS reaches the build key (sse2only=1, asan=1 get their own cache entry);
	without it the key is unchanged."""
	monkeypatch.delenv("MINIMAP2_TEST_MAKE_ARGS", raising=False)
	k0 = minimap2_build.build_key()
	monkeypatch.setenv("MINIMAP2_TEST_MAKE_ARGS", "sse2only=1 'CFLAGS=-O1 -g'")
	assert minimap2_build.make_args() == ["sse2only=1", "CFLAGS=-O1 -g"]
	assert minimap2_build.build_key() != k0
	monkeypatch.setenv("MINIMAP2_TEST_MAKE_ARGS", "")
	assert minimap2_build.build_key() == k0


def test_out_of_tree_build_leaves_tree_clean(mm2_bin, sdust_bin, lite_bin, tree_before):
	"""Building the binary, libminimap2.a, `make extra` and a C driver (and mappy, when Cython is
	present) adds no build product to the source tree."""
	minimap2_build.c_driver("drv_pilot")
	if minimap2_build.mappy_missing() is None:
		minimap2_build.mappy_dir()
	src = minimap2_build.src_dir()
	for b in ([mm2_bin] if os.environ.get("MINIMAP2_BIN") is None else []) + [sdust_bin, lite_bin]:
		assert not os.path.realpath(b).startswith(os.path.realpath(src) + os.sep)
	assert minimap2_build.tree_products(src) == tree_before
