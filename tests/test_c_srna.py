"""Area C: short RNA-seq reads (`-x splice:sr`), junction jumps (`-j`, `--pass1`,
`--jump-min-match`, mm_jump_split in jump.c, mm_idx_jjump_read and
mm_idx_jump_get in index.c) and the two-pass workflow of README.md.

Truth comes from the simulated gene models of the standard transcript genome
(conftest.make_synth_tx) and from hand-built loci (helpers_c)."""

import pytest

import helpers_c as C
import minimap2_helpers as H


@pytest.fixture(scope="module")
def tx(synth_tx, tmp_path_factory):
	"""The transcript genome plus its junctions as BED6 (gff2junc style) and the overhang reads."""
	d = tmp_path_factory.mktemp("srna")
	g = synth_tx.genome
	jbed = H.write_text(str(d / "junc.bed"), H.junctions_bed(g.genes))
	oh = C.overhang_reads(g, g.genes)
	reads = [(n, s) for n, s, _ in oh] + [(n + "rc", H.revcomp(s)) for n, s, _ in oh]
	truth = dict((n, t) for n, _, t in oh)
	truth.update((n + "rc", t) for n, _, t in oh)
	fa = H.write_fasta(str(d / "overhang.fa"), reads)
	plain = H.map_reads(synth_tx.fasta, fa, ["-a"], preset="splice:sr")
	return dict(dir=d, genome=g, fasta=synth_tx.fasta, bed12=synth_tx.bed12, jbed=jbed, oh_fa=fa, oh_truth=truth,
		plain=plain.stdout, plain_ok=_jumped(plain.sam(), truth))


def _jumped(sam, truth):
	"""{read name: True if the primary alignment is exactly the two true blocks, unclipped}."""
	out = {}
	for (q, _), r in sam.primaries().items():
		ctg, blocks, _ = truth[q]
		out[q] = (not r.is_unmapped and r.rname == ctg and r.pos0 == blocks[0][0] and r.left_clip == 0 and r.right_clip == 0
			and r.introns == [(blocks[0][1], blocks[1][0])])
	return out


def test_splice_sr_pairs_against_truth(synth_tx, tmp_path):
	"""README `-ax splice:sr ref.fa r1.fq r2.fq` on stranded paired RNA-seq reads (read 1 sense): every
	mate on its gene's contig and genome strand, at least 90% with exactly the true introns; ts:A on
	the spliced ones gives the read's orientation to the transcript (cookbook.md), so '+' for read 1 and
	'-' for read 2; cs shows each intron as ~gt..ag (or ~ct..ac for a '-' gene); every pair is proper
	(--pairing=weak), and SEQ, CIGAR, NM, MD, cs and the mate fields agree with the genome."""
	g = synth_tx.genome
	pairs = H.simulate_spliced_pairs(g, 150, read_length=100, insert_mean=250, insert_sd=30, seed=301, stranded=True)
	f1, f2 = H.write_pairs(str(tmp_path / "r1.fq"), str(tmp_path / "r2.fq"), pairs)
	sam = H.map_reads(synth_tx.fasta, [f1, f2], ["-a", "--cs"], preset="splice:sr").sam()
	assert H.sam_invariant_errors(sam, g, dict((p.name, [p.r1.seq, p.r2.seq]) for p in pairs)) == []
	assert H.pair_invariant_errors(sam) == []
	assert not any(r.flag & 0x900 for r in sam.records)  # --secondary=no, no supplementary for these reads
	prim = sam.primaries()
	n = n_exact = n_spliced = 0
	for p in pairs:
		for seg, rd, ts in ((1, p.r1, "+"), (2, p.r2, "-")):
			r = prim[(p.name, seg)]
			n += 1
			assert not r.is_unmapped and r.rname == rd.contig and r.strand == rd.strand, (r.line[:120], rd)
			assert r.is_proper_pair and rd.ts_strand == ts
			if r.introns != rd.introns:
				continue
			n_exact += 1
			if r.introns:
				n_spliced += 1
				assert r.tags.get("ts") == ts, r.line[:200]
				sig = "~gt" if rd.gene.strand == "+" else "~ct"
				assert r.tags["cs"].count(sig) == len(r.introns), r.tags["cs"]
	assert n_exact >= 0.9 * n and n_spliced >= 50


def test_splice_sr_paf_names_and_ts(synth_tx, tmp_path):
	"""PAF of paired splice:sr reads with -c: every mate once as primary, named <name>/1 and <name>/2
	(fragment mode), cg:Z carries N for the introns and ts:A sits on exactly the spliced lines."""
	g = synth_tx.genome
	pairs = H.simulate_spliced_pairs(g, 40, read_length=100, seed=302, stranded=True)
	f1, f2 = H.write_pairs(str(tmp_path / "r1.fa"), str(tmp_path / "r2.fa"), pairs)
	paf = H.map_reads(synth_tx.fasta, [f1, f2], ["-c"], preset="splice:sr").paf()
	names = sorted(p.qname for p in paf if p.is_primary)
	assert names == sorted(["%s/%d" % (p.name, s) for p in pairs for s in (1, 2)])
	reads = dict(("%s/1" % p.name, p.r1.seq) for p in pairs)
	reads.update(("%s/2" % p.name, p.r2.seq) for p in pairs)
	assert H.paf_invariant_errors(paf, g, reads) == []
	assert any("N" in p.cigar for p in paf)
	for p in paf:
		assert ("ts" in p.tags) == ("N" in p.cigar), p.line[:200]


def test_splice_sr_single_end(synth_tx, tmp_path):
	"""README `-ax splice:sr ref.fa reads-se.fq`: 100 bp cDNA reads from one file with distinct names
	are single-end (no 0x1 flag) and at least 95% carry exactly their true introns on the true strand."""
	reads = H.simulate_spliced_reads(synth_tx.genome, 120, seed=303, profile="short", length=100)
	fq = H.write_fastq(str(tmp_path / "se.fq"), reads)
	sam = H.map_reads(synth_tx.fasta, fq, ["-a"], preset="splice:sr").sam()
	assert H.sam_invariant_errors(sam, synth_tx.genome, dict((r.name, [r.seq]) for r in reads)) == []
	prim = sam.primaries()
	assert len(prim) == len(reads) and all(r.flag & 0xC1 == 0 for r in sam.records)
	ok = sum(1 for r in reads if prim[(r.name, 0)].introns == r.introns and prim[(r.name, 0)].strand == r.strand)
	assert ok >= 0.95 * len(reads)


def test_jump_annotation_extends_read_ends(tx, index_of):
	"""-j: reads ending 4, 8 or 12 bases into the next exon, or starting that far before an exon end,
	in both orientations. With the annotation as BED12 or as 6-column junction BED, every read is
	extended across its intron: CIGAR = the two true blocks, no clipping. Without -j the short piece
	cannot be placed (at most 5% get the intron). A prebuilt index gives the same output."""
	truth = tx["oh_truth"]
	assert sum(tx["plain_ok"].values()) <= 0.05 * len(truth)
	out = {}
	for name, bed in (("bed12", tx["bed12"]), ("bed6", tx["jbed"])):
		res = H.map_reads(tx["fasta"], tx["oh_fa"], ["-a", "-j", bed], preset="splice:sr")
		got = _jumped(res.sam(), truth)
		assert len(got) == len(truth)
		missed = sorted(q for q, ok in got.items() if not ok)
		assert missed == [], (name, missed[:10])
		out[name] = H.sam_body(res.stdout)
	assert out["bed12"] == out["bed6"]
	mmi = index_of(tx["fasta"], "splice:sr")
	via_idx = H.map_reads(mmi, tx["oh_fa"], ["-a", "-j", tx["bed12"]], preset="splice:sr")
	assert H.sam_body(via_idx.stdout) == out["bed12"]


def test_jump_with_a_mismatch_near_the_junction(tx):
	"""-j on overhang reads (8 or 12 bases into the next exon, or starting that far before an exon end;
	both orientations) with one substitution 2 bases from the junction on the long side: every read
	is still extended across its annotated intron, and NM, the CIGAR's span and cs:Z (the substitution
	and the intron) agree with the genome. (A substitution 4 bases from the junction can be counted
	twice in NM: base behaviour; see Rules.)"""
	other = {"A": "C", "C": "A", "G": "T", "T": "G"}
	reads, truth = [], {}
	for n, seq, t in C.overhang_reads(tx["genome"], tx["genome"].genes, overhangs=(8, 12)):
		k = t[2]
		p = 100 - k - 2 if "_end" in n else k + 1
		m = seq[:p] + other[seq[p]] + seq[p + 1:]
		reads += [(n + "m", m), (n + "mrc", H.revcomp(m))]
		truth[n + "m"] = truth[n + "mrc"] = t
	fa = H.write_fasta(str(tx["dir"] / "overhang_mm.fa"), reads)
	sam = H.map_reads(tx["fasta"], fa, ["-a", "--cs", "-j", tx["jbed"]], preset="splice:sr").sam()
	assert H.sam_invariant_errors(sam, tx["genome"], dict((n, [s]) for n, s in reads)) == []
	got = _jumped(sam, truth)
	assert len(got) == len(reads) and all(got.values()), sorted(q for q, ok in got.items() if not ok)[:10]
	assert all(r.tags["NM"] == 1 and r.tags["cs"].count("*") == 1 for r in sam.records if r.is_primary)


@pytest.mark.parametrize("min_match", [5, 8, 9, 13])
def test_jump_min_match(tx, min_match):
	"""--jump-min-match INT, the minimum matching length to create a jump: a read is extended across
	the annotated intron when its overhang (4, 8 or 12 bases) is at least INT; a shorter overhang
	keeps the intron only if base alignment found it without -j."""
	truth = tx["oh_truth"]
	res = H.map_reads(tx["fasta"], tx["oh_fa"], ["-a", "-j", tx["jbed"], "--jump-min-match", str(min_match)], preset="splice:sr")
	got = _jumped(res.sam(), truth)
	wrong = sorted(q for q, ok in got.items() if ok != (truth[q][2] >= min_match or tx["plain_ok"][q]))
	assert wrong == [], wrong[:10]
	assert any(not ok for ok in got.values()) == (min_match > 4)


@pytest.mark.parametrize("score,used", [(4, False), (5, True), (6, True)])
def test_pass1_score_threshold(tx, tmp_path, score, used):
	"""--pass1 FILE: junctions in --write-junc's 6-column BED; rows with scores lower than 5 are
	ignored. All true introns at one score: with score >= 5 every overhang read is extended; with
	score 4 the output is that of no --pass1 at all."""
	g = tx["genome"]
	rows = "".join("%s\t%d\t%d\tx\t%d\t%s\n" % (gg.contig, s, e, score, gg.strand) for gg in g.genes for s, e in gg.introns)
	p1 = H.write_text(str(tmp_path / "pass1.bed"), rows)
	res = H.map_reads(tx["fasta"], tx["oh_fa"], ["-a", "--pass1", p1], preset="splice:sr")
	if used:
		assert all(_jumped(res.sam(), tx["oh_truth"]).values())
	else:
		assert H.sam_body(res.stdout) == H.sam_body(tx["plain"])


def test_two_pass_workflow(tx, synth_tx, tmp_path):
	"""README's two-pass alignment: `-x splice:sr --write-junc` on paired reads, then `-ax splice:sr
	--pass1=junc.bed`. Pass 1 writes 6-column BED (minimap2.1): every junction is a true intron with
	the gene's strand and score 6 (GT 3 + AG 3), named after its read. In pass 2 an overhang read is
	extended across its intron when pass 1 reported that intron, and otherwise only if base
	alignment finds it alone."""
	g = tx["genome"]
	pairs = H.simulate_spliced_pairs(g, 250, read_length=100, insert_mean=250, insert_sd=30, seed=304)
	f1, f2 = H.write_pairs(str(tmp_path / "r1.fa"), str(tmp_path / "r2.fa"), pairs)
	wj = H.map_reads(tx["fasta"], [f1, f2], ["--write-junc"], preset="splice:sr")
	true_introns = dict(((gg.contig, s, e), gg.strand) for gg in g.genes for s, e in gg.introns)
	names = set(p.name for p in pairs)
	found = set()
	lines = wj.stdout.splitlines()
	assert lines
	for line in lines:
		f = line.split("\t")
		assert len(f) == 6, line
		key = (f[0], int(f[1]), int(f[2]))
		assert key in true_introns and f[5] == true_introns[key] and f[4] == "6" and f[3] in names, line
		found.add(key)
	assert len(found) >= 0.8 * len(true_introns)
	p1 = H.write_text(str(tmp_path / "junc.bed"), wj.stdout)
	got = _jumped(H.map_reads(tx["fasta"], tx["oh_fa"], ["-a", "--pass1", p1], preset="splice:sr").sam(), tx["oh_truth"])
	for q, ok in got.items():
		ctg, blocks, _ = tx["oh_truth"][q]
		assert ok == ((ctg, blocks[0][1], blocks[1][0]) in found or tx["plain_ok"][q]), q


@pytest.fixture(scope="module")
def alt(tmp_path_factory):
	d = tmp_path_factory.mktemp("alt")
	loc = C.alt_acceptor_locus(311)
	loc.write(str(d / "alt.fa"))
	ctg = loc.genome["chr1"]
	reads = []
	for k in (6, 10):
		e1s, e1e = loc.e1
		seq = ctg[e1e - (100 - k):e1e] + ctg[loc.ja[1]:loc.ja[1] + k]
		assert seq[-k:] == ctg[loc.jb[1]:loc.jb[1] + k]  # matches both acceptors
		reads.append(("alt%d" % k, seq))
		reads.append(("alt%drc" % k, H.revcomp(seq)))
	fa = H.write_fasta(str(d / "reads.fa"), reads)

	def bed(name, ivs, score=6):
		return H.write_text(str(d / name), "".join("chr1\t%d\t%d\tx\t%d\t+\n" % (s, e, score) for s, e in ivs))
	return dict(loc=loc, reads=fa, names=[n for n, _ in reads], bed=bed)


def _intron_of(alt, args):
	sam = H.map_reads(alt["loc"].fasta, alt["reads"], ["-a"] + args, preset="splice:sr").sam()
	prim = sam.primaries()
	return set(tuple(prim[(n, 0)].introns) for n in alt["names"])


def test_jump_prefers_annotation_over_pass1(alt):
	"""minimap2.1: when both -j and --pass1 are present, junctions in -j are preferred over those in
	--pass1 when there is ambiguity. A read whose last 6 or 10 bases match the start of two alternative
	3' exons is extended to the junction given by -j, whichever of the two that is; with one file
	alone, to the single junction it lists. Also records base behaviour: with both junctions in the
	same file the read is not extended at all."""
	ja, jb = alt["loc"].ja, alt["loc"].jb
	A, B = alt["bed"]("a.bed", [ja]), alt["bed"]("b.bed", [jb])
	both = alt["bed"]("ab.bed", [ja, jb])
	assert _intron_of(alt, ["-j", A]) == {(ja,)}
	assert _intron_of(alt, ["--pass1", B]) == {(jb,)}
	assert _intron_of(alt, ["-j", A, "--pass1", B]) == {(ja,)}
	assert _intron_of(alt, ["-j", B, "--pass1", A]) == {(jb,)}
	assert _intron_of(alt, ["-j", both]) == {()}
	assert _intron_of(alt, ["--pass1", both]) == {()}
	assert _intron_of(alt, ["-j", both, "--pass1", A]) == {()}


def test_jump_needs_splice_mode(tx):
	"""mm_jump_split runs only for spliced alignment: `-x sr -j anno.bed` writes the same records as
	`-x sr`, and junctions given to -j never appear as introns there."""
	a = H.map_reads(tx["fasta"], tx["oh_fa"], ["-a", "--frag=no"], preset="sr")
	b = H.map_reads(tx["fasta"], tx["oh_fa"], ["-a", "--frag=no", "-j", tx["bed12"]], preset="sr")
	assert H.sam_body(a.stdout) == H.sam_body(b.stdout)
	assert all(not r.introns for r in b.sam().records)


@pytest.mark.pinned
def test_splice_sr_pairs_with_annotation_pinned(synth_tx, tmp_path):
	"""Pins the current output: `-ax splice:sr -j anno.bed12 ref.fa r1.fa r2.fa` (README) on 120
	simulated pairs, SAM without @PG, byte for byte."""
	pairs = H.simulate_spliced_pairs(synth_tx.genome, 120, read_length=100, seed=305, qual=None)
	f1, f2 = H.write_pairs(str(tmp_path / "r1.fa"), str(tmp_path / "r2.fa"), pairs)
	res = H.map_reads(synth_tx.fasta, [f1, f2], ["-a", "-j", synth_tx.bed12], preset="splice:sr")
	H.check_pinned("c/splice_sr_pairs_j.sam.gz", H.strip_pg(res.stdout))
