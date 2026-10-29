"""Whole-program checks: the README's MT-human/MT-orang example, simulated ONT-like,
HiFi-like, spliced and paired-end reads against their truth, one
libminimap2 API case through a C driver, one mappy case and one
paftools.js case."""

import os
import shutil
import subprocess

import pytest

import minimap2_helpers as H


@pytest.fixture(scope="module")
def mt_sam(mt):
	"""README: `minimap2 -a test/MT-human.fa test/MT-orang.fa`, plus --cs. (--cs with --MD writes only
	MD, so tests that need both run minimap2 twice.)"""
	return H.map_reads(mt.human, mt.orang, ["-a", "--cs"])


def test_readme_example_sam(mt, mt_sam, mm2_bin):
	"""Header from the target; one primary alignment of MT_orang on the forward strand covering most of
	it; SEQ is the query; CIGAR spans the query; NM, MD and cs equal their recomputation; AS equals the
	CIGAR's score under the default scoring (-A2 -B4 -O4,24 -E2,1)."""
	sam = mt_sam.sam()
	assert "cs" in sam.records[0].tags
	assert sam.hd.get("VN") == "1.6"
	assert sam.sq == [("MT_human", 16569)]
	version = H.run_mm2(["--version"], check=True).stdout.strip()
	assert len(sam.pg) == 1 and sam.pg[0]["ID"] == "minimap2" and sam.pg[0]["PN"] == "minimap2"
	assert sam.pg[0]["VN"] == version and sam.pg[0]["CL"].startswith("minimap2 ")
	orang = mt.orang_seqs["MT_orang"]
	assert H.sam_invariant_errors(sam, mt.human_seqs, {"MT_orang": [orang]}) == []
	prim = sam.primaries()
	assert list(prim) == [("MT_orang", 0)]
	r = prim[("MT_orang", 0)]
	assert r.rname == "MT_human" and not r.is_reverse and 30 <= r.mapq <= 60
	aligned = len(orang) - r.left_clip - r.right_clip
	assert aligned >= 0.95 * len(orang)
	s = orang[r.left_clip:len(orang) - r.right_clip]
	assert r.tags["AS"] == H.cigar_score(mt.human_seqs["MT_human"], r.pos0, orang, r.cigar)
	assert 0.05 < H.compute_de(mt.human_seqs["MT_human"], r.pos0, orang, r.cigar) < 0.25
	assert r.tags["nn"] == 0 and r.tags["tp"] == "P" and len(s) == H.cigar_query_span(H.strip_clips(r.cigar))
	md = H.map_reads(mt.human, mt.orang, ["-a", "--MD"]).sam()
	assert [x.tags["MD"] for x in md.records] == [H.md_as_minimap2(H.compute_md(mt.human_seqs["MT_human"], x.pos0,
		orang, x.cigar)) for x in md.records]
	assert [x.cigar for x in md.records] == [x.cigar for x in sam.records]


def test_readme_prebuilt_index_gives_the_same_alignments(mt, mt_sam, tmp_path):
	"""README: `minimap2 -x map-ont -d MT-human-ont.mmi test/MT-human.fa` then
	`minimap2 -a MT-human-ont.mmi test/MT-orang.fa`. map-ont is the default preset, so mapping through
	the prebuilt index must give the same SAM records and @SQ lines as indexing on the fly."""
	mmi = str(tmp_path / "MT-human-ont.mmi")
	r = H.run_mm2(["-x", "map-ont", "-d", mmi, mt.human], check=True)
	assert r.stdout == "" and os.path.getsize(mmi) > 0
	with open(mmi, "rb") as f:
		assert f.read(4) == b"MMI\x02"  # MM_IDX_MAGIC
	via_idx = H.map_reads(mmi, mt.orang, ["-a", "--cs"])
	assert H.strip_pg(via_idx.stdout) == H.strip_pg(mt_sam.stdout)


def test_readme_example_paf(mt):
	"""Default output is PAF with 12 columns and the man page's tags; with -c it also carries cg:Z, and
	the coordinates, columns 10/11, NM and de agree with the alignment recomputed from the sequences."""
	orang = mt.orang_seqs["MT_orang"]
	paf = H.map_reads(mt.human, mt.orang).paf()
	assert len(paf) == 1
	p = paf[0]
	assert (p.qname, p.qlen, p.strand, p.tname, p.tlen) == ("MT_orang", len(orang), "+", "MT_human", 16569)
	assert set(p.tags) == {"tp", "cm", "s1", "s2", "dv", "rl"} and p.is_primary
	assert H.paf_invariant_errors(paf, mt.human_seqs, {"MT_orang": orang}) == []
	pc = H.map_reads(mt.human, mt.orang, ["-c", "--cs"]).paf()
	assert len(pc) == 1 and pc[0].cigar is not None
	assert {"NM", "ms", "AS", "nn", "tp", "cm", "s1", "s2", "de", "rl", "cg", "cs"} <= set(pc[0].tags)
	assert H.paf_invariant_errors(pc, mt.human_seqs, {"MT_orang": orang}) == []


@pytest.fixture(scope="module")
def ont_reads(synth):
	return H.simulate_long_reads(synth.genome, 40, profile="ont", seed=101, unique=True)


def test_ont_like_reads_land_at_truth(synth, index_of, ont_reads, tmp_path):
	"""map-ont through a prebuilt index: every read from unique sequence maps primarily to its true
	contig and strand with its unclipped start within 30 bp of the truth; SAM records are consistent
	with the genome (NM and the long-form cs recomputed)."""
	fq = H.write_fastq(str(tmp_path / "ont.fq"), ont_reads)
	sam = H.map_reads(index_of(synth.genome, "map-ont"), fq, ["-a", "--cs=long"], preset="map-ont").sam()
	assert all(r.tags["cs"].startswith("=") for r in sam.records if not r.is_unmapped)
	assert H.sam_invariant_errors(sam, synth.genome, dict((r.name, [r.seq]) for r in ont_reads)) == []
	prim = sam.primaries()
	assert len(prim) == len(ont_reads)
	for r in ont_reads:
		rec = prim[(r.name, 0)]
		err = H.placement_error(rec, r)
		assert err is not None and err <= 30, (rec.line[:150], r)
		assert rec.mapq > 0 and rec.query_length == len(r)


def test_hifi_like_reads_paf(synth, index_of, tmp_path):
	"""map-hifi, PAF with -c --cs: one primary line per read at its true contig and strand within 10 bp
	at both ends; cg, NM, cs, de and columns 10/11 agree with the recomputed alignment."""
	reads = H.simulate_long_reads(synth.genome, 30, profile="hifi", seed=102, unique=True)
	fq = H.write_fastq(str(tmp_path / "hifi.fq"), reads)
	paf = H.map_reads(index_of(synth.genome, "map-hifi"), fq, ["-c", "--cs"], preset="map-hifi").paf()
	assert H.paf_invariant_errors(paf, synth.genome, dict((r.name, r.seq) for r in reads)) == []
	truth = dict((r.name, r) for r in reads)
	pri = [p for p in paf if p.is_primary]
	assert sorted(p.qname for p in pri) == sorted(truth)
	for p in pri:
		err = H.paf_placement_error(p, truth[p.qname])
		assert err is not None and err <= 10, (p.line[:150], truth[p.qname])


def test_spliced_reads_recover_their_introns(synth_tx, tmp_path):
	"""splice:hq on HiFi-like cDNA reads (both orientations): every read maps primarily to its gene's
	contig and the true genome strand; at least 95% recover exactly the true introns; on spliced
	reads ts:A gives the read's orientation relative to the transcript (cookbook.md); cs shows each intron as ~gt..ag on
	the transcript strand."""
	reads = H.simulate_spliced_reads(synth_tx.genome, 60, seed=103, profile="hifi")
	fq = H.write_fastq(str(tmp_path / "tx.fq"), reads)
	sam = H.map_reads(synth_tx.fasta, fq, ["-a", "--cs"], preset="splice:hq").sam()
	assert H.sam_invariant_errors(sam, synth_tx.genome, dict((r.name, [r.seq]) for r in reads)) == []
	prim = sam.primaries()
	n_exact = 0
	for r in reads:
		rec = prim[(r.name, 0)]
		assert not rec.is_unmapped and rec.rname == r.contig and rec.strand == r.strand, (rec.line[:150], r)
		if r.introns:  # cookbook.md: ts:A is inferred from GT-AG, so only spliced reads carry it
			assert rec.tags.get("ts") == r.ts_strand
		n_exact += rec.introns == r.introns
		if rec.introns == r.introns and r.introns:
			want = "~gt" if r.gene.strand == "+" else "~ct"
			assert rec.tags["cs"].count(want) == len(r.introns)
	assert n_exact >= 0.95 * len(reads)


def test_paired_short_reads(synth, index_of, tmp_path):
	"""sr on FR pairs from unique sequence: both mates at the truth (unclipped start within 3 bp), flag
	0x2 on both, TLEN the 5'-to-5' length and within 10 bp of the simulated insert; mate fields agree."""
	pairs = H.simulate_pairs(synth.genome, 200, seed=104, insert_mean=300, insert_sd=30, unique=True)
	f1, f2 = H.write_pairs(str(tmp_path / "r1.fq"), str(tmp_path / "r2.fq"), pairs)
	sam = H.map_reads(index_of(synth.genome, "sr"), [f1, f2], ["-a", "--MD"], preset="sr").sam()
	assert all("MD" in r.tags for r in sam.records if not r.is_unmapped)
	assert H.sam_invariant_errors(sam, synth.genome, dict((p.name, [p.r1.seq, p.r2.seq]) for p in pairs)) == []
	assert H.pair_invariant_errors(sam) == []
	prim = sam.primaries()
	for p in pairs:
		a, b = prim[(p.name, 1)], prim[(p.name, 2)]
		for rec, read in ((a, p.r1), (b, p.r2)):
			err = H.placement_error(rec, read)
			assert err is not None and err <= 3, (rec.line[:150], read)
		assert a.is_proper_pair and b.is_proper_pair
		assert a.tlen == H.tlen_5prime(a, b)
		assert abs(abs(a.tlen) - p.insert) <= 10


def test_c_api_map_truth(synth, tmp_path):
	"""mm_idx_reader_open/read, mm_mapopt_update, mm_tbuf_init and mm_map() through a C driver, as
	example.c uses them: each read's primary hit is on the true contig and strand within 20 bp, and
	every hit's CIGAR spans its query and reference intervals."""
	reads = H.simulate_long_reads(synth.genome, 30, profile="hifi", seed=105, unique=True, qual=None)
	fa = H.write_fasta(str(tmp_path / "r.fa"), reads)
	tsv = str(tmp_path / "truth.tsv")
	with open(tsv, "w") as f:
		for r in reads:
			f.write("%s\t%s\t%d\t%d\t%s\n" % (r.name, r.contig, r.pos, r.end, r.strand))
	H.assert_c_case("drv_pilot", "map_truth", synth.fasta, "map-hifi", fa, tsv, 20)


def test_c_api_idx_str():
	"""mm_idx_str() builds an index from strings; mm_idx_index_name, mm_idx_name2id and mm_idx_getseq read it back."""
	H.assert_c_case("drv_pilot", "idx_str")


@pytest.mark.mappy
def test_mappy_matches_the_command_line(mappy, mt):
	"""python/README.rst's example: Aligner on MT-human.fa, seq(), revcomp(), fastx_read() and map().
	Every Alignment equals the corresponding `minimap2 -c --cs` (and `-c --MD`) PAF line (default options)."""
	a = mappy.Aligner(mt.human)
	assert a and a.seq_names == ["MT_human"] and a.k == 15 and a.w == 10 and a.n_seq == 1
	human = mt.human_seqs["MT_human"]
	s = a.seq("MT_human", 100, 200)
	assert s == human[100:200] and mappy.revcomp(s) == H.revcomp(s)
	recs = list(mappy.fastx_read(mt.orang))
	assert [(n, len(q)) for n, q, _ in recs] == [("MT_orang", 16499)] and recs[0][2] is None
	hits = list(a.map(recs[0][1], cs=True, MD=True))
	paf = H.map_reads(mt.human, mt.orang, ["-c", "--cs"]).paf()
	paf_md = H.map_reads(mt.human, mt.orang, ["-c", "--MD"]).paf()
	assert len(hits) == len(paf) == len(paf_md) >= 1
	for h, p, pm in zip(hits, paf, paf_md):
		assert (h.ctg, h.ctg_len, h.r_st, h.r_en, h.q_st, h.q_en) == (p.tname, p.tlen, p.tstart, p.tend, p.qstart, p.qend)
		assert ("+" if h.strand == 1 else "-") == p.strand and h.mapq == p.mapq and h.is_primary == p.is_primary
		assert (h.mlen, h.blen, h.NM) == (p.nmatch, p.alen, p.tags["NM"])
		assert h.cigar_str == p.cigar and h.cs == p.tags["cs"] and h.MD == pm.tags["MD"]


@pytest.mark.k8
def test_paftools_sam2paf_agrees_with_paf_output(k8, mt, mt_sam):
	"""misc/paftools.js sam2paf on the README example's SAM gives the same 12 columns, cg:Z and cs:Z as
	minimap2's own `-c --cs` PAF (the MT sequences hold no ambiguous base, so columns 10/11 agree too)."""
	r = H.run_paftools(["sam2paf", "-"], stdin=H.map_reads(mt.human, mt.orang, ["-a", "--cs"]).stdout)
	conv = H.parse_paf(r.stdout)
	own = H.map_reads(mt.human, mt.orang, ["-c", "--cs"]).paf()
	assert len(conv) == len(own) >= 1
	for c, o in zip(conv, own):
		assert c.line.split("\t")[:12] == o.line.split("\t")[:12]
		assert c.tags["cg"] == o.tags["cg"] and c.tags["cs"] == o.tags["cs"] and c.tags["NM"] == o.tags["NM"]


@pytest.mark.external
def test_samtools_reads_the_output(mt_sam, tmp_path):
	samtools = shutil.which("samtools")
	if samtools is None:
		pytest.skip("samtools not on PATH")
	p = str(tmp_path / "x.sam")
	with open(p, "w") as f:
		f.write(mt_sam.stdout)
	r = subprocess.run([samtools, "view", "-c", p], stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60)
	assert r.returncode == 0, r.stderr
	assert int(r.stdout) == len(mt_sam.sam().records)


@pytest.mark.external
def test_pysam_reads_the_output(mt_sam, tmp_path):
	pysam = pytest.importorskip("pysam")
	p = str(tmp_path / "x.sam")
	with open(p, "w") as f:
		f.write(mt_sam.stdout)
	with pysam.AlignmentFile(p, "r") as f:
		recs = list(f)
	ours = mt_sam.sam().records
	assert [(r.query_name, r.flag, r.reference_start, r.cigarstring) for r in recs] == [(r.qname, r.flag, r.pos0, r.cigar) for r in ours]
