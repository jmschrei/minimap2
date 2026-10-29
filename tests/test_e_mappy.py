"""mappy, the Python binding, built from the sources under test.

Expected values come from python/README.rst, the simulator's truth, the
genome (cs, MD, NM, mlen and blen recomputed), hand-worked end-extension
arithmetic under minimap2's scoring model, and the command line: every
Aligner argument mirrors a minimap2 option, so mappy's alignments must equal
`minimap2 -c --secondary=yes` with that option (and `--no-hash-name` when no
read name is passed to map()).

python/README.rst documents mate 2 of a pair as having positions on the query;
mappy reports them on its reverse complement, so the pair tests compare mate 2's query
interval and cs only for hits that cover the whole mate.
"""

import gzip
import os
import re
import sys

import pytest

import helpers_e as E
import minimap2_helpers as H

pytestmark = pytest.mark.mappy


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def wd(tmp_path_factory):
	return tmp_path_factory.mktemp("e_mappy")


@pytest.fixture(scope="module")
def hifi(synth, wd):
	"""24 HiFi-like reads drawn from anywhere in the standard genome (repeats included), both strands."""
	reads = H.simulate_long_reads(synth.genome, 24, profile="hifi", seed=501)
	return reads, H.write_fastq(str(wd / "hifi.fq"), reads)


@pytest.fixture(scope="module")
def hifi_unique(synth, wd):
	"""20 HiFi-like reads from unique sequence, with a few N bases (n_rate) to exercise ambiguous columns."""
	reads = H.simulate_long_reads(synth.genome, 20, profile="hifi", seed=502, unique=True, n_rate=0.002)
	return reads, H.write_fastq(str(wd / "hifi_unique.fq"), reads)


@pytest.fixture(scope="module")
def default_aligner(mappy, synth):
	a = mappy.Aligner(synth.fasta)
	assert a
	return a


def _pos(h, read):
	"""Distance between a hit's reference interval and the read's true interval (None if contig/strand differ)."""
	if h.ctg != read.contig or (h.strand > 0) != (read.strand == "+"):
		return None
	return max(abs(h.r_st - read.pos), abs(h.r_en - read.end))


# ---------------------------------------------------------------------------
# Module-level names
# ---------------------------------------------------------------------------


def test_version_matches_the_sources(mappy):
	"""mappy.__version__ equals setup.py's version and the release part of MM_VERSION in minimap.h."""
	setup_ver = re.search(r"version\s*=\s*'([^']+)'", E.source_text("setup.py")).group(1)
	mm_ver = re.search(r'#define MM_VERSION "([^"]+)"', E.source_text("minimap.h")).group(1)
	assert mappy.__version__ == setup_ver == mm_ver.split("-")[0]


def test_revcomp(mappy):
	"""python/README.rst: revcomp recognizes IUB codes, preserves case and complements U to A."""
	s = "ACGTNacgtnRYKMSWBDHVU"
	comp = dict(zip("ACGTNacgtnRYKMSWBDHVU", "TGCANtgcanYRMKSWVHDBA"))
	assert mappy.revcomp(s) == "".join(comp[c] for c in reversed(s))
	assert mappy.revcomp(b"AACG") == "CGTT"
	assert mappy.revcomp("") == ""
	assert mappy.revcomp(mappy.revcomp("GATTACAn")) == "GATTACAn"


def test_verbose_get_and_set(mappy):
	"""verbose() returns the current level; verbose(v) sets it and returns the new level."""
	old = mappy.verbose()
	try:
		assert mappy.verbose(0) == 0
		assert mappy.verbose() == 0
		assert mappy.verbose(3) == 3
		assert mappy.verbose() == 3
	finally:
		mappy.verbose(old)
	assert mappy.verbose() == old


# ---------------------------------------------------------------------------
# fastx_read
# ---------------------------------------------------------------------------


def test_fastx_read_fasta_fastq_gzip_and_comments(mappy, wd):
	"""fastx_read yields (name, seq, qual) per record, qual None for FASTA; multi-line FASTA is joined; gzip is
	read transparently; read_comment=True adds the comment (None when absent)."""
	recs = [("r1", "ACGTACGTAC" * 13, "I" * 130, "some comment"), ("r2", "GGGTTT", "#####!", None), ("r3", "A", "5", "x y")]
	fq = H.write_fastq(str(wd / "fx.fq"), recs)
	fa = H.write_fasta(str(wd / "fx.fa.gz"), recs, line_width=7)
	with gzip.open(fa, "rt") as f:
		assert f.read().count("\n") > 3 * 2  # several lines per record
	assert list(mappy.fastx_read(fq)) == [(n, s, q) for n, s, q, _ in recs]
	assert list(mappy.fastx_read(fq, read_comment=True)) == [(n, s, q, c) for n, s, q, c in recs]
	assert list(mappy.fastx_read(fa)) == [(n, s, None) for n, s, _, _ in recs]
	assert list(mappy.fastx_read(fa, read_comment=True)) == [(n, s, None, c) for n, s, _, c in recs]
	empty = H.write_text(str(wd / "empty.fa"), "")
	assert list(mappy.fastx_read(empty)) == []


def test_fastx_read_stdin(mappy, wd):
	"""python/README.rst and the binding: '-' reads standard input (plain or gzip)."""
	code = "import mappy\nfor r in mappy.fastx_read('-', read_comment=True): print('\\t'.join(map(str, r)))\n"
	text = H.fastx_text([("a", "ACGT", "IIII", "c1"), ("b", "TTGCA", "!!!!!", None)], "fastq")
	for data in (text.encode(), gzip.compress(text.encode(), mtime=0)):
		r = H.run_process([sys.executable, "-c", code], stdin=data, env=E.mappy_env(), check=True)
		assert r.stdout.splitlines() == ["a\tACGT\tIIII\tc1", "b\tTTGCA\t!!!!!\tNone"]


# ---------------------------------------------------------------------------
# Aligner construction, index properties and seq()
# ---------------------------------------------------------------------------


def test_aligner_from_fasta(mappy, synth, default_aligner):
	"""Aligner(fasta) builds an index with the default -k15 -w10 (minimap2.1); n_seq and seq_names follow the
	FASTA; seq() returns the stored bases, clips `end` to the contig and returns None for an unknown name or an
	empty range (python/README.rst)."""
	a, g = default_aligner, synth.genome
	assert bool(a) is True
	assert (a.k, a.w, a.n_seq, a.seq_names) == (15, 10, len(g.names), g.names)
	for n in g.names:
		assert a.seq(n) == g[n].upper()
	c = g["chr1"]
	assert a.seq("chr1", 1234, 1300) == c[1234:1300]
	assert a.seq("chr2", 0, 1) == g["chr2"][0]
	assert a.seq("chr1", len(c) - 7, len(c) + 100) == c[-7:]
	assert a.seq("chr1", 500) == c[500:]
	assert a.seq("chr1", len(c)) is None
	assert a.seq("chr1", 40, 40) is None and a.seq("chr1", 41, 40) is None
	assert a.seq("no_such_contig", 0, 10) is None
	n_run = [f for f in g.features if f.kind == "N"][0]
	assert a.seq(n_run.contig, n_run.start, n_run.end) == "N" * (n_run.end - n_run.start)


def test_missing_index_is_false(mappy, wd):
	"""python/README.rst: `if not a: raise ...` -- an Aligner whose index failed to load is false; map() then
	yields nothing and seq() and seq_names give None."""
	a = mappy.Aligner(str(wd / "does_not_exist.fa"))
	assert not a and bool(a) is False
	assert list(a.map("ACGT" * 100)) == []
	assert a.seq("chr1", 0, 10) is None
	assert a.seq_names is None


# (k, w, homopolymer-compressed) of each preset, from minimap2.1's preset list
PRESET_INDEX = {
	None: (15, 10, 0), "map-ont": (15, 10, 0), "lr:hq": (19, 19, 0), "map-hifi": (19, 19, 0), "map-pb": (19, 10, 1),
	"asm5": (19, 19, 0), "asm10": (19, 19, 0), "asm20": (19, 10, 0), "splice": (15, 5, 0), "splice:hq": (15, 5, 0),
	"sr": (21, 11, 0), "ava-pb": (19, 5, 1), "ava-ont": (15, 5, 0),
}


@pytest.mark.parametrize("preset", sorted(PRESET_INDEX, key=str))
def test_preset_index_parameters(mappy, synth, wd, preset):
	"""Aligner(preset=) applies the preset's -k, -w and -H (minimap2.1): Aligner.k/.w and the header of the
	index written with fn_idx_out (MM_IDX_MAGIC, w, k, bucket bits, n_seq, MM_I_HPC flag)."""
	k, w, hpc = PRESET_INDEX[preset]
	out = str(wd / ("preset-%s.mmi" % str(preset).replace(":", "_")))
	a = mappy.Aligner(synth.fasta, preset=preset, fn_idx_out=out)
	assert a and (a.k, a.w) == (k, w)
	magic, w2, k2, b, n_seq, flag = E.mmi_header(out)
	assert (magic, w2, k2, b, n_seq, flag & E.mm_flags()["MM_I_HPC"]) == (b"MMI\x02", w, k, 14, 4, hpc)


def test_index_file_matches_the_command_line(mappy, synth, wd):
	"""fn_idx_out writes the same .mmi bytes as `minimap2 [-x preset] -d` (the format is fixed), whatever
	n_threads; Aligner(fn_idx_in=.mmi) loads it with its k, w and names and maps exactly like the Aligner built
	from the FASTA."""
	reads = H.simulate_long_reads(synth.genome, 6, profile="hifi", seed=503)
	for preset in (None, "map-pb"):
		tag = str(preset)
		m1 = str(wd / ("mp1-%s.mmi" % tag))
		m4 = str(wd / ("mp4-%s.mmi" % tag))
		cli = str(wd / ("cli-%s.mmi" % tag))
		a1 = mappy.Aligner(synth.fasta, preset=preset, fn_idx_out=m1, n_threads=1)
		a4 = mappy.Aligner(synth.fasta, preset=preset, fn_idx_out=m4, n_threads=4)
		H.run_mm2((["-x", preset] if preset else []) + ["-t", "2", "-d", cli, synth.fasta], check=True)
		with open(m1, "rb") as f1, open(m4, "rb") as f4, open(cli, "rb") as fc:
			b1, b4, bc = f1.read(), f4.read(), fc.read()
		assert b1 == bc and b4 == bc
		loaded = mappy.Aligner(m1, preset=preset)
		assert (loaded.k, loaded.w, loaded.n_seq, loaded.seq_names) == (a1.k, a1.w, a1.n_seq, a1.seq_names)
		assert loaded.seq("chr2", 100, 200) == synth.genome["chr2"][100:200]
		want = E.mappy_hits(a1, reads)
		assert want and E.mappy_hits(loaded, reads) == want == E.mappy_hits(a4, reads)


def test_k_and_w_arguments(mappy, synth, wd, hifi):
	"""Aligner(k=, w=) override the preset's -k/-w: Aligner.k/.w and the written index say so, and the
	alignments equal `minimap2 -x map-hifi -k17 -w7`."""
	reads, fq = hifi
	out = str(wd / "k17w7.mmi")
	a = mappy.Aligner(synth.fasta, preset="map-hifi", k=17, w=7, fn_idx_out=out)
	assert (a.k, a.w) == (17, 7) and E.mmi_header(out)[1:3] == (7, 17)
	want = E.cli_hits(synth.fasta, fq, ["--no-hash-name", "-k17", "-w7"], preset="map-hifi")
	assert E.hit_problems(E.mappy_hits(a, reads), want) == []


def test_seq_argument_builds_a_one_sequence_index(mappy, synth, wd):
	"""python/README.rst: with seq= the index holds that one sequence, named N/A; fn_idx_in and fn_idx_out
	have no effect. Reads cut from it map exactly onto it on either strand."""
	s = synth.genome["chr2"][5000:25000]
	out = str(wd / "never_written.mmi")
	a = mappy.Aligner(str(wd / "no_such_file.fa"), seq=s, fn_idx_out=out)
	assert a and a.n_seq == 1 and a.seq_names == ["N/A"] and (a.k, a.w) == (15, 10)
	assert not os.path.exists(out)
	for st, en, strand in ((1000, 4000, 1), (12000, 15500, -1)):
		q = s[st:en] if strand > 0 else H.revcomp(s[st:en])
		hits = list(a.map(q))
		assert len(hits) == 1
		h = hits[0]
		assert (h.ctg, h.ctg_len, h.r_st, h.r_en, h.q_st, h.q_en, h.strand) == ("N/A", len(s), st, en, 0, en - st, strand)
		assert (h.cigar_str, h.NM, h.mlen, h.blen, h.is_primary) == ("%dM" % (en - st), 0, en - st, en - st, True)



def test_seq_argument_keeps_the_preset_hpc(mappy, tmp_path):
	"""Aligner(seq=, preset=) indexes the sequence with the preset's -k, -w and -H (minimap2.1: map-pb
	uses homopolymer-compressed minimizers, -Hk19). A 6 kbp sequence of short homopolymer runs and a
	read of 300 runs from it with every run lengthened by 1-3 bases share no exact 19-mer, but their HPC
	sequences are identical: with -H the read maps at its origin over at least 95% of its length (the
	ends, where a lengthened run costs a gap, may be clipped), exactly as with an Aligner built from the
	same sequence as a FASTA file; without -H (map-ont) it does not map."""
	rng = H.Rng(77)
	runs = []
	while sum(n for _, n in runs) < 6000:
		c = rng.choice("ACGT")
		if runs and runs[-1][0] == c:
			continue
		runs.append((c, rng.randint(1, 3)))
	s = "".join(c * n for c, n in runs)
	lo = len("".join(c * n for c, n in runs[:1000]))
	read = "".join(c * (n + rng.randint(1, 3) if i % 2 == 0 else n) for i, (c, n) in enumerate(runs[1000:1300]))
	hi = lo + len("".join(c * n for c, n in runs[1000:1300]))
	fa = H.write_fasta(str(tmp_path / "hp.fa"), [("N/A", s)])
	for preset, maps in (("map-pb", True), ("map-ont", False)):
		a = mappy.Aligner(seq=s, preset=preset)
		hits = [(h.r_st, h.r_en, h.q_st, h.q_en, h.strand) for h in a.map(read)]
		if maps:
			r_st, r_en, q_st, q_en, strand = hits[0]
			assert strand == 1 and abs(r_st - lo) <= 12 and abs(r_en - hi) <= 12 and q_en - q_st >= 0.95 * len(read), hits
			assert hits == [(h.r_st, h.r_en, h.q_st, h.q_en, h.strand) for h in mappy.Aligner(fa, preset=preset).map(read)]
		else:
			assert hits == [], hits


# ---------------------------------------------------------------------------
# map() and Alignment fields
# ---------------------------------------------------------------------------


def test_alignments_against_truth_and_recomputed_fields(mappy, synth, default_aligner, hifi_unique):
	"""Each read from unique sequence has exactly one primary alignment, on its true contig and strand within
	10 bp; every alignment's cigar list and cigar_str agree and span [q_st,q_en) and [r_st,r_en); NM, mlen
	and blen equal their recomputation from the genome (README: blen and mlen exclude ambiguous bases, NM
	counts them); ctg_len is the contig length; mapq is in 0..60; read_num is 1; trans_strand is 0."""
	reads, _ = hifi_unique
	g = synth.genome
	assert sum(r.seq.count("N") for r in reads) > 0
	for r in reads:
		hits = list(default_aligner.map(r.seq))
		prim = [h for h in hits if h.is_primary]
		assert len(prim) == 1, (r, hits)
		d = _pos(prim[0], r)
		assert d is not None and d <= 10, (str(prim[0]), r)
		assert prim[0].mapq >= 30
		for h in hits:
			assert h.cigar_str == "".join("%d%s" % (n, "MIDNSHP=XB"[op]) for n, op in h.cigar)
			assert all(isinstance(x, list) and len(x) == 2 for x in h.cigar)
			assert H.cigar_ref_span(h.cigar_str) == h.r_en - h.r_st
			assert H.cigar_query_span(h.cigar_str) == h.q_en - h.q_st
			assert 0 <= h.q_st < h.q_en <= len(r) and 0 <= h.r_st < h.r_en <= h.ctg_len == len(g[h.ctg])
			q = E.ref_oriented(r.seq, h)
			mlen, blen, ambi = H.compute_mlen_blen(g[h.ctg], h.r_st, q, h.cigar_str)
			assert (h.mlen, h.blen) == (mlen, blen)
			assert h.NM == H.compute_nm(g[h.ctg], h.r_st, q, h.cigar_str) == blen - mlen + ambi
			assert 0 <= h.mapq <= 60 and h.read_num == 1 and h.trans_strand == 0
			assert h.strand in (1, -1) and isinstance(h.is_primary, bool)


def test_cs_md_and_ds(mappy, synth, default_aligner, hifi_unique, hifi):
	"""map(cs=True, MD=True, ds=True) returns all three, each computed on its own: cs equals the short-form cs
	recomputed from the genome, MD the recomputed MD (in minimap2's form, without a final 0),
	and ds is cs with its indel bases bracketed (removing the brackets gives cs). Without the arguments the
	three are empty strings."""
	g = synth.genome
	n_rev = n_ds_brackets = 0
	for r in hifi_unique[0] + hifi[0][:8]:
		for h in default_aligner.map(r.seq, cs=True, MD=True, ds=True):
			q = E.ref_oriented(r.seq, h)
			ref = g[h.ctg]
			assert h.cs == H.compute_cs(ref, h.r_st, q, h.cigar_str), (r.name, str(h))
			assert h.MD == H.md_as_minimap2(H.compute_md(ref, h.r_st, q, h.cigar_str)), (r.name, str(h))
			assert re.sub(r"[\[\]]", "", h.ds) == h.cs
			n_rev += h.strand < 0
			n_ds_brackets += "[" in h.ds
		for h in default_aligner.map(r.seq):
			assert (h.cs, h.ds, h.MD) == ("", "", "")
		h = next(default_aligner.map(r.seq, MD=True))
		assert h.cs == "" and h.ds == "" and h.MD != ""
	assert n_rev > 0 and n_ds_brackets > 0


def test_cs_md_ds_match_the_command_line(mappy, synth, default_aligner, hifi):
	"""cs, MD and ds equal the cs:Z, MD:Z and ds:Z of `minimap2 -c --no-hash-name` with --cs, --MD and --ds
	(three runs: the command line writes only one of them per run)."""
	reads, fq = hifi
	want = {}
	for opt, tag in (("--cs", "cs"), ("--MD", "MD"), ("--ds", "ds")):
		for p in H.map_reads(synth.fasta, fq, ["-c", "--no-hash-name", opt]).paf():
			want.setdefault((p.qname, tag), []).append(p.tags.get(tag))
	for r in reads:
		hits = list(default_aligner.map(r.seq, cs=True, MD=True, ds=True))
		assert [h.cs for h in hits] == want.get((r.name, "cs"), [])
		assert [h.MD for h in hits] == want.get((r.name, "MD"), [])
		assert [h.ds for h in hits] == want.get((r.name, "ds"), [])


def test_str_is_a_paf_line_from_column_3(mappy, synth, default_aligner, hifi):
	"""python/README.rst: str(Alignment) is PAF without the query name and length: q_st q_en strand ctg ctg_len
	r_st r_en mlen blen mapq, then tp:A (P or S), ts:A, cg:Z and the cs, ds and MD tags that were asked for.
	Columns 3-12 and cg equal `minimap2 -c --no-hash-name --secondary=yes` PAF."""
	reads, fq = hifi
	paf = {}
	for p in H.map_reads(synth.fasta, fq, ["-c", "--no-hash-name", "--secondary=yes"]).paf():
		paf.setdefault(p.qname, []).append(p)
	n = 0
	for r in reads:
		hits = list(default_aligner.map(r.seq))
		assert len(hits) == len(paf.get(r.name, []))
		for h, p in zip(hits, paf.get(r.name, [])):
			f = str(h).split("\t")
			assert f[:10] == p.line.split("\t")[2:12]
			assert f[10:] == ["tp:A:" + ("P" if h.is_primary else "S"), "ts:A:.", "cg:Z:" + p.tags["cg"]]
			n += 1
		if hits:
			h = next(default_aligner.map(r.seq, cs=True, MD=True, ds=True))
			f = str(h).split("\t")
			assert f[13:] == ["cs:Z:" + h.cs, "ds:Z:" + h.ds, "MD:Z:" + h.MD]
			assert str(next(default_aligner.map(r.seq, MD=True))).split("\t")[13:] == ["MD:Z:" + h.MD]
	assert n >= len(reads)


@pytest.mark.parametrize("preset", [None, "map-hifi", "map-pb", "lr:hq", "asm20", "sr", "splice", "splice:hq"])
def test_preset_alignments_match_the_command_line(mappy, synth, synth_tx, wd, preset, hifi):
	"""Aligner(preset=) maps like `minimap2 -x preset -c --secondary=yes --no-hash-name` (mappy returns
	secondaries even where the preset hides them): every field the PAF line carries, for every alignment,
	in order. Genomic presets use HiFi-like reads (sr: 150 bp reads), splice presets cDNA reads."""
	if preset and preset.startswith("splice"):
		fa = synth_tx.fasta
		reads = H.simulate_spliced_reads(synth_tx.genome, 16, seed=504, profile="hifi")
	elif preset == "sr":
		fa = synth.fasta
		reads = H.simulate_reads(synth.genome, 60, length=150, seed=505, sub_rate=0.01, indel_rate=0.002)
	else:
		fa = synth.fasta
		reads = hifi[0]
	fq = H.write_fastq(str(wd / ("preset-%s.fq" % str(preset).replace(":", "_"))), reads)
	a = mappy.Aligner(fa, preset=preset)
	got = E.mappy_hits(a, reads)
	want = E.cli_hits(fa, fq, ["--no-hash-name"], preset=preset)
	assert len(got) >= 0.9 * len(reads)
	assert E.hit_problems(got, want) == []


def test_name_seeds_the_tie_breaking_hash(mappy, synth, wd):
	"""map(name=) passes the read name to minimap2, whose hash of it breaks ties between equal hits: with names
	mappy equals the command line's default output, without them `--no-hash-name`; on reads cut from exact
	repeat copies the two runs pick different copies for at least one read."""
	g = synth.genome
	reads = []
	for f in g.features:
		if f.kind == "repeat":
			for i in range(4):
				reads.append(("rep_%s_%d_%d" % (f.contig, f.start, i), g.fetch(f.contig, f.start + 10 * i, f.end - 10 * i, "+")))
	fq = H.write_fasta(str(wd / "rep.fa"), reads)
	a = mappy.Aligner(synth.fasta)
	with_name = E.cli_hits(synth.fasta, fq)
	no_name = E.cli_hits(synth.fasta, fq, ["--no-hash-name"])
	assert E.hit_problems(E.mappy_hits(a, reads, use_name=True), with_name) == []
	assert E.hit_problems(E.mappy_hits(a, reads), no_name) == []
	assert sum(1 for n, _ in reads if with_name.get(n) != no_name.get(n)) > 0


def test_trans_strand_of_spliced_reads(mappy, synth_tx, wd):
	"""python/README.rst: trans_strand is +1 or -1 for the transcript strand (relative to the read, as ts:A),
	0 when unknown. splice:hq on cDNA reads in both orientations: spliced primaries carry the simulated
	orientation and the true introns (N operations from r_st) for at least 90% of reads; the value equals the
	command line's ts:A; str() writes it as ts:A."""
	g = synth_tx.genome
	reads = H.simulate_spliced_reads(g, 30, seed=506, profile="hifi")
	fq = H.write_fastq(str(wd / "tx.fq"), reads)
	a = mappy.Aligner(synth_tx.fasta, preset="splice:hq")
	ts = {}
	for p in H.map_reads(synth_tx.fasta, fq, ["-c", "--no-hash-name"], preset="splice:hq").paf():
		ts.setdefault(p.qname, []).append({"+": 1, "-": -1}.get(p.tags.get("ts"), 0))
	n_ok = n_spliced = 0
	for r in reads:
		hits = list(a.map(r.seq))
		assert [h.trans_strand for h in hits] == ts.get(r.name, [])
		for h in hits:
			assert str(h).split("\t")[11] == "ts:A:" + {1: "+", -1: "-", 0: "."}[h.trans_strand]
		if not r.introns:
			continue
		n_spliced += 1
		h = [x for x in hits if x.is_primary][0]
		pos, introns = h.r_st, []
		for n, op in h.cigar:
			if op == 3:
				introns.append((pos, pos + n))
			if op in (0, 2, 3, 7, 8):
				pos += n
		if h.trans_strand == (1 if r.ts_strand == "+" else -1) and introns == r.introns:
			n_ok += 1
	assert n_spliced >= 20 and n_ok >= 0.9 * n_spliced


# ---------------------------------------------------------------------------
# Aligner arguments: filters and band widths
# ---------------------------------------------------------------------------


def test_min_cnt_min_chain_score_min_dp_score(mappy, synth, wd):
	"""A chain has at most len-k+1 minimizers and a chain score below the read length times k, and an
	alignment scores at most A=2 per base: thresholds above those bounds leave nothing. Thresholds that
	matter equal the command line's -n, -m and -s."""
	g = synth.genome
	reads = H.simulate_long_reads(g, 30, profile="ont", seed=507, unique=True, length=(300, 1500))
	fq = H.write_fastq(str(wd / "short_ont.fq"), reads)
	longest = max(len(r) for r in reads)
	for kw in (dict(min_cnt=longest), dict(min_chain_score=15 * longest), dict(min_dp_score=2 * longest + 1)):
		a = mappy.Aligner(synth.fasta, **kw)
		assert a and E.mappy_hits(a, reads) == {}, kw
	base = E.mappy_hits(mappy.Aligner(synth.fasta), reads)
	assert len(base) >= 0.9 * len(reads)
	for kw, opt in ((dict(min_cnt=40), ["-n", "40"]), (dict(min_chain_score=300), ["-m", "300"]),
			(dict(min_dp_score=900), ["-s", "900"])):
		got = E.mappy_hits(mappy.Aligner(synth.fasta, **kw), reads)
		assert got != base, kw
		assert E.hit_problems(got, E.cli_hits(synth.fasta, fq, ["--no-hash-name"] + opt)) == [], kw


def test_bw_short_reads_with_a_deletion(mappy, synth, wd):
	"""bw (-r) bounds the diagonal difference within a chain. 150 bp reads with a 40 bp deletion in the middle,
	-x sr (no long-join rechaining): bw=50 joins the halves into one alignment with a 40D spanning 190 bp of
	reference; bw=10 cannot, so each read gives two alignments without a deletion. Both equal `-x sr -r`."""
	c = synth.genome["chr2"]
	reads = []
	for i in range(6):
		st = 2000 + 3000 * i
		reads.append(("del%d" % i, c[st:st + 75] + c[st + 115:st + 190], st))
	fa = H.write_fasta(str(wd / "del40.fa"), [x[:2] for x in reads])
	for bw in (50, 10):
		a = mappy.Aligner(synth.fasta, preset="sr", bw=bw)
		for name, seq, st in reads:
			hits = list(a.map(seq))
			if bw == 50:
				assert len(hits) == 1, (name, [str(h) for h in hits])
				h = hits[0]
				assert (h.ctg, h.r_st, h.r_en, h.q_st, h.q_en) == ("chr2", st, st + 190, 0, 150)
				assert [x for x in h.cigar if x[1] == 2] == [[40, 2]]
			else:
				assert len(hits) == 2 and all("D" not in h.cigar_str for h in hits), (name, [str(h) for h in hits])
		want = E.cli_hits(synth.fasta, fa, ["--no-hash-name", "-r", str(bw)], preset="sr")
		assert E.hit_problems(E.mappy_hits(a, [x[:2] for x in reads]), want) == []


def test_bw_long_and_long_join(mappy, synth, wd):
	"""bw_long (second value of -r) is the band of rechaining: a 4 kbp read with a 1500 bp deletion aligns in
	one piece with a 1500D by default (bw_long 20000) but in two pieces with bw_long=1000 (-r500,1000). With
	long-join off (extra_flags MM_F_NO_LJOIN, --no-long-join) bw alone decides: bw=100 splits a 300 bp
	deletion, bw=400 keeps it. All equal the command line."""
	c = synth.genome["chr2"]
	r1500 = ("d1500", c[1000:3000] + c[4500:6500])
	r300 = ("d300", c[20000:22000] + c[22300:24300])
	fa = H.write_fasta(str(wd / "longdel.fa"), [r1500, r300])
	cases = (
		(dict(), [], r1500, [("chr2", 1000, 6500, "2000M1500D2000M")]),
		(dict(bw_long=1000), ["-r", "500,1000"], r1500, [("chr2", 1000, 3000), ("chr2", 4500, 6500)]),
		(dict(extra_flags=E.mm_flags()["MM_F_NO_LJOIN"], bw=100), ["--no-long-join", "-r", "100"], r300,
			[("chr2", 20000, 22000), ("chr2", 22300, 24300)]),
		(dict(extra_flags=E.mm_flags()["MM_F_NO_LJOIN"], bw=400), ["--no-long-join", "-r", "400"], r300,
			[("chr2", 20000, 24300, "2000M300D2000M")]),
	)
	for kw, opt, (name, seq), want in cases:
		a = mappy.Aligner(synth.fasta, **kw)
		hits = sorted(a.map(seq), key=lambda h: h.r_st)
		assert len(hits) == len(want), (kw, [str(h) for h in hits])
		for h, w in zip(hits, want):
			assert h.ctg == w[0] and h.is_primary
			if len(w) == 4:  # one alignment across the deletion; the gap may sit a few bases either way
				assert (h.r_st, h.r_en, h.q_st, h.q_en) == (w[1], w[2], 0, len(seq))
				assert [x[0] for x in h.cigar if x[1] == 2] == [int(re.search(r"(\d+)D", w[3]).group(1))]
			else:  # one alignment per half, each within a few bases of its half
				assert abs(h.r_st - w[1]) <= 5 and abs(h.r_en - w[2]) <= 5 and "D" not in h.cigar_str.replace("1D", "")
		got = E.mappy_hits(a, [(name, seq)])
		assert E.hit_problems(got, E.cli_hits(synth.fasta, fa, ["--no-hash-name"] + opt), names=[name]) == [], kw


def test_best_n_limits_secondary_alignments(mappy, synth, wd):
	"""best_n (-N) is the number of secondary alignments kept: a read cut from inside one copy of a
	three-copy exact repeat has its primary plus two secondaries on the other copies by default, and
	min(best_n, 2) secondaries with best_n=1 or 2. Equal to `minimap2 -N`."""
	g = synth.genome
	fam = {}
	for f in g.features:
		if f.kind == "repeat":
			fam.setdefault(f.info[0], []).append(f)
	reads, copies = [], {}
	for fid, fs in sorted(fam.items()):
		assert len(fs) == 3
		f = fs[0]
		name = "fam%d" % fid
		reads.append((name, g.fetch(f.contig, f.start + 20, f.end - 20, "+")))
		copies[name] = set((x.contig, x.start + 20) for x in fs)  # exact copies: trimming 20 bp is strand-symmetric
	fa = H.write_fasta(str(wd / "fam.fa"), reads)
	for best_n, opt in ((None, []), (1, ["-N", "1"]), (2, ["-N", "2"])):
		a = mappy.Aligner(synth.fasta, best_n=best_n) if best_n is not None else mappy.Aligner(synth.fasta)
		for name, seq in reads:
			hits = list(a.map(seq))
			n_sec = sum(1 for h in hits if not h.is_primary)
			assert sum(1 for h in hits if h.is_primary) == 1
			assert n_sec == (2 if best_n is None else min(best_n, 2)), (best_n, [str(h) for h in hits])
			assert all(h.cigar_str == "%dM" % len(seq) and h.NM == 0 for h in hits)
			assert set((h.ctg, h.r_st) for h in hits) <= copies[name]
		want = E.cli_hits(synth.fasta, fa, ["--no-hash-name"] + opt)
		assert E.hit_problems(E.mappy_hits(a, reads), want) == [], best_n


def test_max_chain_skip(mappy, wd):
	"""max_chain_skip (--max-chain-skip) bounds how many already-chained anchors the chaining DP skips before
	giving up on a predecessor search. On reads across a 9-copy tandem repeat it changes the alignments of
	some reads, and mappy equals the command line with the same value."""
	rng = H.Rng(109)
	unit = rng.dna(100)
	ref = rng.dna(5000) + unit * 9 + rng.dna(5000)
	fa = H.write_fasta(str(wd / "tandem.fa"), [("tandem", ref)])
	reads = []
	for i in range(40):
		st, en = rng.randint(3000, 4900), rng.randint(5900, 7900)
		reads.append(("t%d" % i, H.mutate(rng, ref[st:en], 0.05)))
	fq = H.write_fasta(str(wd / "tandem_reads.fa"), reads)
	base = E.mappy_hits(mappy.Aligner(fa), reads)
	for skip in (0, 2):
		got = E.mappy_hits(mappy.Aligner(fa, max_chain_skip=skip), reads)
		assert got != base, skip
		assert E.hit_problems(got, E.cli_hits(fa, fq, ["--no-hash-name", "--max-chain-skip", str(skip)])) == [], skip


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------


def _other(c):
	return {"A": "C", "C": "G", "G": "T", "T": "A"}[c]


@pytest.fixture(scope="module")
def tails(synth):
	"""2 kbp of unique chr2 followed by a short tail that an extension either takes or clips:
	mm: 2 mismatches then 5 matches; n: 6 N then 4 matches; dLmM: an L bp deletion then M matches."""
	c = synth.genome["chr2"]
	st = 10000
	core = c[st:st + 2000]
	t = {
		"mm": (core + _other(c[st + 2000]) + _other(c[st + 2001]) + c[st + 2002:st + 2007], ("X", 2, 5)),
		"n": (core + "N" * 6 + c[st + 2006:st + 2010], ("N", 6, 4)),
	}
	for L, m in ((2, 5), (15, 12), (15, 21)):
		t["d%dm%d" % (L, m)] = (core + c[st + 2000 + L:st + 2000 + L + m], ("D", L, m))
	return t


def _extends(kind, k, m, a, b, q, e, q2, e2, sc_ambi):
	"""Whether extension takes the tail: its m matches (m*a) must outscore what precedes them."""
	if kind == "X":
		cost = k * b
	elif kind == "N":
		cost = k * sc_ambi
	else:
		cost = min(q + k * e, q2 + k * e2)
	return m * a > cost


# (Aligner kwargs, the scoring they amount to (a, b, q, e, q2, e2, sc_ambi), command-line options)
SCORING_CASES = [
	(dict(), (2, 4, 4, 2, 24, 1, 1), []),
	(dict(scoring=(2, 6, 4, 2)), (2, 6, 4, 2, 4, 2, 1), ["-A2", "-B6", "-O4", "-E2"]),
	(dict(scoring=(1, 4, 4, 2)), (1, 4, 4, 2, 4, 2, 1), ["-A1", "-B4", "-O4", "-E2"]),
	(dict(scoring=(2, 4, 8, 2)), (2, 4, 8, 2, 8, 2, 1), ["-A2", "-B4", "-O8", "-E2"]),
	(dict(scoring=(2, 4, 4, 4)), (2, 4, 4, 4, 4, 4, 1), ["-A2", "-B4", "-O4", "-E4"]),
	(dict(scoring=(2, 4, 4, 2, 8, 1)), (2, 4, 4, 2, 8, 1, 1), ["-A2", "-B4", "-O4,8", "-E2,1"]),
	(dict(scoring=(2, 4, 4, 2, 24, 1, 2)), (2, 4, 4, 2, 24, 1, 2), ["-A2", "-B4", "-O4,24", "-E2,1", "--score-N", "2"]),
	(dict(sc_ambi=2), (2, 4, 4, 2, 24, 1, 2), ["--score-N", "2"]),
	(dict(scoring=(2, 4, 4, 2, 24, 1, 2), sc_ambi=1), (2, 4, 4, 2, 24, 1, 1), ["-A2", "-B4", "-O4,24", "-E2,1", "--score-N", "1"]),
]


@pytest.mark.parametrize("kw,sc,opt", SCORING_CASES, ids=[str(c[0]) or "default" for c in SCORING_CASES])
def test_scoring_and_sc_ambi_decide_end_extension(mappy, synth, wd, tails, kw, sc, opt):
	"""python/README.rst: scoring is (match, mismatch, gap open, gap extension[, long-gap open, long-gap
	extension[, ambiguous-base penalty]]); with four values both gap pieces use the same costs; sc_ambi also
	sets the ambiguous-base penalty and wins over a seventh scoring value. Each tail is kept exactly when its
	matches outscore the mismatches, N columns or two-piece gap cost min(q+k*e, q2+k*e2) before them, so the
	read is aligned to its end or clipped right after the 2 kbp core. The alignments equal the command line
	with -A -B -O -E --score-N."""
	a = mappy.Aligner(synth.fasta, **kw)
	recs = []
	for name, (seq, (kind, k, m)) in sorted(tails.items()):
		hits = list(a.map(seq))
		assert len(hits) == 1, (name, [str(h) for h in hits])
		h = hits[0]
		ext = _extends(kind, k, m, *sc)
		assert (h.ctg, h.r_st, h.q_st) == ("chr2", 10000, 0)
		assert h.q_en == (len(seq) if ext else 2000), (name, sc, str(h))
		if ext and kind == "D":
			assert h.cigar_str == "2000M%dD%dM" % (k, m)
		recs.append((name, seq))
	fa = H.write_fasta(str(wd / "tails.fa"), recs)
	want = E.cli_hits(synth.fasta, fa, ["--no-hash-name"] + opt)
	assert E.hit_problems(E.mappy_hits(a, recs), want) == []


# ---------------------------------------------------------------------------
# extra_flags
# ---------------------------------------------------------------------------


def test_extra_flags_strand_and_eqx(mappy, synth, hifi_unique):
	"""extra_flags ORs minimap.h flags into the mapping options, on the Aligner (every call) or on one map()
	call: MM_F_FOR_ONLY keeps forward-strand alignments only and MM_F_REV_ONLY reverse-strand ones (reads
	still land at their truth on the allowed strand); MM_F_EQX writes = and X operations that agree with the
	genome base by base. A flag given to map() does not stick to the Aligner."""
	F = E.mm_flags()
	reads, _ = hifi_unique
	g = synth.genome
	assert set(r.strand for r in reads) == {"+", "-"}
	fwd = mappy.Aligner(synth.fasta, extra_flags=F["MM_F_FOR_ONLY"])
	both = mappy.Aligner(synth.fasta)
	for r in reads:
		for aln, kw, allowed in ((fwd, {}, 1), (both, dict(extra_flags=F["MM_F_REV_ONLY"]), -1)):
			hits = list(aln.map(r.seq, **kw))
			assert all(h.strand == allowed for h in hits)
			if (r.strand == "+") == (allowed > 0):
				prim = [h for h in hits if h.is_primary]
				assert len(prim) == 1 and _pos(prim[0], r) <= 10
		strands = set(h.strand for h in both.map(r.seq))
		assert (1 if r.strand == "+" else -1) in strands
	eqx = mappy.Aligner(synth.fasta, extra_flags=F["MM_F_EQX"])
	n_x = 0
	for r in reads:
		plain = list(both.map(r.seq))
		hits = list(eqx.map(r.seq))
		assert [(h.r_st, h.r_en, h.q_st, h.q_en) for h in hits] == [(h.r_st, h.r_en, h.q_st, h.q_en) for h in plain]
		per_call = list(both.map(r.seq, extra_flags=F["MM_F_EQX"]))
		assert [h.cigar_str for h in per_call] == [h.cigar_str for h in hits]
		for h in hits:
			assert set(op for _, op in h.cigar) <= {1, 2, 7, 8}
			q, ref = E.ref_oriented(r.seq, h), g[h.ctg]
			x, y = h.r_st, 0
			for n, op in h.cigar:
				if op in (7, 8):
					same = [H.nt4(ref[x + i]) == H.nt4(q[y + i]) for i in range(n)]
					assert all(same) if op == 7 else not any(same)
					n_x += op == 8
				x += n if op in (2, 7, 8) else 0
				y += n if op in (1, 7, 8) else 0
	assert n_x > 0


# ---------------------------------------------------------------------------
# Pairs, buffers
# ---------------------------------------------------------------------------


def test_pairs_seq2_and_read_num(mappy, synth, wd):
	"""python/README.rst: map(seq, seq2) aligns an FR pair; read_num tells the mates apart. With -x sr each
	mate of pairs from unique sequence lands at its truth (within 3 bp) on opposite strands, and every
	alignment's contig, reference interval, strand, MAPQ, primary flag, CIGAR, NM, mlen and blen equal the
	command line's /1 and /2 PAF lines (query interval and cs for mate 2 only when the alignment covers the
	whole mate, since mappy reports mate 2's query interval on its reverse complement). Mates are 140-160 bp, so mappy's reverse
	complement of mate 2 meets odd and even lengths."""
	pairs = H.simulate_pairs(synth.genome, 30, seed=508, read_length=(140, 160), insert_mean=300, insert_sd=30, unique=True)
	assert len(set(len(p.r2.seq) % 2 for p in pairs)) == 2
	f1, f2 = H.write_pairs(str(wd / "p1.fq"), str(wd / "p2.fq"), pairs)
	paf = {}
	for p in H.map_reads(synth.fasta, [f1, f2], ["-c", "--secondary=yes", "--no-hash-name", "--cs"], preset="sr").paf():
		paf.setdefault(p.qname, []).append(p)
	a = mappy.Aligner(synth.fasta, preset="sr")
	for p in pairs:
		hits = list(a.map(p.r1.seq, p.r2.seq, cs=True))
		assert set(h.read_num for h in hits) == {1, 2}
		for rn, read in ((1, p.r1), (2, p.r2)):
			mine = [h for h in hits if h.read_num == rn]
			cli = paf["%s/%d" % (p.name, rn)]
			assert len(mine) == len(cli)
			prim = [h for h in mine if h.is_primary]
			assert len(prim) == 1 and _pos(prim[0], read) <= 3, (p, [str(h) for h in mine])
			for h, c in zip(mine, cli):
				k = E.hit_key(h)
				w = E.paf_key(c)
				assert k[:4] + k[6:] == w[:4] + w[6:]
				full = h.q_en - h.q_st == len(read)
				if rn == 1 or full:
					assert (h.q_st, h.q_en, h.cs) == (c.qstart, c.qend, c.tags["cs"])
		s = [h.strand for h in hits if h.is_primary]
		assert sorted(s) == [-1, 1]


def test_max_frag_len_pairs_a_repetitive_mate(mappy, synth, wd):
	"""max_frag_len (-F) is the longest fragment that pairs. Mate 1 from unique sequence 300 bp from an exact
	repeat copy, mate 2 inside the copy: within the sr preset's 800 bp the pair is mapped as a fragment and
	mate 2 lands on its own copy; with max_frag_len=100, given to the Aligner or to map(), it is placed alone
	with MAPQ 0. Equal to `-x sr -F`."""
	g = synth.genome
	pairs = []
	for f in g.features:
		if f.kind != "repeat" or not g.is_unique(f.contig, f.start - 450, f.start - 20):
			continue
		r1 = g.fetch(f.contig, f.start - 400, f.start - 250, "+")
		r2 = g.fetch(f.contig, f.start + 50, f.start + 200, "-")
		pairs.append((f, r1, r2))
	assert len(pairs) >= 2
	a = mappy.Aligner(synth.fasta, preset="sr")
	a100 = mappy.Aligner(synth.fasta, preset="sr", max_frag_len=100)
	for f, r1, r2 in pairs:
		for aln, kw, paired in ((a, {}, True), (a100, {}, False), (a, dict(max_frag_len=100), False),
				(a100, dict(max_frag_len=800), True)):
			m2 = [h for h in aln.map(r1, r2, **kw) if h.read_num == 2 and h.is_primary]
			assert len(m2) == 1
			if paired:
				assert (m2[0].ctg, m2[0].r_st, m2[0].mapq > 0) == (f.contig, f.start + 50, True), (f, str(m2[0]))
			else:
				assert m2[0].mapq == 0, (f, kw, str(m2[0]))
	f1 = H.write_fasta(str(wd / "fl1.fa"), [("fl%d" % i, x[1]) for i, x in enumerate(pairs)])
	f2 = H.write_fasta(str(wd / "fl2.fa"), [("fl%d" % i, x[2]) for i, x in enumerate(pairs)])
	for aln, opt in ((a, []), (a100, ["-F", "100"])):
		cli = {}
		for p in H.map_reads(synth.fasta, [f1, f2], ["-c", "--secondary=yes", "--no-hash-name"] + opt, preset="sr").paf():
			cli.setdefault(p.qname, []).append(p)
		for i, (f, r1, r2) in enumerate(pairs):
			hits = list(aln.map(r1, r2))
			for rn in (1, 2):
				mine = [E.hit_key(h) for h in hits if h.read_num == rn]
				want = [E.paf_key(c) for c in cli["fl%d/%d" % (i, rn)]]
				assert [k[:4] + k[6:] for k in mine] == [k[:4] + k[6:] for k in want], (opt, rn)


def test_thread_buffer_reuse(mappy, synth, default_aligner, hifi):
	"""ThreadBuffer() is a reusable per-thread buffer: passing one buffer to every map() call, single-end and
	paired, gives the same alignments as a fresh buffer per call."""
	reads, _ = hifi
	buf = mappy.ThreadBuffer()
	sr = mappy.Aligner(synth.fasta, preset="sr")
	for r in reads:
		assert [E.hit_key(h) for h in default_aligner.map(r.seq, buf=buf, cs=True)] == \
			[E.hit_key(h) for h in default_aligner.map(r.seq, cs=True)]
		a, b = r.seq[:150], H.revcomp(r.seq[300:450])
		assert [E.hit_key(h) for h in sr.map(a, b, buf=buf)] == [E.hit_key(h) for h in sr.map(a, b)]


def test_threads_with_a_buffer_each(mappy, synth, default_aligner, hifi, hifi_unique):
	"""One Aligner shared by four Python threads, each with its own ThreadBuffer (a per-thread buffer;
	cmappy.h releases the GIL while it maps), gives every read the alignments, cs and MD that one
	thread gives it, single-end and for -x sr pairs. Each thread maps every read twice, starting at
	a different one, so the threads map different reads at the same time."""
	import threading
	reads = [r.seq for r in hifi[0] + hifi_unique[0]]
	sr = mappy.Aligner(synth.fasta, preset="sr")
	pairs = [(s[:150], H.revcomp(s[200:350])) for s in reads]

	def one(buf, i):
		se = [E.hit_key(h) + (h.cs, h.MD) for h in default_aligner.map(reads[i], buf=buf, cs=True, MD=True)]
		pe = [E.hit_key(h) + (h.read_num,) for h in sr.map(pairs[i][0], pairs[i][1], buf=buf)]
		return se, pe
	want = [one(None, i) for i in range(len(reads))]
	assert sum(1 for se, pe in want if se and pe) >= 0.9 * len(reads)
	n_threads, got, errors = 4, {}, []
	start = threading.Barrier(n_threads)

	def work(t):
		try:
			buf = mappy.ThreadBuffer()
			start.wait()
			for rnd in range(2):
				n = len(reads)
				for j in range(n):  # a different starting read per thread, odd threads backwards
					i = (j + 11 * t + 5 * rnd) % n if t % 2 == 0 else (n - 1 - j + 11 * t) % n
					got[(t, rnd, i)] = one(buf, i)
		except Exception as e:  # reported in the main thread
			errors.append(repr(e))
	ths = [threading.Thread(target=work, args=(t,)) for t in range(n_threads)]
	for th in ths:
		th.start()
	for th in ths:
		th.join()
	assert errors == []
	assert len(got) == n_threads * 2 * len(reads)
	for (t, rnd, i), v in got.items():
		assert v == want[i], (t, rnd, i)
