"""PAF output (format.c): the 12 columns and the tags with and without -c, agreement with the SAM
output of the same alignments, --paf-no-hit, -y, cs/MD/ds in PAF, segment-name suffixes in the
fragment mode, and --write-junc.

Expected values come from minimap2.1 (the PAF and tag tables, --write-junc), the README notes
in tests/README.md, and recomputation from the genome.
"""

import re

import pytest

import helpers_d as D
import minimap2_helpers as H


@pytest.fixture(scope="module")
def mixed(synth, tmp_path_factory):
	m = D.mixed_reads(synth.genome)
	d = tmp_path_factory.mktemp("d_paf")
	fq = D.write_mixed(d / "mixed.fq", m)
	return m, fq


@pytest.fixture(scope="module")
def paf_runs(synth, mixed):
	"""Default PAF, -c PAF and -a SAM of the mixed reads."""
	fq = mixed[1]
	return dict(plain=H.map_reads(synth.fasta, fq), c=H.map_reads(synth.fasta, fq, ["-c"]),
		sam=H.map_reads(synth.fasta, fq, ["-a"]))


_DV = re.compile(r"^(0|[01]\.[0-9]{4})$")


def test_paf_without_alignment(synth, mixed, paf_runs):
	"""Default PAF (minimap2.1, OUTPUT FORMAT): 12 columns, then tp, cm, s1, dv and rl, plus s2 on
	primary lines; no base-level tags. Query length and coordinates are in range, each mapped read
	that is not chimeric has one tp:A:P line, unmapped reads are absent, reads come in input order, and dv:f is 0 or a
	value with four decimals."""
	m = mixed[0]
	paf = paf_runs["plain"].paf()
	assert H.paf_invariant_errors(paf, synth.genome, m.seqs) == []
	for p in paf:
		assert p.qlen == len(m.seqs[p.qname])
		want = {"tp", "cm", "s1", "dv", "rl"} | ({"s2"} if p.tags["tp"] in "PI" else set())
		assert set(p.tags) == want, p.line[:200]
		dv = [f for f in p.line.split("\t") if f.startswith("dv:f:")][0][5:]
		assert _DV.match(dv) and 0.0 <= float(dv) <= 1.0
	names = [p.qname for p in paf]
	order = [n for n in m.names() if n in set(names)]
	assert [n for i, n in enumerate(names) if i == 0 or names[i - 1] != n] == order
	assert set(m.names()) - set(names) == set(m.names("random"))
	for n in set(names):
		k = sum(1 for p in paf if p.qname == n and p.is_primary)
		assert k >= 1 if m.kind[n] == "chimeric" else k == 1, n


def test_paf_c_agrees_with_sam(synth, mixed, paf_runs):
	"""With -c, each PAF line describes the same alignment as the SAM record in the same place of
	the -a output (both list a read's hits in one order): target, strand, start (POS - 1), end
	(POS + reference span), query interval (the clips, counted from the read's 5' end), MAPQ,
	cg:Z (the CIGAR without clips), NM, AS, ms, nn, tp, cm, s1, de and rl. -c adds NM, ms, AS, nn,
	de and cg to the tags (tests/README.md); columns 10/11, NM and de equal their recomputation."""
	m = mixed[0]
	paf = paf_runs["c"].paf()
	assert H.paf_invariant_errors(paf, synth.genome, m.seqs) == []
	recs = [r for r in paf_runs["sam"].sam().records if not r.is_unmapped]
	assert len(paf) == len(recs)
	for p, r in zip(paf, recs):
		assert {"NM", "ms", "AS", "nn", "de", "cg", "tp", "cm", "s1", "rl"} <= set(p.tags)
		assert "dv" not in p.tags
		qlen = len(m.seqs[r.qname])
		lc, rc = H.cigar_clips(r.cigar)
		qs, qe = (rc, qlen - lc) if r.is_reverse else (lc, qlen - rc)
		assert (p.qname, p.qlen, p.qstart, p.qend, p.strand) == (r.qname, qlen, qs, qe, r.strand)
		assert (p.tname, p.tstart, p.tend, p.mapq) == (r.rname, r.pos0, r.ref_end0, r.mapq)
		assert p.cigar == H.strip_clips(r.cigar)
		for t in ("NM", "AS", "ms", "nn", "tp", "cm", "s1", "de", "rl"):
			assert p.tags[t] == r.tags[t], (t, p.line[:100], r.line[:100])


def test_paf_no_hit(synth, mixed, paf_runs):
	"""--paf-no-hit (minimap2.1; tests/README.md): an unmapped query gets the line `name qlen 0 0
	* * 0 0 0 0 0 0` plus rl:i, in input order; every other line is unchanged."""
	m = mixed[0]
	r = H.map_reads(synth.fasta, mixed[1], ["--paf-no-hit"])
	got = D.lines(r.stdout)
	base = D.lines(paf_runs["plain"].stdout)
	assert [l for l in got if l.split("\t")[0] not in m.names("random")] == base
	for n in m.names("random"):
		ls = [l for l in got if l.split("\t")[0] == n]
		assert len(ls) == 1
		f = ls[0].split("\t")
		assert f[:12] == [n, str(len(m.seqs[n])), "0", "0", "*", "*", "0", "0", "0", "0", "0", "0"]
		assert len(f) == 13 and re.match(r"^rl:i:[0-9]+$", f[12])
	names = [l.split("\t")[0] for l in got]
	firsts = [n for i, n in enumerate(names) if i == 0 or names[i - 1] != n]
	assert firsts == m.names()
	assert H.parse_paf(r.stdout)[names.index(m.names("random")[0])].is_unmapped


def test_paf_copy_comment(synth, mixed, paf_runs):
	"""-y in PAF: the read's comment is appended to each of its lines, with and without -c;
	lines of reads without a comment are unchanged."""
	m = mixed[0]
	for key, args in (("plain", []), ("c", ["-c"])):
		got = D.lines(H.map_reads(synth.fasta, mixed[1], args + ["-y"]).stdout)
		base = D.lines(paf_runs[key].stdout)
		assert len(got) == len(base)
		for g, b in zip(got, base):
			c = m.comments[b.split("\t")[0]]
			assert g == (b if c is None else b + "\t" + c)


@pytest.mark.parametrize("args,tag", [(["--cs"], "cs"), (["--cs=long"], "cs"), (["-c", "--MD"], "MD")])
def test_paf_cs_md(synth, mixed, paf_runs, args, tag):
	"""--cs and --cs=long in PAF imply alignment (the coordinates are those of -c, and cg:Z is absent
	without -c); --MD with -c. The tag equals its recomputation from the genome. (--MD and --ds
	without -c write nothing in PAF: known-bugs/d-md-ds-need-c-in-paf.md.)"""
	m = mixed[0]
	r = H.map_reads(synth.fasta, mixed[1], args).paf()
	c = paf_runs["c"].paf()
	assert len(r) == len(c)
	for p, q in zip(r, c):
		assert tag in p.tags and ("cg" in p.tags) == ("-c" in args)
		assert p.line.split("\t")[:12] == q.line.split("\t")[:12]
		p.tags["cg"] = q.cigar
	assert H.paf_invariant_errors(r, synth.genome, m.seqs) == []


def test_paf_ds(synth, mixed):
	"""--ds in PAF: ds:Z without its brackets is the cs:Z of the --cs run, and equals its
	recomputation (helpers_d.expected_ds)."""
	m = mixed[0]
	ds = H.map_reads(synth.fasta, mixed[1], ["-c", "--ds"]).paf()
	cs = H.map_reads(synth.fasta, mixed[1], ["-c", "--cs"]).paf()
	assert len(ds) == len(cs)
	for d, c in zip(ds, cs):
		assert D.strip_ds(d.tags["ds"]) == c.tags["cs"]
		s = H.paf_query_part(d, m.seqs[d.qname])
		assert d.tags["ds"] == D.expected_ds(synth.genome[d.tname], d.tstart, s, d.cigar)


def test_fragment_names_in_paf(synth, index_of, tmp_path):
	"""In the fragment mode PAF appends /1 and /2 to the mates' names, also to names that already
	end in /1 or /2, while SAM drops a trailing /1 or /2 (tests/README.md)."""
	pairs = H.simulate_pairs(synth.genome, 6, seed=551, unique=True)
	f1, f2 = H.write_pairs(str(tmp_path / "a1.fq"), str(tmp_path / "a2.fq"), pairs)
	paf = H.map_reads(synth.fasta, [f1, f2], preset="sr").paf()
	assert sorted(set(p.qname for p in paf)) == sorted("%s/%d" % (p.name, k) for p in pairs for k in (1, 2))
	g1, g2 = H.write_pairs(str(tmp_path / "b1.fq"), str(tmp_path / "b2.fq"), pairs, suffix=True)
	paf = H.map_reads(synth.fasta, [g1, g2], preset="sr").paf()
	assert sorted(set(p.qname for p in paf)) == sorted("%s/%d/%d" % (p.name, k, k) for p in pairs for k in (1, 2))
	sam = H.map_reads(synth.fasta, [g1, g2], ["-a"], preset="sr").sam()
	assert sorted(set(r.qname for r in sam.records)) == sorted(p.name for p in pairs)


@pytest.fixture(scope="module")
def junc_data(tmp_path_factory):
	g, signals, dup_gene = D.junction_genome()
	d = tmp_path_factory.mktemp("d_junc")
	fa = g.write_fasta(str(d / "junc.fa"))
	reads = H.simulate_spliced_reads(g, 60, seed=561, profile="hifi", min_end_exon=40)
	fq = H.write_fastq(str(d / "junc.fq"), reads)
	return g, signals, dup_gene, fa, fq, reads


def test_write_junc(junc_data):
	"""--write-junc (minimap2.1): 6-column BED lines (contig, intron start and end, read name,
	donor plus acceptor score, strand), with GT 3, GC 2, AT 1 at donors and AG 3, AC 1 at
	acceptors; alignments with MAPQ below 10 are left out, and only primary and supplementary
	alignments with a transcript strand are written. The expected lines are derived from the SAM
	records of the same run and the genome. Reads of the first gene, which is duplicated, map
	with MAPQ 0 and are absent; introns with every scored signal (and one scored 0) occur; the
	strand is the gene's strand."""
	g, signals, dup_gene, fa, fq, reads = junc_data
	junc = H.map_reads(fa, fq, ["--write-junc"], preset="splice:hq")
	got = D.lines(junc.stdout)
	assert all(len(l.split("\t")) == 6 for l in got)
	sam = H.map_reads(fa, fq, ["-a"], preset="splice:hq").sam()
	assert got == D.expected_junctions(sam, g)
	truth = dict((r.name, r) for r in reads)
	scores = set()
	for l in got:
		ctg, s, e, name, sc, strand = l.split("\t")
		r = truth[name]
		assert r.gene.name != dup_gene
		if r.contig == ctg and (int(s), int(e)) in r.introns:
			assert strand == r.gene.strand
			d, a = signals[(ctg, int(s), int(e))]
			assert int(sc) == D.DONOR_SCORE.get(d, 0) + D.ACCEPTOR_SCORE.get(a, 0)
			scores.add(int(sc))
	assert {6, 5, 4, 2, 0} <= scores, scores
	assert any(r.gene.name == dup_gene for r in reads)
	assert any(rec.mapq < 10 and "ts" in rec.tags and rec.introns for rec in sam.records)


@pytest.mark.pinned
def test_pinned_paf_output(synth, mixed):
	"""Pins the current output: `-c --ds --paf-no-hit -y` on the mixed reads."""
	r = H.map_reads(synth.fasta, mixed[1], ["-c", "--ds", "--paf-no-hit", "-y"])
	H.check_pinned("d/mixed-c-ds-nohit-y.paf.gz", r.stdout)


@pytest.mark.pinned
def test_pinned_write_junc(junc_data):
	"""Pins the current output of --write-junc (splice:hq) on the junction genome."""
	g, signals, dup_gene, fa, fq, reads = junc_data
	H.check_pinned("d/junc-write-junc.bed", H.map_reads(fa, fq, ["--write-junc"], preset="splice:hq").stdout)
