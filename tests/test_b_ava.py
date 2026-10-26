"""Area B: all-vs-all overlapping (-x ava-ont, -X, -D, -P, --dual, --no-long-join), strand
restrictions (--for-only, --rev-only, --qstrand), --max-qlen and --heap-sort, checked against the
simulated reads' positions."""

import os

import pytest

import helpers_b as B
import minimap2_helpers as H


@pytest.fixture(scope="module")
def ava(tmp_path_factory):
	"""30 ONT-like reads of 3-6 kbp from one 30 kbp random contig (so most pairs overlap), the FASTA,
	and the true overlap length of every overlapping pair (keyed by the sorted pair of names)."""
	d = tmp_path_factory.mktemp("b_ava")
	g = H.make_genome(seed=31, contig_lengths=(30000,))
	reads = H.simulate_long_reads(g, 30, profile="ont", seed=8, length=(3000, 6000), name_prefix="r")
	fa = H.write_fasta(str(d / "avar.fa"), reads)
	truth = {}
	for i, a in enumerate(reads):
		for b in reads[i + 1:]:
			o = min(a.end, b.end) - max(a.pos, b.pos)
			if o > 0:
				truth[tuple(sorted((a.name, b.name)))] = (o, a.strand == b.strand)
	return reads, fa, truth


def _pairs(paf):
	out = {}
	for p in paf:
		out.setdefault((p.qname, p.tname), []).append(p)
	return out


def test_ava_ont_reports_true_overlaps(ava):
	"""-x ava-ont on reads against themselves: no self-mapping; each pair at most once, with the query
	name lexicographically smaller (--dual=no); every pair overlapping by at least 1 kbp is reported,
	none that does not overlap; the strand is '+' exactly when the reads come from the same strand; the
	query interval's length is within 20% (+100 bp) of the true overlap."""
	reads, fa, truth = ava
	paf = H.map_reads(fa, fa, [], preset="ava-ont").paf()
	pairs = _pairs(paf)
	assert all(q < t for q, t in pairs), [k for k in pairs if k[0] >= k[1]]
	assert all(len(v) == 1 for v in pairs.values())
	for k, (o, same) in truth.items():
		if o >= 1000:
			assert k in pairs, (k, o)
	for k, v in pairs.items():
		assert k in truth, k
		o, same = truth[k]
		p = v[0]
		assert (p.strand == "+") == same, p.line
		assert abs((p.qend - p.qstart) - o) <= 0.2 * o + 100, (p.line, o)


def test_ava_dual_yes_reports_both_directions(ava):
	"""--dual=yes (minimap2.1: with no, skip pairs whose query name is greater than the target name):
	every pair of -x ava-ont appears in both directions."""
	reads, fa, truth = ava
	one = set(_pairs(H.map_reads(fa, fa, [], preset="ava-ont").paf()))
	both = set(_pairs(H.map_reads(fa, fa, ["--dual=yes"], preset="ava-ont").paf()))
	assert both == one | set((t, q) for q, t in one)


def test_X_equals_its_expansion(ava):
	"""-X is "equivalent to '-DP --dual=no --no-long-join'" (minimap2.1): byte-identical PAF."""
	reads, fa, truth = ava
	a = H.map_reads(fa, fa, ["-X"]).stdout
	assert a and a == H.map_reads(fa, fa, ["-D", "-P", "--dual=no", "--no-long-join"]).stdout


def test_D_drops_self_mappings(ava):
	"""-D ("if query sequence name/length are identical to the target name/length, ignore diagonal
	anchors"): with -P every read maps onto itself end to end on '+', with -D -P no read does, and the
	other pairs keep being reported."""
	reads, fa, truth = ava
	with_self = _pairs(H.map_reads(fa, fa, ["-P"]).paf())
	no_self = _pairs(H.map_reads(fa, fa, ["-P", "-D"]).paf())
	for r in reads:
		s = with_self[(r.name, r.name)]
		assert any(p.strand == "+" and p.qstart < 30 and p.qend > len(r) - 30 and p.tstart == p.qstart for p in s), s
		assert (r.name, r.name) not in no_self
	assert set(k for k in with_self if k[0] != k[1]) == set(no_self)


@pytest.fixture(scope="module")
def stranded(tmp_path_factory):
	"""20 HiFi-like reads (2-3 kbp) from a random 100 kbp genome, about half on each strand."""
	d = tmp_path_factory.mktemp("b_strand")
	g = H.make_genome(seed=21, contig_lengths=(60000, 40000))
	reads = H.simulate_long_reads(g, 20, profile="hifi", seed=7, length=(2000, 3000))
	assert 5 <= sum(r.strand == "+" for r in reads) <= 15
	return g, g.write_fasta(str(d / "g.fa")), reads, H.write_fasta(str(d / "r.fa"), reads)


def test_for_only_rev_only(stranded):
	"""--for-only maps reads only to the forward strand of the reference and --rev-only only to the
	reverse strand (minimap2.1): each read from the allowed strand maps at its locus (within 10 bp,
	-c) and no read from the other strand maps at all; together they are rejected (exit 1)."""
	g, fa, reads, rq = stranded
	for opt, strand in (("--for-only", "+"), ("--rev-only", "-")):
		paf = H.map_reads(fa, rq, ["-c", opt]).paf()
		assert all(p.strand == strand for p in paf)
		names = set(p.qname for p in paf if p.is_primary)
		for r in reads:
			if r.strand == strand:
				p = [x for x in paf if x.qname == r.name and x.is_primary][0]
				assert H.paf_placement_error(p, r) <= 10, (opt, p.line[:150], r)
			else:
				assert r.name not in names
	r = H.map_reads(fa, rq, ["--for-only", "--rev-only"], check=False)
	assert r.returncode == 1 and r.messages and r.stdout == ""


def test_qstrand(stranded):
	"""--qstrand (align on the query's strand): the first 12 PAF columns equal the default output,
	without and with -c; with -c, a '+' hit's cg:Z and cs:Z are unchanged and a '-' hit's describe the query, forward,
	against the reverse complement of the target interval (recomputed from the genome), with the same
	AS:i and NM:i. --qstrand is rejected with -a, -H, --frag=yes and --splice (exit 1)."""
	g, fa, reads, rq = stranded
	seqs = dict((r.name, r.seq) for r in reads)
	plain = H.map_reads(fa, rq, ["-c", "--cs"]).paf()
	qs = H.map_reads(fa, rq, ["-c", "--cs", "--qstrand"]).paf()
	assert len(plain) == len(qs) >= 20
	for a, b in zip(plain, qs):
		assert a.line.split("\t")[:12] == b.line.split("\t")[:12]
		assert (a.tags["AS"], a.tags["NM"]) == (b.tags["AS"], b.tags["NM"])
		if a.strand == "+":
			assert (a.cigar, a.tags["cs"]) == (b.cigar, b.tags["cs"])
		else:
			t = H.revcomp(g[b.tname][b.tstart:b.tend])
			q = seqs[b.qname][b.qstart:b.qend]
			assert H.cigar_query_span(b.cigar) == len(q) and H.cigar_ref_span(b.cigar) == len(t)
			assert H.compute_cs(t, 0, q, b.cigar) == b.tags["cs"] and H.compute_nm(t, 0, q, b.cigar) == b.tags["NM"]
	cols = lambda out: [l.split("\t")[:12] for l in out.splitlines()]
	assert cols(H.map_reads(fa, rq, ["--qstrand"]).stdout) == cols(H.map_reads(fa, rq).stdout)
	for bad in (["-a"], ["-H"], ["--frag=yes"], ["--splice"]):
		r = H.map_reads(fa, rq, ["--qstrand"] + bad, check=False)
		assert r.returncode == 1 and r.messages and r.stdout == "", bad


def test_max_qlen(stranded, tmp_path):
	"""--max-qlen NUM ("filter out query sequences longer than NUM"): a read of length L maps with
	--max-qlen L and is left out with L-1, in PAF (no line) and SAM (flag 4)."""
	g, fa, reads, rq = stranded
	r = reads[0]
	one = H.write_fasta(str(tmp_path / "one.fa"), [r])
	assert len(H.map_reads(fa, one, ["--max-qlen", len(r)]).paf()) >= 1
	assert H.map_reads(fa, one, ["--max-qlen", len(r) - 1]).stdout == ""
	sam = H.map_reads(fa, one, ["-a", "--max-qlen", len(r) - 1]).sam()
	assert [x.flag for x in sam.records] == [4]


@pytest.mark.parametrize("preset,kind", [("map-ont", "ont"), ("map-hifi", "hifi"), ("asm20", "ctg"), ("lr:hq", "hifi")])
def test_heap_sort_same_output(synth, tmp_path, preset, kind):
	"""--heap-sort=yes ("sort anchors with heap merge, instead of radix sort") yields the same anchors,
	hence byte-identical SAM (@PG aside), on reads anywhere on the standard synthetic genome
	(repeats, N runs, low-complexity tracts)."""
	g = synth.genome
	if kind == "ctg":
		reads = H.simulate_contigs(g, 2, length=(10000, 20000), seed=301, inversions=1, inversion_len=(800, 1500))
	else:
		reads = H.simulate_long_reads(g, 15, profile=kind, seed=302, length=(1000, 4000))
	fa = H.write_fasta(str(tmp_path / "r.fa"), reads)
	a = H.map_reads(synth.fasta, fa, ["-a"], preset=preset).stdout
	b = H.map_reads(synth.fasta, fa, ["-a", "--heap-sort=yes"], preset=preset).stdout
	assert H.strip_pg(a) == H.strip_pg(b) and a.count("\n") > 10
	c = H.map_reads(synth.fasta, fa, ["-a", "--heap-sort=no"], preset=preset).stdout
	assert H.strip_pg(a) == H.strip_pg(c)


def test_ava_min_chain_score(tmp_path):
	"""-x ava-ont includes -m100 (minimap2.1): ten dovetail pairs sharing 92-110 bases (read A ends
	with the shared piece, read B starts with it) give overlaps whose chaining score, from brute-force
	k15/w5 minimizers of the reads, falls on both sides of 100; exactly the pairs scoring at least 100
	are reported, with that score as s1:i."""
	rng = H.Rng(71)
	recs = []
	for L in range(92, 112, 2):
		S = rng.dna(L)
		recs += [("ovA%d" % L, rng.dna(400) + S), ("ovB%d" % L, S + rng.dna(400))]
	fa = H.write_fasta(str(tmp_path / "ov.fa"), recs)
	seqs = dict(recs)
	idx = B.MinimizerIndex(seqs, 5, 15)
	want = {}
	for L in range(92, 112, 2):
		a, b = "ovA%d" % L, "ovB%d" % L
		anc = [x for x in B.anchors(seqs[a], idx) if x.contig == b and x.strand == "+"]
		want[(a, b)] = B.chain_dp(anc, 15)[0]
	assert min(want.values()) < 100 <= max(want.values())
	got = dict(((p.qname, p.tname), p.tags["s1"]) for p in H.map_reads(fa, fa, [], preset="ava-ont").paf())
	assert got == dict((k, v) for k, v in want.items() if v >= 100), (got, want)
