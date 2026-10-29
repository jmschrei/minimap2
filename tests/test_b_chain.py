"""Chaining (lchain.c) and the chaining options (-m -n -g -r --no-long-join --rmq --rmq-inner
--max-chain-skip --max-chain-iter --chain-gap-scale --chain-skip-scale), the cm:i, s1:i and dv:f tags,
and the sorting helpers of misc.c.

Chaining scores are derived from brute-force minimizers (helpers_b) and the chaining score of
Li (2018) that minimap2.1 describes for -m ("the approximate number of matching bases minus a
concave gap penalty"), on a genome of random sequence without repeats."""

import math
import os

import pytest

import helpers_b as B
import minimap2_helpers as H


@pytest.fixture(scope="module")
def clean(tmp_path_factory):
	"""A 60 + 40 kbp random genome (no repeats, no N), its FASTA and brute-force minimizer index."""
	d = tmp_path_factory.mktemp("b_chain")
	g = H.make_genome(seed=21, contig_lengths=(60000, 40000))

	class D(object):
		pass
	o = D()
	o.genome, o.fasta, o.index, o.dir = g, g.write_fasta(str(d / "clean.fa")), B.MinimizerIndex(g, 10, 15), str(d)
	return o


# ---------------------------------------------------------------------------
# C driver: lchain.c and misc.c sorts
# ---------------------------------------------------------------------------


def test_collinear_anchors_form_one_chain():
	"""mg_lchain_dp and mg_lchain_rmq: 40 anchors on one diagonal form one chain holding all of them in
	order, scoring the number of query bases their 15-mers cover (no gap cost on a diagonal)."""
	H.assert_c_case("drv_b_chain", "collinear")


def test_chains_never_cross_contigs_strands_or_max_gap():
	"""Anchors on different contigs or strands never share a chain; two runs on one diagonal
	share one only when the gap between them is within max_dist (-g)."""
	H.assert_c_case("drv_b_chain", "separate")


def test_gap_cost_formula():
	"""The score across an indel of dd bases is min(dg, k) minus
	gap_scale*0.01*k*dd + skip_scale*0.01*k*dg + 0.5*log2(dd+1) (truncated), for --chain-gap-scale
	0.8, 1 and 2 and --chain-skip-scale 0 and 0.3, insertions and deletions, dp and rmq chaining;
	runs whose diagonals differ by more than the band (-r) are separate chains."""
	H.assert_c_case("drv_b_chain", "gap_cost")


def test_band_edge():
	"""The band (-r) is inclusive: runs whose diagonals differ by exactly bw form one chain, by bw+1 no
	chain crosses between them (DP chaining and RMQ chaining with --rmq-inner 1000 at bw 100 and 333;
	RMQ chaining with --rmq-inner 0 at bw 100)."""
	H.assert_c_case("drv_b_chain", "band_edge")


def test_min_count_and_score_filters():
	"""-n and -m in the chaining functions: a 12-anchor chain scoring 125 survives min_cnt 12 and
	min_sc 125 and is dropped at 13 or 126; no anchors give no chains."""
	H.assert_c_case("drv_b_chain", "min_cnt_sc")


def test_max_chain_iter():
	"""--max-chain-iter: when each true anchor is followed by a decoy, checking 2 predecessors chains
	all 20 true anchors (score 205), checking 1 does not."""
	H.assert_c_case("drv_b_chain", "max_iter")


def test_random_anchor_sets_against_brute_force():
	"""300 random anchor sets: with the heuristics switched off (--max-chain-skip and
	--max-chain-iter huge) the best chain scores the brute-force optimum; all chains are collinear on
	one contig and strand; the default heuristics, --max-chain-skip 0, --max-chain-iter 1 and RMQ
	chaining never beat it, the two extreme settings lose score on some sets, and --rmq-inner 1000
	helps RMQ chaining more often than it hurts compared with --rmq-inner 0."""
	H.assert_c_case("drv_b_chain", "random_dp", 3, 300)


def test_sort_helpers():
	"""radix_sort_128x (by .x, keeping every element), radix_sort_64 and ks_ksmall_uint32_t (k-th
	smallest) agree with qsort on sizes 0 to 70000 and on keys with few distinct values."""
	H.assert_c_case("drv_b_chain", "sorts", 5)


# ---------------------------------------------------------------------------
# Command line: cm:i, s1:i, dv:f and the chaining options
# ---------------------------------------------------------------------------


def _expected_chain(read, index, **kw):
	anc = [a for a in B.anchors(read.seq, index) if a.contig == read.contig and a.strand == read.strand]
	return B.chain_dp(anc, index.k, **kw)


def _primary(paf):
	out = {}
	for p in paf:
		if p.is_primary:
			assert p.qname not in out, p.line
			out[p.qname] = p
	return out


def _dv(read, index, p, chain):
	"""dv:f as minimap2.1 describes it ("approximate per-base sequence divergence"), estimated from
	minimizer survival: of the n query minimizers present in the reference between the chain's first
	and last anchor (plus one more at each end that has room for one), m are on the chain; then
	dv = 1 - (m/n)^(1/k), printed with 4 decimals, 0 when every one is chained."""
	k = index.k
	mini = sorted(pos for pos, h, _ in B.minimizers(read.seq, index.w, k) if index.count(h) > 0)
	cpos = sorted(a.qpos for a in chain)
	n_tot = mini.index(cpos[-1]) - mini.index(cpos[0]) + 1
	if p.qstart > k and p.tstart > k:
		n_tot += 1
	if p.qlen - p.qstart > k and p.tlen - p.tend > k:
		n_tot += 1
	m = len(chain)
	dv = 0.0 if m >= n_tot else 1.0 - (float(m) / n_tot) ** (1.0 / k)
	return "0" if dv == 0 else "%.4f" % dv


@pytest.fixture(scope="module")
def subs_reads(clean):
	"""24 reads with 2% substitutions on both strands, their FASTA, and per read the brute-force
	chain (score, anchors) on the true contig and strand."""
	reads = H.simulate_reads(clean.genome, 24, length=(800, 2500), seed=3, sub_rate=0.02, indel_rate=0, qual=None, name_prefix="s")
	fa = H.write_fasta(os.path.join(clean.dir, "subs.fa"), reads)
	exp = {}
	for r in reads:
		sc, chain, _ = _expected_chain(r, clean.index)
		assert sc == B.covered_bases(chain)
		exp[r.name] = (sc, chain)
	return reads, fa, exp


@pytest.mark.parametrize("args", [[], ["--rmq"], ["--rmq", "--rmq-inner", "0"], ["--heap-sort=yes"],
	["--max-chain-skip", "1000000", "--max-chain-iter", "1000000"]], ids=["dp", "rmq", "rmq-inner0", "heap", "no-heuristics"])
def test_chain_score_and_count_substitutions(clean, subs_reads, args):
	"""Reads with 2% substitutions (both strands): the primary's s1:i (chaining score) is the number
	of read bases covered by the chained minimizers' k-mers and cm:i the number of minimizers on the
	chain, both from brute-force minimizers; dv:f follows from the chained and skipped minimizers
	(no -c). Holds for DP and RMQ chaining, heap-sorted seeds and with the heuristics off."""
	reads, fa, exp = subs_reads
	prim = _primary(H.map_reads(clean.fasta, fa, args).paf())
	assert sorted(prim) == sorted(r.name for r in reads)
	for r in reads:
		sc, chain = exp[r.name]
		p = prim[r.name]
		assert (p.tname, p.strand) == (r.contig, r.strand), p.line
		assert (p.tags["s1"], p.tags["cm"]) == (sc, len(chain)), (r.name, p.line)
		assert p.line.split("dv:f:")[1].split("\t")[0] == _dv(r, clean.index, p, chain), (r.name, p.line)


@pytest.mark.parametrize("args,gap,skip", [(["--chain-gap-scale", "1"], 1.0, 0.0),
	(["--chain-gap-scale", "2", "--chain-skip-scale", "0.3"], 2.0, 0.3)])
def test_chain_score_with_indels(clean, args, gap, skip):
	"""Reads with substitutions and 1-8 bp indels: s1:i and cm:i equal the brute-force chaining DP
	with the gap cost scaled by --chain-gap-scale and --chain-skip-scale; reads whose chain crosses
	a gap whose cost lies within 0.005 of an integer are skipped (minimap2 approximates log2 to within
	0.005 and computes the linear term in single precision, so truncation could differ there), and
	at least 25 of 30 are checked."""
	reads = H.simulate_reads(clean.genome, 30, length=(800, 2000), seed=4, sub_rate=0.01, indel_rate=0.003, max_indel=8, qual=None, name_prefix="i")
	fa = H.write_fasta(os.path.join(clean.dir, "indel.fa"), reads)
	prim = _primary(H.map_reads(clean.fasta, fa, args).paf())
	checked = 0
	for r in reads:
		anc = [a for a in B.anchors(r.seq, clean.index) if a.contig == r.contig and a.strand == r.strand]
		sc, chain, costs = B.chain_dp(anc, 15, gap_scale=gap, skip_scale=skip)
		if any(c > 0 and min(c - math.floor(c), math.ceil(c) - c) < 0.005 for c in costs):
			continue
		checked += 1
		p = prim[r.name]
		assert (p.tags["s1"], p.tags["cm"]) == (sc, len(chain)), (r.name, p.line)
	assert checked >= 25


def test_min_chain_score_and_count_boundaries(clean, tmp_path):
	"""-m INT discards chains scoring < INT and -n INT chains of < INT minimizers: for each of three
	substitution-only reads, -m s1 and -n cm keep the read mapped and -m s1+1 or -n cm+1 leave it
	unmapped (s1, cm derived from brute-force minimizers)."""
	reads = H.simulate_reads(clean.genome, 3, length=(300, 600), seed=8, sub_rate=0.03, indel_rate=0, qual=None, name_prefix="m")
	for r in reads:
		sc, chain, _ = _expected_chain(r, clean.index)
		fa = H.write_fasta(str(tmp_path / (r.name + ".fa")), [r])
		for args, mapped in ((["-m", sc], True), (["-m", sc + 1], False), (["-n", len(chain)], True), (["-n", len(chain) + 1], False)):
			paf = H.map_reads(clean.fasta, fa, args).paf()
			got = any(p.is_primary and p.tname == r.contig for p in paf)
			assert got == mapped, (r.name, args, paf)


def _gap_read(clean, insert):
	"""A read of 2000 reference bases, `insert` unrelated bases replacing as many reference bases, and
	2000 more reference bases: both flanks on one diagonal. Returns the read and the query distance
	between the last anchor before and the first after the insert (brute force)."""
	ref = clean.genome["chr1"]
	rng = H.Rng(5)
	seq = ref[10000:12000] + rng.dna(insert) + ref[12000 + insert:14000 + insert]
	anc = [a for a in B.anchors(seq, clean.index) if a.contig == "chr1" and a.strand == "+" and a.x - a.y == 10000]
	left = max(a.y for a in anc if a.y < 2000)
	right = min(a.y for a in anc if a.y >= 2000 + insert)
	return seq, right - left


def test_max_gap_g_boundary(clean, tmp_path):
	"""-g NUM ("stop chain elongation if there are no minimizers within NUM bp"): a read with a 6 kbp
	unrelated insert between two 2 kbp reference flanks on one diagonal forms one chain covering both
	flanks when -g is at least the distance G between the anchors facing the insert, and two chains
	(one per flank) at -g G-1 and at the default 5000 (--no-long-join, so the long-join pass with its
	larger band does not reconnect them)."""
	seq, G = _gap_read(clean, 6000)
	fa = H.write_fasta(str(tmp_path / "gap.fa"), [("gap", seq)])
	for g, n in ((G, 1), (G + 500, 1), (G - 1, 2), (5000, 2)):
		paf = H.map_reads(clean.fasta, fa, ["-g", g, "--no-long-join"]).paf()
		assert len(paf) == n and all(p.tname == "chr1" and p.strand == "+" and p.tstart - p.qstart == 10000 for p in paf), (g, paf)
		if n == 1:
			assert (paf[0].qstart, paf[0].qend) == (0, len(seq)) or paf[0].qend - paf[0].qstart > 9900
		else:
			assert sorted((p.qstart < 2000, p.qend > 8000) for p in paf) == [(False, True), (True, False)]


def _deletion_read(clean, d):
	ref = clean.genome["chr1"]
	return ref[20000:22000] + ref[22000 + d:24000 + d]


def test_bandwidth_r_and_long_join(clean, tmp_path):
	"""-r NUM1[,NUM2] and --no-long-join on a read with a 600 bp deletion (the anchors' diagonals
	differ by exactly 600): with -r 600,600 one chain spans it, with -r 599,599 or --no-long-join two
	chains; with the default -r 500,20k the long-join pass (NUM2) reconnects them into one, as does
	-r 599. With -c the joined chain aligns as 2000M600D2000M and the split ones as two 2000M."""
	seq = _deletion_read(clean, 600)
	fa = H.write_fasta(str(tmp_path / "del.fa"), [("del", seq)])
	for args, n in ((["-r", "600,600"], 1), (["-r", "599,599"], 2), (["--no-long-join"], 2), ([], 1), (["-r", "599"], 1),
			(["-r", "1k,20k", "--no-long-join"], 1)):
		paf = H.map_reads(clean.fasta, fa, args).paf()
		assert len(paf) == n and all(p.tname == "chr1" and p.strand == "+" and p.is_primary for p in paf), (args, paf)
	pc = H.map_reads(clean.fasta, fa, ["-c"]).paf()
	assert [(p.qstart, p.qend, p.tstart, p.tend, p.cigar) for p in pc] == [(0, 4000, 20000, 24600, "2000M600D2000M")]
	pc = H.map_reads(clean.fasta, fa, ["-c", "--no-long-join"]).paf()
	assert sorted((p.qstart, p.qend, p.tstart, p.tend, p.cigar) for p in pc) == [(0, 2000, 20000, 22000, "2000M"),
		(2000, 4000, 22600, 24600, "2000M")]
	assert H.paf_invariant_errors(pc, clean.genome, {"del": seq}) == []


@pytest.mark.pinned
def test_chaining_heuristics_pinned(tmp_path):
	"""Pins the current output of mg_lchain_dp and mg_lchain_rmq with their heuristics on, on 150 random
	anchor sets (two contigs, both strands, collinear runs with indels and scattered decoys): DP with the
	defaults, --max-chain-skip 0, --max-chain-iter 1, and --max-chain-iter 3 with --max-chain-skip 2; RMQ
	with --rmq-inner 1000, 40 and 0. Which chains these heuristics form is documented only as
	"heuristics" (test_random_anchor_sets_against_brute_force checks that they never beat the optimum),
	so each result is pinned as one line per set and setting: the number of chains, the best score and a
	hash of every chain's score, count and anchors, independent of the order of the chains. Every chain
	is also checked to be collinear on one contig and strand."""
	out = str(tmp_path / "heuristics.tsv")
	H.assert_c_case("drv_b_chain", "heuristics", 11, 150, out)
	H.check_pinned("b/chain-heuristics.tsv.gz", H.read_text(out))
