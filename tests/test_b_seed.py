"""Area B: seeding (seed.c), the minimizer occurrence threshold (-f, -U, --min-occ-floor, -e,
--q-occ-frac; mm_idx_cal_max_occ, mm_mapopt_update) and ALT contigs (--alt; mm_idx_alt_read).

Expected values come from brute-force minimizers (helpers_b.minimizers, minimap2's documented
(w,k) rule) on a genome with repeat families of 3, 6 and 12 exact copies."""

import os

import pytest

import helpers_b as B
import minimap2_helpers as H


@pytest.fixture(scope="module")
def cn(tmp_path_factory):
	"""The copy-number genome, its brute-force minimizer index (w10 k15), and one query per repeat
	family: the family's first copy with 700 unique bases on each side."""
	d = tmp_path_factory.mktemp("b_seed")
	g, placed = B.copy_number_genome()
	fa = g.write_fasta(str(d / "cn.fa"))
	idx = B.MinimizerIndex(g, 10, 15)
	queries = []
	for fam in range(3):
		_, c, s, e, _ = [x for x in placed if x[0] == fam][0]
		q = g[c][s - 700:e + 700]
		path = H.write_fasta(str(d / ("q%d.fa" % fam)), [("q%d" % fam, q)])
		queries.append((path, q, c, s - 700))

	class D(object):
		pass
	o = D()
	o.genome, o.placed, o.fasta, o.index, o.queries, o.dir = g, placed, fa, idx, queries, str(d)
	return o


def test_query_occurrence_filter():
	"""mm_seed_mz_flt (--q-occ-frac): a query minimizer is dropped when its count in the query
	exceeds both the reference occurrence threshold and FLOAT times the number of query minimizers
	(minimap2.1); boundaries checked on both thresholds; 0 disables the filter."""
	H.assert_c_case("drv_b_seed", "mz_flt")


@pytest.mark.parametrize("max_occ,max_max_occ,dist", [(2, 4095, 0), (5, 4095, 0), (11, 4095, 0), (100, 4095, 0),
	(2, 4095, 500), (2, 4095, 100), (2, 8, 100), (5, 4095, 50)])
def test_collect_matches_against_brute_force(cn, max_occ, max_max_occ, dist):
	"""mm_seed_collect_all/mm_seed_select/mm_collect_matches on queries spanning a 3-, 6- and
	12-copy repeat: the number of anchors (sum of occurrences of the kept minimizers), the number of
	kept query minimizers and, without sampling, the repetitive length (the union of the dropped
	minimizers' k-mers) equal the brute-force values for each occurrence threshold; with -e
	sampling (dist > 0) the least frequent minimizers of each repetitive run are kept, one per dist
	bases, unless above max_max_occ."""
	for path, q, _, _ in cn.queries:
		n_a, n_kept, rep = B.collect_expectation(q, cn.index, max_occ, max_max_occ, dist)
		H.assert_c_case("drv_b_seed", "collect", cn.fasta, 10, 15, path, max_occ, max_max_occ, dist, n_a, n_kept,
			-1 if rep is None else rep)


def test_cal_max_occ(cn):
	"""mm_idx_cal_max_occ (-f FLOAT, "ignore top FLOAT fraction of most frequent minimizers"):
	one more than the count at rank floor((1-f) n) of the sorted per-minimizer counts, for fractions
	landing on counts 1, 3, 6 and 12; f <= 0 disables the threshold (INT32_MAX)."""
	args = []
	for f in (0.0002, 0.001, 0.005, 0.01, 0.05, 0.5, 0.9):
		args += [f, cn.index.mid_occ(f)]
	assert len(set(args[1::2])) >= 4
	H.assert_c_case("drv_b_seed", "cal_max_occ", cn.fasta, 10, 15, *args)


@pytest.mark.parametrize("f,lo,hi", [(0.0002, 10, 1000000), (0.0002, 3, 6), (0.01, 2, 1000), (0.01, 5, 1000),
	(0.005, 1, 7), (0.05, 3, 100)])
def test_mapopt_update_occurrence_threshold(cn, f, lo, hi):
	"""mm_mapopt_update: mid_occ = max{-U INT1, min{-U INT2, -f threshold}} (minimap2.1, -U, for
	INT2 > INT1); a mid_occ set directly (-f INT) is kept; bw_long is raised to bw; the splice strand
	flags imply MM_F_SPLICE."""
	want = max(lo, min(hi, cn.index.mid_occ(f)))
	H.assert_c_case("drv_b_seed", "mapopt_update", cn.fasta, 10, 15, f, lo, hi, want)


def _rl_and_cm(cn, args, q_index=None):
	"""Map each repeat query with map-ont defaults plus args; return {qname: (rl, cm of the primary)}."""
	out = {}
	for i, (path, q, c, s) in enumerate(cn.queries):
		if q_index is not None and i != q_index:
			continue
		paf = H.map_reads(cn.fasta, path, args).paf()
		prim = [p for p in paf if p.is_primary]
		assert len(prim) == 1, paf
		p = prim[0]
		assert p.tname == c and p.strand == "+" and abs(p.tstart - p.qstart - s) <= 10, p.line
		out[p.qname] = (p.tags["rl"], p.tags["cm"])
	return out


def test_f_int_sets_occurrence_threshold(cn):
	"""-f INT: minimizers occurring more than INT times are ignored, so rl:i (the length of query
	regions harbouring repetitive seeds) is the union of the k-mers of the query minimizers occurring
	more than INT times (-e 0: no sampling); with the default threshold (13 here: max{10, -f 0.0002})
	no minimizer of a 12-copy repeat is repetitive."""
	for f, mid in ((2, 2), (5, 5), (11, 11)):
		got = _rl_and_cm(cn, ["-f", str(f), "-e", "0"])
		for path, q, _, _ in cn.queries:
			name = os.path.basename(path)[:-3]
			assert got[name][0] == B.repeat_len(q, cn.index, mid), (f, name, got[name])
	assert cn.index.mid_occ(0.0002) == 13
	assert all(rl == 0 for rl, _ in _rl_and_cm(cn, []).values())


@pytest.mark.parametrize("args,mid", [(["-U", "5,1000", "-f", "0.01"], 5), (["-U", "3,1000", "-f", "0.01"], 4),
	(["-U", "3,6"], 6), (["--min-occ-floor", "5", "-f", "0.01"], 5), (["-U", "2", "-f", "0.05"], 2), (["-U", "2,4", "-f", "0.001"], 4)])
def test_occurrence_bounds_U(cn, args, mid):
	"""-U INT1[,INT2] and --min-occ-floor: the threshold is max{INT1, min{INT2, -f}} (minimap2.1);
	checked through rl:i, which equals the union of the minimizers occurring more than the threshold
	(-e 0). The -f fractions give thresholds 4 (0.01), 2 (0.05) and 13 (0.0002, 0.001)."""
	got = _rl_and_cm(cn, args + ["-e", "0"])
	for path, q, _, _ in cn.queries:
		name = os.path.basename(path)[:-3]
		assert got[name][0] == B.repeat_len(q, cn.index, mid), (args, name, got[name])


def test_e_samples_repetitive_minimizers(cn):
	"""-e INT ("sample a high-frequency minimizer every INT basepairs"): with -f 2 the minimizers of the
	3-, 6- and 12-copy repeats are repetitive; -e 0 drops them all, -e 100 keeps round(span/100) of the
	least frequent per repetitive run. Each kept one adds exactly one anchor to the read's chain (its
	copy on the true diagonal), so cm:i grows by the extra kept minimizers and rl:i shrinks."""
	off = _rl_and_cm(cn, ["-f", "2", "-e", "0"])
	on = _rl_and_cm(cn, ["-f", "2", "-e", "100"])
	grew = 0
	for path, q, _, _ in cn.queries:
		name = os.path.basename(path)[:-3]
		_, k0, _ = B.collect_expectation(q, cn.index, 2, 4095, 0)
		_, k1, _ = B.collect_expectation(q, cn.index, 2, 4095, 100)
		assert on[name][1] - off[name][1] == k1 - k0, (name, on[name], off[name], k0, k1)
		assert on[name][0] <= off[name][0]
		grew += k1 > k0
	assert grew == 3


def test_q_occ_frac_drops_query_repeats(tmp_path):
	"""--q-occ-frac: a read made of 15 copies of a 300 bp segment S separated by unrelated 100 bp
	spacers, against a reference holding S once. Each S minimizer occurs 15 times in the read, above
	the reference threshold (10, the -U floor) and above 1% of the read's minimizers, so by default all
	are discarded and the read is unmapped. --q-occ-frac 0 (disabled), --q-occ-frac 0.2 (15 is below
	20% of the minimizers) and -f 20 (threshold above 15) each keep them and the read maps to S."""
	rng = H.Rng(77)
	S = rng.dna(300)
	ref = rng.dna(5000) + S + rng.dna(5000)
	read = "".join(S + rng.dna(100) for _ in range(15))
	fa = H.write_fasta(str(tmp_path / "ref.fa"), [("ref", ref)])
	fq = H.write_fasta(str(tmp_path / "read.fa"), [("read", read)])
	n_min = len(B.minimizers(read, 10, 15))
	assert 0.01 * n_min < 15 < 0.2 * n_min
	assert H.map_reads(fa, fq).stdout == ""
	for args in (["--q-occ-frac", "0"], ["--q-occ-frac", "0.2"], ["-f", "20"]):
		paf = H.map_reads(fa, fq, args).paf()
		assert paf and all(p.tname == "ref" and 4990 <= p.tstart and p.tend <= 5310 for p in paf), (args, paf)


def test_alt_read_api(tmp_path, synth):
	"""mm_idx_alt_read: names in the list (first word per line, extra columns allowed) that are
	targets become ALT and are counted; unknown names are ignored; a missing file returns < 0."""
	alt = H.write_text(str(tmp_path / "alt.txt"), "chr2\textra column\nno-such-contig\nshort1\n")
	H.assert_c_case("drv_b_seed", "alt_read", synth.fasta, alt, 2, "chr2", "short1")


@pytest.fixture(scope="module")
def many(tmp_path_factory):
	"""helpers_b.many_copy_genome() (a 120-copy family at 2% divergence), its brute-force minimizer index
	and four queries: a copy with 600 unique bases on each side."""
	d = tmp_path_factory.mktemp("b_many")
	g, copies = B.many_copy_genome()
	idx = B.MinimizerIndex(g, 10, 15)
	queries = []
	for i, (c, s, e, _) in enumerate(copies[::30]):
		queries.append(("m%d" % i, g[c][s - 600:e + 600], c, s - 600))
	assert len(queries) == 4
	return g, idx, g.write_fasta(str(d / "many.fa")), H.write_fasta(str(d / "mq.fa"), [(n, q) for n, q, _, _ in queries]), queries, str(d)


@pytest.mark.parametrize("args,thres", [(["-f", "0.0002"], None), (["-f", "0.002"], None), (["-f", "50"], 50), (["-f", "85"], 85),
	(["-U", "88,1000000", "-f", "0.002"], 88), (["-U", "10,70", "-f", "0.0002"], 70)])
def test_many_copy_occurrence_thresholds(many, args, thres):
	"""High-occurrence seeds: on a genome whose 120-copy family has minimizers occurring up to 96 times,
	-f FLOAT (the brute-force count at rank floor((1-FLOAT)n), plus one), -f INT and -U set a threshold
	inside the family's counts; with -e 0 each query's rl:i is the union of the k-mers of its minimizers
	occurring more often (minimap2.1: -f, -U), and the copy plus its unique flanks still maps at its place."""
	g, idx, fa, qfa, queries, _ = many
	if thres is None:
		thres = max(10, idx.mid_occ(float(args[1])))
	fam = [len(v) for v in idx.occ.values() if len(v) > 20]
	assert min(fam) <= thres < max(fam)  # the threshold splits the family's minimizers
	paf = H.map_reads(fa, qfa, args + ["-e", "0"]).paf()
	prim = dict((p.qname, p) for p in paf if p.is_primary)
	n_rep = 0
	for name, q, c, s in queries:
		p = prim[name]
		assert p.tname == c and p.strand == "+" and abs(p.tstart - p.qstart - s) <= 10, p.line
		assert p.tags["rl"] == B.repeat_len(q, idx, thres), (args, name, p.tags["rl"])
		n_rep += p.tags["rl"] > 0
	assert n_rep >= 2


@pytest.mark.parametrize("max_occ,max_max_occ,dist", [(80, 4095, 0), (60, 4095, 50), (85, 90, 30), (40, 4095, 500)])
def test_collect_matches_many_copies(many, max_occ, max_max_occ, dist):
	"""mm_collect_matches with high-occurrence minimizers (counts up to 96): anchors, kept minimizers and
	(without sampling) the repetitive length equal the brute-force values; with -e sampling the least
	frequent repetitive minimizers of each run are kept, one per dist bases, unless above max_max_occ."""
	g, idx, fa, _, queries, d = many
	for name, q, _, _ in queries:
		path = H.write_fasta(str(d) + "/%s.fa" % name, [(name, q)])
		n_a, n_kept, rep = B.collect_expectation(q, idx, max_occ, max_max_occ, dist)
		H.assert_c_case("drv_b_seed", "collect", fa, 10, 15, path, max_occ, max_max_occ, dist, n_a, n_kept, -1 if rep is None else rep)


@pytest.mark.pinned
def test_sampled_minimizers_pinned(cn, tmp_path):
	"""Pins the current output: which high-occurrence minimizers -e sampling keeps when several are
	equally frequent (mm_seed_select; minimap2.1 says only "sample a high-frequency minimizer every INT
	basepairs"). On the queries spanning the 3-, 6- and 12-copy exact repeats every minimizer of a family
	occurs equally often, so the sample is decided by the tie-break; the query positions of the kept
	minimizers are pinned for four threshold and spacing settings, after the counts are checked against
	brute force as in test_collect_matches_against_brute_force."""
	out = []
	for path, q, _, _ in cn.queries:
		for max_occ, max_max_occ, dist in ((2, 4095, 500), (2, 4095, 100), (2, 8, 100), (5, 4095, 50)):
			n_a, n_kept, _ = B.collect_expectation(q, cn.index, max_occ, max_max_occ, dist)
			pos = str(tmp_path / "pos.txt")
			H.assert_c_case("drv_b_seed", "collect", cn.fasta, 10, 15, path, max_occ, max_max_occ, dist, n_a, n_kept, -1, pos)
			out.append("# %s %d %d %d\n" % (os.path.basename(path), max_occ, max_max_occ, dist) + H.read_text(pos))
	H.check_pinned("b/sampled-minimizers.txt.gz", "".join(out))
