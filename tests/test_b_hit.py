"""Area B: hit selection (hit.c): primary and secondary chains (-p -N -M --hard-mask-level
--mask-len -P), mapping quality on unique and repetitive reads, ALT contigs (--alt --alt-drop) and
tie breaking (--no-hash-name).

Chain scores and coordinates are derived from brute-force minimizers (helpers_b)."""

import os

import pytest

import helpers_b as B
import minimap2_helpers as H


def _chain(read_seq, index, contig, strand="+"):
	anc = [a for a in B.anchors(read_seq, index) if a.contig == contig and a.strand == strand]
	sc, ch, _ = B.chain_dp(anc, index.k)
	return sc, ch[0].y + 1 - ch[0].span, ch[-1].y + 1


@pytest.fixture(scope="module")
def overlap(tmp_path_factory):
	"""Read = A (3000 bases of chr1) + C (1000 unrelated bases); chr2 holds A[2000:] + C. Two chains:
	H1 on chr1 over the read's first 3000 bases and H2 on chr2 over its last 2000, overlapping by
	about 1000 bases of the read. Returns the paths, the read and both chains' (score, qs, qe)."""
	d = tmp_path_factory.mktemp("b_overlap")
	rng = H.Rng(41)
	c1, c2 = rng.dna(30000), rng.dna(30000)
	A, C = c1[10000:13000], rng.dna(1000)
	c2 = c2[:20000] + A[2000:] + C + c2[22000:]
	seqs = {"chr1": c1, "chr2": c2}
	fa = H.write_fasta(str(d / "ov.fa"), list(seqs.items()))
	read = A + C
	rq = H.write_fasta(str(d / "ovr.fa"), [("ov", read)])
	idx = B.MinimizerIndex(seqs, 10, 15)
	return fa, rq, read, _chain(read, idx, "chr1"), _chain(read, idx, "chr2")


def _lines(fa, rq, args):
	return sorted((p.tname, p.tags["tp"], p.tags["s1"], p.tags.get("s2")) for p in H.map_reads(fa, rq, args).paf())


def test_mask_level_M(overlap):
	"""-M FLOAT: H2 is secondary to H1 when overlap/shorter - uncovered/longer > FLOAT, where overlap
	is their shared read length, uncovered the part of H2 outside H1 (mm_set_parent; minimap2.1:
	"mark as secondary a chain that overlaps with a better chain by FLOAT or more of the shorter
	chain"), computed here from the brute-force chain coordinates. Just below that value H2 is a
	secondary (tp:A:S, counted in H1's s2:i; -p 0.5 so it is printed), just above it a second
	primary (tp:A:P, H1's s2:i 0)."""
	fa, rq, read, h1, h2 = overlap
	(s1, qs1, qe1), (s2, qs2, qe2) = h1, h2
	ol = qe1 - qs2
	uncov = qe2 - qe1
	frac = float(ol) / min(qe1 - qs1, qe2 - qs2) - float(uncov) / max(qe1 - qs1, qe2 - qs2)
	assert 0.1 < frac < 0.3 and s1 > s2
	sec = [("chr1", "P", s1, s2), ("chr2", "S", s2, None)]
	pri = [("chr1", "P", s1, 0), ("chr2", "P", s2, 0)]
	assert _lines(fa, rq, ["-M", "%.4f" % (frac - 0.005), "-p", "0.5"]) == sec
	assert _lines(fa, rq, ["-M", "%.4f" % (frac + 0.005), "-p", "0.5"]) == pri
	assert _lines(fa, rq, []) == pri  # default -M 0.5


def test_hard_mask_level(overlap):
	"""--hard-mask-level ("honor option -M and disable a heuristic to save unmapped subsequences"):
	the uncovered part no longer counts, so H2 is secondary when overlap/shorter > -M."""
	fa, rq, read, h1, h2 = overlap
	(s1, qs1, qe1), (s2, qs2, qe2) = h1, h2
	frac = float(qe1 - qs2) / min(qe1 - qs1, qe2 - qs2)
	assert 0.4 < frac < 0.6
	assert _lines(fa, rq, ["--hard-mask-level", "-M", "%.4f" % (frac - 0.005), "-p", "0.5"]) == [("chr1", "P", s1, s2), ("chr2", "S", s2, None)]
	assert _lines(fa, rq, ["--hard-mask-level", "-M", "%.4f" % (frac + 0.005), "-p", "0.5"]) == [("chr1", "P", s1, 0), ("chr2", "P", s2, 0)]


def test_mask_len(overlap):
	"""--mask-len NUM ("keep an alignment if dropping it leaves an unaligned region on query longer than
	NUM"): with -M 0.1 H2 is secondary only if the read part it alone covers (the uncovered length,
	from the brute-force chains) is at most NUM."""
	fa, rq, read, h1, h2 = overlap
	(s1, qs1, qe1), (s2, qs2, qe2) = h1, h2
	uncov = qe2 - qe1
	assert _lines(fa, rq, ["-M", "0.1", "-p", "0.5", "--mask-len", uncov]) == [("chr1", "P", s1, s2), ("chr2", "S", s2, None)]
	assert _lines(fa, rq, ["-M", "0.1", "-p", "0.5", "--mask-len", uncov - 1]) == [("chr1", "P", s1, 0), ("chr2", "P", s2, 0)]


def test_all_chains_P(overlap):
	"""-P ("retain all chains and don't attempt to set primary chains; -p and -N have no effect"): both
	chains are printed as tp:A:S with their scores, with or without -p 0.99 -N 1."""
	fa, rq, read, h1, h2 = overlap
	want = [("chr1", "S", h1[0], None), ("chr2", "S", h2[0], None)]
	assert _lines(fa, rq, ["-P"]) == want
	assert _lines(fa, rq, ["-P", "-p", "0.99", "-N", "1"]) == want


@pytest.fixture(scope="module")
def cn(tmp_path_factory):
	d = tmp_path_factory.mktemp("b_cn")
	g, placed = B.copy_number_genome()
	return g, placed, g.write_fasta(str(d / "cn.fa")), str(d)


def _ustart(p):
	"""Reference position of the read's unclipped leftmost base on the forward strand."""
	return p.tstart - (p.qstart if p.strand == "+" else p.qlen - p.qend)


@pytest.mark.parametrize("fam,length", [(0, 350), (1, 250), (2, 180)])
def test_repeat_reads_secondaries_and_N(cn, fam, length):
	"""A read inside one copy of a 3-, 6- or 12-copy exact repeat: the primary and every secondary lie
	on copies of the repeat (at the read's offset within the copy, on the copy's strand), all with
	MAPQ 0; there are min(copies - 1, -N) secondaries for -N 5 (default), 2 and 1; -N 0 warns and
	prints the primary only; -P prints every copy (-c gives the same placements within 10 bp)."""
	g, placed, fa, d = cn
	copies = [x for x in placed if x[0] == fam]
	_, c, s, e, strand = copies[0]
	seq = g[c][s + 20:s + 20 + length]
	rq = H.write_fasta(os.path.join(d, "rep%d.fa" % fam), [("rep%d" % fam, seq)])
	starts = set()
	for _, cc, ss, ee, st in copies:  # where the read lies in each copy, on the reference's forward strand
		off = 20 if st == strand else (ee - ss) - 20 - length
		starts.add((cc, ss + off, "+" if st == strand else "-"))
	for args, n_sec in (([], min(len(copies) - 1, 5)), (["-N", "2"], 2), (["-N", "1"], 1), (["-N", "0"], 0)):
		r = H.map_reads(fa, rq, args)
		paf = r.paf()
		assert [p.tags["tp"] for p in paf].count("P") == 1 and [p.tags["tp"] for p in paf].count("S") == n_sec, (args, paf)
		for p in paf:
			assert p.mapq == 0 and p.qend - p.qstart >= length - 20
			assert any(p.tname == cc and abs(_ustart(p) - ss) <= 2 and p.strand == st for cc, ss, st in starts), p.line
		assert len(set((p.tname, p.tstart) for p in paf)) == len(paf)
		assert bool(r.messages) == (args == ["-N", "0"])
	allc = H.map_reads(fa, rq, ["-P"]).paf()
	assert len(allc) == len(copies) and all(p.tags["tp"] == "S" for p in allc)
	for p in H.map_reads(fa, rq, ["-c"]).paf():
		assert any(p.tname == cc and abs(_ustart(p) - ss) <= 10 for cc, ss, _ in starts), p.line


def test_unique_reads_mapq(cn):
	"""Reads from unique sequence map to their true locus with MAPQ 60 and no secondary."""
	g, placed, fa, d = cn
	reads = []
	rng = H.Rng(3)
	while len(reads) < 20:
		c = g.names[rng.below(2)]
		s = rng.below(g.length(c) - 2000)
		if all(not (pc == c and ps - 2000 < s + 2000 and s < pe + 2000) for _, pc, ps, pe, _ in placed):
			reads.append(H.SimRead("u%d" % len(reads), g[c][s:s + 1500], None, [H.Segment(c, s, s + 1500, "+", 0, 1500, "1500M", [])]))
	rq = H.write_fasta(os.path.join(d, "uniq.fa"), reads)
	paf = H.map_reads(fa, rq).paf()
	assert sorted(p.qname for p in paf) == sorted(r.name for r in reads)
	truth = dict((r.name, r) for r in reads)
	for p in paf:
		assert p.is_primary and p.mapq == 60 and H.paf_placement_error(p, truth[p.qname]) <= 10, p.line


def test_no_hash_name(cn):
	"""--no-hash-name ("produce the same alignment for identical sequences regardless of their sequence
	names"): 30 copies of one read from a 3-copy repeat under different names all get the same
	primary locus with --no-hash-name; without it the name breaks the tie between the equal copies, so
	the primaries spread over more than one copy."""
	g, placed, fa, d = cn
	_, c, s, e, _ = [x for x in placed if x[0] == 0][0]
	seq = g[c][s + 20:s + 370]
	rq = H.write_fasta(os.path.join(d, "same.fa"), [("name%d" % i, seq) for i in range(30)])
	def loci(args):
		return set((p.tname, p.tstart, p.strand) for p in H.map_reads(fa, rq, args).paf() if p.is_primary)
	assert len(loci(["--no-hash-name"])) == 1
	assert len(loci(["--no-hash-name", "-c"])) == 1
	assert len(loci([])) > 1


@pytest.fixture(scope="module")
def divergent(tmp_path_factory):
	"""A 2000 bp unit on chr1 and a copy with 3% substitutions on chr2; read = the chr1 unit. Returns
	paths and the brute-force chain scores of the read on chr1 (S1) and chr2 (S2)."""
	d = tmp_path_factory.mktemp("b_div")
	rng = H.Rng(46)
	c1, c2 = rng.dna(20000), rng.dna(20000)
	U = c1[8000:10000]
	c2 = c2[:5000] + H.mutate(rng, U, 0.03) + c2[7000:]
	seqs = {"chr1": c1, "chr2": c2}
	fa = H.write_fasta(str(d / "div.fa"), list(seqs.items()))
	rq = H.write_fasta(str(d / "divr.fa"), [("u", U)])
	idx = B.MinimizerIndex(seqs, 10, 15)
	return fa, rq, _chain(U, idx, "chr1")[0], _chain(U, idx, "chr2")[0]


def test_pri_ratio_p(divergent):
	"""-p FLOAT ("minimal secondary-to-primary score ratio to output secondary mappings"): the read's
	primary is the exact copy (s1:i = S1) and its s2:i is S2, the diverged copy's chaining score; the
	secondary on the diverged copy is printed when -p is just below S2/S1 and not when just above
	(S1 - S2 exceeds the 2k tolerance minimap2 also applies)."""
	fa, rq, S1, S2 = divergent
	ratio = float(S2) / S1
	assert 0.5 < ratio < 0.95 and S1 - S2 > 30
	lo = _lines(fa, rq, ["-p", "%.4f" % (ratio - 0.002)])
	hi = _lines(fa, rq, ["-p", "%.4f" % (ratio + 0.002)])
	assert lo == [("chr1", "P", S1, S2), ("chr2", "S", S2, None)]
	assert hi == [("chr1", "P", S1, S2)]
	assert _lines(fa, rq, ["-p", "0"]) == lo and _lines(fa, rq, ["-p", "1"]) == hi


@pytest.fixture(scope="module")
def alt(tmp_path_factory):
	"""chr1 (40 kbp) and chr1_alt, a copy of chr1:10000-20000 with 0.5% substitutions; 12 reads
	simulated from chr1_alt (0.2% substitutions)."""
	d = tmp_path_factory.mktemp("b_alt")
	rng = H.Rng(43)
	c1 = rng.dna(40000)
	g = H.Genome([("chr1", c1), ("chr1_alt", H.mutate(rng, c1[10000:20000], 0.005))])
	fa = g.write_fasta(str(d / "alt.fa"))
	lst = H.write_text(str(d / "alt.txt"), "chr1_alt\n")
	reads = H.simulate_reads(g, 12, length=(800, 1500), seed=4, sub_rate=0.002, indel_rate=0, contigs=["chr1_alt"], qual=None, name_prefix="a")
	return g, fa, lst, reads, H.write_fasta(str(d / "altr.fa"), reads), B.MinimizerIndex(g, 10, 15)


def test_alt_contigs(alt):
	"""--alt FILE and --alt-drop FLOAT ("drop ALT hits by FLOAT fraction when ranking and computing
	mapping quality"): for each read (from the ALT copy) the primary is on chr1_alt without --alt or with
	--alt-drop 0, and on the non-ALT contig whenever its chaining score exceeds the ALT score reduced
	by 15% (round(0.85 * score)), which holds for every read here; MAPQ is never lower with --alt."""
	g, fa, lst, reads, rq, idx = alt
	plain = dict((p.qname, p) for p in H.map_reads(fa, rq).paf() if p.is_primary)
	with_alt = dict((p.qname, p) for p in H.map_reads(fa, rq, ["--alt", lst]).paf() if p.is_primary)
	no_drop = dict((p.qname, p) for p in H.map_reads(fa, rq, ["--alt", lst, "--alt-drop", "0"]).paf() if p.is_primary)
	for r in reads:
		s_alt = _chain(r.seq, idx, "chr1_alt", r.strand)[0]
		s_ref = _chain(r.seq, idx, "chr1", r.strand)[0]
		assert s_alt > s_ref > int(0.85 * s_alt + .499), (r.name, s_alt, s_ref)
		assert plain[r.name].tname == "chr1_alt" and no_drop[r.name].tname == "chr1_alt"
		assert with_alt[r.name].tname == "chr1" and with_alt[r.name].tags["s1"] == s_ref
		assert with_alt[r.name].mapq >= plain[r.name].mapq
		assert (no_drop[r.name].line == plain[r.name].line)
