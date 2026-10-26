"""Area B: the mapping parameters and results of the libminimap2 API (minimap.h): mm_map with the
mm_mapopt_t fields and MM_F_* flags set directly (tests/c/drv_b_api.c) must give what the command
line gives for the corresponding options, field for field of mm_reg1_t and mm_extra_t against the
PAF columns and tags; mm_set_opt sets each preset's mode flags; and the fields without an option
(rank_min_len, rank_frac, rmq_rescue_size, rmq_rescue_ratio, anchor_ext_len, max_max_occ) have their
documented effects."""

import math
import os

import pytest

import helpers_b as B
import minimap2_helpers as H

COLS = ("qname qlen qs qe strand tname tlen rs re mlen blen mapq id parent sam_pri inv split split_inv is_alt cnt score "
	"score0 subsc n_sub div rep_len has_p dp_score dp_max dp_max0 dp_max2 n_ambi trans_strand n_cigar capacity cigar").split()


class Hit(object):
	def __init__(self, line):
		f = line.split("\t")
		assert len(f) == len(COLS), line
		for k, v in zip(COLS, f):
			setattr(self, k, v if k in ("qname", "strand", "tname", "cigar") else float(v) if k == "div" else int(v))

	def as_paf(self, qstrand=False):
		"""The PAF fields minimap2 writes for this hit (format.c, minimap2.1's tag table). In the
		query-strand mode the target interval of a '-' hit is on the reverse strand, and PAF converts
		it to the forward strand."""
		tp = ("I" if self.inv else "P") if self.id == self.parent else ("i" if self.inv else "S")
		rs, re = (self.tlen - self.re, self.tlen - self.rs) if qstrand and self.strand == "-" else (self.rs, self.re)
		t = [self.qname, self.qlen, self.qs, self.qe, self.strand, self.tname, self.tlen, rs, re, self.mlen, self.blen,
			self.mapq, tp, self.cnt, self.score, self.subsc if self.id == self.parent else None]
		if self.has_p:
			t += [self.blen - self.mlen + self.n_ambi, self.dp_max0, self.dp_score, self.n_ambi, self.cigar]
		else:
			t += [None] * 5
		return t + [self.split or None, self.rep_len]


def paf_fields(p):
	return [p.qname, p.qlen, p.qstart, p.qend, p.strand, p.tname, p.tlen, p.tstart, p.tend, p.nmatch, p.alen, p.mapq, p.tags["tp"],
		p.tags["cm"], p.tags["s1"], p.tags.get("s2"), p.tags.get("NM"), p.tags.get("ms"), p.tags.get("AS"), p.tags.get("nn"),
		p.tags.get("cg"), p.tags.get("zd"), p.tags["rl"]]


def api_map(target, reads, preset="-", settings=()):
	r = H.assert_c_case("drv_b_api", "map", target, reads, preset, *settings)
	return [Hit(l) for l in r.stdout.splitlines() if l and not l.startswith(("#", "ok", "not ok"))]


@pytest.fixture(scope="module")
def data(synth, tmp_path_factory):
	"""ONT-like reads anywhere on the standard synthetic genome (so some are repetitive or split), plus
	reads with an inversion, a 600 bp deletion and a chimeric read."""
	d = tmp_path_factory.mktemp("b_api")
	g = synth.genome
	reads = H.simulate_long_reads(g, 20, profile="ont", seed=401, length=(1000, 4000), qual=None)
	ref = g["chr1"]
	extra = [("inv", ref[30000:32000] + H.revcomp(ref[32000:32400]) + ref[32400:34400]),
		("del600", ref[20000:22000] + ref[22600:24600])]
	extra += [(c.name, c.seq) for c in H.simulate_chimeric_reads(g, 3, profile="hifi", length=(2000, 4000), seed=402)]
	fa = H.write_fasta(str(d / "api.fa"), [(r.name, r.seq) for r in reads] + extra)
	return fa


def _compare(synth, fa, preset, settings, cli):
	api = [h.as_paf("+QSTRAND" in settings) for h in api_map(synth.fasta, fa, preset, settings)]
	paf = [paf_fields(p) for p in H.map_reads(synth.fasta, fa, (["-x", preset] if preset != "-" else []) + cli).paf()]
	assert api == paf, (settings, cli, [x for x in api if x not in paf][:3], [x for x in paf if x not in api][:3])
	return api


# (mm_mapopt_t settings and flags for the API, the command-line options that set them)
PAIRS = [
	([], []),
	(["+CIGAR"], ["-c"]),
	(["+FOR_ONLY", "+CIGAR"], ["--for-only", "-c"]),
	(["+REV_ONLY"], ["--rev-only"]),
	(["+NO_LJOIN", "+CIGAR"], ["--no-long-join", "-c"]),
	(["+RMQ", "rmq_inner_dist=100"], ["--rmq", "--rmq-inner", "100"]),
	(["+HEAP_SORT", "+NO_HASH_NAME"], ["--heap-sort=yes", "--no-hash-name"]),
	(["+QSTRAND", "+NO_INV", "+CIGAR"], ["--qstrand", "-c"]),
	(["+HARD_MLEVEL", "mask_level=0.3", "pri_ratio=0.5", "best_n=2"], ["--hard-mask-level", "-M", "0.3", "-p", "0.5", "-N", "2"]),
	(["mask_len=500", "mask_level=0.2"], ["--mask-len", "500", "-M", "0.2"]),
	(["+ALL_CHAINS", "+NO_DIAG"], ["-P", "-D"]),
	(["+NO_END_FLT", "+CIGAR", "min_ksw_len=50", "end_bonus=5"], ["--no-end-flt", "-c", "--min-dp-len", "50", "--end-bonus", "5"]),
	(["bw=600", "bw_long=600", "max_gap=7000", "min_cnt=5", "min_chain_score=100"], ["-r", "600,600", "-g", "7000", "-n", "5", "-m", "100"]),
	(["max_chain_skip=0", "max_chain_iter=3", "chain_gap_scale=2", "chain_skip_scale=0.3"],
		["--max-chain-skip", "0", "--max-chain-iter", "3", "--chain-gap-scale", "2", "--chain-skip-scale", "0.3"]),
	(["+CIGAR", "a=3", "b=5", "q=5", "e=3", "q2=30", "e2=1", "transition=2", "sc_ambi=3"],
		["-c", "-A", "3", "-B", "5", "-O", "5,30", "-E", "3,1", "-b", "2", "--score-N", "3"]),
	(["+CIGAR", "zdrop=200", "zdrop_inv=50", "min_dp_max=300", "max_clip_ratio=0.5"], ["-c", "-z", "200,50", "-s", "300", "--max-clip-ratio", "0.5"]),
	(["mid_occ_frac=0.01", "min_mid_occ=2", "max_mid_occ=100", "occ_dist=0", "q_occ_frac=0"], ["-f", "0.01", "-U", "2,100", "-e", "0", "--q-occ-frac", "0"]),
	(["mid_occ=5", "max_qlen=3000"], ["-f", "5", "--max-qlen", "3000"]),
	(["+CIGAR", "max_sw_mat=1000"], ["-c", "--cap-sw-mem", "1000"]),
]


@pytest.mark.parametrize("settings,cli", PAIRS, ids=[" ".join(p[1]) or "default" for p in PAIRS])
def test_api_fields_and_flags_match_command_line(synth, data, settings, cli):
	"""mm_map with the mm_mapopt_t fields / MM_F_* flags set directly reports, for every hit, the
	coordinates, strand, mlen/blen, mapq, primary/secondary/inversion type (id, parent, inv), cnt,
	score, subsc, split, rep_len and, with MM_F_CIGAR, NM (blen - mlen + n_ambi), dp_max0, dp_score,
	n_ambi and the CIGAR that minimap2 prints in PAF for the corresponding options."""
	_compare(synth, data, "-", settings, cli)


@pytest.mark.parametrize("preset", ["map-pb", "map-hifi", "lr:hq", "asm20", "map-iclr"])
def test_api_presets_match_command_line(synth, data, preset):
	"""mm_set_opt(preset) followed by mm_map gives the command line's -x preset -c output."""
	_compare(synth, data, preset, ["+CIGAR"], ["-c"])


def test_api_result_fields_consistent(synth, data):
	"""mm_reg1_t fields without a PAF column: score0 is the chaining score before a split (equal to
	score unless the hit was split, larger when it was); dp_max0 >= dp_max2 for primaries; div is -1
	or in [0, 1]; at most one hit per read has sam_pri, and it is a primary."""
	hits = api_map(synth.fasta, data, "-", ["+CIGAR"])
	assert len(hits) > 20
	for h in hits:
		assert h.score0 == h.score if not h.split else h.score0 >= h.score, (h.qname, h.score, h.score0, h.split)
		if h.id == h.parent and not h.inv:
			assert h.dp_max0 >= h.dp_max2
		assert h.div == -1 or 0 <= h.div <= 1
		if h.sam_pri:
			assert h.id == h.parent
	assert any(h.split for h in hits) and any(h.inv for h in hits)


@pytest.mark.parametrize("preset,flags", [
	("map-ont", []), ("lr", []), ("map-pb", []), ("map10k", []), ("lr:hq", []), ("map-hifi", []), ("map-ccs", []),
	("map-iclr", []), ("map-iclr-prerender", []),
	("ava-ont", ["ALL_CHAINS", "NO_DIAG", "NO_DUAL", "NO_LJOIN"]), ("ava-pb", ["ALL_CHAINS", "NO_DIAG", "NO_DUAL", "NO_LJOIN"]),
	("asm5", ["RMQ"]), ("asm10", ["RMQ"]), ("asm20", ["RMQ"]), ("lr:hqae", ["RMQ"]),
	("sr", ["SR", "FRAG_MODE", "NO_PRINT_2ND", "HEAP_SORT"]), ("splice", ["SPLICE"]),
	("splice:sr", ["SPLICE", "FRAG_MODE", "NO_PRINT_2ND", "HEAP_SORT"])])
def test_preset_mode_flags(preset, flags):
	"""mm_set_opt sets the mode flags each preset implies (minimap2.1: ava-* include -X, i.e.
	MM_F_ALL_CHAINS, MM_F_NO_DIAG, MM_F_NO_DUAL and MM_F_NO_LJOIN; asm* and lr:hqae --rmq; sr and
	splice:sr --heap-sort=yes and --secondary=no; sr --sr --frag=yes; splice --splice) and no other
	among MM_F_NO_DIAG/NO_DUAL/CIGAR/NO_LJOIN/SR/FRAG_MODE/NO_PRINT_2ND/SPLICE/FOR_ONLY/REV_ONLY/
	HEAP_SORT/ALL_CHAINS/NO_END_FLT/HARD_MLEVEL/RMQ/QSTRAND/NO_INV/NO_HASH_NAME/EQX."""
	H.assert_c_case("drv_b_api", "opt_flags", preset, *flags)


def test_rank_min_len_and_rank_frac(cn_api):
	"""rank_min_len and rank_frac (no option): for a read of at least rank_min_len bases with two hits
	whose dp_max are within rank_frac of each other, dp_max is re-estimated (mm_update_dp_max) from the
	mismatches and gaps, so the hits of a 380 bp repeat read with substitutions get a dp_max different
	from dp_max0 with rank_min_len 100, keep dp_max0 with the default 500 (above the read length), and
	keep it with rank_frac above 1."""
	fa, rq = cn_api
	on = api_map(fa, rq, "-", ["+CIGAR", "rank_min_len=100"])
	off = api_map(fa, rq, "-", ["+CIGAR"])
	assert len(on) >= 2 and len(off) == len(on)
	assert any(h.dp_max != h.dp_max0 for h in on)
	assert all(h.dp_max == h.dp_max0 for h in off)
	frac_off = api_map(fa, rq, "-", ["+CIGAR", "rank_min_len=100", "rank_frac=1.01"])
	assert all(h.dp_max == h.dp_max0 for h in frac_off)


@pytest.fixture(scope="module")
def cn_api(tmp_path_factory):
	"""The copy-number genome and a 380 bp read from the first copy of its 3-copy (400 bp) family, with
	3% substitutions."""
	d = tmp_path_factory.mktemp("b_cn_api")
	g, placed = B.copy_number_genome()
	_, c, s, e, _ = [x for x in placed if x[0] == 0][0]
	read = H.mutate(H.Rng(5), g[c][s + 10:s + 390], 0.03)
	return g.write_fasta(str(d / "cn.fa")), H.write_fasta(str(d / "r.fa"), [("rep", read)])


def test_rmq_rescue_size_and_ratio(tmp_path):
	"""rmq_rescue_size and rmq_rescue_ratio (no option): the long-join pass re-chains with the long
	band when the best chain leaves more than rmq_rescue_size query bases uncovered or covers more
	than rmq_rescue_ratio of the query. A read with a 600 bp deletion (two chains of half the read
	each) is joined into one hit by default and stays two hits when neither trigger can fire
	(rmq_rescue_size 1000000, rmq_rescue_ratio 1)."""
	g = H.make_genome(seed=21, contig_lengths=(60000, 40000))
	fa = g.write_fasta(str(tmp_path / "g.fa"))
	rq = H.write_fasta(str(tmp_path / "del.fa"), [("del", g["chr1"][20000:22000] + g["chr1"][22600:24600])])
	assert len(api_map(fa, rq)) == 1
	assert len(api_map(fa, rq, "-", ["rmq_rescue_size=1000000", "rmq_rescue_ratio=1"])) == 2
	assert len(api_map(fa, rq, "-", ["rmq_rescue_size=1000000"])) == 1


def test_anchor_ext_len(tmp_path):
	"""anchor_ext_len (no option): the local score s of --end-seed-pen is computed over the terminal
	anchor plus anchor_ext_len bases on each side (minimap2.1, --end-seed-pen). For a 19 bp first exon
	holding one 15-mer anchor 1521 bp from the next anchor, s is 15 with anchor_ext_len 0 and 19 (the
	whole exon) with 20; with --end-seed-pen 8 the threshold ln(1521) + 8 = 15.3 lies between, so
	the exon is dropped at 0 and kept at 20."""
	rng = H.Rng(1)
	ex1, intron, ex2 = rng.dna(19), "GT" + rng.dna(1496) + "AG", rng.dna(600)
	ex1, ex2 = ex1[:-1] + "C", "T" + ex2[1:]
	ref = rng.dna(3000) + ex1 + intron + ex2 + rng.dna(3000)
	read = ex1 + ex2
	fa = H.write_fasta(str(tmp_path / "sp.fa"), [("ref", ref)])
	rq = H.write_fasta(str(tmp_path / "spr.fa"), [("tx", read)])
	idx = B.MinimizerIndex({"ref": ref}, 5, 15)
	anc = sorted([a for a in B.anchors(read, idx) if a.strand == "+"], key=lambda a: a.x)
	g = anc[1].x - anc[0].x
	assert 15 < math.log(g) + 8 < 19
	for ext, cigar in ((0, "600M"), (20, "19M1500N600M")):
		h = api_map(fa, rq, "splice", ["+CIGAR", "anchor_ext_len=%d" % ext, "anchor_ext_shift=8"])
		assert [x.cigar for x in h] == [cigar], (ext, [x.cigar for x in h])


def test_max_max_occ(tmp_path):
	"""max_max_occ (no option): minimizers occurring more often are dropped even when -e sampling
	would keep them. With -f 2 -e 100 on a read spanning a 12-copy repeat, the chain gains one anchor
	per sampled repetitive minimizer (cnt from brute-force sampling) with max_max_occ 4095, none with
	max_max_occ 11, and cnt matches the brute-force count in both cases."""
	g, placed = B.copy_number_genome()
	idx = B.MinimizerIndex(g, 10, 15)
	_, c, s, e, _ = [x for x in placed if x[0] == 2][0]
	q = g[c][s - 700:e + 700]
	fa = g.write_fasta(str(tmp_path / "cn.fa"))
	rq = H.write_fasta(str(tmp_path / "q.fa"), [("q", q)])
	base = api_map(fa, rq, "-", ["mid_occ=2", "occ_dist=0"])
	hi = api_map(fa, rq, "-", ["mid_occ=2", "occ_dist=100", "max_max_occ=4095"])
	lo = api_map(fa, rq, "-", ["mid_occ=2", "occ_dist=100", "max_max_occ=11"])
	k0 = B.collect_expectation(q, idx, 2, 4095, 0)[1]
	k_hi = B.collect_expectation(q, idx, 2, 4095, 100)[1]
	k_lo = B.collect_expectation(q, idx, 2, 11, 100)[1]
	prim = lambda hs: [h for h in hs if h.id == h.parent][0]
	assert k_hi > k0 and k_lo < k_hi
	assert prim(hi).cnt - prim(base).cnt == k_hi - k0
	assert prim(lo).cnt - prim(base).cnt == k_lo - k0
