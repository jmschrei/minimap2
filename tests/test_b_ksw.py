"""Area B: the ksw2 alignment kernels (ksw2_extz2_sse.c, ksw2_extd2_sse.c, ksw2_exts2_sse.c,
ksw2_ll_sse.c, ksw2_dispatch.c) through tests/c/drv_b_ksw.c, against a brute-force Gotoh DP written
in the driver from the kernels' documented contract (ksw2.h). Every random case is generated from a
fixed seed inside the driver."""

import pytest

import minimap2_helpers as H


@pytest.mark.parametrize("kind", [0, 1], ids=["extz2", "extd2"])
def test_global_alignment_matches_reference(kind):
	"""ksw_extz2_sse (affine, a k-long gap costs O+k*E) and ksw_extd2_sse (two-piece,
	min{O1+k*E1, O2+k*E2}) in global mode: on 600 random pairs (substitutions, indels, N bases,
	random scoring within ksw2's limits) the score equals the brute-force optimum; the CIGAR spans both
	sequences and scores exactly that; KSW_EZ_SCORE_ONLY, KSW_EZ_RIGHT, KSW_EZ_APPROX_MAX and
	KSW_EZ_REV_CIGAR give the same score, and REV_CIGAR is the reversed CIGAR."""
	H.assert_c_case("drv_b_ksw", "global", kind, 11 + kind, 600, 120)


@pytest.mark.parametrize("kind", [0, 1], ids=["extz2", "extd2"])
def test_extension_matches_reference(kind):
	"""KSW_EZ_EXTZ_ONLY (the left and right extensions of every minimap2 alignment) on 600 random
	pairs: max (best score over all cells), mqe (best cell ending the query) and mte (ending the
	target) equal the reference; reach_end is set exactly when mqe + end_bonus > max (--end-bonus);
	the CIGAR ends at the max cell and scores max, or reaches the query end and scores mqe."""
	H.assert_c_case("drv_b_ksw", "extension", kind, 21 + kind, 600, 120)


def test_generic_matrix_transitions():
	"""KSW_EZ_GENERIC_SC with minimap2's transition matrix (-b, set for map-iclr): scores follow the
	full matrix in both kernels."""
	H.assert_c_case("drv_b_ksw", "generic_matrix")


def test_gap_left_and_right_alignment():
	"""A one-base deletion or insertion in a 6-base homopolymer is placed at its leftmost position by
	default and at its rightmost with KSW_EZ_RIGHT (ksw2.h: "right-align gaps")."""
	H.assert_c_case("drv_b_ksw", "gap_alignment")


def test_zdrop_stops_extension():
	"""Z-drop (-z): an exact 120-base match followed by unrelated sequence stops with zdropped set at
	zdrop 50, keeping the match's score as max; with zdrop 100000 the extension equals the reference."""
	H.assert_c_case("drv_b_ksw", "zdrop")


def test_zdrop_threshold_is_strict():
	"""Z-drop (-z; the minimap2 paper: stop when the score falls more than Z + e*|diagonal change| below
	the best cell): on 600 extensions of an exact 20-80 base match into unrelated tails, with Z the
	largest such drop from a brute-force DP (over each anti-diagonal's best cell; e2 for ksw_extd2),
	zdrop Z runs to the end and Z-1 stops. Cases where tied cells would change the drop are skipped."""
	H.assert_c_case("drv_b_ksw", "zdrop_boundary", 71, 600)


def test_band_limits_global_alignment():
	"""The band (-r): a 300-base pair differing by a d-base deletion aligns as 120M dD (180-d)M when
	the band w >= d, and its end is unreachable (zdropped, score KSW_NEG_INF) when w < d."""
	H.assert_c_case("drv_b_ksw", "band_deletion")


def test_empty_sequences():
	"""A zero-length query or target aligns nothing and leaves ksw_extz_t reset."""
	H.assert_c_case("drv_b_ksw", "empty")


def test_local_alignment_ll_i16():
	"""ksw_ll_qinit + ksw_ll_i16 (inversion detection and splice-end filtering): on 1000 random pairs
	the score is the brute-force Smith-Waterman optimum and (qe, te) ends a best local alignment."""
	H.assert_c_case("drv_b_ksw", "ll_i16", 31, 1000, 150)


def test_spliced_kernel_places_intron():
	"""ksw_exts2_sse with the splice preset's scoring (-A1 -B2 -O2,32 -E1): exon1+exon2 against
	exon1 + a 400-base GT..AG intron + exon2 aligns as 60M400N60M with KSW_EZ_SPLICE_FOR, and on the
	reverse complements (CT..AC) with KSW_EZ_SPLICE_REV; a long gap is an intron (N) even without the
	splice-signal flags (minimap2.1, splice mode)."""
	H.assert_c_case("drv_b_ksw", "exts2_intron")


def test_spliced_kernel_long_gap_threshold():
	"""ksw_exts2_sse with the splice preset's scoring: a k-base gap is a deletion while 2 + k < 32 and an
	intron (N) once 2 + k > 32, the cheaper of O1+k*E1 and O2 (minimap2.1: in the splice mode the
	second gap penalty is a flat O2 and long deletions are introns); the score is 120 minus that cost."""
	H.assert_c_case("drv_b_ksw", "exts2_long_gap")


def test_simd_variants_agree():
	"""The SSE2 and SSE4.1 builds of ksw_extz2/ksw_extd2 behind the runtime dispatch return the same
	scores, end cells and CIGARs as the entry point on 300 random cases with bands, Z-drop, end
	bonus and every flag minimap2 uses (a no-op in builds without the dispatch)."""
	H.assert_c_case("drv_b_ksw", "variants", 41, 300, 150)


def test_simd_variants_agree_spliced():
	"""As test_simd_variants_agree for ksw_exts2 (both splice strands, flanking and complex splice
	models, extension and approximate max)."""
	H.assert_c_case("drv_b_ksw", "variants_splice", 51, 300, 150)


@pytest.mark.pinned
def test_narrow_band_edge(tmp_path):
	"""Pins the current output of ksw_extz2/ksw_extd2 on 600 random pairs whose band w is |qlen - tlen|
	plus 0 to 40, so the end cell lies on or near the band edge (global and extension, with end bonus,
	Z-drop, KSW_EZ_RIGHT and KSW_EZ_APPROX_MAX; some with a 5-40 base indel). There the kernels can report
	a score no alignment reaches and a CIGAR leaving the band (known-bugs/b-ksw2-banded-score.md), so the
	pin guards the exact results. Independently: every CIGAR is well formed; a global CIGAR not cut by
	Z-drop covers both sequences and scores at most the unbanded optimum; an extension CIGAR ends at the
	cell the kernel reports (the query end with reach_end, else the max cell)."""
	out = str(tmp_path / "band_edge.tsv")
	H.assert_c_case("drv_b_ksw", "band_edge", 61, 600, 150, out)
	H.check_pinned("b/ksw-band-edge.tsv.gz", H.read_text(out))
