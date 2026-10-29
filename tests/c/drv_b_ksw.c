/* C driver: the ksw2 alignment kernels (ksw2.h) against a brute-force
 * reference written here from their documented contract.
 *
 * ksw_extz2_sse    affine gaps, a k-long gap costs q+k*e
 * ksw_extd2_sse    two-piece affine gaps, min{q+k*e, q2+k*e2}
 * ksw_exts2_sse    spliced alignment (an intron is an N operation)
 * ksw_ll_qinit/ksw_ll_i16  local (Smith-Waterman) score and end positions
 *
 * The reference is a full Gotoh DP with two gap pieces over cells inside the
 * band |t - q| <= w (t: target index, q: query index), the band ksw2 uses. A
 * global alignment covers both sequences; an extension (KSW_EZ_EXTZ_ONLY)
 * starts at the first base of both and may stop anywhere: ksw_extz_t::max is
 * the best score over all cells (at least 0), mqe the best score of a cell
 * that ends the query, mte one that ends the target. The kernels may break
 * ties between equally good alignments as they like, so CIGARs are checked by
 * their span and their score, not compared with a reference CIGAR.
 *
 * Calls the kernels with hand-built ksw_extz_t structs (ksw2.h); the
 * per-variant cases look up ksw_extz2_sse2/ksw_extz2_sse41 etc. as weak
 * symbols, which exist only in the default x86 build with runtime dispatch. */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include "ksw2.h"
#include "minimap2_test.h"

#define REF_NEG (-(1<<28))

typedef void (*extz2_f)(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez);
typedef void (*extd2_f)(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez);
typedef void (*exts2_f)(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t noncan, int zdrop, int end_bonus, int8_t junc_bonus, int8_t junc_pen, int flag, const uint8_t *junc, ksw_extz_t *ez);

/************************
 * Deterministic random *
 ************************/

static uint64_t rng_state;

static uint64_t rng_next(void) /* splitmix64 */
{
	uint64_t z = (rng_state += 0x9E3779B97F4A7C15ULL);
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
	z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
	return z ^ (z >> 31);
}

static int rng_int(int lo, int hi) { return lo + (int)(rng_next() % (uint64_t)(hi - lo + 1)); }

/* a random target and a query derived from it with substitutions, indels and (optionally) N */
static void rand_pair(int max_len, int n_rate1000, uint8_t *t, int *tlen, uint8_t *q, int *qlen)
{
	int i, j, l = rng_int(1, max_len), sub = rng_int(0, 300), indel = rng_int(0, 120);
	for (i = 0; i < l; ++i) t[i] = rng_int(0, 3);
	if (n_rate1000 > 0)
		for (i = 0; i < l; ++i)
			if (rng_int(0, 999) < n_rate1000) t[i] = 4;
	for (i = j = 0; i < l && j < max_len; ++i) {
		int x = rng_int(0, 999);
		if (x < indel / 2) { /* deletion of 1-6 bases */
			i += rng_int(0, 5);
		} else if (x < indel) { /* insertion of 1-6 bases */
			int k, n = rng_int(1, 6);
			for (k = 0; k < n && j < max_len; ++k) q[j++] = rng_int(0, 3);
			--i;
		} else if (x < indel + sub) {
			q[j++] = (t[i] + rng_int(1, 3)) & 3;
		} else q[j++] = t[i];
	}
	if (j == 0) q[j++] = rng_int(0, 3);
	if (n_rate1000 > 0)
		for (i = 0; i < j; ++i)
			if (rng_int(0, 999) < n_rate1000) q[i] = 4;
	*tlen = l, *qlen = j;
}

/* minimap2's scoring matrix: match a, mismatch -b, transitions -ts (if non-zero), anything with N -sc_ambi */
static void gen_mat(int8_t *mat, int a, int b, int ts, int sc_ambi)
{
	int i, j;
	for (i = 0; i < 5; ++i)
		for (j = 0; j < 5; ++j)
			mat[i * 5 + j] = i == 4 || j == 4? -sc_ambi : i == j? a : -b;
	if (ts) mat[0 * 5 + 2] = mat[2 * 5 + 0] = mat[1 * 5 + 3] = mat[3 * 5 + 1] = -ts;
}

/*************
 * Reference *
 *************/

typedef struct {
	int qlen, tlen;
	int32_t *H; /* (qlen+1) x (tlen+1); H[i][j]: best score of query[0..i) against target[0..j) */
} ref_t;

static inline int32_t max2(int32_t a, int32_t b) { return a > b? a : b; }

/* gap cost of a k-long gap: min{q+k*e, q2+k*e2} */
static inline int32_t gap_cost(int k, int q, int e, int q2, int e2)
{
	int32_t c1 = q + k * e, c2 = q2 + k * e2;
	return c1 < c2? c1 : c2;
}

static void ref_dp(ref_t *R, int qlen, const uint8_t *qs, int tlen, const uint8_t *ts, const int8_t *mat, int q, int e, int q2, int e2, int w)
{
	int i, j, W = tlen + 1;
	int32_t *E1, *E2, *F1, *F2;
	R->qlen = qlen, R->tlen = tlen;
	R->H = (int32_t*)malloc((size_t)(qlen + 1) * W * 4);
	E1 = (int32_t*)malloc((size_t)(qlen + 1) * W * 4);
	E2 = (int32_t*)malloc((size_t)(qlen + 1) * W * 4);
	F1 = (int32_t*)malloc((size_t)(qlen + 1) * W * 4);
	F2 = (int32_t*)malloc((size_t)(qlen + 1) * W * 4);
	if (w < 0) w = qlen > tlen? qlen : tlen;
	for (i = 0; i <= qlen; ++i) {
		for (j = 0; j <= tlen; ++j) {
			size_t k = (size_t)i * W + j;
			int32_t h = REF_NEG;
			E1[k] = E2[k] = F1[k] = F2[k] = REF_NEG;
			if (i == 0 && j == 0) { R->H[k] = 0; continue; }
			if (i > 0 && j > 0) {
				int dt = (j - 1) - (i - 1);
				if (dt > w || dt < -w) { R->H[k] = REF_NEG; continue; } /* outside the band */
			}
			if (j > 0) { /* deletion: target base without a query base */
				E1[k] = max2(R->H[k-1] - q - e, E1[k-1] - e);
				E2[k] = max2(R->H[k-1] - q2 - e2, E2[k-1] - e2);
				h = max2(h, max2(E1[k], E2[k]));
			}
			if (i > 0) {
				F1[k] = max2(R->H[k-W] - q - e, F1[k-W] - e);
				F2[k] = max2(R->H[k-W] - q2 - e2, F2[k-W] - e2);
				h = max2(h, max2(F1[k], F2[k]));
			}
			if (i > 0 && j > 0 && R->H[k-W-1] > REF_NEG / 2)
				h = max2(h, R->H[k-W-1] + mat[ts[j-1] * 5 + qs[i-1]]);
			R->H[k] = h < REF_NEG / 2? REF_NEG : h;
		}
	}
	free(E1); free(E2); free(F1); free(F2);
}

static inline int32_t ref_at(const ref_t *R, int i, int j) { return R->H[(size_t)i * (R->tlen + 1) + j]; }

/* extension: max over all cells (the empty alignment scores 0), best cell ending the query/target */
static void ref_ext(const ref_t *R, int32_t *mx, int32_t *mqe, int32_t *mte)
{
	int i, j;
	*mx = 0, *mqe = *mte = REF_NEG;
	for (i = 1; i <= R->qlen; ++i)
		for (j = 1; j <= R->tlen; ++j) {
			int32_t h = ref_at(R, i, j);
			if (h > *mx) *mx = h;
			if (i == R->qlen && h > *mqe) *mqe = h;
			if (j == R->tlen && h > *mte) *mte = h;
		}
}

static void ref_free(ref_t *R) { free(R->H); R->H = 0; }

/* score and spans of a CIGAR under the same model; returns 0 if it is malformed */
static int cigar_eval(int n_cigar, const uint32_t *cigar, const uint8_t *qs, const uint8_t *ts, const int8_t *mat, int q, int e, int q2, int e2,
					  int *qspan, int *tspan, int32_t *score)
{
	int k, i = 0, j = 0;
	int32_t s = 0;
	for (k = 0; k < n_cigar; ++k) {
		int op = cigar[k] & 0xf, len = cigar[k] >> 4, l;
		if (len <= 0) return 0;
		if (k > 0 && (cigar[k-1] & 0xf) == (uint32_t)op) return 0; /* adjacent identical operators */
		if (op == KSW_CIGAR_MATCH) {
			for (l = 0; l < len; ++l) s += mat[ts[j + l] * 5 + qs[i + l]];
			i += len, j += len;
		} else if (op == KSW_CIGAR_INS) {
			s -= gap_cost(len, q, e, q2, e2), i += len;
		} else if (op == KSW_CIGAR_DEL) {
			s -= gap_cost(len, q, e, q2, e2), j += len;
		} else return 0;
	}
	*qspan = i, *tspan = j, *score = s;
	return 1;
}

/***********
 * Helpers *
 ***********/

static void run_kernel(int kind, int variant, int qlen, const uint8_t *qs, int tlen, const uint8_t *ts, const int8_t *mat,
					   int q, int e, int q2, int e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez);

#ifdef __ELF__
extern void ksw_extz2_sse2(void*, int, const uint8_t*, int, const uint8_t*, int8_t, const int8_t*, int8_t, int8_t, int, int, int, int, ksw_extz_t*) __attribute__((weak));
extern void ksw_extz2_sse41(void*, int, const uint8_t*, int, const uint8_t*, int8_t, const int8_t*, int8_t, int8_t, int, int, int, int, ksw_extz_t*) __attribute__((weak));
extern void ksw_extd2_sse2(void*, int, const uint8_t*, int, const uint8_t*, int8_t, const int8_t*, int8_t, int8_t, int8_t, int8_t, int, int, int, int, ksw_extz_t*) __attribute__((weak));
extern void ksw_extd2_sse41(void*, int, const uint8_t*, int, const uint8_t*, int8_t, const int8_t*, int8_t, int8_t, int8_t, int8_t, int, int, int, int, ksw_extz_t*) __attribute__((weak));
extern void ksw_exts2_sse2(void*, int, const uint8_t*, int, const uint8_t*, int8_t, const int8_t*, int8_t, int8_t, int8_t, int8_t, int, int, int8_t, int8_t, int, const uint8_t*, ksw_extz_t*) __attribute__((weak));
extern void ksw_exts2_sse41(void*, int, const uint8_t*, int, const uint8_t*, int8_t, const int8_t*, int8_t, int8_t, int8_t, int8_t, int, int, int8_t, int8_t, int, const uint8_t*, ksw_extz_t*) __attribute__((weak));
#endif

static extz2_f extz2_variant(int v)
{
#ifdef __ELF__
	if (v == 1) return ksw_extz2_sse2;
	if (v == 2) return ksw_extz2_sse41;
#endif
	return v == 0? ksw_extz2_sse : 0;
}

static extd2_f extd2_variant(int v)
{
#ifdef __ELF__
	if (v == 1) return ksw_extd2_sse2;
	if (v == 2) return ksw_extd2_sse41;
#endif
	return v == 0? ksw_extd2_sse : 0;
}

static exts2_f exts2_variant(int v)
{
#ifdef __ELF__
	if (v == 1) return ksw_exts2_sse2;
	if (v == 2) return ksw_exts2_sse41;
#endif
	return v == 0? ksw_exts2_sse : 0;
}

/* kind 0: ksw_extz2 (q2, e2 ignored); kind 1: ksw_extd2. variant 0: the dispatching entry point */
static void run_kernel(int kind, int variant, int qlen, const uint8_t *qs, int tlen, const uint8_t *ts, const int8_t *mat,
					   int q, int e, int q2, int e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
{
	if (kind == 0) extz2_variant(variant)(0, qlen, qs, tlen, ts, 5, mat, q, e, w, zdrop, end_bonus, flag, ez);
	else extd2_variant(variant)(0, qlen, qs, tlen, ts, 5, mat, q, e, q2, e2, w, zdrop, end_bonus, flag, ez);
}

static void enc(const char *s, uint8_t *x)
{
	for (; *s; ++s, ++x) *x = *s == 'A'? 0 : *s == 'C'? 1 : *s == 'G'? 2 : *s == 'T'? 3 : 4;
}

static void ez_free(ksw_extz_t *ez) { free(ez->cigar); memset(ez, 0, sizeof(*ez)); }

static void rand_params(int kind, int *a, int *b, int *q, int *e, int *q2, int *e2)
{
	*a = rng_int(1, 3);
	*q = rng_int(1, 12), *e = rng_int(1, 4);
	if (kind == 0) {
		*q2 = *q, *e2 = *e;
	} else { /* minimap2's constraint (mm_check_opt): e > e2 and q + e < q2 + e2 */
		*e2 = rng_int(0, *e - 1);
		*q2 = rng_int(*q + *e - *e2 + 1, *q + *e - *e2 + 30);
	}
	*b = rng_int(1, 2 * (*q + *e)); /* ksw2 requires the mismatch penalty to be at most 2(q+e) */
}

/*********
 * Cases *
 *********/

/* args: kind(0 extz2, 1 extd2) seed n max_len. Global alignment without a band (w < 0, or a w no
 * cell exceeds; see band_deletion and band_edge for narrow bands): the score
 * equals the reference;
 * the CIGAR covers both sequences and scores exactly that; SCORE_ONLY, KSW_EZ_RIGHT, KSW_EZ_APPROX_MAX
 * (minimap2's gap filling) and KSW_EZ_REV_CIGAR give the same score; REV_CIGAR reverses the CIGAR. */
static int case_global(int argc, char *argv[])
{
	int kind, n, max_len, it, n_banded = 0;
	uint8_t *t, *qq;
	if (argc < 4) return MMT_FAIL("usage: global <kind> <seed> <n> <max_len>"), 1;
	kind = atoi(argv[0]), rng_state = strtoull(argv[1], 0, 10), n = atoi(argv[2]), max_len = atoi(argv[3]);
	t = (uint8_t*)malloc(max_len + 1), qq = (uint8_t*)malloc(max_len + 1);
	for (it = 0; it < n; ++it) {
		int a, b, q, e, q2, e2, tlen, qlen, w, qspan, tspan, sc_ambi = rng_int(1, 3), ts = 0;
		int32_t want, cs;
		int8_t mat[25];
		ref_t R;
		ksw_extz_t ez, ez2, ez3, ez4;
		rand_params(kind, &a, &b, &q, &e, &q2, &e2);
		if (sc_ambi >= b) sc_ambi = b > 1? b - 1 : 1;
		rand_pair(max_len, it % 4 == 3? 30 : 0, t, &tlen, qq, &qlen);
		w = it % 2 == 0? -1 : (qlen > tlen? qlen : tlen) + rng_int(0, 10); /* no cell is outside the band */
		gen_mat(mat, a, b, ts, sc_ambi);
		memset(&ez, 0, sizeof(ez)); memset(&ez2, 0, sizeof(ez2)); memset(&ez3, 0, sizeof(ez3)); memset(&ez4, 0, sizeof(ez4));
		ref_dp(&R, qlen, qq, tlen, t, mat, q, e, q2, e2, w);
		want = ref_at(&R, qlen, tlen);
		run_kernel(kind, 0, qlen, qq, tlen, t, mat, q, e, q2, e2, w, -1, 0, 0, &ez);
		if (w >= 0) ++n_banded;
		MMT_CHECKF(ez.score == want && !ez.zdropped, "case %d (kind %d, qlen %d, tlen %d, w %d, a%d b%d q%d e%d q2%d e2%d): score %d, reference %d",
				   it, kind, qlen, tlen, w, a, b, q, e, q2, e2, ez.score, want);
		MMT_CHECKF(cigar_eval(ez.n_cigar, ez.cigar, qq, t, mat, q, e, q2, e2, &qspan, &tspan, &cs) && qspan == qlen && tspan == tlen && cs == want,
				   "case %d: CIGAR spans %d/%d (want %d/%d), scores %d (want %d)", it, qspan, tspan, qlen, tlen, cs, want);
		run_kernel(kind, 0, qlen, qq, tlen, t, mat, q, e, q2, e2, w, -1, 0, KSW_EZ_SCORE_ONLY, &ez2);
		MMT_CHECKF(ez2.score == want && ez2.n_cigar == 0, "case %d: SCORE_ONLY score %d, n_cigar %d", it, ez2.score, ez2.n_cigar);
		ez_free(&ez2);
		run_kernel(kind, 0, qlen, qq, tlen, t, mat, q, e, q2, e2, w, -1, 0, KSW_EZ_RIGHT, &ez3);
		MMT_CHECKF(ez3.score == want && cigar_eval(ez3.n_cigar, ez3.cigar, qq, t, mat, q, e, q2, e2, &qspan, &tspan, &cs) && qspan == qlen && tspan == tlen && cs == want,
				   "case %d: KSW_EZ_RIGHT score %d / CIGAR score %d, want %d", it, ez3.score, cs, want);
		run_kernel(kind, 0, qlen, qq, tlen, t, mat, q, e, q2, e2, w, -1, 0, KSW_EZ_APPROX_MAX, &ez2);
		MMT_CHECKF(ez2.score == want && cigar_eval(ez2.n_cigar, ez2.cigar, qq, t, mat, q, e, q2, e2, &qspan, &tspan, &cs) && qspan == qlen && tspan == tlen && cs == want,
				   "case %d: KSW_EZ_APPROX_MAX score %d / CIGAR score %d, want %d", it, ez2.score, cs, want);
		ez_free(&ez2);
		run_kernel(kind, 0, qlen, qq, tlen, t, mat, q, e, q2, e2, w, -1, 0, KSW_EZ_REV_CIGAR, &ez4);
		MMT_CHECK_INT(ez4.score, ==, want);
		if (MMT_CHECKF(ez4.n_cigar == ez.n_cigar, "case %d: REV_CIGAR has %d operators, plain %d", it, ez4.n_cigar, ez.n_cigar)) {
			int k, same = 1;
			for (k = 0; k < ez.n_cigar; ++k)
				if (ez4.cigar[k] != ez.cigar[ez.n_cigar - 1 - k]) same = 0;
			MMT_CHECKF(same, "case %d: REV_CIGAR is not the reversed CIGAR", it);
		}
		ref_free(&R); ez_free(&ez); ez_free(&ez3); ez_free(&ez4);
	}
	MMT_CHECKF(n_banded >= n / 4, "too few cases with an explicit band (%d)", n_banded);
	free(t); free(qq);
	return 0;
}

/* args: kind seed n max_len. Extension (KSW_EZ_EXTZ_ONLY, no Z-drop): max, mqe and mte equal the
 * reference; the cells they name score what they claim; without an end bonus large enough the CIGAR
 * ends at (max_q, max_t) and scores max, with one it reaches the query end (reach_end) and scores
 * mqe; reach_end is set exactly when mqe + end_bonus > max. */
static int case_extension(int argc, char *argv[])
{
	int kind, n, max_len, it, n_reach = 0, n_noreach = 0;
	uint8_t *t, *qq;
	if (argc < 4) return MMT_FAIL("usage: extension <kind> <seed> <n> <max_len>"), 1;
	kind = atoi(argv[0]), rng_state = strtoull(argv[1], 0, 10), n = atoi(argv[2]), max_len = atoi(argv[3]);
	t = (uint8_t*)malloc(max_len + 1), qq = (uint8_t*)malloc(max_len + 1);
	for (it = 0; it < n; ++it) {
		int a, b, q, e, q2, e2, tlen, qlen, w, qspan, tspan, end_bonus, flag;
		int32_t mx, mqe, mte, cs;
		int8_t mat[25];
		ref_t R;
		ksw_extz_t ez;
		rand_params(kind, &a, &b, &q, &e, &q2, &e2);
		rand_pair(max_len, 0, t, &tlen, qq, &qlen);
		if (it % 2 == 0) { /* append a random tail to the query so that extension tends to stop early */
			int k, extra = rng_int(1, 30);
			for (k = 0; k < extra && qlen < max_len; ++k) qq[qlen++] = rng_int(0, 3);
		}
		w = it % 2 == 0? -1 : (qlen > tlen? qlen : tlen) + rng_int(0, 10);
		end_bonus = it % 4 == 0? 0 : rng_int(0, 60);
		flag = KSW_EZ_EXTZ_ONLY | (it % 5 == 1? KSW_EZ_RIGHT : 0) | (it % 7 == 2? KSW_EZ_REV_CIGAR : 0);
		gen_mat(mat, a, b, 0, 1);
		ref_dp(&R, qlen, qq, tlen, t, mat, q, e, q2, e2, w);
		ref_ext(&R, &mx, &mqe, &mte);
		memset(&ez, 0, sizeof(ez));
		run_kernel(kind, 0, qlen, qq, tlen, t, mat, q, e, q2, e2, w, -1, end_bonus, flag, &ez);
		MMT_CHECKF((int32_t)ez.max == mx, "case %d (kind %d, qlen %d, tlen %d, w %d): max %d, reference %d", it, kind, qlen, tlen, w, ez.max, mx);
		MMT_CHECKF(ez.mqe == (mqe <= REF_NEG / 2? KSW_NEG_INF : mqe), "case %d: mqe %d, reference %d", it, ez.mqe, mqe);
		MMT_CHECKF(ez.mte == (mte <= REF_NEG / 2? KSW_NEG_INF : mte), "case %d: mte %d, reference %d", it, ez.mte, mte);
		if (mx > 0) MMT_CHECKF(ez.max_q >= 0 && ez.max_t >= 0 && ref_at(&R, ez.max_q + 1, ez.max_t + 1) == mx, "case %d: max cell (%d,%d)", it, ez.max_q, ez.max_t);
		if (ez.mqe > KSW_NEG_INF) MMT_CHECKF(ez.mqe_t >= 0 && ref_at(&R, qlen, ez.mqe_t + 1) == mqe, "case %d: mqe cell t=%d", it, ez.mqe_t);
		if (ez.mte > KSW_NEG_INF) MMT_CHECKF(ez.mte_q >= 0 && ref_at(&R, ez.mte_q + 1, tlen) == mte, "case %d: mte cell q=%d", it, ez.mte_q);
		MMT_CHECKF(ez.reach_end == (ez.mqe > KSW_NEG_INF && ez.mqe + end_bonus > (int)ez.max), "case %d: reach_end %d with mqe %d, bonus %d, max %d",
				   it, ez.reach_end, ez.mqe, end_bonus, ez.max);
		if (ez.n_cigar > 0) {
			uint32_t *cg = ez.cigar, tmp;
			int k;
			if (flag & KSW_EZ_REV_CIGAR)
				for (k = 0; k < ez.n_cigar>>1; ++k)
					tmp = cg[k], cg[k] = cg[ez.n_cigar-1-k], cg[ez.n_cigar-1-k] = tmp;
			if (!MMT_CHECKF(cigar_eval(ez.n_cigar, ez.cigar, qq, t, mat, q, e, q2, e2, &qspan, &tspan, &cs), "case %d: malformed CIGAR", it)) {
			} else if (ez.reach_end) {
				++n_reach;
				MMT_CHECKF(qspan == qlen && tspan == ez.mqe_t + 1 && cs == ez.mqe, "case %d: reach_end CIGAR spans %d/%d scores %d; want %d/%d %d",
						   it, qspan, tspan, cs, qlen, ez.mqe_t + 1, ez.mqe);
			} else {
				++n_noreach;
				MMT_CHECKF(qspan == ez.max_q + 1 && tspan == ez.max_t + 1 && cs == (int32_t)ez.max, "case %d: CIGAR spans %d/%d scores %d; want %d/%d %d",
						   it, qspan, tspan, cs, ez.max_q + 1, ez.max_t + 1, ez.max);
			}
		} else MMT_CHECKF(mx == 0 && !ez.reach_end, "case %d: no CIGAR but max %d", it, mx);
		ref_free(&R); ez_free(&ez);
	}
	MMT_CHECKF(n_reach >= n / 10 && n_noreach >= n / 10, "too few cases reaching (%d) or not reaching (%d) the query end", n_reach, n_noreach);
	free(t); free(qq);
	return 0;
}

/* args: none. The transition matrix (minimap2's -b) needs KSW_EZ_GENERIC_SC: with it, scores follow
 * the full matrix; global scores equal the reference for both kernels. */
static int case_generic_matrix(int argc, char *argv[])
{
	int it, kind;
	uint8_t t[200], qq[200];
	rng_state = 77;
	for (kind = 0; kind < 2; ++kind) {
		for (it = 0; it < 150; ++it) {
			int a, b, q, e, q2, e2, tlen, qlen, qspan, tspan, ts;
			int32_t want, cs;
			int8_t mat[25];
			ref_t R;
			ksw_extz_t ez;
			rand_params(kind, &a, &b, &q, &e, &q2, &e2);
			ts = rng_int(1, b);
			rand_pair(150, 10, t, &tlen, qq, &qlen);
			gen_mat(mat, a, b, ts, 1);
			ref_dp(&R, qlen, qq, tlen, t, mat, q, e, q2, e2, -1);
			want = ref_at(&R, qlen, tlen);
			memset(&ez, 0, sizeof(ez));
			run_kernel(kind, 0, qlen, qq, tlen, t, mat, q, e, q2, e2, -1, -1, 0, KSW_EZ_GENERIC_SC, &ez);
			MMT_CHECKF(ez.score == want, "kind %d case %d: score %d, reference %d (b %d, transition %d)", kind, it, ez.score, want, b, ts);
			MMT_CHECKF(cigar_eval(ez.n_cigar, ez.cigar, qq, t, mat, q, e, q2, e2, &qspan, &tspan, &cs) && qspan == qlen && tspan == tlen && cs == want,
					   "kind %d case %d: CIGAR score %d", kind, it, cs);
			ref_free(&R); ez_free(&ez);
		}
	}
	return 0;
}

/* args: none. A deletion in a homopolymer run is placed at the leftmost position by default and at
 * the rightmost with KSW_EZ_RIGHT (ksw2.h: "right-align gaps"); the same for an insertion. */
static int case_gap_alignment(int argc, char *argv[])
{
	/* target GCTAG AAAAAA CGTCA ; query GCTAG AAAAA CGTCA (one A deleted) */
	static const char *ts = "GCTAGAAAAAACGTCAGT", *qs_del = "GCTAGAAAAACGTCAGT", *qs_ins = "GCTAGAAAAAAACGTCAGT";
	uint8_t t[32], qd[32], qi[32];
	int8_t mat[25];
	int kind, tl = strlen(ts), qdl = strlen(qs_del), qil = strlen(qs_ins);
	enc(ts, t), enc(qs_del, qd), enc(qs_ins, qi);
	gen_mat(mat, 2, 4, 0, 1);
	for (kind = 0; kind < 2; ++kind) {
		ksw_extz_t ez;
		memset(&ez, 0, sizeof(ez));
		run_kernel(kind, 0, qdl, qd, tl, t, mat, 4, 2, 24, 1, -1, -1, 0, 0, &ez);
		MMT_CHECKF(ez.n_cigar == 3 && ez.cigar[0] == (5<<4|KSW_CIGAR_MATCH) && ez.cigar[1] == (1<<4|KSW_CIGAR_DEL) && ez.cigar[2] == (12<<4|KSW_CIGAR_MATCH),
				   "kind %d: left-aligned deletion should be 5M1D12M (n_cigar %d)", kind, ez.n_cigar);
		ez_free(&ez);
		run_kernel(kind, 0, qdl, qd, tl, t, mat, 4, 2, 24, 1, -1, -1, 0, KSW_EZ_RIGHT, &ez);
		MMT_CHECKF(ez.n_cigar == 3 && ez.cigar[0] == (10<<4|KSW_CIGAR_MATCH) && ez.cigar[1] == (1<<4|KSW_CIGAR_DEL) && ez.cigar[2] == (7<<4|KSW_CIGAR_MATCH),
				   "kind %d: right-aligned deletion should be 10M1D7M (n_cigar %d)", kind, ez.n_cigar);
		ez_free(&ez);
		run_kernel(kind, 0, qil, qi, tl, t, mat, 4, 2, 24, 1, -1, -1, 0, 0, &ez);
		MMT_CHECKF(ez.n_cigar == 3 && ez.cigar[0] == (5<<4|KSW_CIGAR_MATCH) && ez.cigar[1] == (1<<4|KSW_CIGAR_INS) && ez.cigar[2] == (13<<4|KSW_CIGAR_MATCH),
				   "kind %d: left-aligned insertion should be 5M1I13M (n_cigar %d)", kind, ez.n_cigar);
		ez_free(&ez);
		run_kernel(kind, 0, qil, qi, tl, t, mat, 4, 2, 24, 1, -1, -1, 0, KSW_EZ_RIGHT, &ez);
		MMT_CHECKF(ez.n_cigar == 3 && ez.cigar[0] == (11<<4|KSW_CIGAR_MATCH) && ez.cigar[1] == (1<<4|KSW_CIGAR_INS) && ez.cigar[2] == (7<<4|KSW_CIGAR_MATCH),
				   "kind %d: right-aligned insertion should be 11M1I7M (n_cigar %d)", kind, ez.n_cigar);
		ez_free(&ez);
	}
	return 0;
}

/* args: none. Z-drop: a 120-base exact match followed by 200 unrelated bases. Without Z-drop (or with
 * a large one) the extension runs to the end of the matrix; with zdrop 50 it stops (zdropped) with max
 * equal to the score of the exact match (240 at -A2) at its last cell; the CIGAR is 120M. */
static int case_zdrop(int argc, char *argv[])
{
	uint8_t t[400], qq[400];
	int8_t mat[25];
	int i, kind;
	rng_state = 99;
	for (i = 0; i < 120; ++i) t[i] = qq[i] = rng_int(0, 3);
	for (i = 120; i < 320; ++i) t[i] = rng_int(0, 3), qq[i] = rng_int(0, 3);
	gen_mat(mat, 2, 4, 0, 1);
	for (kind = 0; kind < 2; ++kind) {
		ksw_extz_t ez;
		int32_t mx, mqe, mte;
		ref_t R;
		ref_dp(&R, 320, qq, 320, t, mat, 4, 2, 24, 1, -1);
		ref_ext(&R, &mx, &mqe, &mte);
		ref_free(&R);
		memset(&ez, 0, sizeof(ez));
		run_kernel(kind, 0, 320, qq, 320, t, mat, 4, 2, 24, 1, -1, 50, 0, KSW_EZ_EXTZ_ONLY, &ez);
		MMT_CHECKF(ez.zdropped == 1, "kind %d: zdrop 50 did not stop the extension", kind);
		MMT_CHECKF(ez.max >= 240 && ez.max_q >= 119 && ez.max_t >= 119 && ez.max_q < 140 && ez.max_t < 140, "kind %d: max %d at (%d,%d)", kind, ez.max, ez.max_q, ez.max_t);
		ez_free(&ez);
		run_kernel(kind, 0, 320, qq, 320, t, mat, 4, 2, 24, 1, -1, 100000, 0, KSW_EZ_EXTZ_ONLY, &ez);
		MMT_CHECKF(ez.zdropped == 0 && (int32_t)ez.max == mx && ez.mqe == mqe && ez.mte == mte, "kind %d: zdrop 100000: zdropped %d max %d (ref %d)", kind, ez.zdropped, ez.max, mx);
		ez_free(&ez);
	}
	return 0;
}

/* args: none. The band: a query equal to a 300-base target except for a d-base deletion after base
 * 120 (d = 10, 35, 60). With w >= d the global alignment is 120M dD (180-d)M with score
 * 2(300-d) - gap(d); with w < d its end lies outside the band |t - q| <= w, so it is not reached
 * (zdropped, score KSW_NEG_INF). The kernels fill whole 16-cell blocks, so cells just outside the
 * band may be scored; nothing is asserted about them. */
static int case_band_deletion(int argc, char *argv[])
{
	uint8_t t[300], qq[300];
	int8_t mat[25];
	int i, kind, w, d;
	rng_state = 5;
	for (i = 0; i < 300; ++i) t[i] = rng_int(0, 3);
	gen_mat(mat, 2, 4, 0, 1);
	for (d = 10; d <= 60; d += 25) {
		int ql = 0;
		for (i = 0; i < 300; ++i)
			if (i < 120 || i >= 120 + d) qq[ql++] = t[i];
		for (kind = 0; kind < 2; ++kind) {
			int q2 = kind? 24 : 4, e2 = kind? 1 : 2;
			for (w = d - 4; w <= d + 4; ++w) {
				ksw_extz_t ez;
				memset(&ez, 0, sizeof(ez));
				run_kernel(kind, 0, ql, qq, 300, t, mat, 4, 2, q2, e2, w, -1, 0, 0, &ez);
				if (w >= d) {
					MMT_CHECKF(!ez.zdropped && ez.score == 2 * (300 - d) - gap_cost(d, 4, 2, q2, e2), "kind %d d %d w %d: score %d", kind, d, w, ez.score);
					MMT_CHECKF(ez.n_cigar == 3 && ez.cigar[0] == (120<<4|KSW_CIGAR_MATCH) && ez.cigar[1] == ((uint32_t)d<<4|KSW_CIGAR_DEL)
							   && ez.cigar[2] == ((uint32_t)(180 - d)<<4|KSW_CIGAR_MATCH), "kind %d d %d w %d: CIGAR is not 120M%dD%dM", kind, d, w, d, 180 - d);
				} else {
					MMT_CHECKF(ez.zdropped && ez.score == KSW_NEG_INF, "kind %d d %d w %d: end outside the band but zdropped %d score %d", kind, d, w, ez.zdropped, ez.score);
				}
				ez_free(&ez);
			}
		}
	}
	return 0;
}

/* args: none. qlen or tlen 0: nothing is aligned and ksw_extz_t is reset */
static int case_empty(int argc, char *argv[])
{
	uint8_t s[4] = {0, 1, 2, 3};
	int8_t mat[25];
	int kind;
	gen_mat(mat, 2, 4, 0, 1);
	for (kind = 0; kind < 2; ++kind) {
		ksw_extz_t ez;
		memset(&ez, 0, sizeof(ez));
		ez.n_cigar = 7, ez.max = 9;
		run_kernel(kind, 0, 0, s, 4, s, mat, 4, 2, 24, 1, -1, -1, 0, 0, &ez);
		MMT_CHECKF(ez.n_cigar == 0 && ez.max == 0 && ez.score == KSW_NEG_INF && ez.max_q == -1 && ez.max_t == -1, "kind %d qlen 0", kind);
		run_kernel(kind, 0, 4, s, 0, s, mat, 4, 2, 24, 1, -1, -1, 0, KSW_EZ_EXTZ_ONLY, &ez);
		MMT_CHECKF(ez.n_cigar == 0 && ez.max == 0 && ez.mqe == KSW_NEG_INF, "kind %d tlen 0", kind);
		ez_free(&ez);
	}
	return 0;
}

/* local alignment score with gap cost gapo + k*gape: brute force */
static int32_t ref_local(int qlen, const uint8_t *qs, int tlen, const uint8_t *ts, const int8_t *mat, int gapo, int gape, int32_t *Hout)
{
	int i, j, W = tlen + 1;
	int32_t best = 0, *E = (int32_t*)calloc((size_t)(qlen + 1) * W, 4), *F = (int32_t*)calloc((size_t)(qlen + 1) * W, 4);
	for (i = 0; i <= qlen; ++i)
		for (j = 0; j <= tlen; ++j) {
			size_t k = (size_t)i * W + j;
			int32_t h = 0;
			E[k] = F[k] = REF_NEG;
			if (i > 0 && j > 0) {
				E[k] = max2(Hout[k-1] - gapo - gape, E[k-1] - gape);
				F[k] = max2(Hout[k-W] - gapo - gape, F[k-W] - gape);
				h = max2(0, Hout[k-W-1] + mat[ts[j-1] * 5 + qs[i-1]]);
				h = max2(h, max2(E[k], F[k]));
			}
			Hout[k] = h;
			if (h > best) best = h;
		}
	free(E); free(F);
	return best;
}

/* args: seed n max_len. ksw_ll_qinit(size 2) + ksw_ll_i16 with a positive gap open penalty (as
 * minimap2 calls it; -O must be positive): the returned score is the best local alignment score;
 * (qe, te) is a cell where a best local alignment ends. */
static int case_ll_i16(int argc, char *argv[])
{
	int n, max_len, it;
	uint8_t *t, *qq;
	if (argc < 3) return MMT_FAIL("usage: ll_i16 <seed> <n> <max_len>"), 1;
	rng_state = strtoull(argv[0], 0, 10), n = atoi(argv[1]), max_len = atoi(argv[2]);
	t = (uint8_t*)malloc(max_len + 64), qq = (uint8_t*)malloc(max_len + 64);
	for (it = 0; it < n; ++it) {
		int a = rng_int(1, 3), b = rng_int(1, 8), gapo = rng_int(1, 12), gape = rng_int(1, 4), tlen, qlen, qe, te, sc, k;
		int32_t want, *Hm;
		int8_t mat[25];
		void *qp;
		rand_pair(max_len, it % 5 == 0? 20 : 0, t, &tlen, qq, &qlen);
		if (it % 2) { /* unrelated flanks around the related part */
			int fl = rng_int(1, 30);
			memmove(qq + fl, qq, qlen);
			for (k = 0; k < fl; ++k) qq[k] = rng_int(0, 3);
			qlen += fl;
		}
		gen_mat(mat, a, b, 0, 1);
		Hm = (int32_t*)malloc((size_t)(qlen + 1) * (tlen + 1) * 4);
		want = ref_local(qlen, qq, tlen, t, mat, gapo, gape, Hm);
		qp = ksw_ll_qinit(0, 2, qlen, qq, 5, mat);
		sc = ksw_ll_i16(qp, tlen, t, gapo, gape, &qe, &te);
		MMT_CHECKF(sc == want, "case %d (qlen %d, tlen %d, a%d b%d o%d e%d): score %d, reference %d", it, qlen, tlen, a, b, gapo, gape, sc, want);
		if (want > 0)
			MMT_CHECKF(qe >= 0 && qe < qlen && te >= 0 && te < tlen && Hm[(size_t)(qe + 1) * (tlen + 1) + te + 1] == want,
					   "case %d: end cell (%d,%d) does not score %d", it, qe, te, want);
		free(qp); free(Hm);
	}
	free(t); free(qq);
	return 0;
}

/* args: none. ksw_exts2_sse: exon1 (60 bp) + exon2 (60 bp) against exon1 + GT..AG intron (400 bp) +
 * exon2 aligns as 60M400N60M with KSW_EZ_SPLICE_FOR; the reverse complement against the reverse
 * complement of the target (an intron reading CT..AC) as 60M400N60M with KSW_EZ_SPLICE_REV. A gap
 * longer than the long-gap threshold is an intron (N) with or without the splice-signal flags
 * (minimap2.1, splice mode: long deletions are taken as introns). The splice preset's scoring:
 * -A1 -B2 -O2,32 -E1. */
static int case_exts2_intron(int argc, char *argv[])
{
	char ex1[64], ex2[64], intr[512], tgt[1024], qry[256], rtgt[1024], rqry[256];
	uint8_t t[1024], qq[256];
	int8_t mat[25];
	int i, tl, ql, v;
	ksw_extz_t ez;
	rng_state = 4242;
	for (i = 0; i < 60; ++i) ex1[i] = "ACGT"[rng_int(0, 3)], ex2[i] = "ACGT"[rng_int(0, 3)];
	ex1[59] = 'C', ex2[0] = 'T', ex1[60] = ex2[60] = 0; /* no G next to the splice sites: their position is unambiguous */
	for (i = 0; i < 400; ++i) intr[i] = "ACGT"[rng_int(0, 3)];
	intr[0] = 'G', intr[1] = 'T', intr[398] = 'A', intr[399] = 'G', intr[400] = 0;
	sprintf(tgt, "%s%s%s", ex1, intr, ex2);
	sprintf(qry, "%s%s", ex1, ex2);
	tl = strlen(tgt), ql = strlen(qry);
	for (i = 0; i < tl; ++i) rtgt[tl - 1 - i] = tgt[i] == 'A'? 'T' : tgt[i] == 'C'? 'G' : tgt[i] == 'G'? 'C' : 'A';
	for (i = 0; i < ql; ++i) rqry[ql - 1 - i] = qry[i] == 'A'? 'T' : qry[i] == 'C'? 'G' : qry[i] == 'G'? 'C' : 'A';
	rtgt[tl] = rqry[ql] = 0;
	gen_mat(mat, 1, 2, 0, 1);
	for (v = 0; v < 3; ++v) {
		exts2_f f = exts2_variant(v);
		if (f == 0) continue;
		enc(tgt, t), enc(qry, qq);
		memset(&ez, 0, sizeof(ez));
		f(0, ql, qq, tl, t, 5, mat, 2, 1, 32, 9, -1, -1, 9, 5, KSW_EZ_SPLICE_FOR, 0, &ez);
		MMT_CHECKF(ez.n_cigar == 3 && ez.cigar[0] == (60<<4|KSW_CIGAR_MATCH) && ez.cigar[1] == (400<<4|KSW_CIGAR_N_SKIP) && ez.cigar[2] == (60<<4|KSW_CIGAR_MATCH),
				   "variant %d: forward intron not 60M400N60M (n_cigar %d)", v, ez.n_cigar);
		MMT_CHECKF(ez.score > 0 && ez.score <= 120 - 32, "variant %d: score %d", v, ez.score);
		ez_free(&ez);
		enc(rtgt, t), enc(rqry, qq);
		f(0, ql, qq, tl, t, 5, mat, 2, 1, 32, 9, -1, -1, 9, 5, KSW_EZ_SPLICE_REV, 0, &ez);
		MMT_CHECKF(ez.n_cigar == 3 && ez.cigar[0] == (60<<4|KSW_CIGAR_MATCH) && ez.cigar[1] == (400<<4|KSW_CIGAR_N_SKIP) && ez.cigar[2] == (60<<4|KSW_CIGAR_MATCH),
				   "variant %d: reverse intron not 60M400N60M (n_cigar %d)", v, ez.n_cigar);
		ez_free(&ez);
		enc(tgt, t), enc(qry, qq);
		f(0, ql, qq, tl, t, 5, mat, 2, 1, 32, 9, -1, -1, 9, 5, 0, 0, &ez); /* no splice signal scoring: the intron is still placed exactly */
		MMT_CHECKF(ez.n_cigar == 3 && ez.cigar[0] == (60<<4|KSW_CIGAR_MATCH) && ez.cigar[1] == (400<<4|KSW_CIGAR_N_SKIP) && ez.cigar[2] == (60<<4|KSW_CIGAR_MATCH),
				   "variant %d: intron without a splice flag not 60M400N60M (n_cigar %d)", v, ez.n_cigar);
		ez_free(&ez);
		enc(tgt, t), enc(qry, qq);
		f(0, ql, qq, 200, t, 5, mat, 2, 1, 32, 9, -1, -1, 9, 5, KSW_EZ_SPLICE_FOR, 0, &ez); /* target cut inside the intron: a long deletion */
		MMT_CHECKF(ez.score < 60, "variant %d: score %d aligning into an incomplete intron", v, ez.score);
		ez_free(&ez);
	}
	return 0;
}

/* args: none. ksw_exts2_sse without splice signals (-A1 -B2 -O2,32 -E1,0, the splice preset): a
 * k-base deletion costs q + k*e = 2 + k as a deletion and q2 = 32 as an intron (minimap2.1, splice
 * mode: long deletions are introns; the second gap extension is not used), so 60M kD 60M is a D for
 * 2 + k < 32 and an N for 2 + k > 32 (k = 30 ties and is skipped). */
static int case_exts2_long_gap(int argc, char *argv[])
{
	uint8_t t[256], qq[128];
	int8_t mat[25];
	int i, k, v;
	rng_state = 515;
	gen_mat(mat, 1, 2, 0, 1);
	for (k = 20; k <= 40; ++k) {
		if (k == 30) continue;
		for (i = 0; i < 120 + k; ++i) t[i] = rng_int(0, 3);
		for (i = 0; i < 120; ++i) qq[i] = i < 60? t[i] : t[i + k];
		t[59] = 1, t[60 + k] = 3, qq[59] = 1, qq[60] = 3; /* exon ends differ from the gap's ends: one placement */
		t[60] = 0, t[59 + k] = 0;
		for (v = 0; v < 3; ++v) {
			exts2_f f = exts2_variant(v);
			ksw_extz_t ez;
			int op = 2 + k < 32? KSW_CIGAR_DEL : KSW_CIGAR_N_SKIP;
			if (f == 0) continue;
			memset(&ez, 0, sizeof(ez));
			f(0, 120, qq, 120 + k, t, 5, mat, 2, 1, 32, 9, -1, -1, 9, 5, 0, 0, &ez);
			MMT_CHECKF(ez.n_cigar == 3 && ez.cigar[0] == (60<<4|KSW_CIGAR_MATCH) && ez.cigar[1] == ((uint32_t)k<<4|op) && ez.cigar[2] == (60<<4|KSW_CIGAR_MATCH),
					   "variant %d k %d: want 60M%d%c60M (n_cigar %d)", v, k, k, op == KSW_CIGAR_DEL? 'D' : 'N', ez.n_cigar);
			MMT_CHECKF(ez.score == 120 - (2 + k < 32? 2 + k : 32), "variant %d k %d: score %d", v, k, ez.score);
			ez_free(&ez);
		}
	}
	return 0;
}

/* args: seed n max_len. The SSE2 and SSE4.1 builds of each kernel (present in the runtime-dispatch
 * build on x86 ELF platforms) return the same scores, positions and CIGARs as the entry point. */
static int case_variants(int argc, char *argv[])
{
	int n, max_len, it, v, n_variants = 0;
	uint8_t *t, *qq;
	if (argc < 3) return MMT_FAIL("usage: variants <seed> <n> <max_len>"), 1;
	rng_state = strtoull(argv[0], 0, 10), n = atoi(argv[1]), max_len = atoi(argv[2]);
	for (v = 1; v <= 2; ++v) n_variants += extz2_variant(v) && extd2_variant(v) && exts2_variant(v);
	if (n_variants == 0) {
		mmt_note("no per-variant kernels in this build (sse2only, NEON or non-ELF); nothing to compare");
		MMT_CHECK(1);
		return 0;
	}
	t = (uint8_t*)malloc(max_len + 1), qq = (uint8_t*)malloc(max_len + 1);
	for (it = 0; it < n; ++it) {
		int a, b, q, e, q2, e2, tlen, qlen, kind = it % 2, w = it % 4 == 0? -1 : rng_int(0, 50), zdrop = it % 2? rng_int(10, 200) : -1;
		int flag = (it % 5 == 0? KSW_EZ_EXTZ_ONLY : 0) | (it % 7 == 0? KSW_EZ_RIGHT : 0) | (it % 6 == 0? KSW_EZ_APPROX_MAX : 0) | (it % 9 == 0? KSW_EZ_GENERIC_SC : 0);
		int8_t mat[25];
		ksw_extz_t ez0, ez1;
		rand_params(kind, &a, &b, &q, &e, &q2, &e2);
		rand_pair(max_len, 10, t, &tlen, qq, &qlen);
		gen_mat(mat, a, b, flag & KSW_EZ_GENERIC_SC? rng_int(1, b) : 0, 1);
		memset(&ez0, 0, sizeof(ez0));
		run_kernel(kind, 0, qlen, qq, tlen, t, mat, q, e, q2, e2, w, zdrop, 5, flag, &ez0);
		for (v = 1; v <= 2; ++v) {
			int k, same;
			if (kind == 0? extz2_variant(v) == 0 : extd2_variant(v) == 0) continue;
			memset(&ez1, 0, sizeof(ez1));
			run_kernel(kind, v, qlen, qq, tlen, t, mat, q, e, q2, e2, w, zdrop, 5, flag, &ez1);
			same = ez0.score == ez1.score && ez0.max == ez1.max && ez0.max_q == ez1.max_q && ez0.max_t == ez1.max_t && ez0.mqe == ez1.mqe
				&& ez0.mqe_t == ez1.mqe_t && ez0.mte == ez1.mte && ez0.mte_q == ez1.mte_q && ez0.zdropped == ez1.zdropped
				&& ez0.reach_end == ez1.reach_end && ez0.n_cigar == ez1.n_cigar;
			for (k = 0; same && k < ez0.n_cigar; ++k)
				if (ez0.cigar[k] != ez1.cigar[k]) same = 0;
			MMT_CHECKF(same, "case %d (kind %d, flag %x): variant %d differs from the entry point", it, kind, flag, v);
			ez_free(&ez1);
		}
		ez_free(&ez0);
	}
	free(t); free(qq);
	return 0;
}

/* args: seed n max_len. As variants, for the spliced kernel. */
static int case_variants_splice(int argc, char *argv[])
{
	int n, max_len, it, v;
	uint8_t *t, *qq;
	if (argc < 3) return MMT_FAIL("usage: variants_splice <seed> <n> <max_len>"), 1;
	rng_state = strtoull(argv[0], 0, 10), n = atoi(argv[1]), max_len = atoi(argv[2]);
	if (!exts2_variant(1) && !exts2_variant(2)) {
		mmt_note("no per-variant kernels in this build; nothing to compare");
		MMT_CHECK(1);
		return 0;
	}
	t = (uint8_t*)malloc(max_len + 1), qq = (uint8_t*)malloc(max_len + 1);
	for (it = 0; it < n; ++it) {
		int a = rng_int(1, 2), b = rng_int(1, 4), q = rng_int(1, 6), e = 1, q2 = q + e + rng_int(1, 30), tlen, qlen, zdrop = it % 2? rng_int(20, 200) : -1;
		int flag = (it % 3 == 0? KSW_EZ_SPLICE_FOR : it % 3 == 1? KSW_EZ_SPLICE_REV : 0) | (it % 4 == 0? KSW_EZ_EXTZ_ONLY : 0)
			| (it % 5 == 0? KSW_EZ_SPLICE_FLANK : 0) | (it % 2 == 0? KSW_EZ_SPLICE_CMPLX : 0) | (it % 7 == 0? KSW_EZ_APPROX_MAX : 0);
		int8_t mat[25];
		ksw_extz_t ez0, ez1;
		rand_pair(max_len, 0, t, &tlen, qq, &qlen);
		gen_mat(mat, a, b, 0, 1);
		memset(&ez0, 0, sizeof(ez0));
		ksw_exts2_sse(0, qlen, qq, tlen, t, 5, mat, q, e, q2, 9, zdrop, 5, 9, 5, flag, 0, &ez0);
		for (v = 1; v <= 2; ++v) {
			exts2_f f = exts2_variant(v);
			int k, same;
			if (!f) continue;
			memset(&ez1, 0, sizeof(ez1));
			f(0, qlen, qq, tlen, t, 5, mat, q, e, q2, 9, zdrop, 5, 9, 5, flag, 0, &ez1);
			same = ez0.score == ez1.score && ez0.max == ez1.max && ez0.max_q == ez1.max_q && ez0.max_t == ez1.max_t && ez0.mqe == ez1.mqe
				&& ez0.mte == ez1.mte && ez0.zdropped == ez1.zdropped && ez0.reach_end == ez1.reach_end && ez0.n_cigar == ez1.n_cigar;
			for (k = 0; same && k < ez0.n_cigar; ++k)
				if (ez0.cigar[k] != ez1.cigar[k]) same = 0;
			MMT_CHECKF(same, "case %d (flag %x): spliced variant %d differs from the entry point", it, flag, v);
			ez_free(&ez1);
		}
		ez_free(&ez0);
	}
	free(t); free(qq);
	return 0;
}

/* the largest |t - q| over the cells a CIGAR passes through, from (0,0) (-1 for a malformed CIGAR) */
static int cigar_max_diag(int n_cigar, const uint32_t *cigar)
{
	int k, i = 0, j = 0, d = 0;
	for (k = 0; k < n_cigar; ++k) {
		int op = cigar[k] & 0xf, len = cigar[k] >> 4;
		if (op == KSW_CIGAR_MATCH) i += len, j += len;
		else if (op == KSW_CIGAR_INS) i += len;
		else if (op == KSW_CIGAR_DEL) j += len;
		else return -1;
		if ((j - i > 0? j - i : i - j) > d) d = j - i > 0? j - i : i - j;
	}
	return d;
}

/* args: seed n max_len out. Narrow bands: w is |qlen - tlen| plus 0 to 40, so the end cell of a global
 * alignment lies on or near the band edge, where the kernels may report a score the banded optimum
 * does not reach. What the kernels guarantee there: the global
 * CIGAR covers both sequences and is a valid alignment, which scores no more than the unbanded
 * optimum; an extension's CIGAR ends at the cell it reports. Every result (score, max and end cells,
 * flags, CIGAR) is written to out, one line per case, for a pin. */
static int case_band_edge(int argc, char *argv[])
{
	int n, max_len, it, n_diff = 0, n_outside = 0, n_edge = 0;
	uint8_t *t, *qq;
	FILE *fp;
	if (argc < 4) return MMT_FAIL("usage: band_edge <seed> <n> <max_len> <out>"), 1;
	rng_state = strtoull(argv[0], 0, 10), n = atoi(argv[1]), max_len = atoi(argv[2]);
	MMT_REQUIRE((fp = fopen(argv[3], "w")) != 0);
	t = (uint8_t*)malloc(max_len + 64), qq = (uint8_t*)malloc(max_len + 64);
	for (it = 0; it < n; ++it) {
		int a, b, q, e, q2, e2, tlen, qlen, w, qspan, tspan, kind = it % 2, ext = it % 3 == 2, k;
		int flag = (ext? KSW_EZ_EXTZ_ONLY : 0) | (it % 5 == 1? KSW_EZ_RIGHT : 0) | (it % 7 == 3? KSW_EZ_APPROX_MAX : 0);
		int zdrop = it % 4 == 3? rng_int(20, 100) : -1, end_bonus = ext? rng_int(0, 20) : 0;
		int32_t want, cs;
		int8_t mat[25];
		ref_t R;
		ksw_extz_t ez;
		rand_params(kind, &a, &b, &q, &e, &q2, &e2);
		rand_pair(max_len, it % 4 == 1? 20 : 0, t, &tlen, qq, &qlen);
		if (it % 3 == 0) { /* a long indel near the middle, so the path drifts towards the band edge */
			int l = rng_int(5, 40), pos;
			if (rng_int(0, 1)) {
				pos = qlen / 2;
				memmove(qq + pos + l, qq + pos, qlen - pos);
				for (k = 0; k < l; ++k) qq[pos + k] = rng_int(0, 3);
				qlen += l;
			} else {
				pos = tlen / 2;
				memmove(t + pos + l, t + pos, tlen - pos);
				for (k = 0; k < l; ++k) t[pos + k] = rng_int(0, 3);
				tlen += l;
			}
		}
		w = (qlen > tlen? qlen - tlen : tlen - qlen) + rng_int(0, 40);
		gen_mat(mat, a, b, 0, 1);
		ref_dp(&R, qlen, qq, tlen, t, mat, q, e, q2, e2, -1);
		want = ref_at(&R, qlen, tlen);
		memset(&ez, 0, sizeof(ez));
		run_kernel(kind, 0, qlen, qq, tlen, t, mat, q, e, q2, e2, w, zdrop, end_bonus, flag, &ez);
		fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t", it, kind, qlen, tlen, w, flag, ez.score, ez.max, ez.max_q, ez.max_t,
				ez.mqe, ez.mqe_t, ez.mte, ez.mte_q, ez.zdropped, ez.reach_end);
		for (k = 0; k < ez.n_cigar; ++k) fprintf(fp, "%d%c", ez.cigar[k]>>4, "MIDN"[ez.cigar[k]&0xf]);
		fputc('\n', fp);
		if (ez.n_cigar > 0) {
			int ok = cigar_eval(ez.n_cigar, ez.cigar, qq, t, mat, q, e, q2, e2, &qspan, &tspan, &cs);
			int d = cigar_max_diag(ez.n_cigar, ez.cigar);
			MMT_CHECKF(ok, "case %d (kind %d, qlen %d, tlen %d, w %d, flag %x): malformed CIGAR", it, kind, qlen, tlen, w, flag);
			if (!ok) { ref_free(&R); ez_free(&ez); continue; }
			if (!ext && !ez.zdropped) {
				MMT_CHECKF(cs <= want, "case %d: the CIGAR scores %d, above the unbanded optimum %d", it, cs, want);
				MMT_CHECKF(qspan == qlen && tspan == tlen, "case %d (kind %d, qlen %d, tlen %d, w %d, flag %x): global CIGAR spans %d/%d", it, kind, qlen, tlen, w, flag, qspan, tspan);
				n_diff += cs != ez.score;
			} else if (ext && ez.reach_end) {
				MMT_CHECKF(qspan == qlen && tspan == ez.mqe_t + 1, "case %d: reach_end CIGAR spans %d/%d, want %d/%d", it, qspan, tspan, qlen, ez.mqe_t + 1);
			} else if (ext) {
				MMT_CHECKF(qspan == ez.max_q + 1 && tspan == ez.max_t + 1, "case %d: extension CIGAR spans %d/%d, want %d/%d", it, qspan, tspan, ez.max_q + 1, ez.max_t + 1);
			}
			n_outside += d > w;
			n_edge += d >= w - 2;
		}
		ref_free(&R); ez_free(&ez);
	}
	mmt_note("band_edge: %d global scores differ from their CIGAR's score; %d CIGARs leave the band; %d come within 2 of its edge", n_diff, n_outside, n_edge);
	fclose(fp);
	free(t); free(qq);
	return 0;
}

/* Z-drop of an unbanded extension by brute force (minimap2 paper, Z-drop; ksw_apply_zdrop): for each
 * anti-diagonal r = q + t, its best cell (H, t); a new best score moves the max cell, otherwise a cell
 * down and right of the max cell drops max - H - e*|(t - max_t) - (q - max_q)|. Returns the largest
 * drop, or INT32_MIN if an anti-diagonal's best cell is tied so that the drop or the max cell depends
 * on which tied cell a kernel picks (not documented). */
static int32_t ref_zdrop(const ref_t *R, int e)
{
	int r, i, j;
	int32_t mx = 0, max_t = -1, max_q = -1, best = INT32_MIN;
	for (r = 0; r <= R->qlen + R->tlen - 2; ++r) {
		int32_t h = INT32_MIN, d = INT32_MIN, n_tie = 0, d_set = 0, amb = 0;
		for (j = 0; j <= r && j < R->tlen; ++j) { /* t = j, q = r - j */
			i = r - j;
			if (i >= R->qlen) continue;
			if (ref_at(R, i + 1, j + 1) > h) h = ref_at(R, i + 1, j + 1), n_tie = 1;
			else if (ref_at(R, i + 1, j + 1) == h) ++n_tie;
		}
		if (h > mx) {
			if (n_tie > 1) return INT32_MIN;
			for (j = 0; j <= r && j < R->tlen; ++j)
				if (r - j < R->qlen && ref_at(R, r - j + 1, j + 1) == h) max_t = j, max_q = r - j;
			mx = h;
			continue;
		}
		for (j = 0; j <= r && j < R->tlen; ++j) { /* every tied best cell must give the same verdict */
			int32_t dj = INT32_MIN;
			i = r - j;
			if (i >= R->qlen || ref_at(R, i + 1, j + 1) != h) continue;
			if (j >= max_t && i >= max_q) {
				int tl = j - max_t, ql = i - max_q;
				dj = mx - h - e * (tl > ql? tl - ql : ql - tl);
			}
			if (d_set && dj != d) amb = 1;
			d = dj, d_set = 1;
		}
		if (amb) return INT32_MIN;
		if (d > best) best = d;
	}
	return best;
}

/* args: seed n. The Z-drop threshold is strict: an exact match of 20-80 bases followed by unrelated
 * tails, extended with KSW_EZ_EXTZ_ONLY and no band; with Z the largest drop along the extension
 * (brute force), zdrop Z runs to the end and zdrop Z-1 stops (zdropped). ksw_extd2 weighs the change
 * of diagonal by its second gap extension e2, ksw_extz2 by e. Cases whose anti-diagonal maxima are tied
 * in a way that changes the drop are skipped. */
static int case_zdrop_boundary(int argc, char *argv[])
{
	int n, it, n_checked = 0;
	uint8_t t[256], qq[256];
	if (argc < 2) return MMT_FAIL("usage: zdrop_boundary <seed> <n>"), 1;
	rng_state = strtoull(argv[0], 0, 10), n = atoi(argv[1]);
	for (it = 0; it < n; ++it) {
		int a, b, q, e, q2, e2, kind = it % 2, i, L = rng_int(20, 80), tl = L + rng_int(20, 100), ql = L + rng_int(20, 100);
		int32_t Z;
		int8_t mat[25];
		ref_t R;
		ksw_extz_t ez;
		rand_params(kind, &a, &b, &q, &e, &q2, &e2);
		for (i = 0; i < L; ++i) t[i] = qq[i] = rng_int(0, 3);
		for (i = L; i < tl; ++i) t[i] = rng_int(0, 3);
		for (i = L; i < ql; ++i) qq[i] = rng_int(0, 3);
		gen_mat(mat, a, b, 0, 1);
		ref_dp(&R, ql, qq, tl, t, mat, q, e, q2, e2, -1);
		Z = ref_zdrop(&R, kind? e2 : e);
		ref_free(&R);
		if (Z == INT32_MIN || Z < 1) continue;
		++n_checked;
		memset(&ez, 0, sizeof(ez));
		run_kernel(kind, 0, ql, qq, tl, t, mat, q, e, q2, e2, -1, Z, 0, KSW_EZ_EXTZ_ONLY, &ez);
		MMT_CHECKF(!ez.zdropped, "case %d (kind %d, L %d, qlen %d, tlen %d, a%d b%d q%d e%d q2%d e2%d): zdrop %d = the largest drop stopped the extension", it, kind, L, ql, tl, a, b, q, e, q2, e2, Z);
		ez_free(&ez);
		run_kernel(kind, 0, ql, qq, tl, t, mat, q, e, q2, e2, -1, Z - 1, 0, KSW_EZ_EXTZ_ONLY, &ez);
		MMT_CHECKF(ez.zdropped, "case %d (kind %d, L %d, qlen %d, tlen %d, a%d b%d q%d e%d q2%d e2%d): zdrop %d did not stop the extension", it, kind, L, ql, tl, a, b, q, e, q2, e2, Z - 1);
		ez_free(&ez);
	}
	mmt_note("zdrop_boundary: %d of %d cases without ties checked", n_checked, n);
	MMT_CHECKF(n_checked >= n / 6, "too few cases without ties (%d of %d)", n_checked, n);
	return 0;
}

static const mmt_case_t cases[] = {
	{ "global", case_global },
	{ "extension", case_extension },
	{ "generic_matrix", case_generic_matrix },
	{ "gap_alignment", case_gap_alignment },
	{ "zdrop", case_zdrop },
	{ "zdrop_boundary", case_zdrop_boundary },
	{ "band_deletion", case_band_deletion },
	{ "band_edge", case_band_edge },
	{ "empty", case_empty },
	{ "ll_i16", case_ll_i16 },
	{ "exts2_intron", case_exts2_intron },
	{ "exts2_long_gap", case_exts2_long_gap },
	{ "variants", case_variants },
	{ "variants_splice", case_variants_splice },
	{ 0, 0 }
};

MMT_MAIN(cases)
