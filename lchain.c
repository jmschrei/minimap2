#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include "mmpriv.h"
#include "kalloc.h"
#include "krmq.h"

// x86 only: on other targets the compiler may fuse the scalar and vector float ops differently
#if defined(__SSE2__) && (defined(__x86_64__) || defined(__i386__)) && defined(FLT_EVAL_METHOD) && FLT_EVAL_METHOD == 0
#define MG_LCHAIN_SSE2
#include <emmintrin.h>
#endif

/* MG_LCHAIN_AVX2_UNIT: this file compiled again with -mavx2 (Makefile, lchain_avx2.o) to hold only mg_lchain_dp_avx2(),
 * mg_lchain_dp() with 8-lane comput_sc8(); MG_LCHAIN_HAVE_AVX2: lchain.o calls it when the CPU and OS support AVX2 */
#if defined(MG_LCHAIN_AVX2_UNIT) && !(defined(MG_LCHAIN_SSE2) && defined(__AVX2__))
#undef MG_LCHAIN_SSE2 // the AVX2 unit is then empty
#endif
#if defined(MG_LCHAIN_AVX2_UNIT) && defined(MG_LCHAIN_SSE2)
#include <immintrin.h>
#define MG_LC_AVX2
#endif
#if !defined(MG_LCHAIN_AVX2_UNIT) && defined(MG_LCHAIN_HAVE_AVX2) && defined(MG_LCHAIN_SSE2)
#include <pthread.h>
#endif

#ifndef MG_LCHAIN_AVX2_UNIT
static int64_t mg_chain_bk_end(int32_t max_drop, const mm128_t *z, const int32_t *f, const int64_t *p, int32_t *t, int64_t k)
{
	int64_t i = z[k].y, end_i = -1, max_i = i;
	int32_t max_s = 0;
	if (i < 0 || t[i] != 0) return i;
	do {
		int32_t s;
		t[i] = 2;
		end_i = i = p[i];
		s = i < 0? z[k].x : (int32_t)z[k].x - f[i];
		if (s > max_s) max_s = s, max_i = i;
		else if (max_s - s > max_drop) break;
	} while (i >= 0 && t[i] == 0);
	for (i = z[k].y; i >= 0 && i != end_i; i = p[i]) // reset modified t[]
		t[i] = 0;
	return max_i;
}

uint64_t *mg_chain_backtrack(void *km, int64_t n, const int32_t *f, const int64_t *p, int32_t *v, int32_t *t, int32_t min_cnt, int32_t min_sc, int32_t max_drop, int32_t *n_u_, int32_t *n_v_)
{
	mm128_t *z;
	uint64_t *u;
	int64_t i, k, n_z, n_v;
	int32_t n_u;

	*n_u_ = *n_v_ = 0;
	for (i = 0, n_z = 0; i < n; ++i) // precompute n_z
		if (f[i] >= min_sc) ++n_z;
	if (n_z == 0) return 0;
	z = Kmalloc(km, mm128_t, n_z);
	for (i = 0, k = 0; i < n; ++i) // populate z[]
		if (f[i] >= min_sc) z[k].x = f[i], z[k++].y = i;
	radix_sort_128x(z, z + n_z);

	memset(t, 0, n * 4);
	for (k = n_z - 1, n_v = n_u = 0; k >= 0; --k) { // populate v[] and u[]
		if (t[z[k].y] == 0) {
			int64_t n_v0 = n_v, end_i;
			int32_t sc;
			end_i = mg_chain_bk_end(max_drop, z, f, p, t, k);
			for (i = z[k].y; i != end_i; i = p[i])
				v[n_v++] = i, t[i] = 1;
			sc = i < 0? z[k].x : (int32_t)z[k].x - f[i];
			if (sc >= min_sc && n_v > n_v0 && n_v - n_v0 >= min_cnt) {
				z[n_z - 1 - n_u].x = (uint64_t)sc << 32 | (n_v - n_v0); // at most one chain per z[k], so this slot is z[k] or already read
				++n_u;
			} else n_v = n_v0;
		}
	}
	u = Kmalloc(km, uint64_t, n_u);
	for (k = 0; k < n_u; ++k)
		u[k] = z[n_z - 1 - k].x;
	kfree(km, z);
	assert(n_v < INT32_MAX);
	*n_u_ = n_u, *n_v_ = n_v;
	return u;
}
#else
uint64_t *mg_chain_backtrack(void *km, int64_t n, const int32_t *f, const int64_t *p, int32_t *v, int32_t *t, int32_t min_cnt, int32_t min_sc, int32_t max_drop, int32_t *n_u_, int32_t *n_v_);
#endif // MG_LCHAIN_AVX2_UNIT

#if !defined(MG_LCHAIN_AVX2_UNIT) || defined(MG_LC_AVX2)
static mm128_t *compact_a(void *km, int32_t n_u, uint64_t *u, int32_t n_v, int32_t *v, mm128_t *a)
{
	mm128_t *b, *w;
	uint64_t *u2;
	int64_t i, j, k;

	// write the result to b[]
	b = Kmalloc(km, mm128_t, n_v);
	for (i = 0, k = 0; i < n_u; ++i) {
		int32_t k0 = k, ni = (int32_t)u[i];
		for (j = 0; j < ni; ++j)
			b[k++] = a[v[k0 + (ni - j - 1)]];
	}
	kfree(km, v);

	// sort u[] and a[] by the target position, such that adjacent chains may be joined
	w = Kmalloc(km, mm128_t, n_u);
	for (i = k = 0; i < n_u; ++i) {
		w[i].x = b[k].x, w[i].y = (uint64_t)k<<32|i;
		k += (int32_t)u[i];
	}
	radix_sort_128x(w, w + n_u);
	u2 = Kmalloc(km, uint64_t, n_u);
	for (i = k = 0; i < n_u; ++i) {
		int32_t j = (int32_t)w[i].y, n = (int32_t)u[j];
		u2[i] = u[j];
		memcpy(&a[k], &b[w[i].y>>32], n * sizeof(mm128_t));
		k += n;
	}
	memcpy(u, u2, n_u * 8);
	memcpy(b, a, k * sizeof(mm128_t)); // write _a_ to _b_ and deallocate _a_ because _a_ is oversized, sometimes a lot
	kfree(km, a); kfree(km, w); kfree(km, u2);
	return b;
}

static inline int32_t comput_sc(const mm128_t *ai, const mm128_t *aj, int32_t max_dist_x, int32_t max_dist_y, int32_t bw, float chn_pen_gap, float chn_pen_skip, int is_cdna, int n_seg)
{
	int32_t dq = (int32_t)ai->y - (int32_t)aj->y, dr, dd, dg, q_span, sc;
	int32_t sidi = (ai->y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT;
	int32_t sidj = (aj->y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT;
	if (dq <= 0 || dq > max_dist_x) return INT32_MIN;
	dr = (int32_t)(ai->x - aj->x);
	if (sidi == sidj && (dr == 0 || dq > max_dist_y)) return INT32_MIN;
	dd = dr > dq? dr - dq : dq - dr;
	if (sidi == sidj && dd > bw) return INT32_MIN;
	if (n_seg > 1 && !is_cdna && sidi == sidj && dr > max_dist_y) return INT32_MIN;
	dg = dr < dq? dr : dq;
	q_span = aj->y>>32&0xff;
	sc = q_span < dg? q_span : dg;
	if (dd || dg > q_span) {
		float lin_pen, log_pen;
		lin_pen = chn_pen_gap * (float)dd + chn_pen_skip * (float)dg;
		log_pen = dd >= 1? mg_log2(dd + 1) : 0.0f; // mg_log2() only works for dd>=2
		if (is_cdna || sidi != sidj) {
			if (sidi != sidj && dr == 0) ++sc; // possibly due to overlapping paired ends; give a minor bonus
			else if (dr > dq || sidi != sidj) sc -= (int)(lin_pen < log_pen? lin_pen : log_pen); // deletion or jump between paired ends
			else sc -= (int)(lin_pen + .5f * log_pen);
		} else sc -= (int)(lin_pen + .5f * log_pen);
	}
	return sc;
}

#ifdef MG_LCHAIN_SSE2
#ifdef MG_LC_AVX2
#define MG_LC_MIN 8  // fewer predecessors than this take the scalar loop
#define MG_LC_W   8  // lanes per vector
#define MG_LC_BLK 32 // predecessors scored per block; a multiple of MG_LC_W, at most 32
typedef __m256i lc_vi_t;
typedef __m256 lc_vf_t;
#define lc_set1_epi32 _mm256_set1_epi32
#define lc_set1_ps    _mm256_set1_ps
#define lc_storeu(p, v) _mm256_storeu_si256((__m256i*)(p), (v))
#else
#define MG_LC_MIN 16
#define MG_LC_W   4
#define MG_LC_BLK 16
typedef __m128i lc_vi_t;
typedef __m128 lc_vf_t;
#define lc_set1_epi32 _mm_set1_epi32
#define lc_set1_ps    _mm_set1_ps
#define lc_storeu(p, v) _mm_storeu_si128((__m128i*)(p), (v))
#endif

typedef struct {
	lc_vi_t xi, yi, sidi, mdx, mdy, bw, seg_chk, cdna, mdxy, di, bwb, bw2b;
	lc_vf_t gap, skip;
} lc_sse_t;

#ifndef MG_LC_AVX2
static const uint8_t lc_rev4[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

/* comput_sc() for the predecessors j-3, j-2, j-1 and j (lanes 0-3) from the low 32 bits of a[].x and of a[].y and the
 * high 32 bits of a[].y; the same int and float ops in the same order. Returns the lanes whose score is not INT32_MIN,
 * bit 3-l for lane l, and the scores in *out when it is not 0. */
static inline uint32_t comput_sc4(const int32_t *xs, const int32_t *ys, const int32_t *hs, int64_t j, const lc_sse_t *c, __m128i *out)
{
	const __m128i one = _mm_set1_epi32(1), zero = _mm_setzero_si128(), ones = _mm_set1_epi32(-1), ff = _mm_set1_epi32(0xff);
	__m128i xj = _mm_loadu_si128((const __m128i*)&xs[j-3]), yj = _mm_loadu_si128((const __m128i*)&ys[j-3]), hj = _mm_loadu_si128((const __m128i*)&hs[j-3]);
	__m128i dq, dr, same, inv, s_inv, drz, lt, dd, dg, qs, qlt, sc, cond, nsame, bonus, use_a, pen, zi, r;
	__m128 lin, z, l2, log_pen;
	dq = _mm_sub_epi32(c->yi, yj);
	dr = _mm_sub_epi32(c->xi, xj);
	same = _mm_cmpeq_epi32(c->sidi, _mm_and_si128(_mm_srli_epi32(hj, MM_SEED_SEG_SHIFT - 32), ff));
	inv = _mm_or_si128(_mm_cmplt_epi32(dq, one), _mm_cmpgt_epi32(dq, c->mdx));
	drz = _mm_cmpeq_epi32(dr, zero);
	lt = _mm_cmplt_epi32(dr, dq);
	dd = _mm_sub_epi32(_mm_xor_si128(_mm_sub_epi32(dr, dq), lt), lt); // dr < dq? dq - dr : dr - dq
	s_inv = _mm_or_si128(_mm_or_si128(drz, _mm_cmpgt_epi32(dq, c->mdy)), _mm_cmpgt_epi32(dd, c->bw));
	s_inv = _mm_or_si128(s_inv, _mm_and_si128(c->seg_chk, _mm_cmpgt_epi32(dr, c->mdy)));
	inv = _mm_or_si128(inv, _mm_and_si128(same, s_inv));
	if (_mm_movemask_epi8(inv) == 0xffff) return 0;
	dg = _mm_or_si128(_mm_and_si128(lt, dr), _mm_andnot_si128(lt, dq));
	qs = _mm_and_si128(hj, ff);
	qlt = _mm_cmpgt_epi32(dg, qs);
	sc = _mm_or_si128(_mm_and_si128(qlt, qs), _mm_andnot_si128(qlt, dg));
	cond = _mm_or_si128(_mm_xor_si128(_mm_cmpeq_epi32(dd, zero), ones), qlt);
	lin = _mm_add_ps(_mm_mul_ps(c->gap, _mm_cvtepi32_ps(dd)), _mm_mul_ps(c->skip, _mm_cvtepi32_ps(dg)));
	zi = _mm_castps_si128(_mm_cvtepi32_ps(_mm_add_epi32(dd, one))); // mg_log2(dd + 1)
	l2 = _mm_cvtepi32_ps(_mm_sub_epi32(_mm_and_si128(_mm_srli_epi32(zi, 23), ff), _mm_set1_epi32(128)));
	z = _mm_castsi128_ps(_mm_add_epi32(_mm_and_si128(zi, _mm_set1_epi32(~(255 << 23))), _mm_set1_epi32(127 << 23)));
	l2 = _mm_add_ps(l2, _mm_sub_ps(_mm_mul_ps(_mm_add_ps(_mm_mul_ps(_mm_set1_ps(-0.34484843f), z), _mm_set1_ps(2.02466578f)), z), _mm_set1_ps(0.67487759f)));
	log_pen = _mm_and_ps(_mm_castsi128_ps(_mm_cmpgt_epi32(dd, zero)), l2);
	nsame = _mm_xor_si128(same, ones);
	bonus = _mm_and_si128(nsame, drz);
	use_a = _mm_and_si128(_mm_or_si128(c->cdna, nsame), _mm_or_si128(_mm_cmpgt_epi32(dr, dq), nsame));
	pen = _mm_or_si128(_mm_and_si128(use_a, _mm_cvttps_epi32(_mm_min_ps(lin, log_pen))),
					   _mm_andnot_si128(use_a, _mm_cvttps_epi32(_mm_add_ps(lin, _mm_mul_ps(_mm_set1_ps(.5f), log_pen)))));
	pen = _mm_or_si128(bonus, _mm_andnot_si128(bonus, pen)); // -1: sc - (-1) is the ++sc bonus
	sc = _mm_sub_epi32(sc, _mm_and_si128(cond, pen));
	r = _mm_or_si128(_mm_and_si128(inv, _mm_set1_epi32(INT32_MIN)), _mm_andnot_si128(inv, sc));
	*out = r;
	return lc_rev4[~_mm_movemask_ps(_mm_castsi128_ps(_mm_cmpeq_epi32(r, _mm_set1_epi32(INT32_MIN)))) & 0xf];
}

/* comput_sc4() when every anchor is in one segment, so sidi == sidj in every lane: a lane is valid only if dq >= 1 and
 * dd <= bw, and a group with none is rejected before the other tests. With ds (x - y per anchor, a[] sorted by x, so
 * dr is in [0, max_dist_x] and dr - dq cannot wrap for a valid lane), a group whose lanes all have |dr - dq| > bw is
 * rejected first, from one load. */
static inline uint32_t comput_sc4_1(const int32_t *xs, const int32_t *ys, const int32_t *hs, const int32_t *ds, int64_t j, const lc_sse_t *c, __m128i *out)
{
	const __m128i one = _mm_set1_epi32(1), zero = _mm_setzero_si128(), ones = _mm_set1_epi32(-1), ff = _mm_set1_epi32(0xff);
	__m128i xj, yj, hj, dq, dr, inv, lt, dd, dg, qs, qlt, sc, cond, use_a, pen, zi, r;
	__m128 lin, z, l2, log_pen;
	if (ds) { // (dr - dq) + bw > 2 * bw, unsigned
		__m128i e = _mm_sub_epi32(c->di, _mm_loadu_si128((const __m128i*)&ds[j-3]));
		if (_mm_movemask_epi8(_mm_cmpgt_epi32(_mm_add_epi32(e, c->bwb), c->bw2b)) == 0xffff) return 0;
	}
	xj = _mm_loadu_si128((const __m128i*)&xs[j-3]), yj = _mm_loadu_si128((const __m128i*)&ys[j-3]);
	dq = _mm_sub_epi32(c->yi, yj);
	dr = _mm_sub_epi32(c->xi, xj);
	lt = _mm_cmplt_epi32(dr, dq);
	dd = _mm_sub_epi32(_mm_xor_si128(_mm_sub_epi32(dr, dq), lt), lt); // dr < dq? dq - dr : dr - dq
	inv = _mm_or_si128(_mm_cmplt_epi32(dq, one), _mm_cmpgt_epi32(dd, c->bw));
	if (_mm_movemask_epi8(inv) == 0xffff) return 0;
	inv = _mm_or_si128(inv, _mm_or_si128(_mm_cmpeq_epi32(dr, zero), _mm_cmpgt_epi32(dq, c->mdxy)));
	hj = _mm_loadu_si128((const __m128i*)&hs[j-3]);
	dg = _mm_or_si128(_mm_and_si128(lt, dr), _mm_andnot_si128(lt, dq));
	qs = _mm_and_si128(hj, ff);
	qlt = _mm_cmpgt_epi32(dg, qs);
	sc = _mm_or_si128(_mm_and_si128(qlt, qs), _mm_andnot_si128(qlt, dg));
	cond = _mm_or_si128(_mm_xor_si128(_mm_cmpeq_epi32(dd, zero), ones), qlt);
	lin = _mm_add_ps(_mm_mul_ps(c->gap, _mm_cvtepi32_ps(dd)), _mm_mul_ps(c->skip, _mm_cvtepi32_ps(dg)));
	zi = _mm_castps_si128(_mm_cvtepi32_ps(_mm_add_epi32(dd, one))); // mg_log2(dd + 1)
	l2 = _mm_cvtepi32_ps(_mm_sub_epi32(_mm_and_si128(_mm_srli_epi32(zi, 23), ff), _mm_set1_epi32(128)));
	z = _mm_castsi128_ps(_mm_add_epi32(_mm_and_si128(zi, _mm_set1_epi32(~(255 << 23))), _mm_set1_epi32(127 << 23)));
	l2 = _mm_add_ps(l2, _mm_sub_ps(_mm_mul_ps(_mm_add_ps(_mm_mul_ps(_mm_set1_ps(-0.34484843f), z), _mm_set1_ps(2.02466578f)), z), _mm_set1_ps(0.67487759f)));
	log_pen = _mm_and_ps(_mm_castsi128_ps(_mm_cmpgt_epi32(dd, zero)), l2);
	use_a = _mm_and_si128(c->cdna, _mm_cmpgt_epi32(dr, dq));
	pen = _mm_or_si128(_mm_and_si128(use_a, _mm_cvttps_epi32(_mm_min_ps(lin, log_pen))),
					   _mm_andnot_si128(use_a, _mm_cvttps_epi32(_mm_add_ps(lin, _mm_mul_ps(_mm_set1_ps(.5f), log_pen)))));
	sc = _mm_sub_epi32(sc, _mm_and_si128(cond, pen));
	r = _mm_or_si128(_mm_and_si128(inv, _mm_set1_epi32(INT32_MIN)), _mm_andnot_si128(inv, sc));
	*out = r;
	return lc_rev4[~_mm_movemask_ps(_mm_castsi128_ps(_mm_cmpeq_epi32(r, _mm_set1_epi32(INT32_MIN)))) & 0xf];
}
#else // MG_LC_AVX2
#define LC_R2(n) n, n + 2*64, n + 1*64, n + 3*64
#define LC_R4(n) LC_R2(n), LC_R2(n + 2*16), LC_R2(n + 1*16), LC_R2(n + 3*16)
#define LC_R6(n) LC_R4(n), LC_R4(n + 2*4), LC_R4(n + 1*4), LC_R4(n + 3*4)
static const uint8_t lc_rev8[256] = { LC_R6(0), LC_R6(2), LC_R6(1), LC_R6(3) }; // bit l -> bit 7-l

/* comput_sc4() on 8 lanes, the predecessors j-7, ..., j; the same ops lane by lane (no FMA: built without -mfma and
 * with -ffp-contract=off). Returns bit 7-l for lane l. */
static inline uint32_t comput_sc8(const int32_t *xs, const int32_t *ys, const int32_t *hs, int64_t j, const lc_sse_t *c, __m256i *out)
{
	const __m256i one = _mm256_set1_epi32(1), zero = _mm256_setzero_si256(), ones = _mm256_set1_epi32(-1), ff = _mm256_set1_epi32(0xff);
	__m256i xj = _mm256_loadu_si256((const __m256i*)&xs[j-7]), yj = _mm256_loadu_si256((const __m256i*)&ys[j-7]), hj = _mm256_loadu_si256((const __m256i*)&hs[j-7]);
	__m256i dq, dr, same, inv, s_inv, drz, lt, dd, dg, qs, qlt, sc, cond, nsame, bonus, use_a, pen, zi, r;
	__m256 lin, z, l2, log_pen;
	dq = _mm256_sub_epi32(c->yi, yj);
	dr = _mm256_sub_epi32(c->xi, xj);
	same = _mm256_cmpeq_epi32(c->sidi, _mm256_and_si256(_mm256_srli_epi32(hj, MM_SEED_SEG_SHIFT - 32), ff));
	inv = _mm256_or_si256(_mm256_cmpgt_epi32(one, dq), _mm256_cmpgt_epi32(dq, c->mdx));
	drz = _mm256_cmpeq_epi32(dr, zero);
	lt = _mm256_cmpgt_epi32(dq, dr);
	dd = _mm256_sub_epi32(_mm256_xor_si256(_mm256_sub_epi32(dr, dq), lt), lt); // dr < dq? dq - dr : dr - dq
	s_inv = _mm256_or_si256(_mm256_or_si256(drz, _mm256_cmpgt_epi32(dq, c->mdy)), _mm256_cmpgt_epi32(dd, c->bw));
	s_inv = _mm256_or_si256(s_inv, _mm256_and_si256(c->seg_chk, _mm256_cmpgt_epi32(dr, c->mdy)));
	inv = _mm256_or_si256(inv, _mm256_and_si256(same, s_inv));
	if (_mm256_movemask_ps(_mm256_castsi256_ps(inv)) == 0xff) return 0;
	dg = _mm256_or_si256(_mm256_and_si256(lt, dr), _mm256_andnot_si256(lt, dq));
	qs = _mm256_and_si256(hj, ff);
	qlt = _mm256_cmpgt_epi32(dg, qs);
	sc = _mm256_or_si256(_mm256_and_si256(qlt, qs), _mm256_andnot_si256(qlt, dg));
	cond = _mm256_or_si256(_mm256_xor_si256(_mm256_cmpeq_epi32(dd, zero), ones), qlt);
	lin = _mm256_add_ps(_mm256_mul_ps(c->gap, _mm256_cvtepi32_ps(dd)), _mm256_mul_ps(c->skip, _mm256_cvtepi32_ps(dg)));
	zi = _mm256_castps_si256(_mm256_cvtepi32_ps(_mm256_add_epi32(dd, one))); // mg_log2(dd + 1)
	l2 = _mm256_cvtepi32_ps(_mm256_sub_epi32(_mm256_and_si256(_mm256_srli_epi32(zi, 23), ff), _mm256_set1_epi32(128)));
	z = _mm256_castsi256_ps(_mm256_add_epi32(_mm256_and_si256(zi, _mm256_set1_epi32(~(255 << 23))), _mm256_set1_epi32(127 << 23)));
	l2 = _mm256_add_ps(l2, _mm256_sub_ps(_mm256_mul_ps(_mm256_add_ps(_mm256_mul_ps(_mm256_set1_ps(-0.34484843f), z), _mm256_set1_ps(2.02466578f)), z), _mm256_set1_ps(0.67487759f)));
	log_pen = _mm256_and_ps(_mm256_castsi256_ps(_mm256_cmpgt_epi32(dd, zero)), l2);
	nsame = _mm256_xor_si256(same, ones);
	bonus = _mm256_and_si256(nsame, drz);
	use_a = _mm256_and_si256(_mm256_or_si256(c->cdna, nsame), _mm256_or_si256(_mm256_cmpgt_epi32(dr, dq), nsame));
	pen = _mm256_or_si256(_mm256_and_si256(use_a, _mm256_cvttps_epi32(_mm256_min_ps(lin, log_pen))),
						  _mm256_andnot_si256(use_a, _mm256_cvttps_epi32(_mm256_add_ps(lin, _mm256_mul_ps(_mm256_set1_ps(.5f), log_pen)))));
	pen = _mm256_or_si256(bonus, _mm256_andnot_si256(bonus, pen)); // -1: sc - (-1) is the ++sc bonus
	sc = _mm256_sub_epi32(sc, _mm256_and_si256(cond, pen));
	r = _mm256_or_si256(_mm256_and_si256(inv, _mm256_set1_epi32(INT32_MIN)), _mm256_andnot_si256(inv, sc));
	*out = r;
	return lc_rev8[~_mm256_movemask_ps(_mm256_castsi256_ps(_mm256_cmpeq_epi32(r, _mm256_set1_epi32(INT32_MIN)))) & 0xff];
}

/* comput_sc8_1()'s first test alone: nonzero unless every lane of j-7, ..., j has |dr - dq| > bw */
static inline uint32_t lc_cand8(const int32_t *ds, int64_t j, const lc_sse_t *c)
{
	__m256i e = _mm256_sub_epi32(c->di, _mm256_loadu_si256((const __m256i*)&ds[j-7]));
	return _mm256_movemask_ps(_mm256_castsi256_ps(_mm256_cmpgt_epi32(_mm256_add_epi32(e, c->bwb), c->bw2b))) ^ 0xff;
}

/* comput_sc4_1() on 8 lanes */
static inline uint32_t comput_sc8_1(const int32_t *xs, const int32_t *ys, const int32_t *hs, const int32_t *ds, int64_t j, const lc_sse_t *c, __m256i *out)
{
	const __m256i one = _mm256_set1_epi32(1), zero = _mm256_setzero_si256(), ones = _mm256_set1_epi32(-1), ff = _mm256_set1_epi32(0xff);
	__m256i xj, yj, hj, dq, dr, inv, lt, dd, dg, qs, qlt, sc, cond, use_a, pen, zi, r;
	__m256 lin, z, l2, log_pen;
	if (ds) { // (dr - dq) + bw > 2 * bw, unsigned
		__m256i e = _mm256_sub_epi32(c->di, _mm256_loadu_si256((const __m256i*)&ds[j-7]));
		if (_mm256_movemask_ps(_mm256_castsi256_ps(_mm256_cmpgt_epi32(_mm256_add_epi32(e, c->bwb), c->bw2b))) == 0xff) return 0;
	}
	xj = _mm256_loadu_si256((const __m256i*)&xs[j-7]), yj = _mm256_loadu_si256((const __m256i*)&ys[j-7]);
	dq = _mm256_sub_epi32(c->yi, yj);
	dr = _mm256_sub_epi32(c->xi, xj);
	lt = _mm256_cmpgt_epi32(dq, dr);
	dd = _mm256_sub_epi32(_mm256_xor_si256(_mm256_sub_epi32(dr, dq), lt), lt); // dr < dq? dq - dr : dr - dq
	inv = _mm256_or_si256(_mm256_cmpgt_epi32(one, dq), _mm256_cmpgt_epi32(dd, c->bw));
	if (_mm256_movemask_ps(_mm256_castsi256_ps(inv)) == 0xff) return 0;
	inv = _mm256_or_si256(inv, _mm256_or_si256(_mm256_cmpeq_epi32(dr, zero), _mm256_cmpgt_epi32(dq, c->mdxy)));
	hj = _mm256_loadu_si256((const __m256i*)&hs[j-7]);
	dg = _mm256_or_si256(_mm256_and_si256(lt, dr), _mm256_andnot_si256(lt, dq));
	qs = _mm256_and_si256(hj, ff);
	qlt = _mm256_cmpgt_epi32(dg, qs);
	sc = _mm256_or_si256(_mm256_and_si256(qlt, qs), _mm256_andnot_si256(qlt, dg));
	cond = _mm256_or_si256(_mm256_xor_si256(_mm256_cmpeq_epi32(dd, zero), ones), qlt);
	lin = _mm256_add_ps(_mm256_mul_ps(c->gap, _mm256_cvtepi32_ps(dd)), _mm256_mul_ps(c->skip, _mm256_cvtepi32_ps(dg)));
	zi = _mm256_castps_si256(_mm256_cvtepi32_ps(_mm256_add_epi32(dd, one))); // mg_log2(dd + 1)
	l2 = _mm256_cvtepi32_ps(_mm256_sub_epi32(_mm256_and_si256(_mm256_srli_epi32(zi, 23), ff), _mm256_set1_epi32(128)));
	z = _mm256_castsi256_ps(_mm256_add_epi32(_mm256_and_si256(zi, _mm256_set1_epi32(~(255 << 23))), _mm256_set1_epi32(127 << 23)));
	l2 = _mm256_add_ps(l2, _mm256_sub_ps(_mm256_mul_ps(_mm256_add_ps(_mm256_mul_ps(_mm256_set1_ps(-0.34484843f), z), _mm256_set1_ps(2.02466578f)), z), _mm256_set1_ps(0.67487759f)));
	log_pen = _mm256_and_ps(_mm256_castsi256_ps(_mm256_cmpgt_epi32(dd, zero)), l2);
	use_a = _mm256_and_si256(c->cdna, _mm256_cmpgt_epi32(dr, dq));
	pen = _mm256_or_si256(_mm256_and_si256(use_a, _mm256_cvttps_epi32(_mm256_min_ps(lin, log_pen))),
						  _mm256_andnot_si256(use_a, _mm256_cvttps_epi32(_mm256_add_ps(lin, _mm256_mul_ps(_mm256_set1_ps(.5f), log_pen)))));
	sc = _mm256_sub_epi32(sc, _mm256_and_si256(cond, pen));
	r = _mm256_or_si256(_mm256_and_si256(inv, _mm256_set1_epi32(INT32_MIN)), _mm256_andnot_si256(inv, sc));
	*out = r;
	return lc_rev8[~_mm256_movemask_ps(_mm256_castsi256_ps(_mm256_cmpeq_epi32(r, _mm256_set1_epi32(INT32_MIN)))) & 0xff];
}
#endif // MG_LC_AVX2
#endif // MG_LCHAIN_SSE2

/* Input:
 *   a[].x: rev<<63 | tid<<32 | tpos
 *   a[].y: flags<<40 | q_span<<32 | q_pos
 * Output:
 *   n_u: #chains
 *   u[]: score<<32 | #anchors (sum of lower 32 bits of u[] is the returned length of a[])
 * input a[] is deallocated on return
 */
#ifdef MG_LC_AVX2
mm128_t *mg_lchain_dp_avx2(int max_dist_x, int max_dist_y, int bw, int max_skip, int max_iter, int min_cnt, int min_sc, float chn_pen_gap, float chn_pen_skip,
						   int is_cdna, int n_seg, int64_t n, mm128_t *a, int *n_u_, uint64_t **_u, void *km)
#else
#if defined(MG_LCHAIN_HAVE_AVX2) && defined(MG_LCHAIN_SSE2)
#define MG_LC_DISPATCH
mm128_t *mg_lchain_dp_avx2(int max_dist_x, int max_dist_y, int bw, int max_skip, int max_iter, int min_cnt, int min_sc, float chn_pen_gap, float chn_pen_skip,
						   int is_cdna, int n_seg, int64_t n, mm128_t *a, int *n_u_, uint64_t **_u, void *km);
int ksw_cpu_avx2(void); // ksw2_dispatch.c
static int mg_lc_avx2;
static pthread_once_t mg_lc_avx2_once = PTHREAD_ONCE_INIT;
static void mg_lc_avx2_init(void) { mg_lc_avx2 = ksw_cpu_avx2(); } // by pthread_once: one chaining thread sets it, the others wait
#endif
mm128_t *mg_lchain_dp(int max_dist_x, int max_dist_y, int bw, int max_skip, int max_iter, int min_cnt, int min_sc, float chn_pen_gap, float chn_pen_skip,
					  int is_cdna, int n_seg, int64_t n, mm128_t *a, int *n_u_, uint64_t **_u, void *km)
#endif
{ // TODO: make sure this works when n has more than 32 bits
	int32_t *f, *t, *v, n_u, n_v, mmax_f = 0, max_drop = bw;
#ifdef MG_LCHAIN_SSE2
	int32_t *xs, *ys, *hs, *ds = 0;
	int one_seg = 0;
	lc_sse_t c;
#endif
	int64_t *p, i, j, max_ii, st = 0;
	uint64_t *u;

#ifdef MG_LC_DISPATCH
	pthread_once(&mg_lc_avx2_once, mg_lc_avx2_init);
	if (mg_lc_avx2)
		return mg_lchain_dp_avx2(max_dist_x, max_dist_y, bw, max_skip, max_iter, min_cnt, min_sc, chn_pen_gap, chn_pen_skip, is_cdna, n_seg, n, a, n_u_, _u, km);
#endif
	if (_u) *_u = 0, *n_u_ = 0;
	if (n == 0 || a == 0) {
		kfree(km, a);
		return 0;
	}
	if (max_dist_x < bw) max_dist_x = bw;
	if (max_dist_y < bw && !is_cdna) max_dist_y = bw;
	if (is_cdna) max_drop = INT32_MAX;
	p = Kmalloc(km, int64_t, n);
	f = Kmalloc(km, int32_t, n);
	v = Kmalloc(km, int32_t, n);
	t = Kcalloc(km, int32_t, n);
#ifdef MG_LCHAIN_SSE2
	xs = Kmalloc(km, int32_t, (n + MG_LC_W - 1) * 3) + MG_LC_W - 1, ys = xs + n + MG_LC_W - 1, hs = ys + n + MG_LC_W - 1; // MG_LC_W - 1 entries before each, read by the last block
	for (j = 1; j < MG_LC_W; ++j) xs[-j] = ys[-j] = hs[-j] = 0;
	if (n_seg == 1 && bw >= 0) {
		uint64_t seg = 0, uns = 0;
		if (bw < max_dist_x) { // otherwise |dr - dq| > bw rejects few groups
			ds = Kmalloc(km, int32_t, n + MG_LC_W - 1) + MG_LC_W - 1;
			for (j = 1; j < MG_LC_W; ++j) ds[-j] = 0;
		}
		if (ds) {
			for (i = 0; i < n; ++i) {
				xs[i] = (int32_t)a[i].x, ys[i] = (int32_t)a[i].y, hs[i] = (int32_t)(a[i].y >> 32), ds[i] = (int32_t)((uint32_t)xs[i] - (uint32_t)ys[i]);
				seg |= (a[i].y ^ a[0].y) & MM_SEED_SEG_MASK, uns |= a[i].x < a[i > 0? i - 1 : 0].x;
			}
		} else {
			for (i = 0; i < n; ++i) {
				xs[i] = (int32_t)a[i].x, ys[i] = (int32_t)a[i].y, hs[i] = (int32_t)(a[i].y >> 32);
				seg |= (a[i].y ^ a[0].y) & MM_SEED_SEG_MASK;
			}
		}
		one_seg = (seg == 0);
		if (ds && (uns || !one_seg)) kfree(km, ds - (MG_LC_W - 1)), ds = 0;
	} else {
		for (i = 0; i < n; ++i)
			xs[i] = (int32_t)a[i].x, ys[i] = (int32_t)a[i].y, hs[i] = (int32_t)(a[i].y >> 32);
	}
	c.mdx = lc_set1_epi32(max_dist_x), c.mdy = lc_set1_epi32(max_dist_y), c.bw = lc_set1_epi32(bw);
	c.seg_chk = lc_set1_epi32(n_seg > 1 && !is_cdna? -1 : 0), c.cdna = lc_set1_epi32(is_cdna? -1 : 0);
	c.gap = lc_set1_ps(chn_pen_gap), c.skip = lc_set1_ps(chn_pen_skip);
	c.mdxy = lc_set1_epi32(max_dist_x < max_dist_y? max_dist_x : max_dist_y);
	c.bwb = lc_set1_epi32((int32_t)((uint32_t)bw ^ 0x80000000u)), c.bw2b = lc_set1_epi32((int32_t)(((uint32_t)bw << 1) ^ 0x80000000u));
	c.di = lc_set1_epi32(0); // set per anchor and read only when ds is
#endif

	// fill the score and backtrack arrays
	for (i = 0, max_ii = -1; i < n; ++i) {
		int64_t max_j = -1, end_j;
		int32_t max_f = a[i].y>>32&0xff, n_skip = 0;
		while (st < i && (a[i].x>>32 != a[st].x>>32 || a[i].x > a[st].x + max_dist_x)) ++st;
		if (i - st > max_iter) st = i - max_iter;
		j = i - 1;
#ifdef MG_LCHAIN_SSE2
		if (i - st >= MG_LC_MIN) { // score MG_LC_BLK predecessors at a time, then scan those whose score is not INT32_MIN as the loop below does
			c.xi = lc_set1_epi32(xs[i]), c.yi = lc_set1_epi32(ys[i]);
			if (ds) c.di = lc_set1_epi32(ds[i]);
			c.sidi = lc_set1_epi32((int32_t)((a[i].y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT));
			while (j >= st) {
				int32_t buf[MG_LC_BLK];
				int64_t j0 = j;
				int k, nb = j - st + 1 < MG_LC_BLK? (int)(j - st + 1) : MG_LC_BLK;
				uint32_t msk = 0;
#ifdef MG_LC_AVX2
				if (ds) { // test every group first, then score those with a lane left: fewer mispredicted branches
					uint32_t gm = 0;
					for (k = 0; k < nb; k += MG_LC_W)
						gm |= (uint32_t)(lc_cand8(ds, j0 - k, &c) != 0) << (k >> 3);
					while (gm) {
						lc_vi_t s4;
						uint32_t m4;
						k = __builtin_ctz(gm) << 3, gm &= gm - 1;
						m4 = comput_sc8_1(xs, ys, hs, 0, j0 - k, &c, &s4);
						if (m4) lc_storeu(&buf[k], s4), msk |= m4 << k;
					}
				} else
#endif
				for (k = 0; k < nb; k += MG_LC_W) {
					lc_vi_t s4;
#ifdef MG_LC_AVX2
					uint32_t m4 = one_seg? comput_sc8_1(xs, ys, hs, ds, j0 - k, &c, &s4) : comput_sc8(xs, ys, hs, j0 - k, &c, &s4);
#else
					uint32_t m4 = one_seg? comput_sc4_1(xs, ys, hs, ds, j0 - k, &c, &s4) : comput_sc4(xs, ys, hs, j0 - k, &c, &s4);
#endif
					if (m4) lc_storeu(&buf[k], s4), msk |= m4 << k;
				}
				msk &= (uint32_t)(((uint64_t)1 << nb) - 1); // drop lanes before st
				while (msk) { // bit k is a[j0-k], its score in buf[k^(MG_LC_W-1)]
					int32_t sc;
					k = __builtin_ctz(msk), msk &= msk - 1;
					j = j0 - k;
					sc = buf[k ^ (MG_LC_W - 1)] + f[j];
					if (sc > max_f) {
						max_f = sc, max_j = j;
						if (n_skip > 0) --n_skip;
					} else if (t[j] == (int32_t)i) {
						if (++n_skip > max_skip)
							goto end_scan;
					}
					if (p[j] >= 0) t[p[j]] = i;
				}
				j = j0 - nb;
			}
		}
#endif
		for (; j >= st; --j) {
			int32_t sc;
			sc = comput_sc(&a[i], &a[j], max_dist_x, max_dist_y, bw, chn_pen_gap, chn_pen_skip, is_cdna, n_seg);
			if (sc == INT32_MIN) continue;
			sc += f[j];
			if (sc > max_f) {
				max_f = sc, max_j = j;
				if (n_skip > 0) --n_skip;
			} else if (t[j] == (int32_t)i) {
				if (++n_skip > max_skip)
					break;
			}
			if (p[j] >= 0) t[p[j]] = i;
		}
#ifdef MG_LCHAIN_SSE2
end_scan:
#endif
		end_j = j;
		if (max_ii < 0 || a[i].x - a[max_ii].x > (int64_t)max_dist_x) {
			int32_t max = INT32_MIN;
			max_ii = -1;
			for (j = i - 1; j >= st; --j)
				if (max < f[j]) max = f[j], max_ii = j;
		}
		if (max_ii >= 0 && max_ii < end_j) {
			int32_t tmp;
			tmp = comput_sc(&a[i], &a[max_ii], max_dist_x, max_dist_y, bw, chn_pen_gap, chn_pen_skip, is_cdna, n_seg);
			if (tmp != INT32_MIN && max_f < tmp + f[max_ii])
				max_f = tmp + f[max_ii], max_j = max_ii;
		}
		f[i] = max_f, p[i] = max_j;
		v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f; // v[] keeps the peak score up to i; f[] is the score ending at i, not always the peak
		if (max_ii < 0 || (a[i].x - a[max_ii].x <= (int64_t)max_dist_x && f[max_ii] < f[i]))
			max_ii = i;
		if (mmax_f < max_f) mmax_f = max_f;
		//fprintf(stderr, "X1\t%ld\t%ld:%d\t%ld\t%ld:%d\t%ld\t%ld\n", (long)i, (long)(a[i].x>>32), (int32_t)a[i].x, (long)max_j, max_j<0?-1L:(long)(a[max_j].x>>32), max_j<0?-1:(int32_t)a[max_j].x, (long)max_f, (long)v[i]);
	}

	u = mg_chain_backtrack(km, n, f, p, v, t, min_cnt, min_sc, max_drop, &n_u, &n_v);
	*n_u_ = n_u, *_u = u; // NB: note that u[] may not be sorted by score here
	kfree(km, p); kfree(km, f); kfree(km, t);
#ifdef MG_LCHAIN_SSE2
	kfree(km, xs - (MG_LC_W - 1));
	if (ds) kfree(km, ds - (MG_LC_W - 1));
#endif
	if (n_u == 0) {
		kfree(km, a); kfree(km, v);
		return 0;
	}
	return compact_a(km, n_u, u, n_v, v, a);
}
#endif // !MG_LCHAIN_AVX2_UNIT || MG_LC_AVX2

#ifndef MG_LCHAIN_AVX2_UNIT
typedef struct lc_elem_s {
	int32_t y;
	int64_t i;
	double pri;
	KRMQ_HEAD(struct lc_elem_s) head;
} lc_elem_t;

#define lc_elem_cmp(a, b) ((a)->y < (b)->y? -1 : (a)->y > (b)->y? 1 : ((a)->i > (b)->i) - ((a)->i < (b)->i))
#define lc_elem_lt2(a, b) ((a)->pri < (b)->pri)
KRMQ_INIT(lc_elem, lc_elem_t, head, lc_elem_cmp, lc_elem_lt2)

KALLOC_POOL_INIT(rmq, lc_elem_t)

/*
 * Ordered set of (y, i) keys for mg_lchain_rmq's inner window. The window is only read in key order (predecessor
 * and descending walk) and by size, so a sorted array cut into chunks replaces a krmq tree. Each chunk keeps its
 * keys in k[b..e) and moves the shorter side on insert or erase; adjacent chunks always hold more than
 * LC_ISET_CK/2 keys between them, which bounds the number of chunks.
 */
#ifndef LC_ISET_CK
#define LC_ISET_CK 64
#endif

typedef struct { int64_t i; int32_t y; } lc_ikey_t;

typedef struct lc_ick_s {
	int32_t b, e;
	struct lc_ick_s *next; // free list
	lc_ikey_t k[LC_ISET_CK];
} lc_ick_t;

typedef struct {
	int32_t n, m; // chunks in use; capacity of c[] and f[]
	unsigned size; // unsigned, as krmq's size, for the cap test
	lc_ick_t **c; // chunks in key order
	lc_ikey_t *f; // f[x] is the first key of c[x]
	lc_ick_t *fl;
	void *km;
} lc_iset_t;

#define lc_ikey_lt(a, b) ((a)->y < (b)->y || ((a)->y == (b)->y && (a)->i < (b)->i))

static lc_ick_t *lc_iset_new_ck(lc_iset_t *s, int32_t x) // an empty chunk at position x
{
	lc_ick_t *c;
	if (s->fl) c = s->fl, s->fl = c->next;
	else c = (lc_ick_t*)kmalloc(s->km, sizeof(lc_ick_t));
	c->b = c->e = 0;
	if (s->n == s->m) {
		s->m = s->m? s->m<<1 : 16;
		s->c = (lc_ick_t**)krealloc(s->km, s->c, s->m * sizeof(lc_ick_t*));
		s->f = (lc_ikey_t*)krealloc(s->km, s->f, s->m * sizeof(lc_ikey_t));
	}
	if (x < s->n) {
		memmove(&s->c[x+1], &s->c[x], (s->n - x) * sizeof(lc_ick_t*));
		memmove(&s->f[x+1], &s->f[x], (s->n - x) * sizeof(lc_ikey_t));
	}
	s->c[x] = c, ++s->n;
	return c;
}

static void lc_iset_del_ck(lc_iset_t *s, int32_t x)
{
	lc_ick_t *c = s->c[x];
	c->next = s->fl, s->fl = c;
	if (x + 1 < s->n) {
		memmove(&s->c[x], &s->c[x+1], (s->n - x - 1) * sizeof(lc_ick_t*));
		memmove(&s->f[x], &s->f[x+1], (s->n - x - 1) * sizeof(lc_ikey_t));
	}
	--s->n;
}

static inline int32_t lc_iset_find_ck(const lc_iset_t *s, const lc_ikey_t *key) // the last chunk whose first key <= key; -1 if none
{
	int32_t lo = 0, hi = s->n;
	while (lo < hi) {
		int32_t mid = (lo + hi) >> 1;
		if (lc_ikey_lt(key, &s->f[mid])) hi = mid;
		else lo = mid + 1;
	}
	return lo - 1;
}

static inline int32_t lc_ick_upper(const lc_ick_t *c, const lc_ikey_t *key) // the first position in [b,e) whose key > key
{
	int32_t lo = c->b, hi = c->e;
	while (lo < hi) {
		int32_t mid = (lo + hi) >> 1;
		if (lc_ikey_lt(key, &c->k[mid])) hi = mid;
		else lo = mid + 1;
	}
	return lo;
}

static void lc_iset_insert(lc_iset_t *s, int32_t y, int64_t i)
{
	lc_ikey_t key;
	lc_ick_t *c;
	int32_t x, p;
	key.i = i, key.y = y;
	if (s->n == 0) {
		c = lc_iset_new_ck(s, 0);
		c->k[0] = key, c->e = 1, s->f[0] = key, s->size = 1;
		return;
	}
	x = lc_iset_find_ck(s, &key);
	if (x < 0) x = 0;
	c = s->c[x];
	p = lc_ikey_lt(&c->k[c->e - 1], &key)? c->e : lc_ick_upper(c, &key);
	if (p > c->b && !lc_ikey_lt(&c->k[p - 1], &key)) return; // already present
	if (p == c->e && c->e == LC_ISET_CK) { // append to a chunk whose right side is full
		if (c->b >= LC_ISET_CK>>1) { // at most half full: shift its keys to the start
			int32_t n = c->e - c->b;
			memmove(c->k, &c->k[c->b], n * sizeof(lc_ikey_t));
			c->b = 0, c->e = p = n;
		} else if (x == s->n - 1) { // the last chunk: start a new one
			c = lc_iset_new_ck(s, x + 1);
			c->k[0] = key, c->e = 1, s->f[x + 1] = key, ++s->size;
			return;
		}
	}
	if (c->e - c->b == LC_ISET_CK) { // full: split in half
		lc_ick_t *d = lc_iset_new_ck(s, x + 1);
		int32_t h = LC_ISET_CK>>1;
		c = s->c[x]; // unchanged; s->c may have moved
		memcpy(d->k, &c->k[h], h * sizeof(lc_ikey_t));
		d->b = 0, d->e = h, c->e = h;
		s->f[x + 1] = d->k[0];
		if (p > h) c = d, p -= h, ++x;
	}
	if (p == c->e && c->e < LC_ISET_CK) {
		c->k[c->e++] = key;
	} else if (p == c->b && c->b > 0) {
		c->k[--c->b] = key;
	} else if (c->e < LC_ISET_CK && (c->b == 0 || c->e - p <= p - c->b)) {
		memmove(&c->k[p + 1], &c->k[p], (c->e - p) * sizeof(lc_ikey_t));
		c->k[p] = key, ++c->e;
	} else {
		memmove(&c->k[c->b - 1], &c->k[c->b], (p - c->b) * sizeof(lc_ikey_t));
		c->k[p - 1] = key, --c->b;
	}
	if (x == 0) s->f[0] = c->k[c->b];
	++s->size;
}

static void lc_iset_merge(lc_iset_t *s, int32_t x) // move c[x+1]'s keys to c[x]
{
	lc_ick_t *c = s->c[x], *d = s->c[x + 1];
	int32_t n = c->e - c->b;
	if (c->b > 0) memmove(c->k, &c->k[c->b], n * sizeof(lc_ikey_t));
	memcpy(&c->k[n], &d->k[d->b], (d->e - d->b) * sizeof(lc_ikey_t));
	c->b = 0, c->e = n + (d->e - d->b);
	lc_iset_del_ck(s, x + 1);
}

static void lc_iset_erase(lc_iset_t *s, int32_t y, int64_t i) // erase the key if present
{
	lc_ikey_t key;
	lc_ick_t *c;
	int32_t x, p;
	if (s->n == 0) return;
	key.i = i, key.y = y;
	c = s->c[0];
	if (c->k[c->b].i == i && c->k[c->b].y == y) { // the smallest key
		x = 0, p = c->b;
	} else {
		x = lc_iset_find_ck(s, &key);
		if (x < 0) return;
		c = s->c[x];
		p = lc_ick_upper(c, &key) - 1; // >= c->b, as c's first key <= key
		if (lc_ikey_lt(&c->k[p], &key)) return; // absent
	}
	--s->size;
	if (c->e - c->b == 1) {
		lc_iset_del_ck(s, x);
		return;
	}
	if (p - c->b < c->e - 1 - p) {
		memmove(&c->k[c->b + 1], &c->k[c->b], (p - c->b) * sizeof(lc_ikey_t));
		++c->b;
	} else {
		memmove(&c->k[p], &c->k[p + 1], (c->e - 1 - p) * sizeof(lc_ikey_t));
		--c->e;
	}
	s->f[x] = c->k[c->b];
	if (x + 1 < s->n && (c->e - c->b) + (s->c[x+1]->e - s->c[x+1]->b) <= LC_ISET_CK>>1)
		lc_iset_merge(s, x);
	else if (x > 0 && (c->e - c->b) + (s->c[x-1]->e - s->c[x-1]->b) <= LC_ISET_CK>>1)
		lc_iset_merge(s, x - 1);
}

static inline int lc_iset_pred(const lc_iset_t *s, int32_t y, int32_t *x_, int32_t *p_) // the last key with key.y < y
{
	int32_t lo = 0, hi = s->n, x;
	const lc_ick_t *c;
	while (lo < hi) {
		int32_t mid = (lo + hi) >> 1;
		if (s->f[mid].y < y) lo = mid + 1;
		else hi = mid;
	}
	if ((x = lo - 1) < 0) return 0;
	c = s->c[x];
	if (c->k[c->e - 1].y < y) {
		*p_ = c->e - 1;
	} else {
		lo = c->b, hi = c->e;
		while (lo < hi) {
			int32_t mid = (lo + hi) >> 1;
			if (c->k[mid].y < y) lo = mid + 1;
			else hi = mid;
		}
		*p_ = lo - 1;
	}
	*x_ = x;
	return 1;
}

static inline int32_t comput_sc_simple(const mm128_t *ai, const mm128_t *aj, float chn_pen_gap, float chn_pen_skip, int32_t *exact, int32_t *width)
{
	int32_t dq = (int32_t)ai->y - (int32_t)aj->y, dr, dd, dg, q_span, sc;
	dr = (int32_t)(ai->x - aj->x);
	*width = dd = dr > dq? dr - dq : dq - dr;
	dg = dr < dq? dr : dq;
	q_span = aj->y>>32&0xff;
	sc = q_span < dg? q_span : dg;
	if (exact) *exact = (dd == 0 && dg <= q_span);
	if (dd || dq > q_span) {
		float lin_pen, log_pen;
		lin_pen = chn_pen_gap * (float)dd + chn_pen_skip * (float)dg;
		log_pen = dd >= 1? mg_log2(dd + 1) : 0.0f; // mg_log2() only works for dd>=2
		sc -= (int)(lin_pen + .5f * log_pen);
	}
	return sc;
}

// The main tree of mg_lchain_rmq as an array min-tree: one leaf per distinct (int32)a[].y of the call, in y order.
// A node keeps its subtree's smallest pri, an anchor holding it (-1 when empty) and whether it is held twice.
typedef struct { double mn; int32_t am, tie; } lc_mnode_t;

static inline int lc_popcnt64(uint64_t x)
{
	x = x - ((x >> 1) & 0x5555555555555555ULL);
	x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
	x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0fULL;
	return (int)((x * 0x0101010101010101ULL) >> 56);
}

static inline int32_t lc_cnt_lt(int64_t y, int32_t ymin, int32_t ymax, int32_t ny, const uint64_t *bits, const int32_t *pre) // distinct y values below y
{
	uint64_t d;
	if (y <= ymin) return 0;
	if (y > ymax) return ny;
	d = (uint64_t)(y - ymin);
	return pre[d>>6] + lc_popcnt64(bits[d>>6] & ((1ULL<<(d&63)) - 1));
}

static inline void lc_mt_fix(lc_mnode_t *T, int64_t k) // recompute the ancestors of node k until one is unchanged
{
	for (k >>= 1; k > 0; k >>= 1) {
		const lc_mnode_t *l = &T[k<<1], *r = l + 1;
		lc_mnode_t z;
		if (r->mn < l->mn) z = *r;
		else {
			z = *l;
			if (l->mn == r->mn && l->am >= 0) z.tie = 1;
		}
		if (z.mn == T[k].mn && z.am == T[k].am && z.tie == T[k].tie) break;
		T[k] = z;
	}
}

// scores the rmq answer j for anchor i, then walks the inner window for a better predecessor; used by both fills
static inline void lc_rmq_score(int64_t i, int64_t j, const mm128_t *a, const int32_t *f, const int64_t *p, int32_t *t, int bw, int max_dist_inner,
								int max_chn_skip, float chn_pen_gap, float chn_pen_skip, const lc_iset_t *si, int32_t *max_f_, int64_t *max_j_)
{
	int32_t sc, exact, width, n_skip = 0, max_f = *max_f_;
	int64_t max_j = *max_j_;
	sc = f[j] + comput_sc_simple(&a[i], &a[j], chn_pen_gap, chn_pen_skip, &exact, &width);
	if (width <= bw && sc > max_f) max_f = sc, max_j = j;
	if (!exact && si->n > 0 && (int32_t)a[i].y > 0) {
		int32_t ck, pos;
		if (lc_iset_pred(si, (int32_t)a[i].y, &ck, &pos)) { // the last key with y < a[i].y; then descending
			const lc_ick_t *c = si->c[ck];
			for (;;) {
				const lc_ikey_t *q = &c->k[pos];
				if (q->y < (int32_t)a[i].y - max_dist_inner) break;
				j = q->i;
				sc = f[j] + comput_sc_simple(&a[i], &a[j], chn_pen_gap, chn_pen_skip, 0, &width);
				if (width <= bw) {
					if (sc > max_f) {
						max_f = sc, max_j = j;
						if (n_skip > 0) --n_skip;
					} else if (t[j] == (int32_t)i) {
						if (++n_skip > max_chn_skip)
							break;
					}
					if (p[j] >= 0) t[p[j]] = i;
				}
				if (--pos < c->b) {
					if (--ck < 0) break;
					c = si->c[ck], pos = c->e - 1;
				}
			}
		}
	}
	*max_f_ = max_f, *max_j_ = max_j;
}

// fills f[], p[] and v[] with the main tree as an array min-tree; returns 0, leaving t[] dirty, when an erase
// comes before its insert or the call is out of its range. A query whose minimum pri is held twice is answered by
// krmq, whose pick among equal pri depends on its shape: a krmq main tree is brought to the same state by
// replaying the inserts and erases made since the last such query, in the order lc_rmq_fill_krmq makes them.
static int lc_rmq_fill_arr(int max_dist, int max_dist_inner, int bw, int max_chn_skip, int cap_rmq_size, float chn_pen_gap, float chn_pen_skip,
						   int64_t n, const mm128_t *a, int32_t *f, int64_t *p, int32_t *t, int32_t *v, void *mem)
{
	int32_t ymin, ymax, ny, c, *rk, *nx, *pre, *ie;
	int64_t i, i0, j, k, range, nw, st = 0, st_inner = 0, ki = 0, ke = 0;
	uint64_t *bits;
	double *pri;
	lc_mnode_t *T;
	lc_elem_t *root = 0;
	lc_iset_t si; // the inner window
	kmp_rmq_t *mp;

	if (n >= INT32_MAX || max_dist < 1) return 0;
	ymin = ymax = (int32_t)a[0].y;
	for (i = 1; i < n; ++i) {
		int32_t y = (int32_t)a[i].y;
		ymin = y < ymin? y : ymin;
		ymax = y > ymax? y : ymax;
	}
	if ((int64_t)ymin - max_dist < INT32_MIN) return 0;
	range = (int64_t)ymax - ymin + 1;
	if (range > 512 * n + 65536) return 0;
	// rank the distinct y values with a bitmap and its prefix counts
	nw = (range + 63) >> 6;
	bits = Kcalloc(mem, uint64_t, nw);
	pre = Kmalloc(mem, int32_t, nw);
	for (i = 0; i < n; ++i) {
		uint32_t d = (uint32_t)(int32_t)a[i].y - (uint32_t)ymin;
		bits[d>>6] |= 1ULL << (d&63);
	}
	for (k = 0, c = 0; k < nw; ++k)
		pre[k] = c, c += lc_popcnt64(bits[k]);
	ny = c;
	rk = Kmalloc(mem, int32_t, n);
	nx = Kmalloc(mem, int32_t, n);
	pri = Kmalloc(mem, double, n);
	ie = Kmalloc(mem, int32_t, n); // ie[st]: i0 when st was erased
	T = Kmalloc(mem, lc_mnode_t, 2 * ny);
	for (i = 0; i < n; ++i) {
		uint32_t d = (uint32_t)(int32_t)a[i].y - (uint32_t)ymin;
		rk[i] = pre[d>>6] + lc_popcnt64(bits[d>>6] & ((1ULL<<(d&63)) - 1));
	}
	for (k = 0; k < ny; ++k) // T[ny..2ny) as the last anchor seen per rank, then nx[] links anchors of one rank in index order
		T[ny + k].am = (int32_t)n;
	for (i = n - 1; i >= 0; --i)
		nx[i] = T[ny + rk[i]].am, T[ny + rk[i]].am = (int32_t)i;
	for (k = 1; k < 2 * ny; ++k)
		T[k].mn = HUGE_VAL, T[k].am = -1, T[k].tie = 0;
	mp = kmp_init_rmq(mem);
	memset(&si, 0, sizeof(si));
	si.km = mem;

	for (i = i0 = 0; i < n; ++i) {
		int64_t max_j = -1, l, r;
		int32_t q_span = a[i].y>>32&0xff, max_f = q_span, yi = (int32_t)a[i].y, bam;
		double best;
		int btie;
		lc_elem_t s, *q;
		// add in-range anchors
		if (i0 < i && a[i0].x != a[i].x) {
			for (j = i0; j < i; ++j) {
				double pj = -(f[j] + 0.5 * chn_pen_gap * ((int32_t)a[j].x + (int32_t)a[j].y));
				if (!(pj < HUGE_VAL)) return 0; // +inf or NaN cannot be told from an empty node
				pri[j] = pj;
				for (k = ny + rk[j]; k > 0; k >>= 1) { // pj joins every subtree on the leaf's path
					if (pj < T[k].mn) T[k].mn = pj, T[k].am = (int32_t)j, T[k].tie = 0;
					else if (pj == T[k].mn && !T[k].tie) T[k].tie = 1;
					else break;
				}
				if (max_dist_inner > 0)
					lc_iset_insert(&si, (int32_t)a[j].y, j);
			}
			i0 = i;
		}
		// get rid of active chains out of range; the tree holds exactly [st, i0)
		while (st < i && (a[i].x>>32 != a[st].x>>32 || a[i].x > a[st].x + max_dist || (unsigned)(i0 - st) > (unsigned)cap_rmq_size)) {
			lc_mnode_t *x = &T[ny + rk[st]];
			if (st >= i0) return 0; // krmq would find st absent and insert it later
			ie[st] = (int32_t)i0;
			if (pri[st] == x->mn) { // st is the leaf's minimum or ties it: recompute the leaf from its other anchors in [st+1, i0)
				x->mn = HUGE_VAL, x->am = -1, x->tie = 0;
				for (k = nx[st]; k < i0; k = nx[k]) {
					if (pri[k] < x->mn) x->mn = pri[k], x->am = (int32_t)k, x->tie = 0;
					else if (pri[k] == x->mn) x->tie = 1;
				}
				lc_mt_fix(T, ny + rk[st]);
			}
			++st;
		}
		if (max_dist_inner > 0)  { // similar to the block above, but applied to the inner tree
			while (st_inner < i && (a[i].x>>32 != a[st_inner].x>>32 || a[i].x > a[st_inner].x + max_dist_inner || si.size > cap_rmq_size)) {
				lc_iset_erase(&si, (int32_t)a[st_inner].y, st_inner);
				++st_inner;
			}
		}
		// RMQ over y in (yi - max_dist, yi), plus anchor 0 at y == yi (krmq's upper key is (yi, 0))
		best = HUGE_VAL, bam = -1, btie = 0;
		for (l = lc_cnt_lt((int64_t)(yi - max_dist) + 1, ymin, ymax, ny, bits, pre) + ny, r = rk[i] + ny; l < r; l >>= 1, r >>= 1) {
			if (l & 1) {
				if (T[l].mn < best) best = T[l].mn, bam = T[l].am, btie = T[l].tie;
				else if (T[l].mn == best && T[l].am >= 0) btie = 1;
				++l;
			}
			if (r & 1) {
				--r;
				if (T[r].mn < best) best = T[r].mn, bam = T[r].am, btie = T[r].tie;
				else if (T[r].mn == best && T[r].am >= 0) btie = 1;
			}
		}
		if (st == 0 && i0 > 0 && (int32_t)a[0].y == yi) {
			if (pri[0] < best) best = pri[0], bam = 0, btie = 0;
			else if (pri[0] == best) btie = 1;
		}
		if (bam >= 0) {
			if (btie) {
				lc_elem_t lo, hi;
				for (; ke < st; ++ke) {
					for (; ki < ie[ke]; ++ki) {
						q = kmp_alloc_rmq(mp);
						q->y = (int32_t)a[ki].y, q->i = ki, q->pri = pri[ki];
						krmq_insert(lc_elem, &root, q, 0);
					}
					s.y = (int32_t)a[ke].y, s.i = ke;
					if (root && (q = krmq_erase(lc_elem, &root, &s, 0)) != 0)
						kmp_free_rmq(mp, q);
				}
				for (; ki < i0; ++ki) {
					q = kmp_alloc_rmq(mp);
					q->y = (int32_t)a[ki].y, q->i = ki, q->pri = pri[ki];
					krmq_insert(lc_elem, &root, q, 0);
				}
				lo.i = INT32_MAX, lo.y = yi - max_dist;
				hi.i = 0, hi.y = yi;
				q = krmq_rmq(lc_elem, root, &lo, &hi);
				assert(q && q->pri == best);
				bam = (int32_t)q->i;
			}
			lc_rmq_score(i, bam, a, f, p, t, bw, max_dist_inner, max_chn_skip, chn_pen_gap, chn_pen_skip, &si, &max_f, &max_j);
		}
		// set max
		assert(max_j < 0 || (a[max_j].x < a[i].x && (int32_t)a[max_j].y < (int32_t)a[i].y));
		f[i] = max_f, p[i] = max_j;
		v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f; // v[] keeps the peak score up to i; f[] is the score ending at i, not always the peak
	}
	return 1;
}

// fills f[], p[] and v[] with a krmq main tree; the fallback when lc_rmq_fill_arr returns 0
static void lc_rmq_fill_krmq(int max_dist, int max_dist_inner, int bw, int max_chn_skip, int cap_rmq_size, float chn_pen_gap, float chn_pen_skip,
							 int64_t n, const mm128_t *a, int32_t *f, int64_t *p, int32_t *t, int32_t *v, void *mem)
{
	int64_t i, i0, st = 0, st_inner = 0;
	lc_elem_t *root = 0;
	lc_iset_t si; // the inner window
	kmp_rmq_t *mp;

	mp = kmp_init_rmq(mem);
	memset(&si, 0, sizeof(si));
	si.km = mem;
	for (i = i0 = 0; i < n; ++i) {
		int64_t max_j = -1;
		int32_t q_span = a[i].y>>32&0xff, max_f = q_span;
		lc_elem_t s, *q, lo, hi;
		// add in-range anchors
		if (i0 < i && a[i0].x != a[i].x) {
			int64_t j;
			for (j = i0; j < i; ++j) {
				q = kmp_alloc_rmq(mp);
				q->y = (int32_t)a[j].y, q->i = j, q->pri = -(f[j] + 0.5 * chn_pen_gap * ((int32_t)a[j].x + (int32_t)a[j].y));
				krmq_insert(lc_elem, &root, q, 0);
				if (max_dist_inner > 0)
					lc_iset_insert(&si, q->y, j);
			}
			i0 = i;
		}
		// get rid of active chains out of range
		while (st < i && (a[i].x>>32 != a[st].x>>32 || a[i].x > a[st].x + max_dist || krmq_size(head, root) > cap_rmq_size)) {
			s.y = (int32_t)a[st].y, s.i = st;
			if (root && (q = krmq_erase(lc_elem, &root, &s, 0)) != 0) // krmq_erase reads *root first; an absent key returns 0
				kmp_free_rmq(mp, q);
			++st;
		}
		if (max_dist_inner > 0)  { // similar to the block above, but applied to the inner tree
			while (st_inner < i && (a[i].x>>32 != a[st_inner].x>>32 || a[i].x > a[st_inner].x + max_dist_inner || si.size > cap_rmq_size)) {
				lc_iset_erase(&si, (int32_t)a[st_inner].y, st_inner);
				++st_inner;
			}
		}
		// RMQ
		lo.i = INT32_MAX, lo.y = (int32_t)a[i].y - max_dist;
		hi.i = 0, hi.y = (int32_t)a[i].y;
		if ((q = krmq_rmq(lc_elem, root, &lo, &hi)) != 0) {
			assert(q->y >= lo.y && q->y <= hi.y);
			lc_rmq_score(i, q->i, a, f, p, t, bw, max_dist_inner, max_chn_skip, chn_pen_gap, chn_pen_skip, &si, &max_f, &max_j);
		}
		// set max
		assert(max_j < 0 || (a[max_j].x < a[i].x && (int32_t)a[max_j].y < (int32_t)a[i].y));
		f[i] = max_f, p[i] = max_j;
		v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f; // v[] keeps the peak score up to i; f[] is the score ending at i, not always the peak
	}
}

mm128_t *mg_lchain_rmq(int max_dist, int max_dist_inner, int bw, int max_chn_skip, int cap_rmq_size, int min_cnt, int min_sc, float chn_pen_gap, float chn_pen_skip,
					   int64_t n, mm128_t *a, int *n_u_, uint64_t **_u, void *km)
{
	int32_t *f,*t, *v, n_u, n_v, max_drop = bw;
	int64_t *p;
	uint64_t *u;
	void *mem_mp = 0;

	if (_u) *_u = 0, *n_u_ = 0;
	if (n == 0 || a == 0) {
		kfree(km, a);
		return 0;
	}
	if (max_dist < bw) max_dist = bw;
	if (max_dist_inner < 0) max_dist_inner = 0;
	if (max_dist_inner > max_dist) max_dist_inner = max_dist;
	p = Kmalloc(km, int64_t, n);
	f = Kmalloc(km, int32_t, n);
	t = Kcalloc(km, int32_t, n);
	v = Kmalloc(km, int32_t, n);
	mem_mp = km_init2(km, 0x10000);

	// fill the score and backtrack arrays
	if (!lc_rmq_fill_arr(max_dist, max_dist_inner, bw, max_chn_skip, cap_rmq_size, chn_pen_gap, chn_pen_skip, n, a, f, p, t, v, mem_mp)) {
		memset(t, 0, n * 4);
		lc_rmq_fill_krmq(max_dist, max_dist_inner, bw, max_chn_skip, cap_rmq_size, chn_pen_gap, chn_pen_skip, n, a, f, p, t, v, mem_mp);
	}
	km_destroy(mem_mp);

	u = mg_chain_backtrack(km, n, f, p, v, t, min_cnt, min_sc, max_drop, &n_u, &n_v);
	*n_u_ = n_u, *_u = u; // NB: note that u[] may not be sorted by score here
	kfree(km, p); kfree(km, f); kfree(km, t);
	if (n_u == 0) {
		kfree(km, a); kfree(km, v);
		return 0;
	}
	return compact_a(km, n_u, u, n_v, v, a);
}
#endif // MG_LCHAIN_AVX2_UNIT
