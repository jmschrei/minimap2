#include <string.h>
#include <assert.h>
#include "ksw2.h"

#if defined(KSW_CPU_DISPATCH) && defined(KSW_AVX512)
// AVX-512BW build: 64 cells per instruction. It computes the same cells as the SSE4.1 build (each row's
// band rounded to 16 cells; a last step of 1-3 16-cell blocks is computed 64 wide and stored 16-48 wide),
// writes the same s[] bytes and breaks ties in max_t as its 4-lane loop does, so every output is identical.
#if !defined(__AVX512BW__) || !defined(__AVX512VL__)
#error "KSW_AVX512 needs -mavx512bw -mavx512vl"
#endif
#include <immintrin.h>

// helpers of the exact-max block
#if defined(__GNUC__) || defined(__clang__)
#define ksw_ctz(x) __builtin_ctz(x)
#else
static inline int ksw_ctz(unsigned x) // x != 0
{
	int i = 0;
	for (; !(x & 1); x >>= 1) ++i;
	return i;
}
#endif

static inline int32_t ksw_ld4(const void *p) // 4 bytes, unaligned
{
	int32_t x;
	memcpy(&x, p, 4);
	return x;
}

void ksw_extz2_avx512(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
{
#define __shl1(cur, prev) _mm512_alignr_epi8((cur), _mm512_alignr_epi64((cur), (prev), 6), 15) /* prev[63], cur[0..62] */
#define __ld64(ptr) _mm512_loadu_si512((const void*)(ptr))
#define __st64(ptr, val) _mm512_storeu_si512((void*)(ptr), (val))
#define __st16(ptr, val) _mm_store_si128((__m128i*)(ptr), _mm512_castsi512_si128(val))
#define __st32(ptr, val) _mm256_storeu_si256((__m256i*)(ptr), _mm512_castsi512_si256(val))
#define __st48(ptr, val) do { __m512i v_ = (val); __st32((ptr), v_); _mm_store_si128((__m128i*)(ptr) + 2, _mm512_extracti32x4_epi32(v_, 2)); } while (0)

	/* the SSE build ORs _mm_cvtsi32_si128(x1), sign extension included, into its first block: x1e_ and v1e_ */
#define __dp_code_block1 \
	z = _mm512_add_epi8(__ld64(&s[t]), qe2_); \
	xt1 = __ld64(&x[t]);                                  /* xt1 <- x[r-1][t..t+63] */ \
	tmp = _mm512_or_si512(__shl1(xt1, x1_), x1e_);    /* tmp <- x[r-1][t-1..t+62] */ \
	x1_ = xt1, xt1 = tmp, x1e_ = zero_; \
	vt1 = __ld64(&v[t]); \
	tmp = _mm512_or_si512(__shl1(vt1, v1_), v1e_); \
	v1_ = vt1, vt1 = tmp, v1e_ = zero_; \
	a = _mm512_add_epi8(xt1, vt1); \
	ut = __ld64(&u[t]); \
	b = _mm512_add_epi8(__ld64(&y[t]), ut);

#define __dp_code_block2(ST) \
	z = _mm512_max_epu8(z, b);                        /* z = max(z, b); this works because both are non-negative */ \
	z = _mm512_min_epu8(z, max_sc_); \
	ST(&u[t], _mm512_sub_epi8(z, vt1)); \
	ST(&v[t], _mm512_sub_epi8(z, ut)); \
	z = _mm512_sub_epi8(z, q_); \
	a = _mm512_sub_epi8(a, z); \
	b = _mm512_sub_epi8(b, z);

#define __dp_score_only(ST) { \
	__m512i z, a, b, xt1, vt1, ut, tmp; \
	__dp_code_block1; \
	z = _mm512_max_epi8(z, a); \
	__dp_code_block2(ST); \
	ST(&x[t], _mm512_max_epi8(a, zero_)); \
	ST(&y[t], _mm512_max_epi8(b, zero_)); \
}

	/* the SSE build's blends on compare results become masked moves; its flag ORs become masked adds, as each flag bit is 0 before */
#define __dp_left(ST) { \
	__m512i d, z, a, b, xt1, vt1, ut, tmp; \
	__dp_code_block1; \
	d = _mm512_maskz_mov_epi8(_mm512_cmpgt_epi8_mask(a, z), flag1_);   /* d = a > z? 1 : 0 */ \
	z = _mm512_max_epi8(z, a); \
	d = _mm512_mask_mov_epi8(d, _mm512_cmpgt_epi8_mask(b, z), flag2_); /* d = b > z? 2 : d */ \
	__dp_code_block2(ST); \
	ST(&x[t], _mm512_max_epi8(a, zero_)); \
	d = _mm512_mask_add_epi8(d, _mm512_cmpgt_epi8_mask(a, zero_), d, flag8_);  /* d |= a > 0? 0x08 : 0 */ \
	ST(&y[t], _mm512_max_epi8(b, zero_)); \
	d = _mm512_mask_add_epi8(d, _mm512_cmpgt_epi8_mask(b, zero_), d, flag16_); \
	ST(&pr[t], d); \
}

#define __dp_right(ST) { \
	__m512i d, z, a, b, xt1, vt1, ut, tmp; \
	__dp_code_block1; \
	d = _mm512_maskz_mov_epi8(_mm512_cmpge_epi8_mask(a, z), flag1_);   /* d = z > a? 0 : 1 */ \
	z = _mm512_max_epi8(z, a); \
	d = _mm512_mask_mov_epi8(d, _mm512_cmpge_epi8_mask(b, z), flag2_); /* d = z > b? d : 2 */ \
	__dp_code_block2(ST); \
	ST(&x[t], _mm512_max_epi8(a, zero_)); \
	d = _mm512_mask_add_epi8(d, _mm512_cmpge_epi8_mask(a, zero_), d, flag8_);  /* d |= 0 > a? 0 : 0x08 */ \
	ST(&y[t], _mm512_max_epi8(b, zero_)); \
	d = _mm512_mask_add_epi8(d, _mm512_cmpge_epi8_mask(b, zero_), d, flag16_); \
	ST(&pr[t], d); \
}

	/* blocks [st_, en_] in steps of 4; a last step of 1-3 blocks loads 64 bytes (mem has the slack) and stores only its blocks */
#define __dp_row(DP) { \
	for (t = st_; t + 3 <= en_; t += 4) DP(__st64); \
	if (t == en_) DP(__st16) \
	else if (t + 1 == en_) DP(__st32) \
	else if (t + 2 == en_) DP(__st48) \
}

	int r, t, qe = q + e, n_col_, *off = 0, *off_end = 0, tlen_, qlen_, last_st, last_en, wl, wr, max_sc, min_sc;
	int with_cigar = !(flag&KSW_EZ_SCORE_ONLY), approx_max = !!(flag&KSW_EZ_APPROX_MAX);
	int32_t *H = 0, H0 = 0, last_H0_t = 0;
	uint8_t *qr, *sf, *mem, *mem2 = 0;
	__m512i q_, qe2_, zero_, flag1_, flag2_, flag8_, flag16_, sc_mch_, sc_mis_, sc_N_, m1_, max_sc_;
	__m128i *u, *v, *x, *y, *s, *p = 0;

	ksw_reset_extz(ez);
	if (m <= 0 || qlen <= 0 || tlen <= 0) return;

	zero_   = _mm512_set1_epi8(0);
	q_      = _mm512_set1_epi8(q);
	qe2_    = _mm512_set1_epi8((q + e) * 2);
	flag1_  = _mm512_set1_epi8(1);
	flag2_  = _mm512_set1_epi8(2);
	flag8_  = _mm512_set1_epi8(0x08);
	flag16_ = _mm512_set1_epi8(0x10);
	sc_mch_ = _mm512_set1_epi8(mat[0]);
	sc_mis_ = _mm512_set1_epi8(mat[1]);
	sc_N_   = mat[m*m-1] == 0? _mm512_set1_epi8(-e) : _mm512_set1_epi8(mat[m*m-1]);
	m1_     = _mm512_set1_epi8(m - 1); // wildcard
	max_sc_ = _mm512_set1_epi8(mat[0] + (q + e) * 2);

	if (w < 0) w = tlen > qlen? tlen : qlen;
	wl = wr = w;
	tlen_ = (tlen + 15) / 16;
	n_col_ = qlen < tlen? qlen : tlen;
	n_col_ = ((n_col_ < w + 1? n_col_ : w + 1) + 15) / 16 + 1;
	qlen_ = (qlen + 15) / 16;
	for (t = 1, max_sc = mat[0], min_sc = mat[1]; t < m * m; ++t) {
		max_sc = max_sc > mat[t]? max_sc : mat[t];
		min_sc = min_sc < mat[t]? min_sc : mat[t];
	}
	if (-min_sc > 2 * (q + e)) return; // otherwise, we won't see any mismatches

	mem = (uint8_t*)kcalloc(km, tlen_ * 6 + qlen_ + 4, 16); // 3 blocks of slack: a row's last step loads up to 48 bytes past its last block
	u = (__m128i*)(((size_t)mem + 15) >> 4 << 4); // 16-byte aligned
	v = u + tlen_, x = v + tlen_, y = x + tlen_, s = y + tlen_, sf = (uint8_t*)(s + tlen_), qr = sf + tlen_ * 16;
	if (!approx_max) {
		H = (int32_t*)kmalloc(km, (tlen_ * 16 + 4) * 4);
		for (t = 0; t < tlen_ * 16 + 4; ++t) H[t] = KSW_NEG_INF;
	}
	if (with_cigar) {
		mem2 = (uint8_t*)kmalloc(km, ((size_t)(qlen + tlen - 1) * n_col_ + 1) * 16);
		p = (__m128i*)(((size_t)mem2 + 15) >> 4 << 4);
		off = (int*)kmalloc(km, (qlen + tlen - 1) * sizeof(int) * 2);
		off_end = off + qlen + tlen - 1;
	}

	for (t = 0; t < qlen; ++t) qr[t] = query[qlen - 1 - t];
	memcpy(sf, target, tlen);

	for (r = 0, last_st = last_en = -1; r < qlen + tlen - 1; ++r) {
		int st = 0, en = tlen - 1, st0, en0, st_, en_;
		int8_t x1, v1;
		uint8_t *qrr = qr + (qlen - 1 - r), *u8 = (uint8_t*)u, *v8 = (uint8_t*)v;
		__m512i x1_, v1_, x1e_, v1e_;
		// find the boundaries
		if (st < r - qlen + 1) st = r - qlen + 1;
		if (en > r) en = r;
		if (st < (r-wr+1)>>1) st = (r-wr+1)>>1; // take the ceil
		if (en > (r+wl)>>1) en = (r+wl)>>1; // take the floor
		if (st > en) {
			ez->zdropped = 1;
			break;
		}
		st0 = st, en0 = en;
		st = st / 16 * 16, en = (en + 16) / 16 * 16 - 1;
		// set boundary conditions
		if (st > 0) {
			if (st - 1 >= last_st && st - 1 <= last_en)
				x1 = ((uint8_t*)x)[st - 1], v1 = v8[st - 1]; // (r-1,s-1) calculated in the last round
			else x1 = v1 = 0; // not calculated; set to zeros
		} else x1 = 0, v1 = r? q : 0;
		if (en >= r) ((uint8_t*)y)[r] = 0, u8[r] = r? q : 0;
		// loop fission: set scores first
		if (!(flag & KSW_EZ_GENERIC_SC)) {
			int en2 = st0 + ((en0 - st0) / 16 + 1) * 16; // the SSE build writes s[st0..en2-1]
			__m512i sq, sr, tmp;
			for (t = st0; t + 64 <= en2; t += 64) {
				sq = _mm512_loadu_si512((const void*)&sf[t]);
				sr = _mm512_loadu_si512((const void*)&qrr[t]);
				tmp = _mm512_mask_blend_epi8(_mm512_cmpeq_epi8_mask(sq, sr), sc_mis_, sc_mch_);
				tmp = _mm512_mask_blend_epi8(_mm512_cmpeq_epi8_mask(sq, m1_) | _mm512_cmpeq_epi8_mask(sr, m1_), tmp, sc_N_);
				_mm512_storeu_si512((void*)((uint8_t*)s + t), tmp);
			}
			if (t < en2) { // 16, 32 or 48 bytes
				__mmask64 tail_ = (__mmask64)(((uint64_t)1 << (en2 - t)) - 1);
				sq = _mm512_maskz_loadu_epi8(tail_, (const void*)&sf[t]);
				sr = _mm512_maskz_loadu_epi8(tail_, (const void*)&qrr[t]);
				tmp = _mm512_mask_blend_epi8(_mm512_cmpeq_epi8_mask(sq, sr), sc_mis_, sc_mch_);
				tmp = _mm512_mask_blend_epi8(_mm512_cmpeq_epi8_mask(sq, m1_) | _mm512_cmpeq_epi8_mask(sr, m1_), tmp, sc_N_);
				_mm512_mask_storeu_epi8((void*)((uint8_t*)s + t), tail_, tmp);
			}
		} else {
			for (t = st0; t <= en0; ++t)
				((uint8_t*)s)[t] = mat[sf[t] * m + qrr[t]];
		}
		// core loop
		x1_ = v1_ = zero_;
		x1e_ = _mm512_inserti32x4(zero_, _mm_cvtsi32_si128(x1), 0);
		v1e_ = _mm512_inserti32x4(zero_, _mm_cvtsi32_si128(v1), 0);
		st_ = st / 16, en_ = en / 16;
		assert(en_ - st_ + 1 <= n_col_);
		if (!with_cigar) { // score only
			__dp_row(__dp_score_only);
		} else if (!(flag&KSW_EZ_RIGHT)) { // gap left-alignment
			__m128i *pr = p + (size_t)r * n_col_ - st_;
			off[r] = st, off_end[r] = en;
			__dp_row(__dp_left);
		} else { // gap right-alignment
			__m128i *pr = p + (size_t)r * n_col_ - st_;
			off[r] = st, off_end[r] = en;
			__dp_row(__dp_right);
		}
		if (!approx_max) { // find the exact max with a 32-bit score array
			int32_t max_H, max_t = 0, need_t = 1;
			// compute H[] and max_H; max_t only on rows where ksw_apply_zdrop can use it
			if (r > 0) {
				int32_t en1 = st0 + (en0 - st0) / 4 * 4, en16 = st0 + (en0 - st0) / 16 * 16, H_en0, ma;
				__m128i max_H_, H1, gt, mx, qe_;
				__m256i max_H8, qe8_;
				__m512i max_H16, qe16_;
				H_en0 = en0 > 0? H[en0-1] + u8[en0] - qe : H[en0] + v8[en0] - qe; // special casing the last element; stored after the tail group, which covers it
				max_H16 = _mm512_set1_epi32(H_en0);
				qe16_  = _mm512_set1_epi32(q + e);
				for (t = st0; t < en16; t += 16) { // H[t]+=v8[t]-qe; lane j holds the SSE build's lane j % 4 at t + j / 4 * 4, each the largest H[t]
					__m512i H16;
					H16 = _mm512_loadu_si512((const void*)&H[t]);
					H16 = _mm512_add_epi32(H16, _mm512_sub_epi32(_mm512_cvtepu8_epi32(_mm_loadu_si128((__m128i*)&v8[t])), qe16_));
					_mm512_storeu_si512((void*)&H[t], H16);
					max_H16 = _mm512_max_epi32(max_H16, H16);
				}
				max_H8 = _mm256_max_epi32(_mm512_castsi512_si256(max_H16), _mm512_extracti64x4_epi64(max_H16, 1));
				qe8_   = _mm512_castsi512_si256(qe16_);
				if (t + 8 <= en1) { // one 8-cell group
					__m256i H8;
					H8 = _mm256_loadu_si256((__m256i*)&H[t]);
					H8 = _mm256_add_epi32(H8, _mm256_sub_epi32(_mm256_cvtepu8_epi32(_mm_loadl_epi64((__m128i*)&v8[t])), qe8_));
					_mm256_storeu_si256((__m256i*)&H[t], H8);
					max_H8 = _mm256_max_epi32(max_H8, H8);
					t += 8;
				}
				// fold lanes 4-7 into 0-3: lane j then holds the largest H[t] over t = st0 + j (mod 4), as in the SSE build
				max_H_ = _mm_max_epi32(_mm256_castsi256_si128(max_H8), _mm256_extracti128_si256(max_H8, 1));
				qe_    = _mm256_castsi256_si128(qe8_);
				if (t < en1) { // one 4-cell group [t, en1)
					H1 = _mm_loadu_si128((__m128i*)&H[t]);
					H1 = _mm_add_epi32(H1, _mm_sub_epi32(_mm_cvtepu8_epi32(_mm_cvtsi32_si128(ksw_ld4(&v8[t]))), qe_));
					_mm_storeu_si128((__m128i*)&H[t], H1);
					max_H_ = _mm_max_epi32(max_H_, H1);
				}
				// the 0-3 cells [en1, en0) as one group; the other lanes add 0, and H has 4 cells of slack
				gt = _mm_cmpgt_epi32(_mm_set1_epi32(en0 - en1), _mm_setr_epi32(0, 1, 2, 3));
				H1 = _mm_loadu_si128((__m128i*)&H[en1]);
				H1 = _mm_add_epi32(H1, _mm_and_si128(_mm_sub_epi32(_mm_cvtepu8_epi32(_mm_cvtsi32_si128(ksw_ld4(&v8[en1]))), qe_), gt));
				_mm_storeu_si128((__m128i*)&H[en1], H1);
				H[en0] = H_en0;
				mx = _mm_set1_epi32(H_en0);
				H1 = _mm_blendv_epi8(mx, H1, gt);
				mx = _mm_max_epi32(max_H_, H1);
				mx = _mm_max_epi32(mx, _mm_shuffle_epi32(mx, 0x4e));
				mx = _mm_max_epi32(mx, _mm_shuffle_epi32(mx, 0xb1));
				max_H = _mm_cvtsi128_si32(mx);
				// otherwise ksw_apply_zdrop returns 0 and changes nothing, whatever max_t, as e >= 0
				need_t = max_H > (int32_t)ez->max || (zdrop >= 0 && (e < 0 || (int32_t)ez->max - max_H > zdrop));
				if (need_t) { // ties go to H[en0], then lanes 0-3 in order (each at its first t), then the tail cells
					if (max_H == H_en0) max_t = en0;
					else if ((ma = ~_mm_movemask_ps(_mm_castsi128_ps(_mm_cmpgt_epi32(mx, max_H_))) & 0xf) != 0)
						for (max_t = st0 + ksw_ctz(ma); H[max_t] != max_H; max_t += 4); // max_H > H[en0], so a cell of this lane holds it
					else max_t = en1 + ksw_ctz(~_mm_movemask_ps(_mm_castsi128_ps(_mm_cmpgt_epi32(mx, H1))) & 0xf);
				}
			} else H[0] = v8[0] - qe - qe, max_H = H[0], max_t = 0; // special casing r==0
			// update ez
			if (en0 == tlen - 1 && H[en0] > ez->mte)
				ez->mte = H[en0], ez->mte_q = r - en0;
			if (r - st0 == qlen - 1 && H[st0] > ez->mqe)
				ez->mqe = H[st0], ez->mqe_t = st0;
			if (need_t && ksw_apply_zdrop(ez, 1, max_H, r, max_t, zdrop, e)) break;
			if (r == qlen + tlen - 2 && en0 == tlen - 1)
				ez->score = H[tlen - 1];
		} else { // find approximate max; Z-drop might be inaccurate, too.
			if (r > 0) {
				if (last_H0_t >= st0 && last_H0_t <= en0 && last_H0_t + 1 >= st0 && last_H0_t + 1 <= en0) {
					int32_t d0 = v8[last_H0_t] - qe;
					int32_t d1 = u8[last_H0_t + 1] - qe;
					if (d0 > d1) H0 += d0;
					else H0 += d1, ++last_H0_t;
				} else if (last_H0_t >= st0 && last_H0_t <= en0) {
					H0 += v8[last_H0_t] - qe;
				} else {
					++last_H0_t, H0 += u8[last_H0_t] - qe;
				}
				if ((flag & KSW_EZ_APPROX_DROP) && ksw_apply_zdrop(ez, 1, H0, r, last_H0_t, zdrop, e)) break;
			} else H0 = v8[0] - qe - qe, last_H0_t = 0;
			if (r == qlen + tlen - 2 && en0 == tlen - 1)
				ez->score = H0;
		}
		last_st = st, last_en = en;
	}
	kfree(km, mem);
	if (!approx_max) kfree(km, H);
	if (with_cigar) { // backtrack
		int rev_cigar = !!(flag & KSW_EZ_REV_CIGAR);
		if (!ez->zdropped && !(flag&KSW_EZ_EXTZ_ONLY)) {
			ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*16, tlen-1, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		} else if (!ez->zdropped && (flag&KSW_EZ_EXTZ_ONLY) && ez->mqe + end_bonus > (int)ez->max) {
			ez->reach_end = 1;
			ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*16, ez->mqe_t, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		} else if (ez->max_t >= 0 && ez->max_q >= 0) {
			ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*16, ez->max_t, ez->max_q, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		}
		kfree(km, mem2); kfree(km, off);
	}
}
#elif defined(KSW_CPU_DISPATCH) && defined(KSW_AVX2)
// AVX2 build: 32 cells per instruction. It computes the same cells as the SSE4.1 build (each row's
// band rounded to 16 cells; an odd last 16-cell block is computed 32 wide and stored 16 wide), writes
// the same s[] bytes and breaks ties in max_t as its 4-lane loop does, so every output is identical.
#ifndef __AVX2__
#error "KSW_AVX2 needs -mavx2"
#endif
#include <immintrin.h>

// helpers of the exact-max block
#if defined(__GNUC__) || defined(__clang__)
#define ksw_ctz(x) __builtin_ctz(x)
#else
static inline int ksw_ctz(unsigned x) // x != 0
{
	int i = 0;
	for (; !(x & 1); x >>= 1) ++i;
	return i;
}
#endif

static inline int32_t ksw_ld4(const void *p) // 4 bytes, unaligned
{
	int32_t x;
	memcpy(&x, p, 4);
	return x;
}

void ksw_extz2_avx2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
{
#define __shl1(cur, prev) _mm256_alignr_epi8((cur), _mm256_permute2x128_si256((cur), (prev), 0x03), 15) /* prev[31], cur[0..30] */
#define __st32(ptr, val) _mm256_storeu_si256((__m256i*)(ptr), (val))
#define __st16(ptr, val) _mm_store_si128((__m128i*)(ptr), _mm256_castsi256_si128(val))

	/* the SSE build ORs _mm_cvtsi32_si128(x1), sign extension included, into its first block: x1e_ and v1e_ */
#define __dp_code_block1 \
	z = _mm256_add_epi8(_mm256_loadu_si256((__m256i*)&s[t]), qe2_); \
	xt1 = _mm256_loadu_si256((__m256i*)&x[t]);        /* xt1 <- x[r-1][t..t+31] */ \
	tmp = _mm256_or_si256(__shl1(xt1, x1_), x1e_);    /* tmp <- x[r-1][t-1..t+30] */ \
	x1_ = xt1, xt1 = tmp, x1e_ = zero_; \
	vt1 = _mm256_loadu_si256((__m256i*)&v[t]); \
	tmp = _mm256_or_si256(__shl1(vt1, v1_), v1e_); \
	v1_ = vt1, vt1 = tmp, v1e_ = zero_; \
	a = _mm256_add_epi8(xt1, vt1); \
	ut = _mm256_loadu_si256((__m256i*)&u[t]); \
	b = _mm256_add_epi8(_mm256_loadu_si256((__m256i*)&y[t]), ut);

#define __dp_code_block2(ST) \
	z = _mm256_max_epu8(z, b);                        /* z = max(z, b); this works because both are non-negative */ \
	z = _mm256_min_epu8(z, max_sc_); \
	ST(&u[t], _mm256_sub_epi8(z, vt1)); \
	ST(&v[t], _mm256_sub_epi8(z, ut)); \
	z = _mm256_sub_epi8(z, q_); \
	a = _mm256_sub_epi8(a, z); \
	b = _mm256_sub_epi8(b, z);

#define __dp_score_only(ST) { \
	__m256i z, a, b, xt1, vt1, ut, tmp; \
	__dp_code_block1; \
	z = _mm256_max_epi8(z, a); \
	__dp_code_block2(ST); \
	ST(&x[t], _mm256_max_epi8(a, zero_)); \
	ST(&y[t], _mm256_max_epi8(b, zero_)); \
}

#define __dp_left(ST) { \
	__m256i d, z, a, b, xt1, vt1, ut, tmp; \
	__dp_code_block1; \
	d = _mm256_and_si256(_mm256_cmpgt_epi8(a, z), flag1_); \
	z = _mm256_max_epi8(z, a); \
	tmp = _mm256_cmpgt_epi8(b, z); \
	d = _mm256_blendv_epi8(d, flag2_, tmp); \
	__dp_code_block2(ST); \
	tmp = _mm256_cmpgt_epi8(a, zero_); \
	ST(&x[t], _mm256_and_si256(tmp, a)); \
	d = _mm256_or_si256(d, _mm256_and_si256(tmp, flag8_)); \
	tmp = _mm256_cmpgt_epi8(b, zero_); \
	ST(&y[t], _mm256_and_si256(tmp, b)); \
	d = _mm256_or_si256(d, _mm256_and_si256(tmp, flag16_)); \
	ST(&pr[t], d); \
}

#define __dp_right(ST) { \
	__m256i d, z, a, b, xt1, vt1, ut, tmp; \
	__dp_code_block1; \
	d = _mm256_andnot_si256(_mm256_cmpgt_epi8(z, a), flag1_); \
	z = _mm256_max_epi8(z, a); \
	tmp = _mm256_cmpgt_epi8(z, b); \
	d = _mm256_blendv_epi8(flag2_, d, tmp); \
	__dp_code_block2(ST); \
	tmp = _mm256_cmpgt_epi8(zero_, a); \
	ST(&x[t], _mm256_andnot_si256(tmp, a)); \
	d = _mm256_or_si256(d, _mm256_andnot_si256(tmp, flag8_)); \
	tmp = _mm256_cmpgt_epi8(zero_, b); \
	ST(&y[t], _mm256_andnot_si256(tmp, b)); \
	d = _mm256_or_si256(d, _mm256_andnot_si256(tmp, flag16_)); \
	ST(&pr[t], d); \
}

	int r, t, qe = q + e, n_col_, *off = 0, *off_end = 0, tlen_, qlen_, last_st, last_en, wl, wr, max_sc, min_sc;
	int with_cigar = !(flag&KSW_EZ_SCORE_ONLY), approx_max = !!(flag&KSW_EZ_APPROX_MAX);
	int32_t *H = 0, H0 = 0, last_H0_t = 0;
	uint8_t *qr, *sf, *mem, *mem2 = 0;
	__m256i q_, qe2_, zero_, flag1_, flag2_, flag8_, flag16_, sc_mch_, sc_mis_, sc_N_, m1_, max_sc_;
	__m128i *u, *v, *x, *y, *s, *p = 0;

	ksw_reset_extz(ez);
	if (m <= 0 || qlen <= 0 || tlen <= 0) return;

	zero_   = _mm256_set1_epi8(0);
	q_      = _mm256_set1_epi8(q);
	qe2_    = _mm256_set1_epi8((q + e) * 2);
	flag1_  = _mm256_set1_epi8(1);
	flag2_  = _mm256_set1_epi8(2);
	flag8_  = _mm256_set1_epi8(0x08);
	flag16_ = _mm256_set1_epi8(0x10);
	sc_mch_ = _mm256_set1_epi8(mat[0]);
	sc_mis_ = _mm256_set1_epi8(mat[1]);
	sc_N_   = mat[m*m-1] == 0? _mm256_set1_epi8(-e) : _mm256_set1_epi8(mat[m*m-1]);
	m1_     = _mm256_set1_epi8(m - 1); // wildcard
	max_sc_ = _mm256_set1_epi8(mat[0] + (q + e) * 2);

	if (w < 0) w = tlen > qlen? tlen : qlen;
	wl = wr = w;
	tlen_ = (tlen + 15) / 16;
	n_col_ = qlen < tlen? qlen : tlen;
	n_col_ = ((n_col_ < w + 1? n_col_ : w + 1) + 15) / 16 + 1;
	qlen_ = (qlen + 15) / 16;
	for (t = 1, max_sc = mat[0], min_sc = mat[1]; t < m * m; ++t) {
		max_sc = max_sc > mat[t]? max_sc : mat[t];
		min_sc = min_sc < mat[t]? min_sc : mat[t];
	}
	if (-min_sc > 2 * (q + e)) return; // otherwise, we won't see any mismatches

	mem = (uint8_t*)kcalloc(km, tlen_ * 6 + qlen_ + 1, 16);
	u = (__m128i*)(((size_t)mem + 15) >> 4 << 4); // 16-byte aligned
	v = u + tlen_, x = v + tlen_, y = x + tlen_, s = y + tlen_, sf = (uint8_t*)(s + tlen_), qr = sf + tlen_ * 16;
	if (!approx_max) {
		H = (int32_t*)kmalloc(km, (tlen_ * 16 + 4) * 4);
		for (t = 0; t < tlen_ * 16 + 4; ++t) H[t] = KSW_NEG_INF;
	}
	if (with_cigar) {
		mem2 = (uint8_t*)kmalloc(km, ((size_t)(qlen + tlen - 1) * n_col_ + 1) * 16);
		p = (__m128i*)(((size_t)mem2 + 15) >> 4 << 4);
		off = (int*)kmalloc(km, (qlen + tlen - 1) * sizeof(int) * 2);
		off_end = off + qlen + tlen - 1;
	}

	for (t = 0; t < qlen; ++t) qr[t] = query[qlen - 1 - t];
	memcpy(sf, target, tlen);

	for (r = 0, last_st = last_en = -1; r < qlen + tlen - 1; ++r) {
		int st = 0, en = tlen - 1, st0, en0, st_, en_;
		int8_t x1, v1;
		uint8_t *qrr = qr + (qlen - 1 - r), *u8 = (uint8_t*)u, *v8 = (uint8_t*)v;
		__m256i x1_, v1_, x1e_, v1e_;
		// find the boundaries
		if (st < r - qlen + 1) st = r - qlen + 1;
		if (en > r) en = r;
		if (st < (r-wr+1)>>1) st = (r-wr+1)>>1; // take the ceil
		if (en > (r+wl)>>1) en = (r+wl)>>1; // take the floor
		if (st > en) {
			ez->zdropped = 1;
			break;
		}
		st0 = st, en0 = en;
		st = st / 16 * 16, en = (en + 16) / 16 * 16 - 1;
		// set boundary conditions
		if (st > 0) {
			if (st - 1 >= last_st && st - 1 <= last_en)
				x1 = ((uint8_t*)x)[st - 1], v1 = v8[st - 1]; // (r-1,s-1) calculated in the last round
			else x1 = v1 = 0; // not calculated; set to zeros
		} else x1 = 0, v1 = r? q : 0;
		if (en >= r) ((uint8_t*)y)[r] = 0, u8[r] = r? q : 0;
		// loop fission: set scores first
		if (!(flag & KSW_EZ_GENERIC_SC)) {
			int en2 = st0 + ((en0 - st0) / 16 + 1) * 16; // the SSE build writes s[st0..en2-1]
			for (t = st0; t + 32 <= en2; t += 32) {
				__m256i sq, sr, tmp, mask;
				sq = _mm256_loadu_si256((__m256i*)&sf[t]);
				sr = _mm256_loadu_si256((__m256i*)&qrr[t]);
				mask = _mm256_or_si256(_mm256_cmpeq_epi8(sq, m1_), _mm256_cmpeq_epi8(sr, m1_));
				tmp = _mm256_cmpeq_epi8(sq, sr);
				tmp = _mm256_blendv_epi8(sc_mis_, sc_mch_, tmp);
				tmp = _mm256_blendv_epi8(tmp,     sc_N_,   mask);
				_mm256_storeu_si256((__m256i*)((uint8_t*)s + t), tmp);
			}
			if (t < en2) {
				__m128i sq, sr, tmp, mask;
				sq = _mm_loadu_si128((__m128i*)&sf[t]);
				sr = _mm_loadu_si128((__m128i*)&qrr[t]);
				mask = _mm_or_si128(_mm_cmpeq_epi8(sq, _mm256_castsi256_si128(m1_)), _mm_cmpeq_epi8(sr, _mm256_castsi256_si128(m1_)));
				tmp = _mm_cmpeq_epi8(sq, sr);
				tmp = _mm_blendv_epi8(_mm256_castsi256_si128(sc_mis_), _mm256_castsi256_si128(sc_mch_), tmp);
				tmp = _mm_blendv_epi8(tmp, _mm256_castsi256_si128(sc_N_), mask);
				_mm_storeu_si128((__m128i*)((uint8_t*)s + t), tmp);
			}
		} else {
			for (t = st0; t <= en0; ++t)
				((uint8_t*)s)[t] = mat[sf[t] * m + qrr[t]];
		}
		// core loop
		x1_ = v1_ = zero_;
		x1e_ = _mm256_inserti128_si256(zero_, _mm_cvtsi32_si128(x1), 0);
		v1e_ = _mm256_inserti128_si256(zero_, _mm_cvtsi32_si128(v1), 0);
		st_ = st / 16, en_ = en / 16;
		assert(en_ - st_ + 1 <= n_col_);
		if (!with_cigar) { // score only
			for (t = st_; t < en_; t += 2) __dp_score_only(__st32);
			if (t == en_) __dp_score_only(__st16);
		} else if (!(flag&KSW_EZ_RIGHT)) { // gap left-alignment
			__m128i *pr = p + (size_t)r * n_col_ - st_;
			off[r] = st, off_end[r] = en;
			for (t = st_; t < en_; t += 2) __dp_left(__st32);
			if (t == en_) __dp_left(__st16);
		} else { // gap right-alignment
			__m128i *pr = p + (size_t)r * n_col_ - st_;
			off[r] = st, off_end[r] = en;
			for (t = st_; t < en_; t += 2) __dp_right(__st32);
			if (t == en_) __dp_right(__st16);
		}
		if (!approx_max) { // find the exact max with a 32-bit score array
			int32_t max_H, max_t = 0, need_t = 1;
			// compute H[] and max_H; max_t only on rows where ksw_apply_zdrop can use it
			if (r > 0) {
				int32_t en1 = st0 + (en0 - st0) / 4 * 4, en8 = st0 + (en0 - st0) / 8 * 8, H_en0, ma;
				__m128i max_H_, H1, gt, mx, qe_;
				__m256i max_H8, qe8_;
				H_en0 = en0 > 0? H[en0-1] + u8[en0] - qe : H[en0] + v8[en0] - qe; // special casing the last element; stored after the tail group, which covers it
				max_H8 = _mm256_set1_epi32(H_en0);
				qe8_   = _mm256_set1_epi32(q + e);
				for (t = st0; t < en8; t += 8) { // H[t]+=v8[t]-qe; lanes 0-3 and 4-7 hold the SSE build's lanes 0-3 at t and t+4, each the largest H[t]
					__m256i H8;
					H8 = _mm256_loadu_si256((__m256i*)&H[t]);
					H8 = _mm256_add_epi32(H8, _mm256_sub_epi32(_mm256_cvtepu8_epi32(_mm_loadl_epi64((__m128i*)&v8[t])), qe8_));
					_mm256_storeu_si256((__m256i*)&H[t], H8);
					max_H8 = _mm256_max_epi32(max_H8, H8);
				}
				// fold lanes 4-7 into 0-3: lane j then holds the largest H[t] over t = st0 + j (mod 4), as in the SSE build
				max_H_ = _mm_max_epi32(_mm256_castsi256_si128(max_H8), _mm256_extracti128_si256(max_H8, 1));
				qe_ = _mm256_castsi256_si128(qe8_);
				if (t < en1) { // one 4-cell group [en8, en1)
					H1 = _mm_loadu_si128((__m128i*)&H[t]);
					H1 = _mm_add_epi32(H1, _mm_sub_epi32(_mm_cvtepu8_epi32(_mm_cvtsi32_si128(ksw_ld4(&v8[t]))), qe_));
					_mm_storeu_si128((__m128i*)&H[t], H1);
					max_H_ = _mm_max_epi32(max_H_, H1);
				}
				// the 0-3 cells [en1, en0) as one group; the other lanes add 0, and H has 4 cells of slack
				gt = _mm_cmpgt_epi32(_mm_set1_epi32(en0 - en1), _mm_setr_epi32(0, 1, 2, 3));
				H1 = _mm_loadu_si128((__m128i*)&H[en1]);
				H1 = _mm_add_epi32(H1, _mm_and_si128(_mm_sub_epi32(_mm_cvtepu8_epi32(_mm_cvtsi32_si128(ksw_ld4(&v8[en1]))), qe_), gt));
				_mm_storeu_si128((__m128i*)&H[en1], H1);
				H[en0] = H_en0;
				mx = _mm_set1_epi32(H_en0);
				H1 = _mm_blendv_epi8(mx, H1, gt);
				mx = _mm_max_epi32(max_H_, H1);
				mx = _mm_max_epi32(mx, _mm_shuffle_epi32(mx, 0x4e));
				mx = _mm_max_epi32(mx, _mm_shuffle_epi32(mx, 0xb1));
				max_H = _mm_cvtsi128_si32(mx);
				// otherwise ksw_apply_zdrop returns 0 and changes nothing, whatever max_t, as e >= 0
				need_t = max_H > (int32_t)ez->max || (zdrop >= 0 && (e < 0 || (int32_t)ez->max - max_H > zdrop));
				if (need_t) { // ties go to H[en0], then lanes 0-3 in order (each at its first t), then the tail cells
					if (max_H == H_en0) max_t = en0;
					else if ((ma = ~_mm_movemask_ps(_mm_castsi128_ps(_mm_cmpgt_epi32(mx, max_H_))) & 0xf) != 0)
						for (max_t = st0 + ksw_ctz(ma); H[max_t] != max_H; max_t += 4); // max_H > H[en0], so a cell of this lane holds it
					else max_t = en1 + ksw_ctz(~_mm_movemask_ps(_mm_castsi128_ps(_mm_cmpgt_epi32(mx, H1))) & 0xf);
				}
			} else H[0] = v8[0] - qe - qe, max_H = H[0], max_t = 0; // special casing r==0
			// update ez
			if (en0 == tlen - 1 && H[en0] > ez->mte)
				ez->mte = H[en0], ez->mte_q = r - en0;
			if (r - st0 == qlen - 1 && H[st0] > ez->mqe)
				ez->mqe = H[st0], ez->mqe_t = st0;
			if (need_t && ksw_apply_zdrop(ez, 1, max_H, r, max_t, zdrop, e)) break;
			if (r == qlen + tlen - 2 && en0 == tlen - 1)
				ez->score = H[tlen - 1];
		} else { // find approximate max; Z-drop might be inaccurate, too.
			if (r > 0) {
				if (last_H0_t >= st0 && last_H0_t <= en0 && last_H0_t + 1 >= st0 && last_H0_t + 1 <= en0) {
					int32_t d0 = v8[last_H0_t] - qe;
					int32_t d1 = u8[last_H0_t + 1] - qe;
					if (d0 > d1) H0 += d0;
					else H0 += d1, ++last_H0_t;
				} else if (last_H0_t >= st0 && last_H0_t <= en0) {
					H0 += v8[last_H0_t] - qe;
				} else {
					++last_H0_t, H0 += u8[last_H0_t] - qe;
				}
				if ((flag & KSW_EZ_APPROX_DROP) && ksw_apply_zdrop(ez, 1, H0, r, last_H0_t, zdrop, e)) break;
			} else H0 = v8[0] - qe - qe, last_H0_t = 0;
			if (r == qlen + tlen - 2 && en0 == tlen - 1)
				ez->score = H0;
		}
		last_st = st, last_en = en;
	}
	kfree(km, mem);
	if (!approx_max) kfree(km, H);
	if (with_cigar) { // backtrack
		int rev_cigar = !!(flag & KSW_EZ_REV_CIGAR);
		if (!ez->zdropped && !(flag&KSW_EZ_EXTZ_ONLY)) {
			ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*16, tlen-1, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		} else if (!ez->zdropped && (flag&KSW_EZ_EXTZ_ONLY) && ez->mqe + end_bonus > (int)ez->max) {
			ez->reach_end = 1;
			ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*16, ez->mqe_t, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		} else if (ez->max_t >= 0 && ez->max_q >= 0) {
			ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*16, ez->max_t, ez->max_q, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		}
		kfree(km, mem2); kfree(km, off);
	}
}
#elif defined(__SSE2__)
#ifdef USE_SIMDE
#include <simde/x86/sse2.h>
#else
#include <emmintrin.h>
#endif

#ifdef KSW_SSE2_ONLY
#undef __SSE4_1__
#endif

#ifdef __SSE4_1__
#ifdef USE_SIMDE
#include <simde/x86/sse4.1.h>
#else
#include <smmintrin.h>
#endif
#endif

// helpers of the exact-max block below
#if defined(__GNUC__) || defined(__clang__)
#define ksw_ctz(x) __builtin_ctz(x)
#else
static inline int ksw_ctz(unsigned x) // x != 0
{
	int i = 0;
	for (; !(x & 1); x >>= 1) ++i;
	return i;
}
#endif

static inline __m128i ksw_load_epu8_epi32(const uint8_t *p) // p[0..3], each zero-extended to 32 bits
{
#ifdef __SSE4_1__
	int32_t x;
	memcpy(&x, p, 4);
	return _mm_cvtepu8_epi32(_mm_cvtsi32_si128(x));
#else
	return _mm_setr_epi32(p[0], p[1], p[2], p[3]);
#endif
}

static inline __m128i ksw_max_epi32(__m128i a, __m128i b)
{
#ifdef __SSE4_1__
	return _mm_max_epi32(a, b);
#else
	__m128i m = _mm_cmpgt_epi32(a, b);
	return _mm_or_si128(_mm_and_si128(m, a), _mm_andnot_si128(m, b));
#endif
}

#ifdef KSW_CPU_DISPATCH
#ifdef __SSE4_1__
void ksw_extz2_sse41(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
#else
void ksw_extz2_sse2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
#endif
#else
void ksw_extz2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
#endif // ~KSW_CPU_DISPATCH
{
#define __dp_code_block1 \
	z = _mm_add_epi8(_mm_load_si128(&s[t]), qe2_); \
	xt1 = _mm_load_si128(&x[t]);                     /* xt1 <- x[r-1][t..t+15] */ \
	tmp = _mm_srli_si128(xt1, 15);                   /* tmp <- x[r-1][t+15] */ \
	xt1 = _mm_or_si128(_mm_slli_si128(xt1, 1), x1_); /* xt1 <- x[r-1][t-1..t+14] */ \
	x1_ = tmp; \
	vt1 = _mm_load_si128(&v[t]);                     /* vt1 <- v[r-1][t..t+15] */ \
	tmp = _mm_srli_si128(vt1, 15);                   /* tmp <- v[r-1][t+15] */ \
	vt1 = _mm_or_si128(_mm_slli_si128(vt1, 1), v1_); /* vt1 <- v[r-1][t-1..t+14] */ \
	v1_ = tmp; \
	a = _mm_add_epi8(xt1, vt1);                      /* a <- x[r-1][t-1..t+14] + v[r-1][t-1..t+14] */ \
	ut = _mm_load_si128(&u[t]);                      /* ut <- u[t..t+15] */ \
	b = _mm_add_epi8(_mm_load_si128(&y[t]), ut);     /* b <- y[r-1][t..t+15] + u[r-1][t..t+15] */

#define __dp_code_block2 \
	z = _mm_max_epu8(z, b);                          /* z = max(z, b); this works because both are non-negative */ \
	z = _mm_min_epu8(z, max_sc_); \
	_mm_store_si128(&u[t], _mm_sub_epi8(z, vt1));    /* u[r][t..t+15] <- z - v[r-1][t-1..t+14] */ \
	_mm_store_si128(&v[t], _mm_sub_epi8(z, ut));     /* v[r][t..t+15] <- z - u[r-1][t..t+15] */ \
	z = _mm_sub_epi8(z, q_); \
	a = _mm_sub_epi8(a, z); \
	b = _mm_sub_epi8(b, z);

	int r, t, qe = q + e, n_col_, *off = 0, *off_end = 0, tlen_, qlen_, last_st, last_en, wl, wr, max_sc, min_sc;
	int with_cigar = !(flag&KSW_EZ_SCORE_ONLY), approx_max = !!(flag&KSW_EZ_APPROX_MAX);
	int32_t *H = 0, H0 = 0, last_H0_t = 0;
	uint8_t *qr, *sf, *mem, *mem2 = 0;
	__m128i q_, qe2_, zero_, flag1_, flag2_, flag8_, flag16_, sc_mch_, sc_mis_, sc_N_, m1_, max_sc_;
	__m128i *u, *v, *x, *y, *s, *p = 0;

	ksw_reset_extz(ez);
	if (m <= 0 || qlen <= 0 || tlen <= 0) return;

	zero_   = _mm_set1_epi8(0);
	q_      = _mm_set1_epi8(q);
	qe2_    = _mm_set1_epi8((q + e) * 2);
	flag1_  = _mm_set1_epi8(1);
	flag2_  = _mm_set1_epi8(2);
	flag8_  = _mm_set1_epi8(0x08);
	flag16_ = _mm_set1_epi8(0x10);
	sc_mch_ = _mm_set1_epi8(mat[0]);
	sc_mis_ = _mm_set1_epi8(mat[1]);
	sc_N_   = mat[m*m-1] == 0? _mm_set1_epi8(-e) : _mm_set1_epi8(mat[m*m-1]);
	m1_     = _mm_set1_epi8(m - 1); // wildcard
	max_sc_ = _mm_set1_epi8(mat[0] + (q + e) * 2);

	if (w < 0) w = tlen > qlen? tlen : qlen;
	wl = wr = w;
	tlen_ = (tlen + 15) / 16;
	n_col_ = qlen < tlen? qlen : tlen;
	n_col_ = ((n_col_ < w + 1? n_col_ : w + 1) + 15) / 16 + 1;
	qlen_ = (qlen + 15) / 16;
	for (t = 1, max_sc = mat[0], min_sc = mat[1]; t < m * m; ++t) {
		max_sc = max_sc > mat[t]? max_sc : mat[t];
		min_sc = min_sc < mat[t]? min_sc : mat[t];
	}
	if (-min_sc > 2 * (q + e)) return; // otherwise, we won't see any mismatches

	mem = (uint8_t*)kcalloc(km, tlen_ * 6 + qlen_ + 1, 16);
	u = (__m128i*)(((size_t)mem + 15) >> 4 << 4); // 16-byte aligned
	v = u + tlen_, x = v + tlen_, y = x + tlen_, s = y + tlen_, sf = (uint8_t*)(s + tlen_), qr = sf + tlen_ * 16;
	if (!approx_max) {
		H = (int32_t*)kmalloc(km, (tlen_ * 16 + 4) * 4);
		for (t = 0; t < tlen_ * 16 + 4; ++t) H[t] = KSW_NEG_INF;
	}
	if (with_cigar) {
		mem2 = (uint8_t*)kmalloc(km, ((size_t)(qlen + tlen - 1) * n_col_ + 1) * 16);
		p = (__m128i*)(((size_t)mem2 + 15) >> 4 << 4);
		off = (int*)kmalloc(km, (qlen + tlen - 1) * sizeof(int) * 2);
		off_end = off + qlen + tlen - 1;
	}

	for (t = 0; t < qlen; ++t) qr[t] = query[qlen - 1 - t];
	memcpy(sf, target, tlen);

	for (r = 0, last_st = last_en = -1; r < qlen + tlen - 1; ++r) {
		int st = 0, en = tlen - 1, st0, en0, st_, en_;
		int8_t x1, v1;
		uint8_t *qrr = qr + (qlen - 1 - r), *u8 = (uint8_t*)u, *v8 = (uint8_t*)v;
		__m128i x1_, v1_;
		// find the boundaries
		if (st < r - qlen + 1) st = r - qlen + 1;
		if (en > r) en = r;
		if (st < (r-wr+1)>>1) st = (r-wr+1)>>1; // take the ceil
		if (en > (r+wl)>>1) en = (r+wl)>>1; // take the floor
		if (st > en) {
			ez->zdropped = 1;
			break;
		}
		st0 = st, en0 = en;
		st = st / 16 * 16, en = (en + 16) / 16 * 16 - 1;
		// set boundary conditions
		if (st > 0) {
			if (st - 1 >= last_st && st - 1 <= last_en)
				x1 = ((uint8_t*)x)[st - 1], v1 = v8[st - 1]; // (r-1,s-1) calculated in the last round
			else x1 = v1 = 0; // not calculated; set to zeros
		} else x1 = 0, v1 = r? q : 0;
		if (en >= r) ((uint8_t*)y)[r] = 0, u8[r] = r? q : 0;
		// loop fission: set scores first
		if (!(flag & KSW_EZ_GENERIC_SC)) {
			for (t = st0; t <= en0; t += 16) {
				__m128i sq, st, tmp, mask;
				sq = _mm_loadu_si128((__m128i*)&sf[t]);
				st = _mm_loadu_si128((__m128i*)&qrr[t]);
				mask = _mm_or_si128(_mm_cmpeq_epi8(sq, m1_), _mm_cmpeq_epi8(st, m1_));
				tmp = _mm_cmpeq_epi8(sq, st);
#ifdef __SSE4_1__
				tmp = _mm_blendv_epi8(sc_mis_, sc_mch_, tmp);
				tmp = _mm_blendv_epi8(tmp,     sc_N_,   mask);
#else
				tmp = _mm_or_si128(_mm_andnot_si128(tmp,  sc_mis_), _mm_and_si128(tmp,  sc_mch_));
				tmp = _mm_or_si128(_mm_andnot_si128(mask, tmp),     _mm_and_si128(mask, sc_N_));
#endif
				_mm_storeu_si128((__m128i*)((uint8_t*)s + t), tmp);
			}
		} else {
			for (t = st0; t <= en0; ++t)
				((uint8_t*)s)[t] = mat[sf[t] * m + qrr[t]];
		}
		// core loop
		x1_ = _mm_cvtsi32_si128(x1);
		v1_ = _mm_cvtsi32_si128(v1);
		st_ = st / 16, en_ = en / 16;
		assert(en_ - st_ + 1 <= n_col_);
		if (!with_cigar) { // score only
			for (t = st_; t <= en_; ++t) {
				__m128i z, a, b, xt1, vt1, ut, tmp;
				__dp_code_block1;
#ifdef __SSE4_1__
				z = _mm_max_epi8(z, a);                          // z = z > a? z : a (signed)
#else // we need to emulate SSE4.1 intrinsics _mm_max_epi8()
				z = _mm_and_si128(z, _mm_cmpgt_epi8(z, zero_));  // z = z > 0? z : 0;
				z = _mm_max_epu8(z, a);                          // z = max(z, a); this works because both are non-negative
#endif
				__dp_code_block2;
#ifdef __SSE4_1__
				_mm_store_si128(&x[t], _mm_max_epi8(a, zero_));
				_mm_store_si128(&y[t], _mm_max_epi8(b, zero_));
#else
				tmp = _mm_cmpgt_epi8(a, zero_);
				_mm_store_si128(&x[t], _mm_and_si128(a, tmp));
				tmp = _mm_cmpgt_epi8(b, zero_);
				_mm_store_si128(&y[t], _mm_and_si128(b, tmp));
#endif
			}
		} else if (!(flag&KSW_EZ_RIGHT)) { // gap left-alignment
			__m128i *pr = p + (size_t)r * n_col_ - st_;
			off[r] = st, off_end[r] = en;
			for (t = st_; t <= en_; ++t) {
				__m128i d, z, a, b, xt1, vt1, ut, tmp;
				__dp_code_block1;
				d = _mm_and_si128(_mm_cmpgt_epi8(a, z), flag1_); // d = a > z? 1 : 0
#ifdef __SSE4_1__
				z = _mm_max_epi8(z, a);                          // z = z > a? z : a (signed)
				tmp = _mm_cmpgt_epi8(b, z);
				d = _mm_blendv_epi8(d, flag2_, tmp);             // d = b > z? 2 : d
#else // we need to emulate SSE4.1 intrinsics _mm_max_epi8() and _mm_blendv_epi8()
				z = _mm_and_si128(z, _mm_cmpgt_epi8(z, zero_));  // z = z > 0? z : 0;
				z = _mm_max_epu8(z, a);                          // z = max(z, a); this works because both are non-negative
				tmp = _mm_cmpgt_epi8(b, z);
				d = _mm_or_si128(_mm_andnot_si128(tmp, d), _mm_and_si128(tmp, flag2_)); // d = b > z? 2 : d; emulating blendv
#endif
				__dp_code_block2;
				tmp = _mm_cmpgt_epi8(a, zero_);
				_mm_store_si128(&x[t], _mm_and_si128(tmp, a));
				d = _mm_or_si128(d, _mm_and_si128(tmp, flag8_));  // d = a > 0? 0x08 : 0
				tmp = _mm_cmpgt_epi8(b, zero_);
				_mm_store_si128(&y[t], _mm_and_si128(tmp, b));
				d = _mm_or_si128(d, _mm_and_si128(tmp, flag16_)); // d = b > 0? 0x10 : 0
				_mm_store_si128(&pr[t], d);
			}
		} else { // gap right-alignment
			__m128i *pr = p + (size_t)r * n_col_ - st_;
			off[r] = st, off_end[r] = en;
			for (t = st_; t <= en_; ++t) {
				__m128i d, z, a, b, xt1, vt1, ut, tmp;
				__dp_code_block1;
				d = _mm_andnot_si128(_mm_cmpgt_epi8(z, a), flag1_); // d = z > a? 0 : 1
#ifdef __SSE4_1__
				z = _mm_max_epi8(z, a);                          // z = z > a? z : a (signed)
				tmp = _mm_cmpgt_epi8(z, b);
				d = _mm_blendv_epi8(flag2_, d, tmp);             // d = z > b? d : 2
#else // we need to emulate SSE4.1 intrinsics _mm_max_epi8() and _mm_blendv_epi8()
				z = _mm_and_si128(z, _mm_cmpgt_epi8(z, zero_));  // z = z > 0? z : 0;
				z = _mm_max_epu8(z, a);                          // z = max(z, a); this works because both are non-negative
				tmp = _mm_cmpgt_epi8(z, b);
				d = _mm_or_si128(_mm_andnot_si128(tmp, flag2_), _mm_and_si128(tmp, d)); // d = z > b? d : 2; emulating blendv
#endif
				__dp_code_block2;
				tmp = _mm_cmpgt_epi8(zero_, a);
				_mm_store_si128(&x[t], _mm_andnot_si128(tmp, a));
				d = _mm_or_si128(d, _mm_andnot_si128(tmp, flag8_));  // d = 0 > a? 0 : 0x08
				tmp = _mm_cmpgt_epi8(zero_, b);
				_mm_store_si128(&y[t], _mm_andnot_si128(tmp, b));
				d = _mm_or_si128(d, _mm_andnot_si128(tmp, flag16_)); // d = 0 > b? 0 : 0x10
				_mm_store_si128(&pr[t], d);
			}
		}
		if (!approx_max) { // find the exact max with a 32-bit score array
			int32_t max_H, max_t = 0, need_t = 1;
			// compute H[] and max_H; max_t only on rows where ksw_apply_zdrop can use it
			if (r > 0) {
				int32_t en1 = st0 + (en0 - st0) / 4 * 4, H_en0, ma;
				__m128i max_H_, H1, gt, mx, qe_;
				H_en0 = en0 > 0? H[en0-1] + u8[en0] - qe : H[en0] + v8[en0] - qe; // special casing the last element; stored after the tail group, which covers it
				max_H_ = _mm_set1_epi32(H_en0);
				qe_    = _mm_set1_epi32(q + e);
				for (t = st0; t < en1; t += 4) { // H[t]+=v8[t] - qe; per lane, the largest H[t]
					H1 = _mm_loadu_si128((__m128i*)&H[t]);
					H1 = _mm_add_epi32(H1, ksw_load_epu8_epi32(&v8[t]));
					H1 = _mm_sub_epi32(H1, qe_);
					_mm_storeu_si128((__m128i*)&H[t], H1);
					max_H_ = ksw_max_epi32(max_H_, H1);
				}
				// the 0-3 cells [en1, en0) as one group; the other lanes add 0, and H has 4 cells of slack
				gt = _mm_cmpgt_epi32(_mm_set1_epi32(en0 - en1), _mm_setr_epi32(0, 1, 2, 3));
				H1 = _mm_loadu_si128((__m128i*)&H[en1]);
				H1 = _mm_add_epi32(H1, _mm_and_si128(_mm_sub_epi32(ksw_load_epu8_epi32(&v8[en1]), qe_), gt));
				_mm_storeu_si128((__m128i*)&H[en1], H1);
				H[en0] = H_en0;
				mx = _mm_set1_epi32(H_en0);
				H1 = _mm_or_si128(_mm_and_si128(gt, H1), _mm_andnot_si128(gt, mx));
				mx = ksw_max_epi32(max_H_, H1);
				mx = ksw_max_epi32(mx, _mm_shuffle_epi32(mx, 0x4e));
				mx = ksw_max_epi32(mx, _mm_shuffle_epi32(mx, 0xb1));
				max_H = _mm_cvtsi128_si32(mx);
				// otherwise ksw_apply_zdrop returns 0 and changes nothing, whatever max_t, as e >= 0
				need_t = max_H > (int32_t)ez->max || (zdrop >= 0 && (e < 0 || (int32_t)ez->max - max_H > zdrop));
				if (need_t) { // ties go to H[en0], then lanes 0-3 in order (each at its first t), then the tail cells
					if (max_H == H_en0) max_t = en0;
					else if ((ma = ~_mm_movemask_ps(_mm_castsi128_ps(_mm_cmpgt_epi32(mx, max_H_))) & 0xf) != 0)
						for (max_t = st0 + ksw_ctz(ma); H[max_t] != max_H; max_t += 4); // max_H > H[en0], so a cell of this lane holds it
					else max_t = en1 + ksw_ctz(~_mm_movemask_ps(_mm_castsi128_ps(_mm_cmpgt_epi32(mx, H1))) & 0xf);
				}
			} else H[0] = v8[0] - qe - qe, max_H = H[0], max_t = 0; // special casing r==0
			// update ez
			if (en0 == tlen - 1 && H[en0] > ez->mte)
				ez->mte = H[en0], ez->mte_q = r - en0;
			if (r - st0 == qlen - 1 && H[st0] > ez->mqe)
				ez->mqe = H[st0], ez->mqe_t = st0;
			if (need_t && ksw_apply_zdrop(ez, 1, max_H, r, max_t, zdrop, e)) break;
			if (r == qlen + tlen - 2 && en0 == tlen - 1)
				ez->score = H[tlen - 1];
		} else { // find approximate max; Z-drop might be inaccurate, too.
			if (r > 0) {
				if (last_H0_t >= st0 && last_H0_t <= en0 && last_H0_t + 1 >= st0 && last_H0_t + 1 <= en0) {
					int32_t d0 = v8[last_H0_t] - qe;
					int32_t d1 = u8[last_H0_t + 1] - qe;
					if (d0 > d1) H0 += d0;
					else H0 += d1, ++last_H0_t;
				} else if (last_H0_t >= st0 && last_H0_t <= en0) {
					H0 += v8[last_H0_t] - qe;
				} else {
					++last_H0_t, H0 += u8[last_H0_t] - qe;
				}
				if ((flag & KSW_EZ_APPROX_DROP) && ksw_apply_zdrop(ez, 1, H0, r, last_H0_t, zdrop, e)) break;
			} else H0 = v8[0] - qe - qe, last_H0_t = 0;
			if (r == qlen + tlen - 2 && en0 == tlen - 1)
				ez->score = H0;
		}
		last_st = st, last_en = en;
		//for (t = st0; t <= en0; ++t) printf("(%d,%d)\t(%d,%d,%d,%d)\t%d\n", r, t, ((int8_t*)u)[t], ((int8_t*)v)[t], ((int8_t*)x)[t], ((int8_t*)y)[t], H[t]); // for debugging
	}
	kfree(km, mem);
	if (!approx_max) kfree(km, H);
	if (with_cigar) { // backtrack
		int rev_cigar = !!(flag & KSW_EZ_REV_CIGAR);
		if (!ez->zdropped && !(flag&KSW_EZ_EXTZ_ONLY)) {
			ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*16, tlen-1, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		} else if (!ez->zdropped && (flag&KSW_EZ_EXTZ_ONLY) && ez->mqe + end_bonus > (int)ez->max) {
			ez->reach_end = 1;
			ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*16, ez->mqe_t, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		} else if (ez->max_t >= 0 && ez->max_q >= 0) {
			ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*16, ez->max_t, ez->max_q, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		}
		kfree(km, mem2); kfree(km, off);
	}
}
#endif // __SSE2__
