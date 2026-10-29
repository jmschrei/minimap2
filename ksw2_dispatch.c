#ifdef KSW_CPU_DISPATCH
#include <stdlib.h>
#include "ksw2.h"
#ifdef _MSC_VER
#include <immintrin.h>
#endif

#define SIMD_SSE     0x1
#define SIMD_SSE2    0x2
#define SIMD_SSE3    0x4
#define SIMD_SSSE3   0x8
#define SIMD_SSE4_1  0x10
#define SIMD_SSE4_2  0x20
#define SIMD_AVX     0x40
#define SIMD_AVX2    0x80
#define SIMD_AVX512F 0x100
#define SIMD_AVX512BW 0x200 // with AVX512F and AVX512VL, and the OS saves the opmask and ZMM state

#ifndef _MSC_VER
// adapted from https://github.com/01org/linux-sgx/blob/master/common/inc/internal/linux/cpuid_gnu.h
void __cpuidex(int cpuid[4], int func_id, int subfunc_id)
{
#if defined(__x86_64__)
	__asm__ volatile ("cpuid"
			: "=a" (cpuid[0]), "=b" (cpuid[1]), "=c" (cpuid[2]), "=d" (cpuid[3])
			: "0" (func_id), "2" (subfunc_id));
#else // on 32bit, ebx can NOT be used as PIC code
	__asm__ volatile ("xchgl %%ebx, %1; cpuid; xchgl %%ebx, %1"
			: "=a" (cpuid[0]), "=r" (cpuid[1]), "=c" (cpuid[2]), "=d" (cpuid[3])
			: "0" (func_id), "2" (subfunc_id));
#endif
}
#endif

// XCR0: bits 1 and 2 are set when the OS saves the XMM and YMM registers, bits 5-7 the opmask and ZMM ones; call only when CPUID reports OSXSAVE
static unsigned ksw_xgetbv0(void)
{
#ifdef _MSC_VER
	return (unsigned)_xgetbv(0);
#else
	unsigned eax, edx;
	__asm__ volatile ("xgetbv" : "=a" (eax), "=d" (edx) : "c" (0));
	return eax;
#endif
}

static int ksw_simd = -1;

static int x86_simd(void)
{
	int flag = 0, cpuid[4], max_id, os_ymm, os_zmm;
	__cpuidex(cpuid, 0, 0);
	max_id = cpuid[0];
	if (max_id == 0) return 0;
	__cpuidex(cpuid, 1, 0);
	if (cpuid[3]>>25&1) flag |= SIMD_SSE;
	if (cpuid[3]>>26&1) flag |= SIMD_SSE2;
	if (cpuid[2]>>0 &1) flag |= SIMD_SSE3;
	if (cpuid[2]>>9 &1) flag |= SIMD_SSSE3;
	if (cpuid[2]>>19&1) flag |= SIMD_SSE4_1;
	if (cpuid[2]>>20&1) flag |= SIMD_SSE4_2;
	if (cpuid[2]>>28&1) flag |= SIMD_AVX;
	os_ymm = (cpuid[2]>>27&1) && (cpuid[2]>>28&1) && (ksw_xgetbv0() & 6) == 6; // OSXSAVE, AVX, and the OS saves YMM state
	os_zmm = os_ymm && (ksw_xgetbv0() & 0xe0) == 0xe0; // and the opmask and ZMM state
	if (max_id >= 7) {
		__cpuidex(cpuid, 7, 0);
		if ((cpuid[1]>>5 &1) && os_ymm) flag |= SIMD_AVX2;
		if (cpuid[1]>>16&1) flag |= SIMD_AVX512F;
		if ((cpuid[1]>>16&1) && (cpuid[1]>>30&1) && (cpuid[1]>>31&1) && os_zmm) flag |= SIMD_AVX512BW; // F, BW and VL
	}
	return flag;
}

int ksw_cpu_avx2(void) // also read by lchain.c's mg_lchain_dp() dispatch
{
	if (ksw_simd < 0) ksw_simd = x86_simd();
	return (ksw_simd & SIMD_AVX2) != 0;
}

void ksw_extz2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
{
	extern void ksw_extz2_sse2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez);
	extern void ksw_extz2_sse41(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez);
#ifdef KSW_HAVE_AVX2
	extern void ksw_extz2_avx2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez);
#endif
#ifdef KSW_HAVE_AVX512
	extern void ksw_extz2_avx512(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez);
#endif
	if (ksw_simd < 0) ksw_simd = x86_simd();
#ifdef KSW_HAVE_AVX512
	if (ksw_simd & SIMD_AVX512BW)
		ksw_extz2_avx512(km, qlen, query, tlen, target, m, mat, q, e, w, zdrop, end_bonus, flag, ez);
	else
#endif
#ifdef KSW_HAVE_AVX2
	if (ksw_simd & SIMD_AVX2)
		ksw_extz2_avx2(km, qlen, query, tlen, target, m, mat, q, e, w, zdrop, end_bonus, flag, ez);
	else
#endif
	if (ksw_simd & SIMD_SSE4_1)
		ksw_extz2_sse41(km, qlen, query, tlen, target, m, mat, q, e, w, zdrop, end_bonus, flag, ez);
	else if (ksw_simd & SIMD_SSE2)
		ksw_extz2_sse2(km, qlen, query, tlen, target, m, mat, q, e, w, zdrop, end_bonus, flag, ez);
	else abort();
}

void ksw_extd2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
{
	extern void ksw_extd2_sse2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez);
	extern void ksw_extd2_sse41(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez);
#ifdef KSW_HAVE_AVX2
	extern void ksw_extd2_avx2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez);
#endif
#ifdef KSW_HAVE_AVX512
	extern void ksw_extd2_avx512(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez);
#endif
	if (ksw_simd < 0) ksw_simd = x86_simd();
#ifdef KSW_HAVE_AVX512
	if (ksw_simd & SIMD_AVX512BW)
		ksw_extd2_avx512(km, qlen, query, tlen, target, m, mat, q, e, q2, e2, w, zdrop, end_bonus, flag, ez);
	else
#endif
#ifdef KSW_HAVE_AVX2
	if (ksw_simd & SIMD_AVX2)
		ksw_extd2_avx2(km, qlen, query, tlen, target, m, mat, q, e, q2, e2, w, zdrop, end_bonus, flag, ez);
	else
#endif
	if (ksw_simd & SIMD_SSE4_1)
		ksw_extd2_sse41(km, qlen, query, tlen, target, m, mat, q, e, q2, e2, w, zdrop, end_bonus, flag, ez);
	else if (ksw_simd & SIMD_SSE2)
		ksw_extd2_sse2(km, qlen, query, tlen, target, m, mat, q, e, q2, e2, w, zdrop, end_bonus, flag, ez);
	else abort();
}

void ksw_exts2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t noncan, int zdrop, int end_bonus, int8_t junc_bonus, int8_t junc_pen, int flag, const uint8_t *junc, ksw_extz_t *ez)
{
	extern void ksw_exts2_sse2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t noncan, int zdrop, int end_bonus, int8_t junc_bonus, int8_t junc_pen, int flag, const uint8_t *junc, ksw_extz_t *ez);
	extern void ksw_exts2_sse41(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t noncan, int zdrop, int end_bonus, int8_t junc_bonus, int8_t junc_pen, int flag, const uint8_t *junc, ksw_extz_t *ez);
#ifdef KSW_HAVE_AVX2
	extern void ksw_exts2_avx2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t noncan, int zdrop, int end_bonus, int8_t junc_bonus, int8_t junc_pen, int flag, const uint8_t *junc, ksw_extz_t *ez);
#endif
#ifdef KSW_HAVE_AVX512
	extern void ksw_exts2_avx512(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t noncan, int zdrop, int end_bonus, int8_t junc_bonus, int8_t junc_pen, int flag, const uint8_t *junc, ksw_extz_t *ez);
#endif
	if (ksw_simd < 0) ksw_simd = x86_simd();
#ifdef KSW_HAVE_AVX512
	if (ksw_simd & SIMD_AVX512BW)
		ksw_exts2_avx512(km, qlen, query, tlen, target, m, mat, q, e, q2, noncan, zdrop, end_bonus, junc_bonus, junc_pen, flag, junc, ez);
	else
#endif
#ifdef KSW_HAVE_AVX2
	if (ksw_simd & SIMD_AVX2)
		ksw_exts2_avx2(km, qlen, query, tlen, target, m, mat, q, e, q2, noncan, zdrop, end_bonus, junc_bonus, junc_pen, flag, junc, ez);
	else
#endif
	if (ksw_simd & SIMD_SSE4_1)
		ksw_exts2_sse41(km, qlen, query, tlen, target, m, mat, q, e, q2, noncan, zdrop, end_bonus, junc_bonus, junc_pen, flag, junc, ez);
	else if (ksw_simd & SIMD_SSE2)
		ksw_exts2_sse2(km, qlen, query, tlen, target, m, mat, q, e, q2, noncan, zdrop, end_bonus, junc_bonus, junc_pen, flag, junc, ez);
	else abort();
}
#endif
