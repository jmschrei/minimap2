#include <stdlib.h>
#include <assert.h>
#if defined(WIN32) || defined(_WIN32)
#include <io.h> // for open(2)
#define idx_ftell(fp) _ftelli64(fp)
#define idx_fseek(fp, off, whence) _fseeki64((fp), (off), (whence))
#else
#include <unistd.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/mman.h> // for madvise(2)
#ifndef MM_NO_PREAD
#define MM_IDX_PREAD // an index loaded on several threads is read with pread(2), each thread reading its own buckets
#endif
#define idx_ftell(fp) ((int64_t)ftello(fp))
#define idx_fseek(fp, off, whence) fseeko((fp), (off_t)(off), (whence))
#endif
#include <fcntl.h>
#include <stdio.h>
#define __STDC_LIMIT_MACROS
#include "kthread.h"
#include "bseq.h"
#include "minimap.h"
#include "mmpriv.h"
#include "ksw2.h"
#include "kvec.h"
#include "khash.h"

#if defined(__SSE2__) && (defined(__x86_64__) || defined(__i386__))
#define MM_SEQ4_SSE2 // little-endian: byte k of mi->S holds bases 2k (low nibble) and 2k+1 (high nibble)
#include <emmintrin.h>
#endif

#if defined(MADV_HUGEPAGE) && !defined(MM_NO_HUGEPAGE)
#define MM_IDX_HUGEPAGE // an index loaded on several threads lives in 2 MB-aligned blocks advised for huge pages
#endif

#define idx_hash(a) ((a)>>1)
#define idx_eq(a, b) ((a)>>1 == (b)>>1)
KHASH_INIT(idx, uint64_t, uint64_t, 1, idx_hash, idx_eq)
typedef khash_t(idx) idxhash_t;

KHASH_MAP_INIT_STR(str, uint32_t)

#define kroundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, ++(x))

typedef struct mm_idx_bucket_s {
	mm128_v a;   // (minimizer, position) array
	int32_t n;   // size of the _p_ array
	uint64_t *p; // position array for minimizers appearing >1 times
	void *h;     // hash table indexing _p_ and minimizers appearing once
} mm_idx_bucket_t;

typedef struct {
	int32_t st, en, cnt;
	int32_t score:30, strand:2;
} mm_idx_intv1_t;

typedef struct mm_idx_intv_s {
	int32_t n, m;
	mm_idx_intv1_t *a;
} mm_idx_intv_t;

typedef struct mm_idx_jjump_s {
	int32_t n, m;
	mm_idx_jjump1_t *a;
} mm_idx_jjump_t;

typedef struct { // mm_idx_t::mem: the blocks holding p[] and the hash tables of an index loaded on several threads
	size_t n, m;
	void **a;
	uint8_t *cur;     // the block being filled: l of cur_m bytes used
	size_t l, cur_m;
} mm_idx_mem_t;

static void mm_idx_mem_destroy(mm_idx_mem_t *m)
{
	size_t i;
	if (m == 0) return;
	for (i = 0; i < m->n; ++i) free(m->a[i]);
	free(m->a); free(m);
}

mm_idx_t *mm_idx_init(int w, int k, int b, int flag)
{
	mm_idx_t *mi;
	if (k*2 < b) b = k * 2;
	if (w < 1) w = 1;
	mi = (mm_idx_t*)calloc(1, sizeof(mm_idx_t));
	mi->w = w, mi->k = k, mi->b = b, mi->flag = flag;
	mi->B = (mm_idx_bucket_t*)calloc(1<<b, sizeof(mm_idx_bucket_t));
	if (!(mm_dbg_flag & 1)) mi->km = km_init();
	return mi;
}

void mm_idx_destroy(mm_idx_t *mi)
{
	uint32_t i;
	if (mi == 0) return;
	if (mi->h) kh_destroy(str, (khash_t(str)*)mi->h);
	if (mi->B) {
		for (i = 0; i < 1U<<mi->b; ++i) {
			free(mi->B[i].a.a);
			if (mi->mem) continue; // p[] and the hash table are in mi->mem's blocks
			free(mi->B[i].p);
			kh_destroy(idx, (idxhash_t*)mi->B[i].h);
		}
	}
	mm_idx_mem_destroy((mm_idx_mem_t*)mi->mem);
	if (mi->spsc) free(mi->spsc);
	if (mi->I) {
		for (i = 0; i < mi->n_seq; ++i)
			free(mi->I[i].a);
		free(mi->I);
	}
	if (mi->J) {
		for (i = 0; i < mi->n_seq; ++i)
			free(mi->J[i].a);
		free(mi->J);
	}
	if (!mi->km) {
		for (i = 0; i < mi->n_seq; ++i)
			free(mi->seq[i].name);
		free(mi->seq);
	} else km_destroy(mi->km);
	free(mi->B); free(mi->S); free(mi);
}

const uint64_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, int *n)
{
	int mask = (1<<mi->b) - 1;
	khint_t k;
	mm_idx_bucket_t *b = &mi->B[minier&mask];
	idxhash_t *h = (idxhash_t*)b->h;
	*n = 0;
	if (h == 0) return 0;
	k = kh_get(idx, h, minier>>mi->b<<1);
	if (k == kh_end(h)) return 0;
	if (kh_key(h, k)&1) { // special casing when there is only one k-mer
		*n = 1;
		return &kh_val(h, k);
	} else {
		*n = (uint32_t)kh_val(h, k);
		return &b->p[kh_val(h, k)>>32];
	}
}

// prefetch the hash slot mm_idx_get(mi, minier) probes first: its flag word, key and value
void mm_idx_prefetch(const mm_idx_t *mi, uint64_t minier)
{
	const mm_idx_bucket_t *b = &mi->B[minier & ((1U<<mi->b) - 1)];
	const idxhash_t *h = (const idxhash_t*)b->h;
	khint_t i;
	if (h == 0 || h->n_buckets == 0) return;
	i = (khint_t)(minier >> mi->b) & (h->n_buckets - 1);
	mm_prefetch(&h->flags[i>>4]);
	mm_prefetch(&h->keys[i]);
	mm_prefetch(&h->vals[i]);
}

typedef struct { // one thread's sums in mm_idx_stat_mt()
	int64_t n1;
	uint64_t sum;
	uint8_t pad[48]; // a cache line per thread
} idx_stat_acc_t;

typedef struct {
	const mm_idx_t *mi;
	idx_stat_acc_t *acc;
} idx_stat_walk_t;

static void idx_stat_bucket(void *data, long i, int tid) // kt_for() callback
{
	idx_stat_walk_t *w = (idx_stat_walk_t*)data;
	idxhash_t *h = (idxhash_t*)w->mi->B[i].h;
	const uint64_t *keys, *vals;
	int64_t n1 = 0;
	uint64_t sum = 0;
	khint_t k;
	if (h == 0) return;
	keys = h->keys, vals = h->vals;
	if (kh_end(h) < 16) {
		for (k = 0; k < kh_end(h); ++k) {
			if (!kh_exist(h, k)) continue;
			if (keys[k]&1) ++n1;
			else sum += (uint32_t)vals[k];
		}
	} else {
		for (k = 0; k < kh_end(h); k += 16) { // the 16 slots of one flag word, without a branch per slot, as in mm_idx_cal_max_occ
			khint32_t fw = h->flags[k>>4];
			int t;
			for (t = 0; t < 16; ++t) {
				uint64_t e = (fw >> (t<<1) & 3) == 0, s = keys[k + t] & 1;
				n1 += e & s;
				sum += (uint32_t)vals[k + t] & -(e & (s ^ 1));
			}
		}
	}
	w->acc[tid].n1 += n1, w->acc[tid].sum += sum;
}

void mm_idx_stat(const mm_idx_t *mi)
{
	mm_idx_stat_mt(mi, 1);
}

void mm_idx_stat_mt(const mm_idx_t *mi, int n_threads)
{
	int64_t n = 0, n1 = 0;
	uint32_t i;
	uint64_t sum = 0, len = 0;
	idx_stat_walk_t w;
	int t;
	if (n_threads < 1) n_threads = 1;
	fprintf(stderr, "[M::%s] kmer size: %d; skip: %d; is_hpc: %d; #seq: %d\n", "mm_idx_stat", mi->k, mi->w, mi->flag&MM_I_HPC, mi->n_seq);
	for (i = 0; i < mi->n_seq; ++i)
		len += mi->seq[i].len;
	for (i = 0; i < 1U<<mi->b; ++i)
		if (mi->B[i].h) n += kh_size((idxhash_t*)mi->B[i].h);
	w.mi = mi;
	w.acc = (idx_stat_acc_t*)calloc(n_threads, sizeof(idx_stat_acc_t));
	kt_for(n_threads, idx_stat_bucket, &w, 1L<<mi->b);
	for (t = 0; t < n_threads; ++t) // integer sums: the same for any split of the buckets
		n1 += w.acc[t].n1, sum += w.acc[t].sum;
	free(w.acc);
	sum += n1; // a singleton occurs once
	fprintf(stderr, "[M::%s::%.3f*%.2f] distinct minimizers: %ld (%.2f%% are singletons); average occurrences: %.3lf; average spacing: %.3lf; total length: %ld\n",
			"mm_idx_stat", realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), (long)n, 100.0*n1/n, (double)sum / n, (double)len / sum, (long)len);
}

int mm_idx_index_name(mm_idx_t *mi)
{
	khash_t(str) *h;
	uint32_t i;
	int has_dup = 0, absent;
	if (mi->h) return 0;
	h = kh_init(str);
	for (i = 0; i < mi->n_seq; ++i) {
		khint_t k;
		k = kh_put(str, h, mi->seq[i].name, &absent);
		if (absent) kh_val(h, k) = i;
		else has_dup = 1;
	}
	mi->h = h;
	if (has_dup && mm_verbose >= 2)
		fprintf(stderr, "[WARNING] some database sequences have identical sequence names\n");
	return has_dup;
}

int mm_idx_name2id(const mm_idx_t *mi, const char *name)
{
	khash_t(str) *h = (khash_t(str)*)mi->h;
	khint_t k;
	if (h == 0) return -2;
	k = kh_get(str, h, name);
	return k == kh_end(h)? -1 : kh_val(h, k);
}

static inline uint8_t seq4_comp(uint32_t c) { return c < 4? 3 - c : c; }

#ifdef MM_SEQ4_SSE2
static inline __m128i seq4_comp16(__m128i c) // seq4_comp() on 16 bytes
{
	__m128i lt = _mm_cmplt_epi8(c, _mm_set1_epi8(4));
	return _mm_or_si128(_mm_and_si128(lt, _mm_sub_epi8(_mm_set1_epi8(3), c)), _mm_andnot_si128(lt, c));
}

static inline __m128i seq4_rev16(__m128i x) // reverse the 16 bytes
{
	x = _mm_shuffle_epi32(x, 0x1b);
	x = _mm_shufflelo_epi16(x, 0xb1);
	x = _mm_shufflehi_epi16(x, 0xb1);
	return _mm_or_si128(_mm_slli_epi16(x, 8), _mm_srli_epi16(x, 8));
}
#endif

// seq[0..en-i) = bases [i, en) of S, a whole 32-bit word (8 bases) or 16 bytes (32 bases) at a time once i is word-aligned
static void seq4_decode(const uint32_t *S, uint64_t i, uint64_t en, uint8_t *p)
{
	for (; i < en && (i & 7); ++i) *p++ = mm_seq4_get(S, i);
#ifdef MM_SEQ4_SSE2
	for (; i + 32 <= en; i += 32, p += 32) {
		const __m128i m = _mm_set1_epi8(0xf);
		__m128i v = _mm_loadu_si128((const __m128i*)&S[i>>3]);
		__m128i lo = _mm_and_si128(v, m), hi = _mm_and_si128(_mm_srli_epi16(v, 4), m);
		_mm_storeu_si128((__m128i*)p, _mm_unpacklo_epi8(lo, hi));
		_mm_storeu_si128((__m128i*)(p + 16), _mm_unpackhi_epi8(lo, hi));
	}
#endif
	for (; i + 8 <= en; i += 8, p += 8) {
		uint32_t w = S[i>>3];
		p[0] = w & 0xf, p[1] = w >> 4 & 0xf, p[2] = w >> 8 & 0xf, p[3] = w >> 12 & 0xf;
		p[4] = w >> 16 & 0xf, p[5] = w >> 20 & 0xf, p[6] = w >> 24 & 0xf, p[7] = w >> 28;
	}
	for (; i < en; ++i) *p++ = mm_seq4_get(S, i);
}

// seq[0..en-i) = reverse complement of bases [i, en) of S, decoded as in seq4_decode()
static void seq4_decode_rev(const uint32_t *S, uint64_t i, uint64_t en, uint8_t *seq)
{
	uint8_t *p;
	if (i >= en) return;
	p = seq + (en - i); // base i goes to p[-1]
	for (; i < en && (i & 7); ++i) *--p = seq4_comp(mm_seq4_get(S, i));
#ifdef MM_SEQ4_SSE2
	for (; i + 32 <= en; i += 32) {
		const __m128i m = _mm_set1_epi8(0xf);
		__m128i v = _mm_loadu_si128((const __m128i*)&S[i>>3]);
		__m128i lo = _mm_and_si128(v, m), hi = _mm_and_si128(_mm_srli_epi16(v, 4), m);
		__m128i a = seq4_comp16(_mm_unpacklo_epi8(lo, hi)), b = seq4_comp16(_mm_unpackhi_epi8(lo, hi));
		p -= 32;
		_mm_storeu_si128((__m128i*)(p + 16), seq4_rev16(a));
		_mm_storeu_si128((__m128i*)p, seq4_rev16(b));
	}
#endif
	for (; i + 8 <= en; i += 8) {
		uint32_t w = S[i>>3];
		p -= 8;
		p[7] = seq4_comp(w & 0xf), p[6] = seq4_comp(w >> 4 & 0xf), p[5] = seq4_comp(w >> 8 & 0xf), p[4] = seq4_comp(w >> 12 & 0xf);
		p[3] = seq4_comp(w >> 16 & 0xf), p[2] = seq4_comp(w >> 20 & 0xf), p[1] = seq4_comp(w >> 24 & 0xf), p[0] = seq4_comp(w >> 28);
	}
	for (; i < en; ++i) *--p = seq4_comp(mm_seq4_get(S, i));
}

int mm_idx_getseq(const mm_idx_t *mi, uint32_t rid, uint32_t st, uint32_t en, uint8_t *seq)
{
	uint64_t st1, en1;
	if (rid >= mi->n_seq || st >= mi->seq[rid].len) return -1;
	if (en > mi->seq[rid].len) en = mi->seq[rid].len;
	st1 = mi->seq[rid].offset + st;
	en1 = mi->seq[rid].offset + en;
	seq4_decode(mi->S, st1, en1, seq);
	return en - st;
}

int mm_idx_getseq_rev(const mm_idx_t *mi, uint32_t rid, uint32_t st, uint32_t en, uint8_t *seq)
{
	uint64_t st1, en1;
	const mm_idx_seq_t *s;
	if (rid >= mi->n_seq || st >= mi->seq[rid].len) return -1;
	s = &mi->seq[rid];
	if (en > s->len) en = s->len;
	st1 = s->offset + (s->len - en);
	en1 = s->offset + (s->len - st);
	seq4_decode_rev(mi->S, st1, en1, seq);
	return en - st;
}

int mm_idx_getseq2(const mm_idx_t *mi, int is_rev, uint32_t rid, uint32_t st, uint32_t en, uint8_t *seq)
{
	if (is_rev) return mm_idx_getseq_rev(mi, rid, st, en, seq);
	else return mm_idx_getseq(mi, rid, st, en, seq);
}

#define MM_OCC_HIST_N 0x10000 // counts below this are tallied in a histogram; larger ones are kept for selection
#define MM_OCC_BUF_N  4096

typedef kvec_t(uint32_t) occ_v;

static void occ_tally(uint64_t *hist, occ_v *large, const uint32_t *c, size_t n)
{
	size_t j;
	for (j = 0; j < n; ++j) {
		if (c[j] < MM_OCC_HIST_N) ++hist[c[j]];
		else kv_push(uint32_t, 0, *large, c[j]);
	}
}

typedef struct { // one thread's tally in mm_idx_cal_max_occ_mt()
	uint64_t *hist;
	occ_v large;
	size_t n1, n_buf;
	uint32_t buf[MM_OCC_BUF_N];
} occ_acc_t;

typedef struct {
	const mm_idx_t *mi;
	occ_acc_t *acc;
} occ_walk_t;

static void occ_bucket(void *data, long i, int tid) // kt_for() callback
{
	occ_walk_t *w = (occ_walk_t*)data;
	occ_acc_t *a = &w->acc[tid];
	idxhash_t *h = (idxhash_t*)w->mi->B[i].h;
	const uint64_t *keys, *vals;
	uint32_t *buf = a->buf;
	size_t n1 = a->n1, n_buf = a->n_buf;
	khint_t k;
	if (h == 0) return;
	keys = h->keys, vals = h->vals;
	if (kh_end(h) < 16) {
		for (k = 0; k < kh_end(h); ++k) {
			if (!kh_exist(h, k)) continue;
			if (keys[k]&1) ++n1;
			else buf[n_buf++] = (uint32_t)vals[k];
		}
	} else {
		for (k = 0; k < kh_end(h); k += 16) { // the 16 slots of one flag word; a branch per slot would mispredict on occupancy
			khint32_t fw = h->flags[k>>4];
			int t;
			for (t = 0; t < 16; ++t) {
				uint32_t e = (fw >> (t<<1) & 3) == 0, s = keys[k + t] & 1;
				n1 += e & s;
				buf[n_buf] = (uint32_t)vals[k + t];
				n_buf += e & (s ^ 1);
			}
			if (n_buf > MM_OCC_BUF_N - 16) occ_tally(a->hist, &a->large, buf, n_buf), n_buf = 0;
		}
	}
	if (n_buf > MM_OCC_BUF_N - 16) occ_tally(a->hist, &a->large, buf, n_buf), n_buf = 0;
	a->n1 = n1, a->n_buf = n_buf;
}

int32_t mm_idx_cal_max_occ(const mm_idx_t *mi, float f)
{
	return mm_idx_cal_max_occ_mt(mi, f, 1);
}

int32_t mm_idx_cal_max_occ_mt(const mm_idx_t *mi, float f, int n_threads)
{
	int i, t;
	size_t n = 0, n1 = 0, kk, cum;
	uint32_t c, thres;
	uint64_t *hist;
	occ_v large = {0,0,0};
	occ_walk_t w;
	if (f <= 0.) return INT32_MAX;
	for (i = 0; i < 1<<mi->b; ++i)
		if (mi->B[i].h) n += kh_size((idxhash_t*)mi->B[i].h);
	if (n == 0) return INT32_MAX;
	if (n_threads < 1) n_threads = 1;
	w.mi = mi;
	w.acc = (occ_acc_t*)calloc(n_threads, sizeof(occ_acc_t));
	for (t = 0; t < n_threads; ++t)
		w.acc[t].hist = (uint64_t*)calloc(MM_OCC_HIST_N, 8);
	kt_for(n_threads, occ_bucket, &w, 1L<<mi->b);
	hist = w.acc[0].hist;
	for (t = 0; t < n_threads; ++t) { // counts and the multiset of large counts do not depend on which thread saw a bucket
		occ_acc_t *a = &w.acc[t];
		occ_tally(a->hist, &a->large, a->buf, a->n_buf);
		n1 += a->n1;
		if (t > 0) {
			for (c = 0; c < MM_OCC_HIST_N; ++c) hist[c] += a->hist[c];
			free(a->hist);
		}
		if (a->large.n > 0) {
			kv_resize(uint32_t, 0, large, large.n + a->large.n);
			memcpy(large.a + large.n, a->large.a, a->large.n * sizeof(uint32_t));
			large.n += a->large.n;
		}
		kfree(0, a->large.a);
	}
	free(w.acc);
	hist[1] += n1;
	kk = (uint32_t)((1. - f) * n);
	if (kk >= n) kk = n - 1; // only for f > 1 or NaN, where the conversion above is undefined
	for (c = 0, cum = 0; c < MM_OCC_HIST_N && cum + hist[c] <= kk; ++c)
		cum += hist[c];
	thres = c < MM_OCC_HIST_N? c + 1 : ks_ksmall_uint32_t(large.n, large.a, kk - cum) + 1;
	free(hist);
	kfree(0, large.a);
	return thres;
}

/*********************************
 * Sort and generate hash tables *
 *********************************/

static void worker_post(void *g, long i, int tid)
{
	int n, n_keys;
	size_t j, start_a, start_p;
	idxhash_t *h;
	mm_idx_t *mi = (mm_idx_t*)g;
	mm_idx_bucket_t *b = &mi->B[i];
	if (b->a.n == 0) return;

	// sort by minimizer
	radix_sort_128x(b->a.a, b->a.a + b->a.n);

	// count and preallocate
	for (j = 1, n = 1, n_keys = 0, b->n = 0; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x>>8 != b->a.a[j-1].x>>8) {
			++n_keys;
			if (n > 1) b->n += n;
			n = 1;
		} else ++n;
	}
	h = kh_init(idx);
	kh_resize(idx, h, n_keys);
	b->p = (uint64_t*)calloc(b->n, 8);

	// create the hash table
	for (j = 1, n = 1, start_a = start_p = 0; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x>>8 != b->a.a[j-1].x>>8) {
			khint_t itr;
			int absent;
			mm128_t *p = &b->a.a[j-1];
			itr = kh_put(idx, h, p->x>>8>>mi->b<<1, &absent);
			assert(absent && j == start_a + n);
			if (n == 1) {
				kh_key(h, itr) |= 1;
				kh_val(h, itr) = p->y;
			} else {
				int k;
				for (k = 0; k < n; ++k)
					b->p[start_p + k] = b->a.a[start_a + k].y;
				radix_sort_64(&b->p[start_p], &b->p[start_p + n]); // sort by position; needed as in-place radix_sort_128x() is not stable
				kh_val(h, itr) = (uint64_t)start_p<<32 | n;
				start_p += n;
			}
			start_a = j, n = 1;
		} else ++n;
	}
	b->h = h;
	assert(b->n == (int32_t)start_p);

	// deallocate and clear b->a
	kfree(0, b->a.a);
	b->a.n = b->a.m = 0, b->a.a = 0;
}
 
static void mm_idx_post(mm_idx_t *mi, int n_threads)
{
	kt_for(n_threads, worker_post, mi, 1<<mi->b);
}

/******************
 * Generate index *
 ******************/

#include <string.h>
#include <zlib.h>
#include "bseq.h"

typedef struct {
	int mini_batch_size;
	uint64_t batch_size, sum_len;
	mm_bseq_file_t *fp;
	mm_idx_t *mi;
} pipeline_t;

typedef struct {
    int n_seq;
	mm_bseq1_t *seq;
	mm128_v a;
} step_t;

static void mm_idx_add(mm_idx_t *mi, int n, const mm128_t *a)
{
	int i, mask = (1<<mi->b) - 1;
	for (i = 0; i < n; ++i) {
		mm128_v *p = &mi->B[a[i].x>>8&mask].a;
		kv_push(mm128_t, 0, *p, a[i]);
	}
}

static void *worker_pipeline(void *shared, int step, void *in)
{
	int i;
    pipeline_t *p = (pipeline_t*)shared;
    if (step == 0) { // step 0: read sequences
        step_t *s;
		if (p->sum_len > p->batch_size) return 0;
        s = (step_t*)calloc(1, sizeof(step_t));
		s->seq = mm_bseq_read(p->fp, p->mini_batch_size, 0, &s->n_seq); // read a mini-batch
		if (s->seq) {
			uint32_t old_m, m;
			assert((uint64_t)p->mi->n_seq + s->n_seq <= UINT32_MAX); // to prevent integer overflow
			// make room for p->mi->seq
			old_m = p->mi->n_seq, m = p->mi->n_seq + s->n_seq;
			kroundup32(m); kroundup32(old_m);
			if (old_m != m)
				p->mi->seq = (mm_idx_seq_t*)krealloc(p->mi->km, p->mi->seq, m * sizeof(mm_idx_seq_t));
			// make room for p->mi->S
			if (!(p->mi->flag & MM_I_NO_SEQ)) {
				uint64_t sum_len, old_max_len, max_len;
				for (i = 0, sum_len = 0; i < s->n_seq; ++i) sum_len += s->seq[i].l_seq;
				old_max_len = (p->sum_len + 7) / 8;
				max_len = (p->sum_len + sum_len + 7) / 8;
				kroundup64(old_max_len); kroundup64(max_len);
				if (old_max_len != max_len) {
					p->mi->S = (uint32_t*)realloc(p->mi->S, max_len * 4);
					memset(&p->mi->S[old_max_len], 0, 4 * (max_len - old_max_len));
				}
			}
			// populate p->mi->seq
			for (i = 0; i < s->n_seq; ++i) {
				mm_idx_seq_t *seq = &p->mi->seq[p->mi->n_seq];
				uint32_t j;
				if (!(p->mi->flag & MM_I_NO_NAME)) {
					seq->name = (char*)kmalloc(p->mi->km, strlen(s->seq[i].name) + 1);
					strcpy(seq->name, s->seq[i].name);
				} else seq->name = 0;
				seq->len = s->seq[i].l_seq;
				seq->offset = p->sum_len;
				seq->is_alt = 0;
				// copy the sequence
				if (!(p->mi->flag & MM_I_NO_SEQ)) {
					for (j = 0; j < seq->len; ++j) { // TODO: this is not the fastest way, but let's first see if speed matters here
						uint64_t o = p->sum_len + j;
						int c = seq_nt4_table[(uint8_t)s->seq[i].seq[j]];
						mm_seq4_set(p->mi->S, o, c);
					}
				}
				// update p->sum_len and p->mi->n_seq
				p->sum_len += seq->len;
				s->seq[i].rid = p->mi->n_seq++;
			}
			return s;
		} else free(s);
    } else if (step == 1) { // step 1: compute sketch
        step_t *s = (step_t*)in;
		for (i = 0; i < s->n_seq; ++i) {
			mm_bseq1_t *t = &s->seq[i];
			if (t->l_seq > 0)
				mm_sketch(0, t->seq, t->l_seq, p->mi->w, p->mi->k, t->rid, p->mi->flag&MM_I_HPC, &s->a);
			else if (mm_verbose >= 2)
				fprintf(stderr, "[WARNING] the length database sequence '%s' is 0\n", t->name);
			free(t->seq); free(t->name);
		}
		free(s->seq); s->seq = 0;
		return s;
    } else if (step == 2) { // dispatch sketch to buckets
        step_t *s = (step_t*)in;
		mm_idx_add(p->mi, s->a.n, s->a.a);
		kfree(0, s->a.a); free(s);
	}
    return 0;
}

mm_idx_t *mm_idx_gen(mm_bseq_file_t *fp, int w, int k, int b, int flag, int mini_batch_size, int n_threads, uint64_t batch_size)
{
	pipeline_t pl;
	if (fp == 0 || mm_bseq_eof(fp)) return 0;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.mini_batch_size = (uint64_t)mini_batch_size < batch_size? mini_batch_size : batch_size;
	pl.batch_size = batch_size;
	pl.fp = fp;
	pl.mi = mm_idx_init(w, k, b, flag);

	kt_pipeline(n_threads < 3? n_threads : 3, worker_pipeline, &pl, 3);
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] collected minimizers\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));

	mm_idx_post(pl.mi, n_threads);
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] sorted minimizers\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));

	return pl.mi;
}

mm_idx_t *mm_idx_build(const char *fn, int w, int k, int flag, int n_threads) // a simpler interface; deprecated
{
	mm_bseq_file_t *fp;
	mm_idx_t *mi;
	fp = mm_bseq_open(fn);
	if (fp == 0) return 0;
	mi = mm_idx_gen(fp, w, k, 14, flag, 1<<18, n_threads, UINT64_MAX);
	mm_bseq_close(fp);
	return mi;
}

mm_idx_t *mm_idx_str(int w, int k, int is_hpc, int bucket_bits, int n, const char **seq, const char **name)
{
	uint64_t sum_len = 0;
	mm128_v a = {0,0,0};
	mm_idx_t *mi;
	khash_t(str) *h;
	int i, flag = 0;

	if (n <= 0) return 0;
	for (i = 0; i < n; ++i) // get the total length
		sum_len += strlen(seq[i]);
	if (is_hpc) flag |= MM_I_HPC;
	if (name == 0) flag |= MM_I_NO_NAME;
	if (bucket_bits < 0) bucket_bits = 14;
	mi = mm_idx_init(w, k, bucket_bits, flag);
	mi->n_seq = n;
	mi->seq = (mm_idx_seq_t*)kcalloc(mi->km, n, sizeof(mm_idx_seq_t)); // ->seq is allocated from km
	mi->S = (uint32_t*)calloc((sum_len + 7) / 8, 4);
	mi->h = h = kh_init(str);
	for (i = 0, sum_len = 0; i < n; ++i) {
		const char *s = seq[i];
		mm_idx_seq_t *p = &mi->seq[i];
		uint32_t j;
		if (name && name[i]) {
			int absent;
			p->name = (char*)kmalloc(mi->km, strlen(name[i]) + 1);
			strcpy(p->name, name[i]);
			kh_put(str, h, p->name, &absent);
			assert(absent);
		}
		p->offset = sum_len;
		p->len = strlen(s);
		p->is_alt = 0;
		for (j = 0; j < p->len; ++j) {
			int c = seq_nt4_table[(uint8_t)s[j]];
			uint64_t o = sum_len + j;
			mm_seq4_set(mi->S, o, c);
		}
		sum_len += p->len;
		if (p->len > 0) {
			a.n = 0;
			mm_sketch(0, s, p->len, w, k, i, is_hpc, &a);
			mm_idx_add(mi, a.n, a.a);
		}
	}
	free(a.a);
	mm_idx_post(mi, 1);
	return mi;
}

/*************
 * index I/O *
 *************/

void mm_idx_dump(FILE *fp, const mm_idx_t *mi)
{
	uint64_t sum_len = 0;
	uint32_t x[5], i;

	x[0] = mi->w, x[1] = mi->k, x[2] = mi->b, x[3] = mi->n_seq, x[4] = mi->flag;
	fwrite(MM_IDX_MAGIC, 1, 4, fp);
	fwrite(x, 4, 5, fp);
	for (i = 0; i < mi->n_seq; ++i) {
		if (mi->seq[i].name) {
			uint8_t l = strlen(mi->seq[i].name);
			fwrite(&l, 1, 1, fp);
			fwrite(mi->seq[i].name, 1, l, fp);
		} else {
			uint8_t l = 0;
			fwrite(&l, 1, 1, fp);
		}
		fwrite(&mi->seq[i].len, 4, 1, fp);
		sum_len += mi->seq[i].len;
	}
	for (i = 0; i < 1<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		khint_t k;
		idxhash_t *h = (idxhash_t*)b->h;
		uint32_t size = h? h->size : 0;
		fwrite(&b->n, 4, 1, fp);
		fwrite(b->p, 8, b->n, fp);
		fwrite(&size, 4, 1, fp);
		if (size == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			uint64_t x[2];
			if (!kh_exist(h, k)) continue;
			x[0] = kh_key(h, k), x[1] = kh_val(h, k);
			fwrite(x, 8, 2, fp);
		}
	}
	if (!(mi->flag & MM_I_NO_SEQ))
		fwrite(mi->S, 4, (sum_len + 7) / 8, fp);
	fflush(fp);
}

/* Loading on several threads: one thread reads the buckets in file order, in
 * chunks, while kt_for builds each bucket's hash table with the kh_resize and
 * kh_put calls of the serial loop, so every table and its slot order is the
 * same as a load on one thread. */

#define MM_IDX_LOAD_CHUNK (1<<25) /* bytes of p[] and hash entries per chunk */

#ifdef MM_IDX_HUGEPAGE
/* Each bucket's p[] and hash table are carved from 256 MB blocks, 2 MB-aligned
 * and advised for huge pages, so the load takes a page fault per 2 MB and
 * lookups miss the TLB less; mm_idx_destroy frees the blocks. */

#define MM_IDX_HUGE  ((size_t)1<<21)
#define MM_IDX_BLOCK ((size_t)1<<28)

static void *mm_idx_huge_alloc(size_t size) /* freed with free() */
{
	void *p;
	if (posix_memalign(&p, MM_IDX_HUGE, size) != 0) {
		fprintf(stderr, "[E::%s] failed to allocate %lu bytes\n", __func__, (unsigned long)size);
		abort();
	}
	madvise(p, size, MADV_HUGEPAGE); /* when refused, the block keeps 4 KB pages */
	return p;
}

static uint8_t *mm_idx_mem_get(mm_idx_mem_t *m, size_t l)
{
	uint8_t *p;
	if (m->cur == 0 || m->l + l > m->cur_m) {
		size_t size = l > MM_IDX_BLOCK? (l + MM_IDX_HUGE - 1) & ~(MM_IDX_HUGE - 1) : MM_IDX_BLOCK;
		if (m->n == m->m) {
			m->m = m->m? m->m<<1 : 16;
			m->a = (void**)realloc(m->a, m->m * sizeof(void*));
		}
		m->a[m->n++] = m->cur = (uint8_t*)mm_idx_huge_alloc(size);
		m->l = 0, m->cur_m = size;
	}
	p = m->cur + m->l;
	m->l += l;
	return p;
}

static inline size_t mm_idx_hp_round(size_t x) { return (x + 63) & ~(size_t)63; }

/* n_buckets of a table after kh_resize(size) and size kh_put calls: kh_put
 * doubles it once when size passes the upper bound */
static khint_t mm_idx_hp_n_buckets(uint32_t size, int *doubled)
{
	khint_t n = size;
	kroundup32(n);
	if (n < 4) n = 4;
	*doubled = (size > (khint_t)(n * __ac_HASH_UPPER + 0.5));
	return *doubled? n<<1 : n;
}

/* a bucket's region: p[], then its table's header, keys, vals and flags */
static size_t mm_idx_hp_bytes(int32_t n, uint32_t size)
{
	size_t l = mm_idx_hp_round((size_t)n * 8);
	if (size > 0) {
		int doubled;
		khint_t nb = mm_idx_hp_n_buckets(size, &doubled);
		l += mm_idx_hp_round(sizeof(idxhash_t)) + 2 * mm_idx_hp_round((size_t)nb * 8) + mm_idx_hp_round(__ac_fsize(nb) * 4);
	}
	return l;
}

/* idxhp is idx with kmalloc, krealloc and kfree taken from the bucket's
 * region (the macros reach the build through khash's argument h), so
 * kh_resize and kh_put lay each table out exactly as for idx */
typedef struct mm_idx_hpbuild_s mm_idx_hpbuild_t;
static void *mm_idx_hp_malloc(mm_idx_hpbuild_t *t, size_t size);
static void *mm_idx_hp_realloc(mm_idx_hpbuild_t *t, void *p, size_t size);
static void mm_idx_hp_free(mm_idx_hpbuild_t *t, void *p);
#define kmalloc(km, size) mm_idx_hp_malloc((mm_idx_hpbuild_t*)h, (size))
#define krealloc(km, p, size) mm_idx_hp_realloc((mm_idx_hpbuild_t*)h, (p), (size))
#define kfree(km, p) mm_idx_hp_free((mm_idx_hpbuild_t*)h, (p))
KHASH_INIT(idxhp, uint64_t, uint64_t, 1, idx_hash, idx_eq)
#undef kmalloc
#undef krealloc
#undef kfree

struct mm_idx_hpbuild_s {
	khash_t(idxhp) h;   /* first, so a pointer to h is one to the build */
	uint64_t *kv[2];    /* keys and vals, kv_cap bytes each, in the order asked for */
	size_t kv_cap, flags_cap;
	int n_kv, n_flags;  /* flags arrays still to come; only the last is the table's */
	uint32_t *flags;
};

static void *mm_idx_hp_malloc(mm_idx_hpbuild_t *t, size_t size) /* a flags array */
{
	if (--t->n_flags > 0) return malloc(size); /* replaced when kh_put doubles the table */
	assert(t->n_flags == 0 && size <= t->flags_cap);
	return t->flags;
}

static void *mm_idx_hp_realloc(mm_idx_hpbuild_t *t, void *p, size_t size) /* keys or vals, grown in place */
{
	assert(size <= t->kv_cap);
	if (p) return p;
	assert(t->n_kv < 2);
	return t->kv[t->n_kv++];
}

static void mm_idx_hp_free(mm_idx_hpbuild_t *t, void *p)
{
	if (p != t->flags) free(p);
}
#endif

typedef struct {
	mm_idx_t *mi;
	uint32_t st, en;  /* buckets [st, en) */
	size_t l, m;      /* bytes used and allocated in buf */
	uint8_t *buf;     /* per bucket: its p[] (b->n*8 bytes), then its hash entries (size*16 bytes); or l bytes of mi->S */
	uint8_t *dst;     /* where in mi->S the l bytes go, or 0 for a chunk of buckets */
	size_t *off;      /* per bucket: where its p[] starts in buf */
	uint32_t *size;   /* per bucket: number of hash entries */
	uint8_t **place;  /* per bucket: its region in mi->mem's blocks (MM_IDX_HUGEPAGE) */
} idx_load_chunk_t;

typedef struct {
	FILE *fp;
	mm_idx_t *mi;
	size_t seq_len, seq_read; /* bytes of mi->S, and bytes read */
	int n_workers;
	const mm_prefetch_t *pf; /* a query read running beside the load: one worker fewer until it ends */
	uint32_t i;       /* next bucket to read */
	int64_t k;        /* chunks read so far */
	int64_t off, end; /* file offset of the next byte to read, and the file's length */
	int bad;          /* set when a field lies past the end of the file, a bucket's n is negative or a read comes up short */
	idx_load_chunk_t c[2]; /* kt_pipeline with two threads holds at most two chunks */
} idx_load_t;

static void idx_load_reserve(idx_load_chunk_t *c, size_t l)
{
	if (c->l + l > c->m) {
		c->m = c->l + l > c->m * 2? c->l + l : c->m * 2;
		c->buf = (uint8_t*)realloc(c->buf, c->m);
	}
}

#define MM_IDX_LOAD_COPY (1<<20) /* bytes of mi->S per copy job */

static void idx_load_copy(void *data, long i, int tid)
{
	idx_load_chunk_t *c = (idx_load_chunk_t*)data;
	size_t st = (size_t)i * MM_IDX_LOAD_COPY, l = c->l - st < MM_IDX_LOAD_COPY? c->l - st : MM_IDX_LOAD_COPY;
	memcpy(c->dst + st, c->buf + st, l);
}

/* builds bucket b's hash table from its size entries x (key, value pairs as
 * in the file); b->n and b->p are set, and with MM_IDX_HUGEPAGE the table goes
 * in b's region, right after p[] */
static void idx_load_table(mm_idx_bucket_t *b, uint32_t size, const uint64_t *x)
{
	uint32_t j;
	idxhash_t *h;
#ifdef MM_IDX_HUGEPAGE
	uint8_t *m = (uint8_t*)b->p + mm_idx_hp_round((size_t)b->n * 8);
	mm_idx_hpbuild_t t;
	int doubled;
	if (size == 0) return;
	memset(&t, 0, sizeof(t));
	h = (idxhash_t*)m, m += mm_idx_hp_round(sizeof(idxhash_t));
	t.kv_cap = (size_t)mm_idx_hp_n_buckets(size, &doubled) * 8;
	t.kv[0] = (uint64_t*)m, m += mm_idx_hp_round(t.kv_cap);
	t.kv[1] = (uint64_t*)m, m += mm_idx_hp_round(t.kv_cap);
	t.flags = (uint32_t*)m, t.flags_cap = __ac_fsize(t.kv_cap / 8) * 4;
	t.n_flags = doubled? 2 : 1;
	kh_resize(idxhp, &t.h, size);
	for (j = 0; j < size; ++j) {
		khint_t k;
		int absent;
		k = kh_put(idxhp, &t.h, x[j<<1], &absent);
		assert(absent);
		kh_val(&t.h, k) = x[j<<1|1];
	}
	assert(t.h.flags == t.flags);
	h->n_buckets = t.h.n_buckets, h->size = t.h.size, h->n_occupied = t.h.n_occupied, h->upper_bound = t.h.upper_bound;
	h->flags = t.h.flags, h->keys = t.h.keys, h->vals = t.h.vals;
#else
	if (size == 0) return;
	h = kh_init(idx);
	kh_resize(idx, h, size);
	for (j = 0; j < size; ++j) {
		khint_t k;
		int absent;
		k = kh_put(idx, h, x[j<<1], &absent);
		assert(absent);
		kh_val(h, k) = x[j<<1|1];
	}
#endif
	b->h = h;
}

static void idx_load_alloc(mm_idx_t *mi, size_t seq_len) /* mi->S, and mi->mem with MM_IDX_HUGEPAGE */
{
	if (!(mi->flag & MM_I_NO_SEQ)) {
#ifdef MM_IDX_HUGEPAGE
		mi->S = (uint32_t*)mm_idx_huge_alloc(seq_len);
#else
		mi->S = (uint32_t*)malloc(seq_len);
#endif
	}
#ifdef MM_IDX_HUGEPAGE
	mi->mem = calloc(1, sizeof(mm_idx_mem_t));
#endif
}

static void idx_load_undo(mm_idx_t *mi) /* frees what idx_load_alloc() and a parallel load gave mi */
{
	uint32_t i;
	for (i = 0; i < 1U<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		if (mi->mem == 0) {
			free(b->p);
			kh_destroy(idx, (idxhash_t*)b->h);
		}
		b->n = 0, b->p = 0, b->h = 0;
	}
	mm_idx_mem_destroy((mm_idx_mem_t*)mi->mem);
	free(mi->S);
	mi->mem = 0, mi->S = 0;
}

static void idx_load_worker(void *data, long i, int tid)
{
	idx_load_chunk_t *c = (idx_load_chunk_t*)data;
	mm_idx_bucket_t *b = &c->mi->B[c->st + i];
	const uint8_t *q = c->buf + c->off[i];
#ifdef MM_IDX_HUGEPAGE
	b->p = (uint64_t*)c->place[i];
#else
	b->p = (uint64_t*)malloc((size_t)b->n * 8);
#endif
	if (b->n > 0) memcpy(b->p, q, (size_t)b->n * 8);
	idx_load_table(b, c->size[i], (const uint64_t*)(q + (size_t)b->n * 8));
}

/* 0 when all n items were read */
static int idx_fread(void *buf, size_t size, size_t n, FILE *fp)
{
	return fread(buf, size, n, fp) == n? 0 : -1;
}

/* Each field's bytes are checked against the file's length before anything
 * is allocated for it. A field past the end, a negative n or a short read
 * sets s->bad and ends the pipeline; idx_load_mt() then frees what was built
 * and hands the part to the serial load. */
static void *idx_load_pipeline(void *shared, int step, void *in)
{
	idx_load_t *s = (idx_load_t*)shared;
	mm_idx_t *mi = s->mi;
	if (step == 0) {
		uint32_t n_b = 1U<<mi->b;
		idx_load_chunk_t *c;
		if (s->bad) return 0;
		if (s->i == n_b) { /* every bucket read: read mi->S in chunks, which the workers copy into place */
			size_t l;
			if (s->seq_read == s->seq_len) return 0;
			if (s->off + (int64_t)(s->seq_len - s->seq_read) > s->end) {
				s->bad = 1;
				return 0;
			}
			c = &s->c[s->k++ & 1];
			c->st = c->en = n_b, c->l = 0;
			idx_load_reserve(c, s->seq_len - s->seq_read < MM_IDX_LOAD_CHUNK? s->seq_len - s->seq_read : MM_IDX_LOAD_CHUNK);
			l = c->m < s->seq_len - s->seq_read? c->m : s->seq_len - s->seq_read;
			if (idx_fread(c->buf, 1, l, s->fp) != 0) {
				s->bad = 1;
				return 0;
			}
			c->l = l, c->dst = (uint8_t*)mi->S + s->seq_read;
			s->seq_read += l, s->off += l;
			return c;
		}
		c = &s->c[s->k++ & 1];
		c->st = s->i, c->l = 0, c->dst = 0;
		if (c->off == 0) {
			c->off = (size_t*)malloc(n_b * sizeof(size_t));
			c->size = (uint32_t*)malloc(n_b * sizeof(uint32_t));
#ifdef MM_IDX_HUGEPAGE
			c->place = (uint8_t**)malloc(n_b * sizeof(uint8_t*));
#endif
		}
		while (s->i < n_b && c->l < MM_IDX_LOAD_CHUNK) {
			mm_idx_bucket_t *b = &mi->B[s->i];
			uint32_t j = s->i - c->st, size;
			int32_t n;
			if (idx_fread(&n, 4, 1, s->fp) != 0 || n < 0 || s->off + 8 + (int64_t)n * 8 > s->end) {
				s->bad = 1;
				return 0;
			}
			idx_load_reserve(c, (size_t)n * 8);
			c->off[j] = c->l;
			if (idx_fread(c->buf + c->l, 8, n, s->fp) != 0 || idx_fread(&size, 4, 1, s->fp) != 0) {
				s->bad = 1;
				return 0;
			}
			c->l += (size_t)n * 8, s->off += 8 + (int64_t)n * 8;
			if (s->off + (int64_t)size * 16 > s->end) {
				s->bad = 1;
				return 0;
			}
			idx_load_reserve(c, (size_t)size * 16);
			if (idx_fread(c->buf + c->l, 16, size, s->fp) != 0) {
				s->bad = 1;
				return 0;
			}
			c->l += (size_t)size * 16, s->off += (int64_t)size * 16;
			b->n = n, c->size[j] = size;
#ifdef MM_IDX_HUGEPAGE
			c->place[j] = mm_idx_mem_get((mm_idx_mem_t*)mi->mem, mm_idx_hp_bytes(n, size));
#endif
			++s->i;
		}
		c->en = s->i;
		return c;
	} else if (step == 1) {
		idx_load_chunk_t *c = (idx_load_chunk_t*)in;
		int n_workers = s->n_workers - mm_prefetch_busy(s->pf);
		if (n_workers < 1) n_workers = 1;
		if (c->dst) kt_for(n_workers, idx_load_copy, c, (c->l + MM_IDX_LOAD_COPY - 1) / MM_IDX_LOAD_COPY);
		else kt_for(n_workers, idx_load_worker, c, c->en - c->st);
	}
	return 0;
}

#ifdef MM_IDX_PREAD
/* Loading with pread: one thread walks the bucket headers (each bucket is a
 * 4-byte n, n*8 bytes of p[], a 4-byte size and size*16 bytes of hash
 * entries), then kt_for reads each bucket's p[] into place and its entries
 * into a per-thread buffer, and builds its table as idx_load_worker does. */

#define MM_IDX_LOAD_SEQ (1<<22) /* bytes of mi->S per job */

typedef struct {
	mm_idx_t *mi;
	int fd;
	uint32_t n_b;
	int64_t *off;     /* per bucket: file offset of its p[] */
	int32_t *n;       /* per bucket: length of p[] */
	uint32_t *size;   /* per bucket: number of hash entries, which start 4 bytes after p[] */
	uint8_t **place;  /* per bucket: its region in mi->mem's blocks (MM_IDX_HUGEPAGE) */
	int64_t s_off;    /* file offset of mi->S */
	size_t seq_len;
	uint8_t **buf;    /* per thread: one bucket's hash entries */
	size_t *buf_m;
} idx_pload_t;

static int idx_pread(int fd, void *buf, size_t l, int64_t off) /* 0 when all l bytes were read */
{
	uint8_t *p = (uint8_t*)buf;
	while (l > 0) {
		ssize_t r = pread(fd, p, l < (size_t)1<<30? l : (size_t)1<<30, (off_t)off);
		if (r < 0 && errno == EINTR) continue;
		if (r <= 0) return -1;
		p += r, l -= r, off += r;
	}
	return 0;
}

/* the offset where the part ends, or -1 when a field is negative or the file
 * ends first; the caller then hands the part to the serial load */
static int64_t idx_pload_walk(idx_pload_t *s, int64_t off, int64_t file_size)
{
	uint32_t i;
	for (i = 0; i < s->n_b; ++i) {
		int32_t n;
		uint32_t size;
		if (off + 4 > file_size || idx_pread(s->fd, &n, 4, off) < 0 || n < 0) return -1;
		s->n[i] = n, s->off[i] = off + 4;
		off += 4 + (int64_t)n * 8;
		if (off + 4 > file_size || idx_pread(s->fd, &size, 4, off) < 0) return -1;
		s->size[i] = size;
		off += 4 + (int64_t)size * 16;
	}
	s->s_off = off;
	off += s->seq_len;
	return off > file_size? -1 : off;
}

static void idx_pload_worker(void *data, long i, int tid)
{
	idx_pload_t *s = (idx_pload_t*)data;
	mm_idx_bucket_t *b;
	size_t l;
	if (i >= s->n_b) { /* a slice of mi->S */
		size_t st = (size_t)(i - s->n_b) * MM_IDX_LOAD_SEQ;
		l = s->seq_len - st < MM_IDX_LOAD_SEQ? s->seq_len - st : MM_IDX_LOAD_SEQ;
		idx_pread(s->fd, (uint8_t*)s->mi->S + st, l, s->s_off + st);
		return;
	}
	b = &s->mi->B[i];
#ifdef MM_IDX_HUGEPAGE
	b->p = (uint64_t*)s->place[i];
#else
	b->p = (uint64_t*)malloc((size_t)b->n * 8);
#endif
	if (b->n > 0) idx_pread(s->fd, b->p, (size_t)b->n * 8, s->off[i]);
	if (s->size[i] == 0) return;
	l = (size_t)s->size[i] * 16;
	if (l > s->buf_m[tid]) {
		s->buf_m[tid] = l;
		free(s->buf[tid]);
		s->buf[tid] = (uint8_t*)malloc(l);
	}
	idx_pread(s->fd, s->buf[tid], l, s->off[i] + (int64_t)b->n * 8 + 4);
	idx_load_table(b, s->size[i], (const uint64_t*)s->buf[tid]);
}

/* loads the buckets and mi->S of the part at fp's position and leaves fp at
 * the part's end (1); with fp and mi untouched, 0 when fp is not a regular
 * file and -1 when the walk finds the part does not fit in it */
static int idx_pload(FILE *fp, mm_idx_t *mi, int n_threads, size_t seq_len)
{
	idx_pload_t s;
	struct stat st;
	int64_t off, end = -1;
	uint32_t i;
	memset(&s, 0, sizeof(s));
	s.fd = fileno(fp);
	if (s.fd < 0 || fstat(s.fd, &st) != 0 || !S_ISREG(st.st_mode)) return 0;
	if ((off = ftello(fp)) < 0) return 0;
	s.mi = mi, s.n_b = 1U<<mi->b, s.seq_len = seq_len;
	s.off = (int64_t*)malloc(s.n_b * sizeof(int64_t));
	s.n = (int32_t*)malloc(s.n_b * sizeof(int32_t));
	s.size = (uint32_t*)malloc(s.n_b * sizeof(uint32_t));
	end = idx_pload_walk(&s, off, st.st_size);
	if (end >= 0) {
		idx_load_alloc(mi, seq_len);
#ifdef MM_IDX_HUGEPAGE
		s.place = (uint8_t**)malloc(s.n_b * sizeof(uint8_t*));
#endif
		for (i = 0; i < s.n_b; ++i) {
			mi->B[i].n = s.n[i];
#ifdef MM_IDX_HUGEPAGE
			s.place[i] = mm_idx_mem_get((mm_idx_mem_t*)mi->mem, mm_idx_hp_bytes(s.n[i], s.size[i]));
#endif
		}
		s.buf = (uint8_t**)calloc(n_threads, sizeof(uint8_t*));
		s.buf_m = (size_t*)calloc(n_threads, sizeof(size_t));
		kt_for(n_threads, idx_pload_worker, &s, (long)s.n_b + (long)((seq_len + MM_IDX_LOAD_SEQ - 1) / MM_IDX_LOAD_SEQ));
		for (i = 0; i < (uint32_t)n_threads; ++i) free(s.buf[i]);
		free(s.buf); free(s.buf_m); free(s.place);
		fseeko(fp, (off_t)end, SEEK_SET);
	}
	free(s.off); free(s.n); free(s.size);
	return end >= 0? 1 : -1;
}
#endif

/* loads the buckets and mi->S of the part at fp's position on n_threads (1);
 * 0, with fp at the first bucket and B[], S and mem empty, when a field lies
 * past the end of the file, a bucket's n is negative or a read comes up short:
 * the serial load then reads the part, so a damaged file gives the same result
 * at every thread count */
static int idx_load_mt(FILE *fp, mm_idx_t *mi, int n_threads, const mm_prefetch_t *pf, uint64_t sum_len)
{
	idx_load_t s;
	uint32_t i;
	int64_t st;
	size_t seq_len = mi->flag & MM_I_NO_SEQ? 0 : (sum_len + 7) / 8 * 4;
#ifdef MM_IDX_PREAD
	int ret = idx_pload(fp, mi, n_threads - mm_prefetch_busy(pf), seq_len); /* one thread fewer while a query read runs beside the load */
	if (ret != 0) return ret > 0;
#endif
	memset(&s, 0, sizeof(s));
	s.fp = fp, s.mi = mi, s.n_workers = n_threads - 1, s.pf = pf, s.seq_len = seq_len;
	if ((st = s.off = idx_ftell(fp)) < 0) return 0;
	if (idx_fseek(fp, 0, SEEK_END) != 0 || (s.end = idx_ftell(fp)) < 0) {
		idx_fseek(fp, st, SEEK_SET);
		return 0;
	}
	idx_fseek(fp, st, SEEK_SET);
	idx_load_alloc(mi, seq_len);
	for (i = 0; i < 2; ++i) s.c[i].mi = mi;
	kt_pipeline(2, idx_load_pipeline, &s, 2);
	for (i = 0; i < 2; ++i) {
		free(s.c[i].buf); free(s.c[i].off); free(s.c[i].size); free(s.c[i].place);
	}
	if (s.bad) {
		idx_load_undo(mi);
		idx_fseek(fp, st, SEEK_SET);
		return 0;
	}
	return 1;
}

static mm_idx_t *mm_idx_load_core(FILE *fp, int n_threads, const mm_prefetch_t *pf)
{
	char magic[4];
	uint32_t x[5], i, size = 0;
	uint64_t sum_len = 0, kv[2] = {0, 0};
	mm_idx_t *mi;

	if (fread(magic, 1, 4, fp) != 4) return 0;
	if (strncmp(magic, MM_IDX_MAGIC, 4) != 0) return 0;
	if (fread(x, 4, 5, fp) != 5) return 0;
	mi = mm_idx_init(x[0], x[1], x[2], x[4]);
	mi->n_seq = x[3];
	mi->seq = (mm_idx_seq_t*)kcalloc(mi->km, mi->n_seq, sizeof(mm_idx_seq_t));
	for (i = 0; i < mi->n_seq; ++i) {
		uint8_t l;
		mm_idx_seq_t *s = &mi->seq[i];
		fread(&l, 1, 1, fp);
		if (l) {
			s->name = (char*)kmalloc(mi->km, l + 1);
			fread(s->name, 1, l, fp);
			s->name[l] = 0;
		}
		fread(&s->len, 4, 1, fp);
		s->offset = sum_len;
		s->is_alt = 0;
		sum_len += s->len;
	}
	if (n_threads > 1 && idx_load_mt(fp, mi, n_threads, pf, sum_len)) return mi;
	/* Serial load. size and kv[] live across buckets: after a short read they
	 * keep the last values read, which decides what a damaged file loads or
	 * how it fails; a file that ends before the first size is an error. */
	for (i = 0; i < 1<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		uint32_t j;
		khint_t k;
		idxhash_t *h;
		fread(&b->n, 4, 1, fp);
		b->p = (uint64_t*)malloc(b->n * 8);
		fread(b->p, 8, b->n, fp);
		if (fread(&size, 4, 1, fp) != 1 && i == 0) {
			fprintf(stderr, "[E::%s] the index file is truncated\n", __func__);
			abort();
		}
		if (size == 0) continue;
		b->h = h = kh_init(idx);
		kh_resize(idx, h, size);
		for (j = 0; j < size; ++j) {
			int absent;
			fread(kv, 8, 2, fp);
			k = kh_put(idx, h, kv[0], &absent);
			assert(absent);
			kh_val(h, k) = kv[1];
		}
	}
	if (!(mi->flag & MM_I_NO_SEQ)) {
		mi->S = (uint32_t*)malloc((sum_len + 7) / 8 * 4);
		fread(mi->S, 4, (sum_len + 7) / 8, fp);
	}
	return mi;
}

mm_idx_t *mm_idx_load(FILE *fp)
{
	return mm_idx_load_core(fp, 1, 0);
}

int64_t mm_idx_is_idx(const char *fn)
{
	int fd, is_idx = 0;
	int64_t ret, off_end;
	char magic[4];

	if (strcmp(fn, "-") == 0) return 0; // read from pipe; not an index
	fd = open(fn, O_RDONLY);
	if (fd < 0) return -1; // error
#ifdef WIN32
	if ((off_end = _lseeki64(fd, 0, SEEK_END)) >= 4) {
		_lseeki64(fd, 0, SEEK_SET);
#else
	if ((off_end = lseek(fd, 0, SEEK_END)) >= 4) {
		lseek(fd, 0, SEEK_SET);
#endif // WIN32
		ret = read(fd, magic, 4);
		if (ret == 4 && strncmp(magic, MM_IDX_MAGIC, 4) == 0)
			is_idx = 1;
	}
	close(fd);
	return is_idx? off_end : 0;
}

mm_idx_reader_t *mm_idx_reader_open(const char *fn, const mm_idxopt_t *opt, const char *fn_out)
{
	int64_t is_idx;
	mm_idx_reader_t *r;
	is_idx = mm_idx_is_idx(fn);
	if (is_idx < 0) return 0; // failed to open the index
	r = (mm_idx_reader_t*)calloc(1, sizeof(mm_idx_reader_t));
	r->is_idx = is_idx;
	if (opt) r->opt = *opt;
	else mm_idxopt_init(&r->opt);
	if (r->is_idx) {
		r->fp.idx = fopen(fn, "rb");
		r->idx_size = is_idx;
	} else r->fp.seq = mm_bseq_open(fn);
	if (fn_out) r->fp_out = fopen(fn_out, "wb");
	return r;
}

void mm_idx_reader_close(mm_idx_reader_t *r)
{
	if (r->is_idx) fclose(r->fp.idx);
	else mm_bseq_close(r->fp.seq);
	if (r->fp_out) fclose(r->fp_out);
	free(r);
}

mm_idx_t *mm_idx_reader_read(mm_idx_reader_t *r, int n_threads)
{
	return mm_idx_reader_read_pf(r, n_threads, 0);
}

mm_idx_t *mm_idx_reader_read_pf(mm_idx_reader_t *r, int n_threads, const mm_prefetch_t *pf)
{
	mm_idx_t *mi;
	if (r->is_idx) {
		mi = mm_idx_load_core(r->fp.idx, n_threads, pf);
		if (mi && mm_verbose >= 2 && (mi->k != r->opt.k || mi->w != r->opt.w || (mi->flag&MM_I_HPC) != (r->opt.flag&MM_I_HPC)))
			fprintf(stderr, "[WARNING]\033[1;31m Indexing parameters (-k, -w or -H) overridden by parameters used in the prebuilt index.\033[0m\n");
	} else
		mi = mm_idx_gen(r->fp.seq, r->opt.w, r->opt.k, r->opt.bucket_bits, r->opt.flag, r->opt.mini_batch_size, n_threads, r->opt.batch_size);
	if (mi) {
		if (r->fp_out) mm_idx_dump(r->fp_out, mi);
		mi->index = r->n_parts++;
	}
	return mi;
}

int mm_idx_reader_eof(const mm_idx_reader_t *r) // TODO: in extremely rare cases, mm_bseq_eof() might not work
{
	return r->is_idx? (feof(r->fp.idx) || ftell(r->fp.idx) == r->idx_size) : mm_bseq_eof(r->fp.seq);
}

#include <ctype.h>
#include <zlib.h>
#include "ksort.h"
#include "kseq.h"
KSTREAM_DECLARE(gzFile, gzread)

int mm_idx_alt_read(mm_idx_t *mi, const char *fn)
{
	int n_alt = 0;
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return -1;
	ks = ks_init(fp);
	if (mi->h == 0) mm_idx_index_name(mi);
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
		char *p;
		int id;
		for (p = str.s; *p && !isspace(*p); ++p) { }
		*p = 0;
		id = mm_idx_name2id(mi, str.s);
		if (id >= 0) mi->seq[id].is_alt = 1, ++n_alt;
	}
	mi->n_alt = n_alt;
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s] found %d ALT contigs\n", __func__, n_alt);
	return n_alt;
}

/***************
 * BED reading *
 ***************/

#define sort_key_bed(a) ((a).st)
KRADIX_SORT_INIT(bed, mm_idx_intv1_t, sort_key_bed, 4)

#define sort_key_end(a) ((a).en)
KRADIX_SORT_INIT(end, mm_idx_intv1_t, sort_key_end, 4)

static mm_idx_intv_t *mm_idx_bed_read_core(const mm_idx_t *mi, const char *fn, int read_junc, int min_sc)
{
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	mm_idx_intv_t *I;

	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	I = CALLOC(mm_idx_intv_t, mi->n_seq);
	ks = ks_init(fp);
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
		mm_idx_intv_t *r;
		mm_idx_intv1_t t = {-1,-1,-1,-1,0};
		char *p, *q, *bl, *bs;
		int32_t i, id = -1, n_blk = 0;
		for (p = q = str.s, i = 0;; ++p) {
			if (*p == 0 || *p == '\t') {
				int32_t c = *p;
				*p = 0;
				if (i == 0) { // chr
					id = mm_idx_name2id(mi, q);
					if (id < 0) break; // unknown name; TODO: throw a warning
				} else if (i == 1) { // start
					t.st = atol(q); // TODO: watch out integer overflow!
					if (t.st < 0) break;
				} else if (i == 2) { // end
					t.en = atol(q);
					if (t.en < 0) break;
				} else if (i == 4) { // BED score
					t.score = *q >= '0' && *q <= '9'? atol(q) : -1;
				} else if (i == 5) { // strand
					t.strand = *q == '+'? 1 : *q == '-'? -1 : 0;
				} else if (i == 9) {
					if (!isdigit(*q)) break;
					n_blk = atol(q);
				} else if (i == 10) {
					bl = q;
				} else if (i == 11) {
					bs = q;
					break;
				}
				if (c == 0) break;
				++i, q = p + 1;
			}
		}
		if (id < 0 || t.st < 0 || t.st >= t.en) continue; // contig ID not found, or other problems
		if (min_sc > 0 && t.score < min_sc) continue;
		r = &I[id];
		if (i >= 11 && read_junc) { // BED12
			int32_t st, sz, en;
			st = strtol(bs, &bs, 10); ++bs;
			sz = strtol(bl, &bl, 10); ++bl;
			en = t.st + st + sz;
			for (i = 1; i < n_blk; ++i) {
				mm_idx_intv1_t s = t;
				if (r->n == r->m) {
					r->m = r->m? r->m + (r->m>>1) : 16;
					r->a = (mm_idx_intv1_t*)realloc(r->a, sizeof(*r->a) * r->m);
				}
				st = strtol(bs, &bs, 10); ++bs;
				sz = strtol(bl, &bl, 10); ++bl;
				s.st = en, s.en = t.st + st;
				en = t.st + st + sz;
				if (s.en > s.st) r->a[r->n++] = s;
			}
		} else {
			if (r->n == r->m) {
				r->m = r->m? r->m + (r->m>>1) : 16;
				r->a = (mm_idx_intv1_t*)realloc(r->a, sizeof(*r->a) * r->m);
			}
			r->a[r->n++] = t;
		}
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);
	return I;
}

static mm_idx_intv_t *mm_idx_bed_read_merge(const mm_idx_t *mi, const char *fn, int read_junc, int min_sc)
{
	long n = 0, n0 = 0;
	int32_t i;
	mm_idx_intv_t *I;
	I = mm_idx_bed_read_core(mi, fn, read_junc, min_sc);
	if (I == 0) return 0;
	for (i = 0; i < mi->n_seq; ++i) {
		int32_t j, j0, k;
		mm_idx_intv_t *intv = &I[i];
		n0 += intv->n;
		radix_sort_bed(intv->a, intv->a + intv->n); // sort by st
		for (j = 1, j0 = 0; j <= intv->n; ++j) { // sort by st and then by end
			if (j == intv->n || intv->a[j].st != intv->a[j0].st) {
				radix_sort_end(intv->a + j0, intv->a + j);
				j0 = j;
			}
		}
		for (j = 1, j0 = 0, k = 0; j <= intv->n; ++j) { // merge intervals with the same (st, en)
			if (j == intv->n || intv->a[j].st != intv->a[j0].st || intv->a[j].en != intv->a[j0].en) {
				intv->a[k] = intv->a[j0];
				intv->a[k++].cnt = j - j0;
				j0 = j;
			}
		}
		intv->a = REALLOC(mm_idx_intv1_t, intv->a, k);
		intv->n = intv->m = k;
		n += k;
	}
	if (mm_verbose >= 3)
		fprintf(stderr, "[%s] read %ld introns, %ld of which are non-redundant\n", __func__, n0, n);
	return I;
}

int mm_idx_bed_read(mm_idx_t *mi, const char *fn, int read_junc)
{
	if (mi->h == 0) mm_idx_index_name(mi);
	mi->I = mm_idx_bed_read_merge(mi, fn, read_junc, -1);
	return 0;
}

int mm_idx_bed_junc(const mm_idx_t *mi, int32_t ctg, int32_t st, int32_t en, uint8_t *s)
{
	int32_t i, left, right;
	mm_idx_intv_t *r;
	memset(s, 0, en - st);
	if (mi->I == 0 || ctg < 0 || ctg >= mi->n_seq) return -1;
	r = &mi->I[ctg];
	left = 0, right = r->n;
	while (right > left) {
		int32_t mid = left + ((right - left) >> 1);
		if (r->a[mid].st >= st) right = mid;
		else left = mid + 1;
	}
	for (i = left; i < r->n; ++i) {
		if (st <= r->a[i].st && en >= r->a[i].en && r->a[i].strand != 0) {
			if (r->a[i].strand > 0) {
				s[r->a[i].st - st] |= 1, s[r->a[i].en - 1 - st] |= 2;
			} else {
				s[r->a[i].st - st] |= 8, s[r->a[i].en - 1 - st] |= 4;
			}
		}
	}
	return left;
}

/*********************************
 * Reading junctions for jumping *
 *********************************/

#define sort_key_jj(a) ((a).off)
KRADIX_SORT_INIT(jj, mm_idx_jjump1_t, sort_key_jj, 4)

#define sort_key_jj2(a) ((a).off2)
KRADIX_SORT_INIT(jj2, mm_idx_jjump1_t, sort_key_jj2, 4)

static void sort_jjump(mm_idx_jjump_t *jj2)
{
	int32_t j0, j, k;
	if (jj2 == 0 || jj2->n == 0) return;
	radix_sort_jj(jj2->a, jj2->a + jj2->n);
	for (j0 = 0, j = 1; j <= jj2->n; ++j) {
		if (j == jj2->n || jj2->a[j0].off != jj2->a[j].off) {
			radix_sort_jj2(jj2->a + j0, jj2->a + j);
			j0 = j;
		}
	}
	// the actual merge
	for (j0 = 0, j = 1, k = 0; j <= jj2->n; ++j) {
		if (j == jj2->n || jj2->a[j0].off != jj2->a[j].off || jj2->a[j0].off2 != jj2->a[j].off2) {
			int32_t t, cnt = 0;
			uint16_t flag = 0;
			for (t = j0; t < j; ++t) cnt += jj2->a[t].cnt, flag |= jj2->a[t].flag;
			jj2->a[k] = jj2->a[j0];
			jj2->a[k].cnt = cnt;
			jj2->a[k++].flag = flag;
			j0 = j;
		}
	}
	jj2->n = k;
	jj2->a = REALLOC(mm_idx_jjump1_t, jj2->a, k);
}

static mm_idx_jjump_t *mm_idx_bed2jjump(const mm_idx_t *mi, const mm_idx_intv_t *I, uint16_t flag)
{
	int32_t i;
	mm_idx_jjump_t *J;
	J = CALLOC(mm_idx_jjump_t, mi->n_seq);
	for (i = 0; i < mi->n_seq; ++i) {
		int32_t j, k;
		const mm_idx_intv_t *intv = &I[i];
		mm_idx_jjump_t *jj = &J[i];
		jj->n = intv->n * 2;
		jj->a = CALLOC(mm_idx_jjump1_t, jj->n);
		for (j = k = 0; j < intv->n; ++j) {
			jj->a[k].off = intv->a[j].st, jj->a[k].off2 = intv->a[j].en, jj->a[k].cnt = intv->a[j].cnt, jj->a[k].strand = intv->a[j].strand, jj->a[k++].flag = flag;
			jj->a[k].off = intv->a[j].en, jj->a[k].off2 = intv->a[j].st, jj->a[k].cnt = intv->a[j].cnt, jj->a[k].strand = intv->a[j].strand, jj->a[k++].flag = flag;
		}
		sort_jjump(jj);
	}
	return J;
}

static mm_idx_jjump_t *mm_idx_jjump_merge(const mm_idx_t *mi, const mm_idx_jjump_t *J0, const mm_idx_jjump_t *J1)
{
	int32_t i;
	mm_idx_jjump_t *J2;
	J2 = CALLOC(mm_idx_jjump_t, mi->n_seq);
	for (i = 0; i < mi->n_seq; ++i) {
		int32_t j, k;
		const mm_idx_jjump_t *jj0 = &J0[i], *jj1 = &J1[i];
		mm_idx_jjump_t *jj2 = &J2[i];
		jj2->n = jj0->n + jj1->n;
		jj2->a = CALLOC(mm_idx_jjump1_t, jj2->n);
		for (j = k = 0; j < jj0->n; ++j) jj2->a[k++] = jj0->a[j];
		for (j = 0; j < jj1->n; ++j) jj2->a[k++] = jj1->a[j];
		sort_jjump(jj2);
	}
	return J2;
}

int mm_idx_jjump_read(mm_idx_t *mi, const char *fn, int flag, int min_sc)
{
	int32_t i, j, n_anno = 0, n_misc = 0;
	mm_idx_intv_t *I;
	mm_idx_jjump_t *J;
	if (mi->h == 0) mm_idx_index_name(mi);
	I = mm_idx_bed_read_merge(mi, fn, 1, min_sc);
	J = mm_idx_bed2jjump(mi, I, flag);
	for (i = 0; i < mi->n_seq; ++i) free(I[i].a);
	free(I);
	if (mi->J) {
		mm_idx_jjump_t *J2;
		J2 = mm_idx_jjump_merge(mi, mi->J, J);
		for (i = 0; i < mi->n_seq; ++i) {
			free(mi->J[i].a); free(J[i].a);
		}
		free(mi->J); free(J);
		mi->J = J2;
	} else mi->J = J;
	for (i = 0; i < mi->n_seq; ++i) {
		for (j = 0; j < mi->J[i].n; ++j)
			if (mi->J[i].a[j].flag & MM_JUNC_ANNO) ++n_anno;
			else ++n_misc;
	}
	if (mm_verbose >= 3)
		fprintf(stderr, "[%s] there are %d annotated and %d other splice positions in the index\n", __func__, n_anno, n_misc);
	return 0;
}

static int32_t mm_idx_jump_get_core(int32_t n, const mm_idx_jjump1_t *a, int32_t x) // similar to mm_idx_find_intv()
{
	int32_t s = 0, e = n;
	if (n == 0) return -1;
	if (x < a[0].off) return -1;
	while (s < e) {
		int32_t mid = s + (e - s) / 2;
		if (x >= a[mid].off && (mid + 1 >= n || x < a[mid+1].off)) return mid;
		else if (x < a[mid].off) e = mid;
		else s = mid + 1;
	}
	assert(0);
}

const mm_idx_jjump1_t *mm_idx_jump_get(const mm_idx_t *db, int32_t cid, int32_t st, int32_t en, int32_t *n)
{
	mm_idx_jjump_t *s;
	int32_t l, r;
	*n = 0;
	if (cid >= db->n_seq || cid < 0 || db->J == 0) return 0;
	if (en < 0 || en > db->seq[cid].len) en = db->seq[cid].len;
	s = &db->J[cid];
	if (s->n == 0) return 0;
	l = mm_idx_jump_get_core(s->n, s->a, st);
	r = mm_idx_jump_get_core(s->n, s->a, en);
	*n = r - l;
	return &s->a[l + 1];
}

/****************
 * splice score *
 ****************/

typedef struct mm_idx_spsc_s {
	uint32_t n, m;
	uint64_t *a; // pos<<56 | score<<1 | acceptor
} mm_idx_spsc_t;

int32_t mm_idx_spsc_read2(mm_idx_t *idx, const char *fn, int32_t max_sc, float scale)
{
	gzFile fp;
	kstring_t str = {0,0,0};
	kstream_t *ks;
	int32_t dret, j;
	int64_t n_read = 0;

	fp = fn && strcmp(fn, "-") != 0? gzopen(fn, "rb") : gzdopen(0, "rb");
	if (fp == 0) return -1;
	if (idx->h == 0) mm_idx_index_name(idx);
	if (max_sc > 63) max_sc = 63;
	idx->spsc = Kcalloc(0, mm_idx_spsc_t, idx->n_seq * 2);
	ks = ks_init(fp);
	while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
		mm_idx_spsc_t *s;
		char *p, *q, *name = 0;
		int32_t i, type = -1, strand = 0, cid = -1, score = -1;
		int64_t pos = -1;
		for (i = 0, p = q = str.s;; ++p) {
			if (*p == '\t' || *p == 0) {
				int c = *p;
				*p = 0;
				if (i == 0) {
					name = q;
				} else if (i == 1) {
					pos = atol(q);
				} else if (i == 2) {
					strand = *q == '+'? 1 : '-'? -1 : 0;
				} else if (i == 3) {
					type = *q == 'D'? 0 : *q == 'A'? 1 : -1;
				} else if (i == 4) {
					score = atoi(q);
					break;
				}
				if (c == 0) break;
				q = p + 1, ++i;
			}
		}
		if (i < 4) continue; // not enough fields
		if (scale > 0.0f && scale < 1.0f)
			score = score > 0.0f? (int)(score * scale + .499) : (int)(score * scale - .499);
		if (score > max_sc) score = max_sc;
		if (score < -max_sc) score = -max_sc;
		cid = mm_idx_name2id(idx, name);
		if (cid < 0 || type < 0 || strand == 0 || pos < 0) continue; // FIXME: give a warning!
		s = &idx->spsc[cid << 1 | (strand > 0? 0 : 1)];
		Kgrow(0, uint64_t, s->a, s->n, s->m);
		if (pos > 0 && pos < idx->seq[cid].len) { // ignore scores at the ends
			s->a[s->n++] = (uint64_t)pos << 8 | (score + KSW_SPSC_OFFSET) << 1 | type;
			++n_read;
		}
	}
	ks_destroy(ks);
	gzclose(fp);
	for (j = 0; j < idx->n_seq * 2; ++j) {
		mm_idx_spsc_t *s = &idx->spsc[j];
		if (s->n > 0)
			radix_sort_64(s->a, s->a + s->n);
	}
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s] read %ld splice scores\n", __func__, (long)n_read);
	return 0;
}

int32_t mm_idx_spsc_read(mm_idx_t *idx, const char *fn, int32_t max_sc)
{
	return mm_idx_spsc_read2(idx, fn, max_sc, 1.0f);
}

static int32_t mm_idx_find_intv(int32_t n, const uint64_t *a, int64_t x)
{
	int32_t s = 0, e = n;
	if (n == 0) return -1;
	if (x < a[0]>>8) return -1;
	while (s < e) {
		int32_t mid = s + (e - s) / 2;
		if (x >= a[mid]>>8 && (mid + 1 >= n || x < a[mid+1]>>8)) return mid;
		else if (x < a[mid]>>8) e = mid;
		else s = mid + 1;
	}
	assert(0);
}

int64_t mm_idx_spsc_get(const mm_idx_t *db, int32_t cid, int64_t st, int64_t en, int32_t rev, uint8_t *sc)
{
	const mm_idx_spsc_t *s;
	if (cid >= db->n_seq || cid < 0 || db->spsc == 0) return -1;
	if (en < 0 || en > db->seq[cid].len) en = db->seq[cid].len;
	memset(sc, 0xff, en - st);
	s = &db->spsc[cid << 1 | (!!rev)];
	if (s->n > 0) {
		int32_t j, l, r;
		l = mm_idx_find_intv(s->n, s->a, st);
		r = mm_idx_find_intv(s->n, s->a, en);
		for (j = l + 1; j <= r; ++j) {
			int64_t x = (s->a[j]>>8) - st;
			uint8_t score = s->a[j] & 0xff;
			assert(x <= en - st);
			if (x == en - st) continue;
			if (sc[x] == 0xff || sc[x] < score) sc[x] = score;
		}
	}
	return en - st;
}
