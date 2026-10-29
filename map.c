#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <pthread.h>
#include "kthread.h"
#include "kvec.h"
#include "kalloc.h"
#include "sdust.h"
#include "mmpriv.h"
#include "bseq.h"
#include "khash.h"

mm_tbuf_t *mm_tbuf_init(void)
{
	mm_tbuf_t *b;
	b = (mm_tbuf_t*)calloc(1, sizeof(mm_tbuf_t));
	if (!(mm_dbg_flag & 1)) b->km = km_init();
	return b;
}

void mm_tbuf_destroy(mm_tbuf_t *b)
{
	if (b == 0) return;
	km_destroy(b->km);
	free(b);
}

void *mm_tbuf_get_km(mm_tbuf_t *b)
{
	return b->km;
}

static int mm_dust_minier(void *km, int n, mm128_t *a, int l_seq, const char *seq, int sdust_thres)
{
	int n_dreg, j, k, u = 0;
	const uint64_t *dreg;
	sdust_buf_t *sdb;
	if (sdust_thres <= 0) return n;
	sdb = sdust_buf_init(km);
	dreg = sdust_core((const uint8_t*)seq, l_seq, sdust_thres, 64, &n_dreg, sdb);
	for (j = k = 0; j < n; ++j) { // squeeze out minimizers that significantly overlap with LCRs
		int32_t qpos = (uint32_t)a[j].y>>1, span = a[j].x&0xff;
		int32_t s = qpos - (span - 1), e = s + span;
		while (u < n_dreg && (int32_t)dreg[u] <= s) ++u;
		if (u < n_dreg && (int32_t)(dreg[u]>>32) < e) {
			int v, l = 0;
			for (v = u; v < n_dreg && (int32_t)(dreg[v]>>32) < e; ++v) { // iterate over LCRs overlapping this minimizer
				int ss = s > (int32_t)(dreg[v]>>32)? s : dreg[v]>>32;
				int ee = e < (int32_t)dreg[v]? e : (uint32_t)dreg[v];
				l += ee - ss;
			}
			if (l <= span>>1) a[k++] = a[j]; // keep the minimizer if less than half of it falls in masked region
		} else a[k++] = a[j];
	}
	sdust_buf_destroy(sdb);
	return k; // the new size
}

static void collect_minimizers(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, mm128_v *mv)
{
	int i, n, sum = 0;
	mv->n = 0;
	for (i = n = 0; i < n_segs; ++i) {
		size_t j;
		mm_sketch(km, seqs[i], qlens[i], mi->w, mi->k, i, mi->flag&MM_I_HPC, mv);
		for (j = n; j < mv->n; ++j)
			mv->a[j].y += sum << 1;
		if (opt->sdust_thres > 0) // mask low-complexity minimizers
			mv->n = n + mm_dust_minier(km, mv->n - n, mv->a + n, qlens[i], seqs[i], opt->sdust_thres);
		sum += qlens[i], n = mv->n;
	}
}

#include "ksort.h"
#define heap_lt(a, b) ((a).x > (b).x)
KSORT_INIT(heap, mm128_t, heap_lt)

static inline int skip_seed(int flag, uint64_t r, const mm_seed_t *q, const char *qname, int qlen, const mm_idx_t *mi, int *is_self)
{
	*is_self = 0;
	if (qname && (flag & (MM_F_NO_DIAG|MM_F_NO_DUAL))) {
		const mm_idx_seq_t *s = &mi->seq[r>>32];
		int cmp;
		cmp = strcmp(qname, s->name);
		if ((flag&MM_F_NO_DIAG) && cmp == 0 && (int)s->len == qlen) {
			if ((uint32_t)r>>1 == (q->q_pos>>1)) return 1; // avoid the diagnonal anchors
			if ((r&1) == (q->q_pos&1)) *is_self = 1; // this flag is used to avoid spurious extension on self chain
		}
		if ((flag&MM_F_NO_DUAL) && cmp > 0) // all-vs-all mode: map once
			return 1;
	}
	if (flag & (MM_F_FOR_ONLY|MM_F_REV_ONLY)) {
		if ((r&1) == (q->q_pos&1)) { // forward strand
			if (flag & MM_F_REV_ONLY) return 1;
		} else {
			if (flag & MM_F_FOR_ONLY) return 1;
		}
	}
	return 0;
}

static mm128_t *collect_seed_hits_heap(void *km, const mm_mapopt_t *opt, int max_occ, const mm_idx_t *mi, const char *qname, const mm128_v *mv, int qlen, int64_t *n_a, int *rep_len,
								  int *n_mini_pos, uint64_t **mini_pos)
{
	int i, n_m, heap_size = 0;
	int64_t j, n_for = 0, n_rev = 0;
	mm_seed_t *m;
	mm128_t *a, *heap;

	m = mm_collect_matches(km, &n_m, qlen, max_occ, opt->max_max_occ, opt->occ_dist, mi, mv, n_a, rep_len, n_mini_pos, mini_pos);

	heap = (mm128_t*)kmalloc(km, n_m * sizeof(mm128_t));
	a = (mm128_t*)kmalloc(km, *n_a * sizeof(mm128_t));

	for (i = 0; i < n_m && i < MM_SEED_PF_DIST; ++i)
		mm_prefetch(m[i].cr);
	for (i = 0, heap_size = 0; i < n_m; ++i) {
		if (i + MM_SEED_PF_DIST < n_m) mm_prefetch(m[i + MM_SEED_PF_DIST].cr);
		if (m[i].n > 0) {
			heap[heap_size].x = m[i].cr[0];
			heap[heap_size].y = (uint64_t)i<<32;
			++heap_size;
		}
	}
	ks_heapmake_heap(heap_size, heap);
	while (heap_size > 0) {
		mm_seed_t *q = &m[heap->y>>32];
		mm128_t *p;
		uint64_t r = heap->x;
		int32_t is_self, rpos = (uint32_t)r >> 1;
		if (!skip_seed(opt->flag, r, q, qname, qlen, mi, &is_self)) {
			if ((r&1) == (q->q_pos&1)) { // forward strand
				p = &a[n_for++];
				p->x = (r&0xffffffff00000000ULL) | rpos;
				p->y = (uint64_t)q->q_span << 32 | q->q_pos >> 1;
			} else { // reverse strand
				p = &a[(*n_a) - (++n_rev)];
				p->x = 1ULL<<63 | (r&0xffffffff00000000ULL) | rpos;
				p->y = (uint64_t)q->q_span << 32 | (qlen - ((q->q_pos>>1) + 1 - q->q_span) - 1);
			}
			p->y |= (uint64_t)q->seg_id << MM_SEED_SEG_SHIFT;
			if (q->is_tandem) p->y |= MM_SEED_TANDEM;
			if (is_self) p->y |= MM_SEED_SELF;
		}
		// update the heap
		if ((uint32_t)heap->y < q->n - 1) {
			++heap[0].y;
			heap[0].x = m[heap[0].y>>32].cr[(uint32_t)heap[0].y];
		} else {
			heap[0] = heap[heap_size - 1];
			--heap_size;
		}
		ks_heapdown_heap(0, heap_size, heap);
	}
	kfree(km, m);
	kfree(km, heap);

	// reverse anchors on the reverse strand, as they are in the descending order
	for (j = 0; j < n_rev>>1; ++j) {
		mm128_t t = a[(*n_a) - 1 - j];
		a[(*n_a) - 1 - j] = a[(*n_a) - (n_rev - j)];
		a[(*n_a) - (n_rev - j)] = t;
	}
	if (*n_a > n_for + n_rev) {
		memmove(a + n_for, a + (*n_a) - n_rev, n_rev * sizeof(mm128_t));
		*n_a = n_for + n_rev;
	}
	return a;
}

// whether two kept seeds on one query strand parity share a position run. One k-mer ends at each reference position,
// so without such a pair (and outside --qstrand) the anchors of collect_seed_hits have distinct keys.
static int seeds_share_run(void *km, int n_m, const mm_seed_t *m)
{
	size_t i, h, mask, size = 16;
	uint64_t *tab;
	if (n_m < 2) return 0;
	while (size < (size_t)n_m * 2) size <<= 1;
	mask = size - 1;
	tab = (uint64_t*)kcalloc(km, size, sizeof(uint64_t));
	for (i = 0; i < (size_t)n_m; ++i) {
		uint64_t x = (uint64_t)(uintptr_t)m[i].cr | (m[i].q_pos & 1); // cr is non-null and 8-byte aligned
		for (h = (size_t)((x * 0x9E3779B97F4A7C15ULL) >> 32) & mask; tab[h] && tab[h] != x; h = (h + 1) & mask) {}
		if (tab[h] == x) break;
		tab[h] = x;
	}
	kfree(km, tab);
	return i < (size_t)n_m;
}

static inline int mm_bit_len(uint64_t v)
{
	int n = 0;
	for (; v; v >>= 1) ++n;
	return n;
}

// sort a[0..n) by x with an LSD radix sort; equals radix_sort_128x() when all keys are distinct. The digits cover only
// the bits that differ between keys: bit 63, then the low bits of bits 32-62, then the low bits of bits 0-31.
static void radix_sort_128x_distinct(void *km, mm128_t *a, int64_t n)
{
	uint64_t x_or = 0, x_and = ~0ULL, v, pm, rm, sm, k;
	uint32_t *cnt, dm;
	int64_t i;
	int pb, rb, nb, np, d, p, dmax;
	mm128_t *b, *src, *dst, *t;
	if (n <= 64 || n > (int64_t)UINT32_MAX) {
		radix_sort_128x(a, a + n);
		return;
	}
	for (i = 0; i < n; ++i) x_or |= a[i].x, x_and &= a[i].x;
	v = x_or ^ x_and;
	if (v == 0) return;
	sm = v >> 63, rb = mm_bit_len(v << 1 >> 33), pb = mm_bit_len((uint32_t)v);
	nb = (int)sm + rb + pb;
	dmax = n < 256? 8 : n < 4096? 10 : 13;
	np = (nb + dmax - 1) / dmax, d = (nb + np - 1) / np, dm = (1U << d) - 1;
	pm = (1ULL << pb) - 1, rm = (1ULL << rb) - 1;
	b = (mm128_t*)kmalloc(km, n * sizeof(mm128_t) + ((size_t)np << d) * sizeof(uint32_t));
	cnt = (uint32_t*)(b + n);
	memset(cnt, 0, ((size_t)np << d) * sizeof(uint32_t));
#define LSD_KEY(x) (((x) & pm) | ((x) >> 32 & rm) << pb | ((x) >> 63 & sm) << (pb + rb))
	for (i = 0; i < n; ++i) {
		k = LSD_KEY(a[i].x);
		for (p = 0; p < np; ++p)
			++cnt[(size_t)p << d | (k >> p * d & dm)];
	}
	for (p = 0; p < np; ++p) {
		uint32_t j, s = 0, c, *q = cnt + ((size_t)p << d);
		for (j = 0; j <= dm; ++j)
			c = q[j], q[j] = s, s += c;
	}
	for (p = 0, src = a, dst = b; p < np; ++p) {
		uint32_t *q = cnt + ((size_t)p << d);
		int sh = p * d;
		for (i = 0; i < n; ++i) {
			k = LSD_KEY(src[i].x);
			dst[q[k >> sh & dm]++] = src[i];
		}
		t = src, src = dst, dst = t;
	}
#undef LSD_KEY
	if (src != a) memcpy(a, src, n * sizeof(mm128_t));
	kfree(km, b);
}

static mm128_t *collect_seed_hits(void *km, const mm_mapopt_t *opt, int max_occ, const mm_idx_t *mi, const char *qname, const mm128_v *mv, int qlen, int64_t *n_a, int *rep_len,
								  int *n_mini_pos, uint64_t **mini_pos, int *distinct)
{
	int i, n_m;
	mm_seed_t *m;
	mm128_t *a;
	m = mm_collect_matches(km, &n_m, qlen, max_occ, opt->max_max_occ, opt->occ_dist, mi, mv, n_a, rep_len, n_mini_pos, mini_pos);
	*distinct = *n_a > 64 && !(opt->flag & MM_F_QSTRAND) && !seeds_share_run(km, n_m, m);
	a = (mm128_t*)kmalloc(km, *n_a * sizeof(mm128_t));
	for (i = 0; i < n_m && i < MM_SEED_PF_DIST; ++i)
		mm_prefetch(m[i].cr);
	for (i = 0, *n_a = 0; i < n_m; ++i) {
		mm_seed_t *q = &m[i];
		const uint64_t *r = q->cr;
		uint32_t k;
		if (i + MM_SEED_PF_DIST < n_m) mm_prefetch(m[i + MM_SEED_PF_DIST].cr);
		for (k = 0; k < q->n; ++k) {
			int32_t is_self, rpos = (uint32_t)r[k] >> 1;
			mm128_t *p;
			if (skip_seed(opt->flag, r[k], q, qname, qlen, mi, &is_self)) continue;
			p = &a[(*n_a)++];
			if ((r[k]&1) == (q->q_pos&1)) { // forward strand
				p->x = (r[k]&0xffffffff00000000ULL) | rpos;
				p->y = (uint64_t)q->q_span << 32 | q->q_pos >> 1;
			} else if (!(opt->flag & MM_F_QSTRAND)) { // reverse strand and not in the query-strand mode
				p->x = 1ULL<<63 | (r[k]&0xffffffff00000000ULL) | rpos;
				p->y = (uint64_t)q->q_span << 32 | (qlen - ((q->q_pos>>1) + 1 - q->q_span) - 1);
			} else { // reverse strand; query-strand
				int32_t len = mi->seq[r[k]>>32].len;
				p->x = 1ULL<<63 | (r[k]&0xffffffff00000000ULL) | (len - (rpos + 1 - q->q_span) - 1); // coordinate only accurate for non-HPC seeds
				p->y = (uint64_t)q->q_span << 32 | q->q_pos >> 1;
			}
			p->y |= (uint64_t)q->seg_id << MM_SEED_SEG_SHIFT;
			if (q->is_tandem) p->y |= MM_SEED_TANDEM;
			if (is_self) p->y |= MM_SEED_SELF;
		}
	}
	kfree(km, m);
	if (*distinct) radix_sort_128x_distinct(km, a, *n_a);
	else radix_sort_128x(a, a + (*n_a));
	return a;
}

static void chain_post(const mm_mapopt_t *opt, int max_chain_gap_ref, const mm_idx_t *mi, void *km, int qlen, int n_segs, const int *qlens, int *n_regs, mm_reg1_t *regs, mm128_t *a)
{
	if (!(opt->flag & MM_F_ALL_CHAINS)) { // don't choose primary mapping(s)
		mm_set_parent(km, opt->mask_level, opt->mask_len, *n_regs, regs, opt->a * 2 + opt->b, opt->flag&MM_F_HARD_MLEVEL, opt->alt_drop);
		if (n_segs <= 1) mm_select_sub(km, opt->pri_ratio, mi->k*2, opt->best_n, 1, opt->max_gap * 0.8, n_regs, regs);
		else mm_select_sub_multi(km, opt->pri_ratio, 0.2f, 0.7f, max_chain_gap_ref, mi->k*2, opt->best_n, n_segs, qlens, n_regs, regs);
	}
}

static mm_reg1_t *align_regs(const mm_mapopt_t *opt, const mm_idx_t *mi, void *km, int qlen, const char *seq, int *n_regs, mm_reg1_t *regs, mm128_t *a)
{
	if (!(opt->flag & MM_F_CIGAR)) return regs;
	regs = mm_align_skeleton(km, opt, mi, qlen, seq, n_regs, regs, a); // this calls mm_filter_regs()
	if (!(opt->flag & MM_F_ALL_CHAINS)) { // don't choose primary mapping(s)
		mm_set_parent(km, opt->mask_level, opt->mask_len, *n_regs, regs, opt->a * 2 + opt->b, opt->flag&MM_F_HARD_MLEVEL, opt->alt_drop);
		mm_select_sub(km, opt->pri_ratio, mi->k*2, opt->best_n, 0, opt->max_gap * 0.8, n_regs, regs);
		mm_set_sam_pri(*n_regs, regs);
	}
	return regs;
}

void mm_map_frag_core(const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, int *n_regs, mm_reg1_t **regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *qname)
{
	int i, j, rep_len, qlen_sum, n_regs0, n_mini_pos, distinct = 0;
	int max_chain_gap_qry, max_chain_gap_ref, is_splice = !!(opt->flag & MM_F_SPLICE), is_sr = !!(opt->flag & MM_F_SR), is_sr_rna = !!(opt->flag & MM_F_SR_RNA);
	uint32_t hash;
	int64_t n_a;
	uint64_t *u, *mini_pos;
	mm128_t *a;
	mm128_v mv = {0,0,0};
	mm_reg1_t *regs0;
	km_stat_t kmst;
	float chn_pen_gap, chn_pen_skip;

	for (i = 0, qlen_sum = 0; i < n_segs; ++i)
		qlen_sum += qlens[i], n_regs[i] = 0, regs[i] = 0;

	if (qlen_sum == 0 || n_segs <= 0 || n_segs > MM_MAX_SEG) return;
	if (opt->max_qlen > 0 && qlen_sum > opt->max_qlen) return;

	hash  = qname && !(opt->flag & MM_F_NO_HASH_NAME)? __ac_X31_hash_string(qname) : 0;
	hash ^= __ac_Wang_hash(qlen_sum) + __ac_Wang_hash(opt->seed);
	hash  = __ac_Wang_hash(hash);

	collect_minimizers(b->km, opt, mi, n_segs, qlens, seqs, &mv);
	if (opt->q_occ_frac > 0.0f) mm_seed_mz_flt(b->km, &mv, opt->mid_occ, opt->q_occ_frac);
	if (opt->flag & MM_F_HEAP_SORT) a = collect_seed_hits_heap(b->km, opt, opt->mid_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
	else a = collect_seed_hits(b->km, opt, opt->mid_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos, &distinct);

	if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
		fprintf(stderr, "RS\t%d\n", rep_len);
		for (i = 0; i < n_a; ++i)
			fprintf(stderr, "SD\t%s\t%d\t%c\t%d\t%d\t%d\n", mi->seq[a[i].x<<1>>33].name, (int32_t)a[i].x, "+-"[a[i].x>>63], (int32_t)a[i].y, (int32_t)(a[i].y>>32&0xff),
					i == 0? 0 : ((int32_t)a[i].y - (int32_t)a[i-1].y) - ((int32_t)a[i].x - (int32_t)a[i-1].x));
	}

	// set max chaining gap on the query and the reference sequence
	if (is_sr)
		max_chain_gap_qry = qlen_sum > opt->max_gap? qlen_sum : opt->max_gap;
	else max_chain_gap_qry = opt->max_gap;
	if (opt->max_gap_ref > 0) {
		max_chain_gap_ref = opt->max_gap_ref; // always honor mm_mapopt_t::max_gap_ref if set
	} else if (opt->max_frag_len > 0) {
		max_chain_gap_ref = opt->max_frag_len - qlen_sum;
		if (max_chain_gap_ref < opt->max_gap) max_chain_gap_ref = opt->max_gap;
	} else max_chain_gap_ref = opt->max_gap;

	chn_pen_gap  = opt->chain_gap_scale * 0.01 * mi->k;
	chn_pen_skip = opt->chain_skip_scale * 0.01 * mi->k;
	if (opt->flag & MM_F_RMQ) {
		a = mg_lchain_rmq(opt->max_gap, opt->rmq_inner_dist, opt->bw, opt->max_chain_skip, opt->rmq_size_cap, opt->min_cnt, opt->min_chain_score,
						  chn_pen_gap, chn_pen_skip, n_a, a, &n_regs0, &u, b->km);
	} else {
		a = mg_lchain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->max_chain_iter, opt->min_cnt, opt->min_chain_score,
						 chn_pen_gap, chn_pen_skip, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km);
	}

	if (opt->bw_long > opt->bw && (opt->flag & (MM_F_SPLICE|MM_F_SR|MM_F_NO_LJOIN)) == 0 && n_segs == 1 && n_regs0 > 1) { // re-chain/long-join for long sequences
		int32_t st = (int32_t)a[0].y, en = (int32_t)a[(int32_t)u[0] - 1].y;
		if (qlen_sum - (en - st) > opt->rmq_rescue_size || en - st > qlen_sum * opt->rmq_rescue_ratio) {
			int32_t i;
			for (i = 0, n_a = 0; i < n_regs0; ++i) n_a += (int32_t)u[i];
			kfree(b->km, u);
			if (distinct) radix_sort_128x_distinct(b->km, a, n_a); // a subset of the first collection's anchors
			else radix_sort_128x(a, a + n_a);
			a = mg_lchain_rmq(opt->max_gap, opt->rmq_inner_dist, opt->bw_long, opt->max_chain_skip, opt->rmq_size_cap, opt->min_cnt, opt->min_chain_score,
							  chn_pen_gap, chn_pen_skip, n_a, a, &n_regs0, &u, b->km);
		}
	} else if (opt->max_occ > opt->mid_occ && rep_len > 0 && !(opt->flag & MM_F_RMQ)) { // re-chain, mostly for short reads
		int rechain = 0;
		if (n_regs0 > 0) { // test if the best chain has all the segments
			int n_chained_segs = 1, max = 0, max_i = -1, max_off = -1, off = 0;
			for (i = 0; i < n_regs0; ++i) { // find the best chain
				if (max < (int)(u[i]>>32)) max = u[i]>>32, max_i = i, max_off = off;
				off += (uint32_t)u[i];
			}
			for (i = 1; i < (int32_t)u[max_i]; ++i) // count the number of segments in the best chain
				if ((a[max_off+i].y&MM_SEED_SEG_MASK) != (a[max_off+i-1].y&MM_SEED_SEG_MASK))
					++n_chained_segs;
			if (n_chained_segs < n_segs)
				rechain = 1;
		} else rechain = 1;
		if (rechain) { // redo chaining with a higher max_occ threshold
			kfree(b->km, a);
			kfree(b->km, u);
			kfree(b->km, mini_pos);
			if (opt->flag & MM_F_HEAP_SORT) a = collect_seed_hits_heap(b->km, opt, opt->max_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
			else a = collect_seed_hits(b->km, opt, opt->max_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos, &distinct);
			a = mg_lchain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->max_chain_iter, opt->min_cnt, opt->min_chain_score,
							 chn_pen_gap, chn_pen_skip, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km);
		}
	}
	b->frag_gap = max_chain_gap_ref;
	b->rep_len = rep_len;

	regs0 = mm_gen_regs(b->km, hash, qlen_sum, n_regs0, u, a, !!(opt->flag&MM_F_QSTRAND));
	if (mi->n_alt) {
		mm_mark_alt(mi, n_regs0, regs0);
		mm_hit_sort(b->km, &n_regs0, regs0, opt->alt_drop); // this step can be merged into mm_gen_regs(); will do if this shows up in profile
	}

	if (mm_dbg_flag & (MM_DBG_PRINT_SEED|MM_DBG_PRINT_CHAIN))
		for (j = 0; j < n_regs0; ++j)
			for (i = regs0[j].as; i < regs0[j].as + regs0[j].cnt; ++i)
				fprintf(stderr, "CN\t%d\t%s\t%d\t%c\t%d\t%d\t%d\n", j, mi->seq[a[i].x<<1>>33].name, (int32_t)a[i].x, "+-"[a[i].x>>63], (int32_t)a[i].y, (int32_t)(a[i].y>>32&0xff),
						i == regs0[j].as? 0 : ((int32_t)a[i].y - (int32_t)a[i-1].y) - ((int32_t)a[i].x - (int32_t)a[i-1].x));

	chain_post(opt, max_chain_gap_ref, mi, b->km, qlen_sum, n_segs, qlens, &n_regs0, regs0, a);
	if (!is_sr && !(opt->flag&MM_F_QSTRAND)) {
		mm_est_err(mi, qlen_sum, n_regs0, regs0, a, n_mini_pos, mini_pos);
		n_regs0 = mm_filter_strand_retained(n_regs0, regs0);
	}

	if (n_segs == 1) { // uni-segment
		regs0 = align_regs(opt, mi, b->km, qlens[0], seqs[0], &n_regs0, regs0, a);
		regs0 = (mm_reg1_t*)realloc(regs0, sizeof(*regs0) * n_regs0);
		mm_set_mapq2(b->km, n_regs0, regs0, opt->min_chain_score, opt->a, rep_len, is_sr || is_sr_rna, is_splice);
		n_regs[0] = n_regs0, regs[0] = regs0;
	} else { // multi-segment
		mm_seg_t *seg;
		seg = mm_seg_gen(b->km, hash, n_segs, qlens, n_regs0, regs0, n_regs, regs, a); // split fragment chain to separate segment chains
		free(regs0);
		for (i = 0; i < n_segs; ++i) {
			mm_set_parent(b->km, opt->mask_level, opt->mask_len, n_regs[i], regs[i], opt->a * 2 + opt->b, opt->flag&MM_F_HARD_MLEVEL, opt->alt_drop); // update mm_reg1_t::parent
			regs[i] = align_regs(opt, mi, b->km, qlens[i], seqs[i], &n_regs[i], regs[i], seg[i].a);
			mm_set_mapq2(b->km, n_regs[i], regs[i], opt->min_chain_score, opt->a, rep_len, is_sr || is_sr_rna, is_splice);
		}
		mm_seg_free(b->km, n_segs, seg);
		if (n_segs == 2 && opt->pe_ori >= 0 && (opt->flag&MM_F_CIGAR))
			mm_pair(b->km, max_chain_gap_ref, opt->pe_bonus, opt->a * 2 + opt->b, opt->a, qlens, n_regs, regs); // pairing
	}

	kfree(b->km, mv.a);
	kfree(b->km, a);
	kfree(b->km, u);
	kfree(b->km, mini_pos);

	if (mi->J && n_segs == 1 && is_splice)
		for (i = 0; i < n_regs0; ++i)
			mm_jump_split(b->km, mi, opt, qlens[0], (const uint8_t*)seqs[0], &regs0[i], 0);

	if (b->km) {
		km_stat(b->km, &kmst);
		if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
			fprintf(stderr, "QM\t%s\t%d\tcap=%ld,nCore=%ld,largest=%ld\n", qname, qlen_sum, kmst.capacity, kmst.n_cores, kmst.largest);
		assert(kmst.n_blocks == kmst.n_cores); // otherwise, there is a memory leak
		if (kmst.largest > 1U<<28 || (opt->cap_kalloc > 0 && kmst.capacity > opt->cap_kalloc)) {
			if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
				fprintf(stderr, "[W::%s] reset thread-local memory after read %s\n", __func__, qname);
			km_destroy(b->km);
			b->km = km_init();
		}
	}
}

void mm_map_frag(const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, int *n_regs, mm_reg1_t **regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *qname)
{
	if ((opt->flag & MM_F_WEAK_PAIRING) && n_segs == 2 && opt->pe_ori >= 0 && (opt->flag&MM_F_CIGAR)) {
		int i;
		for (i = 0; i < n_segs; ++i)
			mm_map_frag_core(mi, 1, &qlens[i], &seqs[i], &n_regs[i], &regs[i], b, opt, qname);
		mm_pair(b->km, opt->max_gap_ref, opt->pe_bonus, opt->a * 2 + opt->b, opt->a, qlens, n_regs, regs);
	} else {
		mm_map_frag_core(mi, n_segs, qlens, seqs, n_regs, regs, b, opt, qname);
	}
}

mm_reg1_t *mm_map(const mm_idx_t *mi, int qlen, const char *seq, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *qname)
{
	mm_reg1_t *regs;
	mm_map_frag(mi, 1, &qlen, &seq, n_regs, &regs, b, opt, qname);
	return regs;
}

/**************************
 * Multi-threaded mapping *
 **************************/

typedef struct {
	int n_processed, n_threads, n_fp;
	int64_t mini_batch_size;
	const mm_mapopt_t *opt;
	mm_bseq_file_t **fp;
	const mm_idx_t *mi;
	kstring_t str;

	int n_parts;
	uint32_t *rid_shift;
	FILE *fp_split, **fp_parts;

	int has_first; // step 0 returns first once instead of reading (a batch read by mm_prefetch_start)
	void *first;
} pipeline_t;

typedef struct {
	const pipeline_t *p;
    int n_seq, n_frag;
	mm_bseq1_t *seq;
	int *n_reg, *seg_off, *n_seg, *rep_len, *frag_gap;
	mm_reg1_t **reg;
	mm_tbuf_t **buf;
	uint64_t *order; // step 1 at -t > 1: kt_for's i-th job is fragment (uint32_t)order[i]; 0 for input order
} step_t;

static void worker_for(void *_data, long i, int tid) // kt_for() callback
{
    step_t *s = (step_t*)_data;
	int qlens[MM_MAX_SEG], j, off, pe_ori = s->p->opt->pe_ori;
	const char *qseqs[MM_MAX_SEG];
	double t = 0.0;
	mm_tbuf_t *b = s->buf[tid];
	if (s->order) i = (uint32_t)s->order[i];
	off = s->seg_off[i];
	assert(s->n_seg[i] <= MM_MAX_SEG);
	if (mm_dbg_flag & MM_DBG_PRINT_QNAME) {
		fprintf(stderr, "QR\t%s\t%d\t%d\n", s->seq[off].name, tid, s->seq[off].l_seq);
		t = realtime();
	}
	for (j = 0; j < s->n_seg[i]; ++j) {
		if (s->n_seg[i] == 2 && ((j == 0 && (pe_ori>>1&1)) || (j == 1 && (pe_ori&1))))
			mm_revcomp_bseq(&s->seq[off + j]);
		qlens[j] = s->seq[off + j].l_seq;
		qseqs[j] = s->seq[off + j].seq;
	}
	if (s->p->opt->flag & MM_F_INDEPEND_SEG) {
		for (j = 0; j < s->n_seg[i]; ++j) {
			mm_map_frag(s->p->mi, 1, &qlens[j], &qseqs[j], &s->n_reg[off+j], &s->reg[off+j], b, s->p->opt, s->seq[off+j].name);
			s->rep_len[off + j] = b->rep_len;
			s->frag_gap[off + j] = b->frag_gap;
		}
	} else {
		mm_map_frag(s->p->mi, s->n_seg[i], qlens, qseqs, &s->n_reg[off], &s->reg[off], b, s->p->opt, s->seq[off].name);
		for (j = 0; j < s->n_seg[i]; ++j) {
			s->rep_len[off + j] = b->rep_len;
			s->frag_gap[off + j] = b->frag_gap;
		}
	}
	for (j = 0; j < s->n_seg[i]; ++j) // flip the query strand and coordinate to the original read strand
		if (s->n_seg[i] == 2 && ((j == 0 && (pe_ori>>1&1)) || (j == 1 && (pe_ori&1)))) {
			int k, t;
			mm_revcomp_bseq(&s->seq[off + j]);
			for (k = 0; k < s->n_reg[off + j]; ++k) {
				mm_reg1_t *r = &s->reg[off + j][k];
				t = r->qs;
				r->qs = qlens[j] - r->qe;
				r->qe = qlens[j] - t;
				r->rev = !r->rev;
				if (r->p) {
					if (r->p->trans_strand == 1) r->p->trans_strand = 2;
					else if (r->p->trans_strand == 2) r->p->trans_strand = 1;
				}
			}
		}
	if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
		fprintf(stderr, "QT\t%s\t%d\t%.6f\n", s->seq[off].name, tid, realtime() - t);
}

static void merge_hits(step_t *s)
{
	int f, i, k0, k, max_seg = 0, *n_reg_part, *rep_len_part, *frag_gap_part, *qlens;
	void *km;
	FILE **fp = s->p->fp_parts;
	const mm_mapopt_t *opt = s->p->opt;

	km = km_init();
	for (f = 0; f < s->n_frag; ++f)
		max_seg = max_seg > s->n_seg[f]? max_seg : s->n_seg[f];
	qlens = CALLOC(int, max_seg + s->p->n_parts * 3);
	n_reg_part = qlens + max_seg;
	rep_len_part = n_reg_part + s->p->n_parts;
	frag_gap_part = rep_len_part + s->p->n_parts;
	for (f = 0, k = k0 = 0; f < s->n_frag; ++f) {
		k0 = k;
		for (i = 0; i < s->n_seg[f]; ++i, ++k) {
			int j, l, t, rep_len = 0;
			qlens[i] = s->seq[k].l_seq;
			for (j = 0, s->n_reg[k] = 0; j < s->p->n_parts; ++j) {
				mm_err_fread(&n_reg_part[j],    sizeof(int), 1, fp[j]);
				mm_err_fread(&rep_len_part[j],  sizeof(int), 1, fp[j]);
				mm_err_fread(&frag_gap_part[j], sizeof(int), 1, fp[j]);
				s->n_reg[k] += n_reg_part[j];
				if (rep_len < rep_len_part[j])
					rep_len = rep_len_part[j];
			}
			s->reg[k] = CALLOC(mm_reg1_t, s->n_reg[k]);
			for (j = 0, l = 0; j < s->p->n_parts; ++j) {
				for (t = 0; t < n_reg_part[j]; ++t, ++l) {
					mm_reg1_t *r = &s->reg[k][l];
					uint32_t capacity;
					mm_err_fread(r, sizeof(mm_reg1_t), 1, fp[j]);
					r->rid += s->p->rid_shift[j];
					if (opt->flag & MM_F_CIGAR) {
						mm_err_fread(&capacity, 4, 1, fp[j]);
						r->p = (mm_extra_t*)calloc(capacity, 4);
						r->p->capacity = capacity;
						mm_err_fread(r->p, r->p->capacity, 4, fp[j]);
					}
				}
			}
			if (!(opt->flag&MM_F_SR) && s->seq[k].l_seq >= opt->rank_min_len)
				mm_update_dp_max(s->seq[k].l_seq, s->n_reg[k], s->reg[k], opt->rank_frac, opt->a, opt->b);
			for (j = 0; j < s->n_reg[k]; ++j) {
				mm_reg1_t *r = &s->reg[k][j];
				if (r->p) r->p->dp_max2 = 0; // reset ->dp_max2 as mm_set_parent() doesn't clear it; necessary with mm_update_dp_max()
				r->subsc = 0; // this may not be necessary
				r->n_sub = 0; // n_sub will be an underestimate as we don't see all the chains now, but it can't be accurate anyway
			}
			mm_hit_sort(km, &s->n_reg[k], s->reg[k], opt->alt_drop);
			mm_set_parent(km, opt->mask_level, opt->mask_len, s->n_reg[k], s->reg[k], opt->a * 2 + opt->b, opt->flag&MM_F_HARD_MLEVEL, opt->alt_drop);
			if (!(opt->flag & MM_F_ALL_CHAINS)) {
				mm_select_sub(km, opt->pri_ratio, s->p->mi->k*2, opt->best_n, 0, opt->max_gap * 0.8, &s->n_reg[k], s->reg[k]);
				mm_set_sam_pri(s->n_reg[k], s->reg[k]);
			}
			mm_set_mapq2(km, s->n_reg[k], s->reg[k], opt->min_chain_score, opt->a, rep_len, !!(opt->flag & (MM_F_SR|MM_F_SR_RNA)), !!(opt->flag & MM_F_SPLICE));
		}
		if (s->n_seg[f] == 2 && opt->pe_ori >= 0 && (opt->flag&MM_F_CIGAR))
			mm_pair(km, frag_gap_part[0], opt->pe_bonus, opt->a * 2 + opt->b, opt->a, qlens, &s->n_reg[k0], &s->reg[k0]);
	}
	free(qlens);
	km_destroy(km);
}

static step_t *read_batch(pipeline_t *p)
{
	int i, j;
	int with_qual = (!!(p->opt->flag & MM_F_OUT_SAM) && !(p->opt->flag & MM_F_NO_QUAL));
	int with_comment = !!(p->opt->flag & MM_F_COPY_COMMENT);
	int frag_mode = (p->n_fp > 1 || !!(p->opt->flag & MM_F_FRAG_MODE));
	step_t *s;
	s = (step_t*)calloc(1, sizeof(step_t));
	if (p->n_fp > 1) s->seq = mm_bseq_read_frag2(p->n_fp, p->fp, p->mini_batch_size, with_qual, with_comment, &s->n_seq);
	else s->seq = mm_bseq_read3(p->fp[0], p->mini_batch_size, with_qual, with_comment, frag_mode, &s->n_seq);
	if (s->seq) {
		s->p = p;
		for (i = 0; i < s->n_seq; ++i)
			s->seq[i].rid = p->n_processed++;
		s->buf = (mm_tbuf_t**)calloc(p->n_threads, sizeof(mm_tbuf_t*));
		for (i = 0; i < p->n_threads; ++i)
			s->buf[i] = mm_tbuf_init();
		s->n_reg = (int*)calloc(5 * s->n_seq, sizeof(int));
		s->seg_off = s->n_reg + s->n_seq; // seg_off, n_seg, rep_len and frag_gap are allocated together with n_reg
		s->n_seg = s->seg_off + s->n_seq;
		s->rep_len = s->n_seg + s->n_seq;
		s->frag_gap = s->rep_len + s->n_seq;
		s->reg = (mm_reg1_t**)calloc(s->n_seq, sizeof(mm_reg1_t*));
		for (i = 1, j = 0; i <= s->n_seq; ++i)
			if (i == s->n_seq || !frag_mode || !mm_qname_same(s->seq[i-1].name, s->seq[i].name)) {
				s->n_seg[s->n_frag] = i - j;
				s->seg_off[s->n_frag++] = j;
				j = i;
			}
		return s;
	}
	free(s);
	return 0;
}

// Fragments longest first (total bases; ties in input order), so a long read late in a batch does not run alone at its end.
// Only for a batch with a long tail: one whose longest fragment, at 64 times the batch's mean cost per base (per-read cost
// spreads that far on ONT and HiFi reads), would take over half of one thread's share. 0 for input order.
static uint64_t *frag_order(const step_t *s, int n_threads)
{
	int i, j, sorted = 1;
	uint64_t *a, max_len = 0, tot = 0;
	if (s->n_frag < 2) return 0;
	a = (uint64_t*)malloc(s->n_frag * sizeof(uint64_t));
	for (i = 0; i < s->n_frag; ++i) {
		uint64_t len = 0;
		for (j = 0; j < s->n_seg[i]; ++j)
			len += s->seq[s->seg_off[i] + j].l_seq;
		if (len > UINT32_MAX) len = UINT32_MAX;
		tot += len;
		if (len > max_len) max_len = len;
		a[i] = (uint64_t)(UINT32_MAX - len) << 32 | (uint32_t)i;
		if (i > 0 && a[i] < a[i - 1]) sorted = 0;
	}
	if (sorted || max_len * n_threads * 128 <= tot) {
		free(a);
		return 0;
	}
	radix_sort_64(a, a + s->n_frag);
	return a;
}

static void *worker_pipeline(void *shared, int step, void *in)
{
	int i, j, k;
    pipeline_t *p = (pipeline_t*)shared;
    if (step == 0) { // step 0: read sequences
		if (p->has_first) {
			p->has_first = 0;
			return p->first;
		}
		return read_batch(p);
    } else if (step == 1) { // step 1: map
		step_t *s = (step_t*)in;
		if (p->n_parts > 0) merge_hits(s);
		else {
			if (p->n_threads > 1) s->order = frag_order(s, p->n_threads);
			kt_for(p->n_threads, worker_for, s, s->n_frag);
			free(s->order);
			s->order = 0;
		}
		return in;
    } else if (step == 2) { // step 2: output
		void *km = 0;
        step_t *s = (step_t*)in;
		const mm_idx_t *mi = p->mi;
		for (i = 0; i < p->n_threads; ++i) mm_tbuf_destroy(s->buf[i]);
		free(s->buf);
		if ((p->opt->flag & MM_F_OUT_CS) && !(mm_dbg_flag & MM_DBG_NO_KALLOC)) km = km_init();
		for (k = 0; k < s->n_frag; ++k) {
			int seg_st = s->seg_off[k], seg_en = s->seg_off[k] + s->n_seg[k];
			for (i = seg_st; i < seg_en; ++i) {
				mm_bseq1_t *t = &s->seq[i];
				if (p->opt->split_prefix && p->n_parts == 0) { // then write to temporary files
					mm_err_fwrite(&s->n_reg[i],    sizeof(int), 1, p->fp_split);
					mm_err_fwrite(&s->rep_len[i],  sizeof(int), 1, p->fp_split);
					mm_err_fwrite(&s->frag_gap[i], sizeof(int), 1, p->fp_split);
					for (j = 0; j < s->n_reg[i]; ++j) {
						mm_reg1_t *r = &s->reg[i][j];
						mm_err_fwrite(r, sizeof(mm_reg1_t), 1, p->fp_split);
						if (p->opt->flag & MM_F_CIGAR) {
							mm_err_fwrite(&r->p->capacity, 4, 1, p->fp_split);
							mm_err_fwrite(r->p, r->p->capacity, 4, p->fp_split);
						}
					}
				} else if (p->opt->flag & MM_F_OUT_JUNC) { // extra logic for --write-junc
					for (j = 0; j < s->n_reg[i]; ++j) {
						const mm_reg1_t *r = &s->reg[i][j];
						if (r->id != r->parent || r->mapq < 10) continue;
						mm_write_junc(&p->str, mi, t, r);
						if (p->str.l > 0) mm_err_puts(p->str.s);
					}
				} else if (s->n_reg[i] > 0) { // the query has at least one hit
					for (j = 0; j < s->n_reg[i]; ++j) {
						const mm_reg1_t *r = &s->reg[i][j];
						assert(!r->sam_pri || r->id == r->parent);
						if ((p->opt->flag & MM_F_NO_PRINT_2ND) && r->id != r->parent)
							continue;
						if (p->opt->flag & MM_F_OUT_SAM)
							mm_write_sam3(&p->str, mi, t, i - seg_st, j, s->n_seg[k], &s->n_reg[seg_st], (const mm_reg1_t*const*)&s->reg[seg_st], km, p->opt->flag, s->rep_len[i]);
						else
							mm_write_paf4(&p->str, mi, t, r, km, p->opt->flag, s->rep_len[i], s->n_seg[k], i - seg_st);
						mm_err_puts(p->str.s);
					}
				} else if ((p->opt->flag & MM_F_PAF_NO_HIT) || ((p->opt->flag & MM_F_OUT_SAM) && !(p->opt->flag & MM_F_SAM_HIT_ONLY))) { // output an empty hit, if requested
					if (p->opt->flag & MM_F_OUT_SAM)
						mm_write_sam3(&p->str, mi, t, i - seg_st, -1, s->n_seg[k], &s->n_reg[seg_st], (const mm_reg1_t*const*)&s->reg[seg_st], km, p->opt->flag, s->rep_len[i]);
					else
						mm_write_paf4(&p->str, mi, t, 0, 0, p->opt->flag, s->rep_len[i], s->n_seg[k], i - seg_st);
					mm_err_puts(p->str.s);
				}
			}
			for (i = seg_st; i < seg_en; ++i) {
				for (j = 0; j < s->n_reg[i]; ++j) free(s->reg[i][j].p);
				free(s->reg[i]);
				free(s->seq[i].seq); free(s->seq[i].name);
				if (s->seq[i].qual) free(s->seq[i].qual);
				if (s->seq[i].comment) free(s->seq[i].comment);
			}
		}
		free(s->reg); free(s->n_reg); free(s->seq); // seg_off, n_seg, rep_len and frag_gap were allocated with reg; no memory leak here
		km_destroy(km);
		if (mm_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] mapped %d sequences\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), s->n_seq);
		free(s);
	}
    return 0;
}

static mm_bseq_file_t **open_bseqs(int n, const char **fn, int quiet)
{
	mm_bseq_file_t **fp;
	int i, j;
	fp = (mm_bseq_file_t**)calloc(n, sizeof(mm_bseq_file_t*));
	for (i = 0; i < n; ++i) {
		if ((fp[i] = mm_bseq_open(fn[i])) == 0) {
			if (mm_verbose >= 1 && !quiet)
				fprintf(stderr, "ERROR: failed to open file '%s': %s\n", fn[i], strerror(errno));
			for (j = 0; j < i; ++j)
				mm_bseq_close(fp[j]);
			free(fp);
			return 0;
		}
	}
	return fp;
}

struct mm_prefetch_s {
	pthread_t tid;
	const char **fn;
	mm_mapopt_t opt; // a copy: the caller may update its own (mm_mapopt_update) while the batch is read
	pipeline_t pl;
	step_t *s;
	int done; // set when the read has ended; see mm_prefetch_busy()
	int left; // set by the first of the ended read and mm_prefetch_release(); the second frees pf
};

static void prefetch_drop(mm_prefetch_t *pf);

static void *prefetch_worker(void *data)
{
	mm_prefetch_t *pf = (mm_prefetch_t*)data;
	pf->pl.fp = open_bseqs(pf->pl.n_fp, pf->fn, 1); // quiet: on a failure mm_map_file_frag_pf opens the files again and reports the error there
	if (pf->pl.fp) pf->s = read_batch(&pf->pl);
	__atomic_store_n(&pf->done, 1, __ATOMIC_RELEASE);
	if (__atomic_exchange_n(&pf->left, 1, __ATOMIC_ACQ_REL)) { // released while reading
		if (pf->pl.fp) prefetch_drop(pf);
		free(pf);
	}
	return 0;
}

int mm_prefetch_busy(const mm_prefetch_t *pf)
{
	return pf && !__atomic_load_n(&pf->done, __ATOMIC_ACQUIRE);
}

mm_prefetch_t *mm_prefetch_start(int n_segs, const char **fn, const mm_mapopt_t *opt, int n_threads)
{
	mm_prefetch_t *pf;
	if (n_segs < 1) return 0;
	pf = (mm_prefetch_t*)calloc(1, sizeof(mm_prefetch_t));
	pf->fn = fn;
	pf->opt = *opt;
	pf->pl.n_fp = n_segs;
	pf->pl.opt = &pf->opt;
	pf->pl.n_threads = n_threads > 1? n_threads : 1;
	pf->pl.mini_batch_size = opt->mini_batch_size;
	if (pthread_create(&pf->tid, 0, prefetch_worker, pf) != 0) {
		free(pf);
		return 0;
	}
	return pf;
}

static void prefetch_drop(mm_prefetch_t *pf) // frees the batch and closes the files of a prefetch whose read has ended and opened them
{
	step_t *s = pf->s;
	int i;
	if (s) {
		for (i = 0; i < pf->pl.n_threads; ++i) mm_tbuf_destroy(s->buf[i]);
		for (i = 0; i < s->n_seq; ++i) {
			free(s->seq[i].seq); free(s->seq[i].name);
			free(s->seq[i].qual); free(s->seq[i].comment);
		}
		free(s->buf); free(s->reg); free(s->n_reg); free(s->seq); free(s);
	}
	for (i = 0; i < pf->pl.n_fp; ++i)
		mm_bseq_close(pf->pl.fp[i]);
	free(pf->pl.fp);
}

// Joins the prefetch thread; returns its files and batch through pl when they are what step 0 would read for these arguments
static void prefetch_take(mm_prefetch_t *pf, pipeline_t *pl, const char **fn, const mm_mapopt_t *opt)
{
	const int64_t step0_flags = MM_F_OUT_SAM | MM_F_NO_QUAL | MM_F_COPY_COMMENT | MM_F_FRAG_MODE;
	int i, same;
	pthread_join(pf->tid, 0);
	if (pf->pl.fp == 0) {
		free(pf);
		return;
	}
	same = (pf->pl.n_fp == pl->n_fp && pf->pl.n_threads == pl->n_threads && pf->opt.mini_batch_size == opt->mini_batch_size
			&& (pf->opt.flag & step0_flags) == (opt->flag & step0_flags));
	for (i = 0; same && i < pl->n_fp; ++i)
		same = (strcmp(pf->fn[i], fn[i]) == 0);
	if (same) {
		pl->fp = pf->pl.fp;
		pl->n_processed = pf->pl.n_processed;
		pl->has_first = 1, pl->first = pf->s;
		if (pf->s) pf->s->p = pl;
	} else prefetch_drop(pf); // the files are read again from the start
	free(pf);
}

void mm_prefetch_destroy(mm_prefetch_t *pf)
{
	if (pf == 0) return;
	pthread_join(pf->tid, 0);
	if (pf->pl.fp) prefetch_drop(pf);
	free(pf);
}

void mm_prefetch_release(mm_prefetch_t *pf)
{
	if (pf == 0) return;
	pthread_detach(pf->tid);
	if (__atomic_exchange_n(&pf->left, 1, __ATOMIC_ACQ_REL)) { // the read has ended
		if (pf->pl.fp) prefetch_drop(pf);
		free(pf);
	}
}

int mm_map_file_frag_pf(const mm_idx_t *idx, int n_segs, const char **fn, const mm_mapopt_t *opt, int n_threads, mm_prefetch_t *pf)
{
	int i, pl_threads;
	pipeline_t pl;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.n_fp = n_segs;
	pl.n_threads = n_threads > 1? n_threads : 1;
	if (pf) prefetch_take(pf, &pl, fn, opt);
	if (n_segs < 1) return -1;
	if (pl.fp == 0) pl.fp = open_bseqs(pl.n_fp, fn, 0);
	if (pl.fp == 0) return -1;
	pl.opt = opt, pl.mi = idx;
	pl.mini_batch_size = opt->mini_batch_size;
	if (opt->split_prefix)
		pl.fp_split = mm_split_init(opt->split_prefix, idx);
	pl_threads = n_threads == 1? 1 : (opt->flag&MM_F_2_IO_THREADS)? 3 : 2;
	kt_pipeline(pl_threads, worker_pipeline, &pl, 3);

	free(pl.str.s);
	if (pl.fp_split) fclose(pl.fp_split);
	for (i = 0; i < pl.n_fp; ++i)
		mm_bseq_close(pl.fp[i]);
	free(pl.fp);
	return 0;
}

int mm_map_file_frag(const mm_idx_t *idx, int n_segs, const char **fn, const mm_mapopt_t *opt, int n_threads)
{
	return mm_map_file_frag_pf(idx, n_segs, fn, opt, n_threads, 0);
}

int mm_map_file(const mm_idx_t *idx, const char *fn, const mm_mapopt_t *opt, int n_threads)
{
	return mm_map_file_frag(idx, 1, &fn, opt, n_threads);
}

int mm_split_merge(int n_segs, const char **fn, const mm_mapopt_t *opt, int n_split_idx)
{
	int i;
	pipeline_t pl;
	mm_idx_t *mi;
	if (n_segs < 1 || n_split_idx < 1) return -1;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.n_fp = n_segs;
	pl.fp = open_bseqs(pl.n_fp, fn, 0);
	if (pl.fp == 0) return -1;
	pl.opt = opt;
	pl.mini_batch_size = opt->mini_batch_size;

	pl.n_parts = n_split_idx;
	pl.fp_parts  = CALLOC(FILE*, pl.n_parts);
	pl.rid_shift = CALLOC(uint32_t, pl.n_parts);
	pl.mi = mi = mm_split_merge_prep(opt->split_prefix, n_split_idx, pl.fp_parts, pl.rid_shift);
	if (pl.mi == 0) {
		free(pl.fp_parts);
		free(pl.rid_shift);
		return -1;
	}
	for (i = n_split_idx - 1; i > 0; --i)
		pl.rid_shift[i] = pl.rid_shift[i - 1];
	for (pl.rid_shift[0] = 0, i = 1; i < n_split_idx; ++i)
		pl.rid_shift[i] += pl.rid_shift[i - 1];
	if (opt->flag & MM_F_OUT_SAM)
		for (i = 0; i < (int32_t)pl.mi->n_seq; ++i)
			printf("@SQ\tSN:%s\tLN:%d\n", pl.mi->seq[i].name, pl.mi->seq[i].len);

	kt_pipeline(2, worker_pipeline, &pl, 3);

	free(pl.str.s);
	mm_idx_destroy(mi);
	free(pl.rid_shift);
	for (i = 0; i < n_split_idx; ++i)
		fclose(pl.fp_parts[i]);
	free(pl.fp_parts);
	for (i = 0; i < pl.n_fp; ++i)
		mm_bseq_close(pl.fp[i]);
	free(pl.fp);
	mm_split_rm_tmp(opt->split_prefix, n_split_idx);
	return 0;
}
