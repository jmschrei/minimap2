/* Area B driver: chaining (lchain.c) and the sorting helpers of misc.c.
 *
 * Builds anchor arrays (mm128_t, the layout mmpriv.h documents for
 * mg_lchain_dp/mg_lchain_rmq) by hand and checks the chains against the
 * chaining score of Li (2018), as minimap2.1 describes it for -m: "matching
 * bases minus a concave gap penalty". Between anchors j and i with query
 * distance dq, target distance dr, dd = |dr - dq| and dg = min(dq, dr), the
 * score gains min(dg, k) and loses chn_pen_gap*dd + chn_pen_skip*dg +
 * 0.5*log2(dd+1) (truncated to an integer; 0 when dd == 0 and dg <= k). The
 * crafted cases use gaps whose cost is not near an integer, so exact log2
 * gives the same truncation as minimap2's approximation.
 *
 * Anchor layout: x = strand<<63 | contig<<32 | target end position,
 * y = span<<32 | query end position. Internal function signatures from
 * mmpriv.h; a speed change that alters them breaks this driver, not the
 * documented behaviour it checks. */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <limits.h>
#include "mmpriv.h"
#include "kalloc.h"
#include "minimap2_test.h"

static uint64_t rng_state;

static uint64_t rng_next(void)
{
	uint64_t z = (rng_state += 0x9E3779B97F4A7C15ULL);
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
	z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
	return z ^ (z >> 31);
}

static int rng_int(int lo, int hi) { return lo + (int)(rng_next() % (uint64_t)(hi - lo + 1)); }

static mm128_t anchor(int rev, int cid, int32_t tend, int32_t qend, int span)
{
	mm128_t a;
	a.x = (uint64_t)rev << 63 | (uint64_t)cid << 32 | (uint32_t)tend;
	a.y = (uint64_t)span << 32 | (uint32_t)qend;
	return a;
}

typedef struct {
	int n_u;
	uint64_t *u;
	mm128_t *a; /* anchors of the chains, chain after chain */
} chains_t;

static int rmq_inner = 1000; /* --rmq-inner for mg_lchain_rmq */

/* method 0: mg_lchain_dp; 1: mg_lchain_rmq. The input copy is consumed by the chaining function. */
static void run_chain(void *km, int method, const mm128_t *in, int64_t n, int max_dist, int bw, int max_skip, int max_iter, int min_cnt, int min_sc,
					  float gap_scale, float skip_scale, int k, chains_t *c)
{
	mm128_t *a = (mm128_t*)kmalloc(km, (n > 0? n : 1) * sizeof(mm128_t));
	memcpy(a, in, n * sizeof(mm128_t));
	if (method == 0)
		c->a = mg_lchain_dp(max_dist, max_dist, bw, max_skip, max_iter, min_cnt, min_sc, gap_scale * 0.01 * k, skip_scale * 0.01 * k, 0, 1, n, a, &c->n_u, &c->u, km);
	else
		c->a = mg_lchain_rmq(max_dist, rmq_inner, bw, max_skip, 100000, min_cnt, min_sc, gap_scale * 0.01 * k, skip_scale * 0.01 * k, n, a, &c->n_u, &c->u, km);
}

static void free_chains(void *km, chains_t *c)
{
	kfree(km, c->a); kfree(km, c->u);
	memset(c, 0, sizeof(*c));
}

static int chain_score(const chains_t *c, int i) { return (int)(c->u[i] >> 32); }
static int chain_cnt(const chains_t *c, int i) { return (int32_t)c->u[i]; }

/* the gap cost of the documented formula */
static int gap_cost(int dq, int dr, int span, float gap_scale, float skip_scale, int k)
{
	int dd = dr > dq? dr - dq : dq - dr, dg = dr < dq? dr : dq;
	double c;
	if (dd == 0 && dg <= span) return 0;
	c = gap_scale * 0.01 * k * dd + skip_scale * 0.01 * k * dg + (dd > 0? 0.5 * log2(dd + 1.0) : 0.0);
	return (int)c;
}

/* args: none. 40 anchors on one diagonal, 7 bases apart, span 15: one chain of 40 anchors whose score
 * is the 15 + 39*7 = 288 query bases they cover (no gap cost on one diagonal); the same with
 * anchors 20 apart (gaps longer than the span: 40*15 = 600). Both chaining methods. */
static int case_collinear(int argc, char *argv[])
{
	void *km = km_init();
	int method, step;
	for (method = 0; method < 2; ++method) {
		for (step = 7; step <= 20; step += 13) {
			mm128_t in[40];
			chains_t c;
			int i, want = step < 15? 15 + 39 * step : 40 * 15;
			for (i = 0; i < 40; ++i) in[i] = anchor(0, 0, 1000 + i * step, 100 + i * step, 15);
			run_chain(km, method, in, 40, 5000, 500, 25, 5000, 3, 40, 0.8f, 0.0f, 15, &c);
			MMT_REQUIRE(c.n_u == 1);
			MMT_CHECKF(chain_score(&c, 0) == want && chain_cnt(&c, 0) == 40, "method %d step %d: score %d cnt %d, want %d 40", method, step, chain_score(&c, 0), chain_cnt(&c, 0), want);
			for (i = 0; i < 40; ++i)
				MMT_CHECKF(c.a[i].x == in[i].x && c.a[i].y == in[i].y, "method %d step %d: anchor %d out of order", method, step, i);
			free_chains(km, &c);
		}
	}
	km_destroy(km);
	return 0;
}

/* args: none. Anchors on two contigs, on the two strands of one contig, and on one diagonal but
 * further apart than max_dist (-g) never share a chain; within max_dist they do. */
static int case_separate(int argc, char *argv[])
{
	void *km = km_init();
	int method, i;
	for (method = 0; method < 2; ++method) {
		mm128_t in[60];
		chains_t c;
		int n = 0, j, s0 = 0, s1 = 0;
		for (i = 0; i < 20; ++i) in[n++] = anchor(0, 0, 1000 + i * 10, 100 + i * 10, 15);
		for (i = 0; i < 20; ++i) in[n++] = anchor(0, 1, 1000 + i * 10, 100 + i * 10, 15); /* another contig */
		for (i = 0; i < 20; ++i) in[n++] = anchor(1, 0, 1000 + i * 10, 100 + i * 10, 15); /* reverse strand */
		radix_sort_128x(in, in + n);
		run_chain(km, method, in, n, 5000, 500, 25, 5000, 3, 40, 0.8f, 0.0f, 15, &c);
		MMT_CHECKF(c.n_u == 3, "method %d: %d chains from three groups", method, c.n_u);
		for (j = 0; j < c.n_u; ++j) MMT_CHECKF(chain_cnt(&c, j) == 20 && chain_score(&c, j) == 15 + 19 * 10, "method %d chain %d: cnt %d score %d", method, j, chain_cnt(&c, j), chain_score(&c, j));
		for (j = 0, i = 0; j < c.n_u; ++j) { /* every chain keeps one contig and strand */
			int l, same = 1;
			for (l = 1; l < chain_cnt(&c, j); ++l)
				if (c.a[i + l].x >> 32 != c.a[i].x >> 32) same = 0;
			MMT_CHECKF(same, "method %d chain %d mixes contigs or strands", method, j);
			i += chain_cnt(&c, j);
		}
		free_chains(km, &c);
		/* two runs of 10 anchors on one diagonal, separated by a gap of G bases on both sequences */
		for (j = 0; j < 2; ++j) {
			int G = j == 0? 4990 : 5010; /* max_dist 5000: dq = dr = G + 10 between the runs' facing anchors */
			n = 0;
			for (i = 0; i < 10; ++i) in[n++] = anchor(0, 0, 1000 + i * 10, 100 + i * 10, 15);
			for (i = 0; i < 10; ++i) in[n++] = anchor(0, 0, 1090 + G + i * 10, 190 + G + i * 10, 15);
			run_chain(km, method, in, n, 5000 + 10 * (j == 0), 500, 25, 5000, 3, 40, 0.8f, 0.0f, 15, &c);
			if (j == 0) s0 = c.n_u; else s1 = c.n_u;
			free_chains(km, &c);
		}
		MMT_CHECKF(s0 == 1 && s1 == 2, "method %d: a gap within max_dist gave %d chains, beyond it %d", method, s0, s1);
	}
	km_destroy(km);
	return 0;
}

/* args: none. Two runs of 20 anchors (10 apart, span 15) whose diagonals differ by dd (an indel of dd
 * bases): one chain scoring (15+19*10) + min(dg,15) - cost(dd) + 19*10 for dd <= bw (dg = 10 between the
 * runs; steps inside a run have dd 0 and dg 10 <= span, so no cost), at --chain-gap-scale
 * 0.8, 1 and 2 and --chain-skip-scale 0 and 0.3; two chains when dd > bw. Deletion (dr > dq) and
 * insertion (dq > dr) both; dp and rmq. */
static int case_gap_cost(int argc, char *argv[])
{
	static const int dds[] = { 1, 3, 9, 17, 40, 117, 333 };
	static const float gscale[] = { 0.8f, 1.0f, 2.0f }, sscale[] = { 0.0f, 0.3f };
	void *km = km_init();
	int method, di, gi, si, ins;
	for (method = 0; method < 2; ++method)
	for (di = 0; di < 7; ++di)
	for (gi = 0; gi < 3; ++gi)
	for (si = 0; si < 2; ++si)
	for (ins = 0; ins < 2; ++ins) {
		mm128_t in[40];
		chains_t c;
		int i, n = 0, dd = dds[di], dq, dr, want, bw = 200;
		double raw;
		for (i = 0; i < 20; ++i) in[n++] = anchor(0, 0, 1000 + i * 10, 100 + i * 10, 15);
		for (i = 0; i < 20; ++i) in[n++] = anchor(0, 0, 1200 + (ins? 0 : dd) + i * 10, 300 + (ins? dd : 0) + i * 10, 15);
		dq = 300 + (ins? dd : 0) - 290, dr = 1200 + (ins? 0 : dd) - 1190;
		raw = gscale[gi] * 0.01 * 15 * dd + sscale[si] * 0.01 * 15 * (dq < dr? dq : dr) + 0.5 * log2(dd + 1.0);
		if (raw - floor(raw) < 0.02 || raw - floor(raw) > 0.98) continue; /* truncation too close to call */
		want = (15 + 19 * 10) + (dq < dr? dq : dr) - gap_cost(dq, dr, 15, gscale[gi], sscale[si], 15) + 19 * 10; /* min(dg, span) = dg = 10 */
		run_chain(km, method, in, n, 5000, bw, 25, 5000, 3, 40, gscale[gi], sscale[si], 15, &c);
		if (dd <= bw) {
			MMT_CHECKF(c.n_u == 1 && chain_score(&c, 0) == want && chain_cnt(&c, 0) == 40, "method %d dd %d gap %.1f skip %.1f %s: %d chains, score %d cnt %d; want 1 chain of %d",
					   method, dd, gscale[gi], sscale[si], ins? "ins" : "del", c.n_u, c.n_u? chain_score(&c, 0) : -1, c.n_u? chain_cnt(&c, 0) : -1, want);
		} else {
			MMT_CHECKF(c.n_u == 2, "method %d dd %d > bw %d: %d chains", method, dd, bw, c.n_u);
		}
		free_chains(km, &c);
	}
	km_destroy(km);
	return 0;
}

/* args: none. The band edge (-r): two runs of 20 anchors whose diagonals differ by exactly bw form
 * one chain of all 40, and by bw+1 no chain holds anchors of both runs, for DP and RMQ chaining
 * (with the inner DP of --rmq-inner, bw 100 and 333; without it, bw 100). Without the inner DP the
 * second run may not chain at all, since RMQ offers each anchor only the best-scoring predecessor. */
static int case_band_edge(int argc, char *argv[])
{
	void *km = km_init();
	int method, bi, extra;
	static const int bws[] = { 100, 333 };
	for (method = 0; method < 3; ++method) /* 0: DP; 1: RMQ with --rmq-inner 1000; 2: RMQ with --rmq-inner 0 */
	for (bi = 0; bi < 2; ++bi)
	for (extra = 0; extra < 2; ++extra) {
		mm128_t in[40];
		chains_t c;
		int i, n = 0, dd = bws[bi] + extra;
		if (method == 2 && bi == 1) continue; /* RMQ alone at dd 333: the gap cost lets later anchors pick other predecessors */
		for (i = 0; i < 20; ++i) in[n++] = anchor(0, 0, 1000 + i * 10, 100 + i * 10, 15);
		for (i = 0; i < 20; ++i) in[n++] = anchor(0, 0, 1200 + dd + i * 10, 300 + i * 10, 15);
		rmq_inner = method == 2? 0 : 1000;
		run_chain(km, method > 0, in, n, 5000, bws[bi], 25, 5000, 3, 40, 0.8f, 0.0f, 15, &c);
		if (!extra) {
			MMT_CHECKF(c.n_u == 1 && chain_cnt(&c, 0) == 40, "method %d bw %d dd %d: %d chains, the first of %d anchors", method, bws[bi], dd, c.n_u, c.n_u? chain_cnt(&c, 0) : 0);
		} else {
			int j, k = 0, mixed = 0;
			for (j = 0; j < c.n_u; k += chain_cnt(&c, j), ++j)
				if ((int32_t)c.a[k].y < 300 && (int32_t)c.a[k + chain_cnt(&c, j) - 1].y >= 300) mixed = 1;
			MMT_CHECKF(c.n_u >= 1 && !mixed, "method %d bw %d dd %d: a chain crosses the band (%d chains)", method, bws[bi], dd, c.n_u);
		}
		free_chains(km, &c);
	}
	rmq_inner = 1000;
	km_destroy(km);
	return 0;
}

/* args: none. -n and -m: a chain of 12 anchors scoring 15+11*10 = 125 is kept with min_cnt 12 and
 * min_sc 125, and dropped with min_cnt 13 or min_sc 126. */
static int case_min_cnt_sc(int argc, char *argv[])
{
	void *km = km_init();
	int method, i;
	mm128_t in[12];
	for (i = 0; i < 12; ++i) in[i] = anchor(0, 0, 500 + i * 10, 50 + i * 10, 15);
	for (method = 0; method < 2; ++method) {
		chains_t c;
		run_chain(km, method, in, 12, 5000, 500, 25, 5000, 12, 125, 0.8f, 0.0f, 15, &c);
		MMT_CHECKF(c.n_u == 1 && chain_score(&c, 0) == 125, "method %d: min_cnt 12 min_sc 125 gave %d chains", method, c.n_u);
		free_chains(km, &c);
		run_chain(km, method, in, 12, 5000, 500, 25, 5000, 13, 125, 0.8f, 0.0f, 15, &c);
		MMT_CHECKF(c.n_u == 0, "method %d: min_cnt 13 kept %d chains", method, c.n_u);
		free_chains(km, &c);
		run_chain(km, method, in, 12, 5000, 500, 25, 5000, 12, 126, 0.8f, 0.0f, 15, &c);
		MMT_CHECKF(c.n_u == 0, "method %d: min_sc 126 kept %d chains", method, c.n_u);
		free_chains(km, &c);
		run_chain(km, method, in, 0, 5000, 500, 25, 5000, 1, 1, 0.8f, 0.0f, 15, &c);
		MMT_CHECKF(c.n_u == 0 && c.a == 0, "method %d: no anchors gave %d chains", method, c.n_u);
		free_chains(km, &c);
	}
	km_destroy(km);
	return 0;
}

/* args: none. --max-chain-iter: 20 true anchors on one diagonal (10 apart), each followed in target
 * order by a decoy anchor that cannot chain with them (its query position is 10000 further). With
 * max_iter 2 each true anchor sees the previous one past its decoy and the 20 chain (score
 * 15 + 19*10 = 205); with max_iter 1 each sees only its decoy, so the true anchors lose links and no
 * chain holds all 20 or scores 205. The decoys chain among themselves either way. */
static int case_max_iter(int argc, char *argv[])
{
	void *km = km_init();
	mm128_t in[40];
	chains_t c;
	int i, j, n = 0, it;
	for (i = 0; i < 20; ++i) {
		in[n++] = anchor(0, 0, 1000 + i * 10, 100 + i * 10, 15);
		in[n++] = anchor(0, 0, 1005 + i * 10, 10100 + i * 37, 15);
	}
	for (it = 1; it <= 3; ++it) {
		int best_cnt = 0, best_sc = 0;
		run_chain(km, 0, in, n, 20000, 500, 25, it, 3, 40, 0.8f, 0.0f, 15, &c);
		for (j = 0, i = 0; j < c.n_u; i += chain_cnt(&c, j), ++j)
			if ((int32_t)c.a[i].y < 10000 && chain_cnt(&c, j) > best_cnt) best_cnt = chain_cnt(&c, j), best_sc = chain_score(&c, j);
		if (it >= 2) MMT_CHECKF(best_cnt == 20 && best_sc == 205, "max_iter %d: best true chain has %d anchors, score %d", it, best_cnt, best_sc);
		else MMT_CHECKF(best_cnt < 20 && best_sc < 205, "max_iter 1: a chain of %d true anchors scoring %d formed", best_cnt, best_sc);
		free_chains(km, &c);
	}
	km_destroy(km);
	return 0;
}

/* brute-force chaining DP of the documented formula over anchors sorted by x; returns the best f */
static int ref_best_chain(int n, const mm128_t *a, int max_dist, int bw, float gap_scale, int k)
{
	int i, j, best = 0, *f = (int*)malloc(n * sizeof(int));
	for (i = 0; i < n; ++i) {
		int span = a[i].y >> 32 & 0xff;
		f[i] = span;
		for (j = i - 1; j >= 0; --j) {
			int32_t dq = (int32_t)a[i].y - (int32_t)a[j].y, dr = (int32_t)a[i].x - (int32_t)a[j].x, dd, dg, sc, sj = a[j].y >> 32 & 0xff;
			if (a[i].x >> 32 != a[j].x >> 32) continue;
			if (dq <= 0 || dr <= 0 || dq > max_dist || dr > max_dist) continue;
			dd = dr > dq? dr - dq : dq - dr;
			if (dd > bw) continue;
			dg = dr < dq? dr : dq;
			sc = (dg < sj? dg : sj);
			if (dd || dg > sj) {
				float lin = gap_scale * 0.01f * k * (float)dd, lg = dd >= 1? mg_log2(dd + 1) : 0.0f;
				sc -= (int)(lin + .5f * lg); /* mg_log2: minimap2's log2 for x >= 2 (mmpriv.h) */
			}
			if (f[j] + sc > f[i]) f[i] = f[j] + sc;
		}
		if (f[i] > best) best = f[i];
	}
	free(f);
	return best;
}

/* args: seed n_sets. Random anchor sets (two contigs, both strands, true chains with indels plus
 * random decoys): with the heuristics off (max_skip and max_iter INT_MAX, "set INT to a large
 * number to switch off this heuristic") the best chain from mg_lchain_dp scores the brute-force
 * optimum; every chain's anchors keep one contig and strand and increase in both coordinates
 * (dp and rmq); the default heuristics, --max-chain-skip 0, --max-chain-iter 1 and RMQ chaining never
 * beat the optimum, and the two extreme settings lose score on some sets; RMQ chaining with full DP
 * for anchors within 1000 (--rmq-inner) finds better best chains than with --rmq-inner 0 more often
 * than worse ones. */
static int case_random_dp(int argc, char *argv[])
{
	void *km = km_init();
	int s, n_sets, h, n_lose[3] = {0, 0, 0}, n_inner_helps = 0, n_inner_hurts = 0, rmq_best[2];
	if (argc < 2) return MMT_FAIL("usage: random_dp <seed> <n_sets>"), 1;
	rng_state = strtoull(argv[0], 0, 10), n_sets = atoi(argv[1]);
	for (s = 0; s < n_sets; ++s) {
		mm128_t *in;
		int n = 0, m = rng_int(20, 400), i, j, want, got = 0, got2 = 0, k = 15;
		chains_t c;
		in = (mm128_t*)malloc(m * sizeof(mm128_t));
		while (n < m) {
			int rev = rng_int(0, 1), cid = rng_int(0, 1), len = rng_int(3, 40), x = rng_int(0, 20000), y = rng_int(0, 5000);
			for (i = 0; i < len && n < m; ++i) {
				int dx = rng_int(3, 25), dy = dx;
				if (rng_int(0, 9) == 0) dx += rng_int(1, 30); /* deletion */
				else if (rng_int(0, 9) == 0) dy += rng_int(1, 30); /* insertion */
				x += dx, y += dy;
				in[n++] = anchor(rev, cid, x, y, k);
			}
			if (rng_int(0, 2) == 0 && n < m) { /* fields drawn last-first in sequence: argument evaluation order is unspecified */
				int ty = rng_int(0, 8000), tx = rng_int(0, 30000), tcid = rng_int(0, 1), trev = rng_int(0, 1);
				in[n++] = anchor(trev, tcid, tx, ty, k);
			}
		}
		radix_sort_128x(in, in + n);
		want = ref_best_chain(n, in, 5000, 500, 0.8f, k);
		run_chain(km, 0, in, n, 5000, 500, INT_MAX, INT_MAX, 1, 1, 0.8f, 0.0f, k, &c);
		for (j = 0, i = 0; j < c.n_u; i += chain_cnt(&c, j), ++j) {
			int l, ok = 1;
			if (chain_score(&c, j) > got) got = chain_score(&c, j);
			for (l = 1; l < chain_cnt(&c, j); ++l)
				if (c.a[i+l].x >> 32 != c.a[i].x >> 32 || (int32_t)c.a[i+l].x <= (int32_t)c.a[i+l-1].x || (int32_t)c.a[i+l].y <= (int32_t)c.a[i+l-1].y) ok = 0;
			MMT_CHECKF(ok, "set %d chain %d: anchors not collinear", s, j);
		}
		MMT_CHECKF(got == want, "set %d (%d anchors): best chain %d, brute force %d", s, n, got, want);
		free_chains(km, &c);
		for (h = 0; h < 2; ++h) { /* RMQ chaining (--rmq-inner 1000, then 0): valid chains, never above the optimum */
			rmq_inner = h == 0? 1000 : 0;
			run_chain(km, 1, in, n, 5000, 500, 25, 5000, 1, 1, 0.8f, 0.0f, k, &c);
			for (j = 0, i = 0, got2 = 0; j < c.n_u; i += chain_cnt(&c, j), ++j) {
				int l, ok = 1;
				if (chain_score(&c, j) > got2) got2 = chain_score(&c, j);
				for (l = 1; l < chain_cnt(&c, j); ++l)
					if (c.a[i+l].x >> 32 != c.a[i].x >> 32 || (int32_t)c.a[i+l].x <= (int32_t)c.a[i+l-1].x || (int32_t)c.a[i+l].y <= (int32_t)c.a[i+l-1].y) ok = 0;
				MMT_CHECKF(ok, "set %d rmq chain %d: anchors not collinear", s, j);
			}
			MMT_CHECKF(got2 <= want && got2 > 0, "set %d: rmq best chain %d (inner %d), optimum %d", s, got2, rmq_inner, want);
			rmq_best[h] = got2;
			free_chains(km, &c);
		}
		n_inner_helps += rmq_best[0] > rmq_best[1];
		n_inner_hurts += rmq_best[0] < rmq_best[1];
		for (h = 0; h < 3; ++h) { /* the heuristics never beat the optimum; count how often they lose */
			static const int skip[3] = { 25, 0, 25 }, iter[3] = { 5000, 5000, 1 };
			got2 = 0;
			run_chain(km, 0, in, n, 5000, 500, skip[h], iter[h], 1, 1, 0.8f, 0.0f, k, &c);
			for (j = 0; j < c.n_u; ++j) if (chain_score(&c, j) > got2) got2 = chain_score(&c, j);
			MMT_CHECKF(got2 <= want, "set %d: heuristic chain %d (max_skip %d, max_iter %d) beats the optimum %d", s, got2, skip[h], iter[h], want);
			n_lose[h] += got2 < want;
			free_chains(km, &c);
		}
		free(in);
	}
	mmt_note("sets where the best chain lost score: max_skip 25 %d, max_skip 0 %d, max_iter 1 %d", n_lose[0], n_lose[1], n_lose[2]);
	mmt_note("rmq-inner 1000 vs 0: better in %d sets, worse in %d", n_inner_helps, n_inner_hurts);
	MMT_CHECKF(n_sets < 100 || n_inner_helps > n_inner_hurts, "--rmq-inner 1000 improved the best chain in %d sets and worsened it in %d", n_inner_helps, n_inner_hurts);
	MMT_CHECKF(n_sets < 100 || (n_lose[1] > 0 && n_lose[2] > 0), "--max-chain-skip 0 (%d sets) or --max-chain-iter 1 (%d sets) never lost score", n_lose[1], n_lose[2]);
	km_destroy(km);
	return 0;
}

static uint64_t fnv64(uint64_t h, uint64_t x)
{
	int k;
	for (k = 0; k < 8; ++k, x >>= 8) h = (h ^ (x & 0xff)) * 0x100000001b3ULL;
	return h;
}

static int cmp_chain_key(const void *a, const void *b)
{
	const uint64_t *x = (const uint64_t*)a, *y = (const uint64_t*)b;
	return x[0] < y[0]? -1 : x[0] > y[0]? 1 : x[1] < y[1]? -1 : x[1] > y[1];
}

/* one line for a chaining result, independent of the order of the chains: the number of chains and a
 * hash of the sorted (score << 32 | count, hash of the chain's anchors) pairs */
static void put_chains(FILE *fp, int set, int cfg, const chains_t *c)
{
	uint64_t *key = (uint64_t*)malloc((c->n_u + 1) * 16), h = 0xcbf29ce484222325ULL;
	int i, j, k = 0;
	for (j = 0; j < c->n_u; ++j) {
		uint64_t ha = 0xcbf29ce484222325ULL;
		for (i = 0; i < chain_cnt(c, j); ++i)
			ha = fnv64(fnv64(ha, c->a[k + i].x), c->a[k + i].y);
		key[2*j] = c->u[j], key[2*j+1] = ha;
		k += chain_cnt(c, j);
	}
	qsort(key, c->n_u, 16, cmp_chain_key);
	for (j = 0; j < c->n_u; ++j) h = fnv64(fnv64(h, key[2*j]), key[2*j+1]);
	fprintf(fp, "%d\t%d\t%d\t%d\t%016llx\n", set, cfg, c->n_u, c->n_u? chain_score(c, 0) : 0, (unsigned long long)h);
	free(key);
}

/* args: seed n_sets out. Random anchor sets as in random_dp, chained with the heuristics on: DP with
 * the defaults (--max-chain-skip 25, --max-chain-iter 5000), --max-chain-skip 0, --max-chain-iter 1 and
 * 3 (--max-chain-skip 2), and RMQ with --rmq-inner 1000, 40 and 0. Their results are not documented
 * beyond "heuristics", so each one's chains (scores, counts and anchors) are written as one line for a
 * pin; every chain is still checked to be collinear on one contig and strand. */
static int case_heuristics(int argc, char *argv[])
{
	static const struct { int rmq, max_skip, max_iter, inner; } cfg[] = {
		{ 0, 25, 5000, 0 }, { 0, 0, 5000, 0 }, { 0, 25, 1, 0 }, { 0, 2, 3, 0 }, { 1, 25, 0, 1000 }, { 1, 25, 0, 40 }, { 1, 25, 0, 0 }
	};
	void *km = km_init();
	int s, n_sets, h;
	FILE *fp;
	if (argc < 3) return MMT_FAIL("usage: heuristics <seed> <n_sets> <out>"), 1;
	rng_state = strtoull(argv[0], 0, 10), n_sets = atoi(argv[1]);
	MMT_REQUIRE((fp = fopen(argv[2], "w")) != 0);
	for (s = 0; s < n_sets; ++s) {
		mm128_t *in;
		int n = 0, m = rng_int(20, 300), i, j;
		in = (mm128_t*)malloc(m * sizeof(mm128_t));
		while (n < m) {
			int rev = rng_int(0, 1), cid = rng_int(0, 1), len = rng_int(3, 40), x = rng_int(0, 20000), y = rng_int(0, 5000);
			for (i = 0; i < len && n < m; ++i) {
				int dx = rng_int(3, 25), dy = dx;
				if (rng_int(0, 6) == 0) dx += rng_int(1, 30);
				else if (rng_int(0, 6) == 0) dy += rng_int(1, 30);
				x += dx, y += dy;
				in[n++] = anchor(rev, cid, x, y, 15);
			}
			if (rng_int(0, 1) == 0 && n < m) { /* fields drawn last-first in sequence: argument evaluation order is unspecified */
				int ty = rng_int(0, 8000), tx = rng_int(0, 30000), tcid = rng_int(0, 1), trev = rng_int(0, 1);
				in[n++] = anchor(trev, tcid, tx, ty, 15);
			}
		}
		radix_sort_128x(in, in + n);
		for (h = 0; h < (int)(sizeof(cfg) / sizeof(cfg[0])); ++h) {
			chains_t c;
			if (cfg[h].rmq) {
				rmq_inner = cfg[h].inner;
				run_chain(km, 1, in, n, 5000, 500, cfg[h].max_skip, 5000, 1, 1, 0.8f, 0.0f, 15, &c);
			} else run_chain(km, 0, in, n, 5000, 500, cfg[h].max_skip, cfg[h].max_iter, 1, 1, 0.8f, 0.0f, 15, &c);
			for (j = 0, i = 0; j < c.n_u; i += chain_cnt(&c, j), ++j) {
				int l, ok = 1;
				for (l = 1; l < chain_cnt(&c, j); ++l)
					if (c.a[i+l].x >> 32 != c.a[i].x >> 32 || (int32_t)c.a[i+l].x <= (int32_t)c.a[i+l-1].x || (int32_t)c.a[i+l].y <= (int32_t)c.a[i+l-1].y) ok = 0;
				MMT_CHECKF(ok, "set %d config %d chain %d: anchors not collinear", s, h, j);
			}
			put_chains(fp, s, h, &c);
			free_chains(km, &c);
		}
		free(in);
	}
	rmq_inner = 1000;
	fclose(fp);
	km_destroy(km);
	return 0;
}

static int cmp_u64(const void *a, const void *b)
{
	uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b;
	return x < y? -1 : x > y? 1 : 0;
}

static int cmp_128(const void *a, const void *b)
{
	const mm128_t *x = (const mm128_t*)a, *y = (const mm128_t*)b;
	return x->x < y->x? -1 : x->x > y->x? 1 : x->y < y->y? -1 : x->y > y->y;
}

static int cmp_u32(const void *a, const void *b)
{
	uint32_t x = *(const uint32_t*)a, y = *(const uint32_t*)b;
	return x < y? -1 : x > y;
}

/* args: seed. radix_sort_128x sorts by .x and keeps the multiset of (x, y) pairs; radix_sort_64 sorts;
 * ks_ksmall_uint32_t returns the k-th smallest (0-based). Sizes 0 to 100000, keys with few distinct
 * values and keys spanning 64 bits. */
static int case_sorts(int argc, char *argv[])
{
	static const int sizes[] = { 0, 1, 2, 5, 31, 64, 65, 200, 1000, 70000 };
	int si, mode;
	rng_state = argc > 0? strtoull(argv[0], 0, 10) : 1;
	for (si = 0; si < 10; ++si)
	for (mode = 0; mode < 3; ++mode) {
		int n = sizes[si], i, ok = 1;
		mm128_t *a = (mm128_t*)malloc((n + 1) * sizeof(mm128_t)), *b = (mm128_t*)malloc((n + 1) * sizeof(mm128_t));
		uint64_t *u = (uint64_t*)malloc((n + 1) * 8), *v = (uint64_t*)malloc((n + 1) * 8);
		uint32_t *w = (uint32_t*)malloc((n + 1) * 4), *ws = (uint32_t*)malloc((n + 1) * 4);
		for (i = 0; i < n; ++i) {
			uint64_t x = mode == 0? rng_next() : mode == 1? rng_next() % 7 : rng_next();
			if (mode == 2) x >>= rng_int(0, 63); /* a separate statement: operand order is unspecified */
			a[i].x = x, a[i].y = rng_next() % 1000;
			u[i] = mode == 1? rng_next() % 5 : rng_next();
			w[i] = mode == 1? (uint32_t)(rng_next() % 9) : (uint32_t)rng_next();
		}
		memcpy(b, a, n * sizeof(mm128_t));
		memcpy(v, u, n * 8);
		radix_sort_128x(a, a + n);
		for (i = 1; i < n; ++i) if (a[i-1].x > a[i].x) ok = 0;
		MMT_CHECKF(ok, "radix_sort_128x n %d mode %d: not sorted by x", n, mode);
		qsort(a, n, sizeof(mm128_t), cmp_128);
		qsort(b, n, sizeof(mm128_t), cmp_128);
		MMT_CHECKF(n == 0 || memcmp(a, b, n * sizeof(mm128_t)) == 0, "radix_sort_128x n %d mode %d: elements changed", n, mode);
		radix_sort_64(u, u + n);
		qsort(v, n, 8, cmp_u64);
		MMT_CHECKF(n == 0 || memcmp(u, v, n * 8) == 0, "radix_sort_64 n %d mode %d differs from qsort", n, mode);
		if (n > 0) {
			int k;
			memcpy(ws, w, n * 4);
			qsort(ws, n, 4, cmp_u32);
			for (k = 0; k < n; k += (n > 50? n / 7 + 1 : 1)) {
				uint32_t *t = (uint32_t*)malloc(n * 4);
				memcpy(t, w, n * 4);
				MMT_CHECKF(ks_ksmall_uint32_t(n, t, k) == ws[k], "ks_ksmall n %d mode %d k %d", n, mode, k);
				free(t);
			}
		}
		free(a); free(b); free(u); free(v); free(w); free(ws);
	}
	return 0;
}

static const mmt_case_t cases[] = {
	{ "collinear", case_collinear },
	{ "separate", case_separate },
	{ "gap_cost", case_gap_cost },
	{ "band_edge", case_band_edge },
	{ "min_cnt_sc", case_min_cnt_sc },
	{ "max_iter", case_max_iter },
	{ "random_dp", case_random_dp },
	{ "heuristics", case_heuristics },
	{ "sorts", case_sorts },
	{ 0, 0 }
};

MMT_MAIN(cases)
