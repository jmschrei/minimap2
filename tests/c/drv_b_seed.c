/* C driver: seeding (seed.c), the occurrence threshold (mm_idx_cal_max_occ,
 * mm_mapopt_update) and ALT contigs (mm_idx_alt_read).
 *
 * Expected values come from the Python side (tests/helpers_b.py), which
 * computes minimizers by brute force from minimap2's documented (w,k) rule and
 * passes counts here as arguments. Calls internal functions declared in
 * mmpriv.h (mm_sketch, mm_collect_matches, mm_seed_mz_flt,
 * mm_idx_cal_max_occ) and builds mm128_v minimizer arrays by hand. */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include "mmpriv.h"
#include "kalloc.h"
#include "minimap2_test.h"

static mm_idx_t *load_index(const char *fn, int w, int k)
{
	mm_idxopt_t iopt;
	mm_mapopt_t mopt;
	mm_idx_reader_t *r;
	mm_idx_t *mi;
	mm_set_opt(0, &iopt, &mopt);
	iopt.w = w, iopt.k = k;
	r = mm_idx_reader_open(fn, &iopt, 0);
	if (r == 0) return 0;
	mi = mm_idx_reader_read(r, 1);
	mm_idx_reader_close(r);
	return mi;
}

/* args: none. mm_seed_mz_flt (--q-occ-frac): of 300 query minimizers, value A occurs 40 times, B 12
 * times and C 5 times, the rest once. A minimizer goes when its count exceeds both q_occ_max and
 * frac * 300: with q_occ_max 10 and frac 0.01 A and B go and C stays; frac 0.1 (30) keeps B;
 * q_occ_max 11 drops B, 12 keeps it; 39 drops A, 40 keeps it; frac 0 or q_occ_max 0 disables the
 * filter. The survivors keep their order. */
static int case_mz_flt(int argc, char *argv[])
{
	static const struct { int q_occ_max; float frac; int keep_a, keep_b; } T[] = {
		{ 10, 0.01f, 0, 0 }, { 10, 0.1f, 0, 1 }, { 40, 0.01f, 1, 1 }, { 39, 0.01f, 0, 1 }, { 10, 0.0f, 1, 1 }, { 0, 0.01f, 1, 1 }, { 11, 0.01f, 0, 0 }, { 12, 0.01f, 0, 1 }
	};
	void *km = km_init();
	int t;
	for (t = 0; t < (int)(sizeof(T) / sizeof(T[0])); ++t) {
		mm128_v mv = {0, 0, 0};
		size_t i, n_a = 0, n_b = 0, n_c = 0, n_c0 = 0, n_other = 0;
		int ok = 1;
		uint64_t last = 0;
		mv.m = 300, mv.a = (mm128_t*)kmalloc(km, mv.m * sizeof(mm128_t));
		for (i = 0; i < 300; ++i) {
			uint64_t v = i % 7 == 0 && n_a < 40? 0xAAAA : i % 11 == 1 && n_b < 12? 0xBBBB : i % 50 == 3? 0xCCCC : 0x10000 + i;
			if (v == 0xAAAA) ++n_a; else if (v == 0xBBBB) ++n_b; else if (v == 0xCCCC) ++n_c0;
			mv.a[mv.n].x = v << 8 | 15, mv.a[mv.n].y = (uint64_t)i << 1;
			++mv.n;
		}
		MMT_REQUIRE(n_a == 40 && n_b == 12);
		mm_seed_mz_flt(km, &mv, T[t].q_occ_max, T[t].frac);
		n_a = n_b = 0;
		for (i = 0; i < mv.n; ++i) {
			uint64_t v = mv.a[i].x >> 8;
			if (v == 0xAAAA) ++n_a; else if (v == 0xBBBB) ++n_b; else if (v == 0xCCCC) ++n_c; else ++n_other;
			if (i > 0 && mv.a[i].y <= last) ok = 0;
			last = mv.a[i].y;
		}
		MMT_CHECKF(n_a == (T[t].keep_a? 40U : 0U) && n_b == (T[t].keep_b? 12U : 0U) && n_c == n_c0 && n_other == 300 - 40 - 12 - n_c0,
				   "test %d (q_occ_max %d, frac %g): kept %zu A, %zu B, %zu C, %zu others", t, T[t].q_occ_max, T[t].frac, n_a, n_b, n_c, n_other);
		MMT_CHECKF(ok, "test %d: order changed", t);
		kfree(km, mv.a);
	}
	km_destroy(km);
	return 0;
}

/* args: target.fa w k query.fa max_occ max_max_occ dist exp_n_a exp_n_mini_pos exp_rep_len(-1: unchecked) [out]
 * mm_collect_matches on the query's minimizers: the anchors (sum of the occurrences of the kept
 * seeds), the kept query minimizers and the repetitive query length agree with brute force. With out,
 * the query positions of the kept minimizers are written to it, one per line. */
static int case_collect(int argc, char *argv[])
{
	mm_idx_t *mi;
	mmt_seqs_t s;
	mm128_v mv = {0, 0, 0};
	mm_seed_t *m;
	void *km = km_init();
	int n_m, rep_len, n_mini_pos, i, w, k, max_occ, max_max_occ, dist, want_rep;
	int64_t n_a, want_n_a;
	uint64_t *mini_pos;
	if (argc < 10) return MMT_FAIL("usage: collect <target.fa> <w> <k> <query.fa> <max_occ> <max_max_occ> <dist> <n_a> <n_mini_pos> <rep_len>"), 1;
	w = atoi(argv[1]), k = atoi(argv[2]), max_occ = atoi(argv[4]), max_max_occ = atoi(argv[5]), dist = atoi(argv[6]);
	want_n_a = atoll(argv[7]), want_rep = atoi(argv[9]);
	mi = load_index(argv[0], w, k);
	MMT_REQUIRE(mi != 0);
	MMT_REQUIRE(mmt_read_seqs(argv[3], &s) == 1);
	mm_sketch(km, s.seq[0], s.len[0], mi->w, mi->k, 0, 0, &mv);
	m = mm_collect_matches(km, &n_m, s.len[0], max_occ, max_max_occ, dist, mi, &mv, &n_a, &rep_len, &n_mini_pos, &mini_pos);
	MMT_CHECKF(n_a == want_n_a, "anchors %lld, brute force %lld", (long long)n_a, (long long)want_n_a);
	MMT_CHECKF(n_mini_pos == atoi(argv[8]) && n_m == n_mini_pos, "kept minimizers %d (n_m %d), brute force %s", n_mini_pos, n_m, argv[8]);
	if (want_rep >= 0) MMT_CHECKF(rep_len == want_rep, "rep_len %d, brute force %d", rep_len, want_rep);
	for (i = 1; i < n_mini_pos; ++i)
		MMT_CHECKF((uint32_t)mini_pos[i] > (uint32_t)mini_pos[i-1] && (mini_pos[i] >> 32) == (uint64_t)k, "mini_pos %d not increasing or span not k", i);
	for (i = 0; i < n_m; ++i)
		MMT_CHECKF(m[i].n > 0 && (int)m[i].n <= (dist > 0 && max_max_occ > max_occ? max_max_occ : max_occ), "seed %d has %u occurrences", i, m[i].n);
	if (argc > 10) {
		FILE *fp = fopen(argv[10], "w");
		MMT_REQUIRE(fp != 0);
		for (i = 0; i < n_mini_pos; ++i) fprintf(fp, "%u\n", (uint32_t)mini_pos[i]);
		fclose(fp);
	}
	kfree(km, m); kfree(km, mini_pos); kfree(km, mv.a);
	km_destroy(km);
	mmt_free_seqs(&s);
	mm_idx_destroy(mi);
	return 0;
}

/* args: target.fa w k f expected [f expected ...]. mm_idx_cal_max_occ: one more than the occurrence
 * count at rank floor((1-f) n) of the n distinct minimizers' counts (-f FLOAT); f <= 0 gives INT32_MAX. */
static int case_cal_max_occ(int argc, char *argv[])
{
	mm_idx_t *mi;
	int i;
	if (argc < 3) return MMT_FAIL("usage: cal_max_occ <target.fa> <w> <k> [<f> <expected>]..."), 1;
	mi = load_index(argv[0], atoi(argv[1]), atoi(argv[2]));
	MMT_REQUIRE(mi != 0);
	for (i = 3; i + 1 < argc; i += 2)
		MMT_CHECKF(mm_idx_cal_max_occ(mi, atof(argv[i])) == atoi(argv[i+1]), "f %s: %d, expected %s", argv[i], mm_idx_cal_max_occ(mi, atof(argv[i])), argv[i+1]);
	MMT_CHECK_INT(mm_idx_cal_max_occ(mi, 0.0f), ==, INT32_MAX);
	MMT_CHECK_INT(mm_idx_cal_max_occ(mi, -1.0f), ==, INT32_MAX);
	mm_idx_destroy(mi);
	return 0;
}

/* args: target.fa w k f min_mid_occ max_mid_occ expected_mid_occ. mm_mapopt_update (minimap.h):
 * mid_occ = max{-U INT1, min{-U INT2, the -f threshold}}; a positive mid_occ (-f INT) is kept; the
 * long-join band is raised to the chaining band; the splice strand flags imply MM_F_SPLICE. */
static int case_mapopt_update(int argc, char *argv[])
{
	mm_idxopt_t iopt;
	mm_mapopt_t mopt;
	mm_idx_t *mi;
	if (argc < 7) return MMT_FAIL("usage: mapopt_update <target.fa> <w> <k> <f> <min> <max> <expected>"), 1;
	mi = load_index(argv[0], atoi(argv[1]), atoi(argv[2]));
	MMT_REQUIRE(mi != 0);
	mm_verbose = 1;
	mm_set_opt(0, &iopt, &mopt);
	mopt.mid_occ_frac = atof(argv[3]), mopt.min_mid_occ = atoi(argv[4]), mopt.max_mid_occ = atoi(argv[5]);
	mm_mapopt_update(&mopt, mi);
	MMT_CHECKF(mopt.mid_occ == atoi(argv[6]), "mid_occ %d, expected %s (f %s, -U %s,%s)", mopt.mid_occ, argv[6], argv[3], argv[4], argv[5]);
	mm_set_opt(0, &iopt, &mopt);
	mopt.mid_occ = 7, mopt.bw = 900, mopt.bw_long = 300; /* -f 7 -r 900,300 */
	mm_mapopt_update(&mopt, mi);
	MMT_CHECK_INT(mopt.mid_occ, ==, 7);
	MMT_CHECK_INT(mopt.bw_long, ==, 900);
	mm_set_opt(0, &iopt, &mopt);
	mopt.flag |= MM_F_SPLICE_FOR;
	mm_mapopt_update(&mopt, mi);
	MMT_CHECK(mopt.flag & MM_F_SPLICE);
	mm_set_opt(0, &iopt, &mopt);
	mm_mapopt_update(&mopt, mi);
	MMT_CHECK(!(mopt.flag & MM_F_SPLICE));
	mm_idx_destroy(mi);
	return 0;
}

/* args: target.fa alt_list n_expected name... mm_idx_alt_read (minimap.h): every listed name (first
 * word of a line) that is a target sequence becomes ALT; unknown names are ignored; returns and
 * stores the number found. */
static int case_alt_read(int argc, char *argv[])
{
	mm_idx_t *mi;
	int i, j, n;
	if (argc < 3) return MMT_FAIL("usage: alt_read <target.fa> <alt.txt> <n> [names...]"), 1;
	mi = load_index(argv[0], 10, 15);
	MMT_REQUIRE(mi != 0);
	mm_verbose = 1;
	n = mm_idx_alt_read(mi, argv[1]);
	MMT_CHECK_INT(n, ==, atoi(argv[2]));
	MMT_CHECK_INT(mi->n_alt, ==, atoi(argv[2]));
	for (i = 0; i < (int)mi->n_seq; ++i) {
		int listed = 0;
		for (j = 3; j < argc; ++j)
			if (strcmp(argv[j], mi->seq[i].name) == 0) listed = 1;
		MMT_CHECKF(mi->seq[i].is_alt == listed, "%s: is_alt %d, listed %d", mi->seq[i].name, mi->seq[i].is_alt, listed);
	}
	MMT_CHECK(mm_idx_alt_read(mi, "/nonexistent/alt-list") < 0);
	mm_idx_destroy(mi);
	return 0;
}

static const mmt_case_t cases[] = {
	{ "mz_flt", case_mz_flt },
	{ "collect", case_collect },
	{ "cal_max_occ", case_cal_max_occ },
	{ "mapopt_update", case_mapopt_update },
	{ "alt_read", case_alt_read },
	{ 0, 0 }
};

MMT_MAIN(cases)
