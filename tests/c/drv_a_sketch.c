/* C driver: mm_sketch() (sketch.c) and the SDUST functions (sdust.c).
 *
 * Internal API: mm_sketch() is declared in mmpriv.h and returns mm128_t
 * arrays (x = hash<<8 | span, y = rid<<32 | lastPos<<1 | strand, as
 * documented above mm_sketch()); the sdust functions are declared in sdust.h.
 * The cases write what the functions return to a file, and the Python tests
 * compare it with a reference minimizer scan and a reference SDUST.
 *
 *   dump <seqs.fa> <w> <k> <hpc> <out>      minimizers of each sequence
 *   km_append <seqs.fa> <w> <k> <hpc>       kalloc vs malloc; appending to a non-empty vector
 *   sdust <seqs.fa> <T> <W> <out>           sdust() and sdust_core() with a reused buffer
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "mmpriv.h"
#include "kalloc.h"
#include "sdust.h"
#include "minimap2_test.h"

static int case_dump(int argc, char *argv[])
{
	mmt_seqs_t s;
	mm128_v a = {0,0,0};
	int32_t i, w, k, hpc;
	FILE *fp;
	if (argc < 5) return MMT_FAIL("usage: dump <seqs.fa> <w> <k> <hpc> <out>"), 1;
	w = atoi(argv[1]), k = atoi(argv[2]), hpc = atoi(argv[3]);
	MMT_REQUIRE(mmt_read_seqs(argv[0], &s) > 0);
	MMT_REQUIRE((fp = fopen(argv[4], "w")) != 0);
	for (i = 0; i < s.n; ++i) {
		size_t j;
		uint32_t last = 0;
		if (s.len[i] == 0) continue; // mm_sketch() requires len > 0
		a.n = 0;
		mm_sketch(0, s.seq[i], s.len[i], w, k, i, hpc, &a);
		for (j = 0; j < a.n; ++j) {
			uint64_t h = a.a[j].x >> 8;
			uint32_t span = a.a[j].x & 0xff, rid = a.a[j].y >> 32, pos = (uint32_t)a.a[j].y >> 1, strand = a.a[j].y & 1;
			MMT_CHECKF(rid == (uint32_t)i, "seq %d minimizer %ld: rid %u", i, (long)j, rid);
			MMT_CHECKF(pos < (uint32_t)s.len[i], "seq %d minimizer %ld: position %u beyond length %d", i, (long)j, pos, s.len[i]);
			MMT_CHECKF(j == 0 || pos >= last, "seq %d minimizer %ld: position %u after %u", i, (long)j, pos, last);
			MMT_CHECKF(h < (1ULL << 2 * k), "seq %d minimizer %ld: hash %llu has more than 2k bits", i, (long)j, (unsigned long long)h);
			MMT_CHECKF(hpc? span >= (uint32_t)k : span == (uint32_t)k, "seq %d minimizer %ld: span %u with k=%d hpc=%d", i, (long)j, span, k, hpc);
			last = pos;
			fprintf(fp, "%d\t%u\t%u\t%llu\t%u\n", i, pos, strand, (unsigned long long)h, span);
		}
	}
	fclose(fp);
	free(a.a);
	mmt_free_seqs(&s);
	return 0;
}

static int same_v(const mm128_v *a, size_t off, const mm128_v *b)
{
	size_t j;
	for (j = 0; j < b->n; ++j)
		if (a->a[off + j].x != b->a[j].x || a->a[off + j].y != b->a[j].y) return 0;
	return 1;
}

static int case_km_append(int argc, char *argv[])
{
	mmt_seqs_t s;
	int32_t i, w, k, hpc;
	void *km;
	if (argc < 4) return MMT_FAIL("usage: km_append <seqs.fa> <w> <k> <hpc>"), 1;
	w = atoi(argv[1]), k = atoi(argv[2]), hpc = atoi(argv[3]);
	MMT_REQUIRE(mmt_read_seqs(argv[0], &s) > 0);
	km = km_init();
	for (i = 0; i < s.n; ++i) {
		mm128_v a = {0,0,0}, b = {0,0,0}, c = {0,0,0};
		if (s.len[i] == 0) continue;
		mm_sketch(0, s.seq[i], s.len[i], w, k, i, hpc, &a);
		mm_sketch(km, s.seq[i], s.len[i], w, k, i, hpc, &b); // kalloc pool instead of malloc
		MMT_CHECK_INT(b.n, ==, a.n);
		MMT_CHECKF(b.n == a.n && same_v(&a, 0, &b), "seq %d: kalloc and malloc sketches differ", i);
		mm_sketch(0, s.seq[i], s.len[i], w, k, i + 1, hpc, &c); // a different rid in y's high 32 bits
		mm_sketch(0, s.seq[i], s.len[i], w, k, i, hpc, &c); // appended after the first sketch
		MMT_CHECK_INT(c.n, ==, 2 * a.n);
		if (c.n == 2 * a.n) {
			size_t j;
			int ok = 1;
			for (j = 0; j < a.n; ++j)
				if (c.a[j].x != a.a[j].x || c.a[j].y != (a.a[j].y + (1ULL << 32))) ok = 0;
			MMT_CHECKF(ok, "seq %d: sketch with rid %d differs from rid %d other than in y>>32", i, i + 1, i);
			MMT_CHECKF(same_v(&c, a.n, &a), "seq %d: appended sketch differs", i);
		}
		free(a.a); kfree(km, b.a); free(c.a);
	}
	km_destroy(km);
	mmt_free_seqs(&s);
	return 0;
}

static int case_sdust(int argc, char *argv[])
{
	mmt_seqs_t s;
	int32_t i, T, W;
	FILE *fp;
	void *km;
	sdust_buf_t *buf;
	if (argc < 4) return MMT_FAIL("usage: sdust <seqs.fa> <T> <W> <out>"), 1;
	T = atoi(argv[1]), W = atoi(argv[2]);
	MMT_REQUIRE(mmt_read_seqs(argv[0], &s) > 0);
	MMT_REQUIRE((fp = fopen(argv[3], "w")) != 0);
	km = km_init();
	buf = sdust_buf_init(km); // one buffer reused for every sequence
	for (i = 0; i < s.n; ++i) {
		int n1, n2, j;
		uint64_t *r1;
		const uint64_t *r2;
		r1 = sdust(0, (const uint8_t*)s.seq[i], -1, T, W, &n1); // l_seq < 0: up to the NUL
		r2 = sdust_core((const uint8_t*)s.seq[i], s.len[i], T, W, &n2, buf);
		MMT_CHECK_INT(n2, ==, n1);
		MMT_CHECKF(n1 != n2 || n1 == 0 || memcmp(r1, r2, n1 * 8) == 0, "seq %d: sdust() and sdust_core() differ", i);
		for (j = 0; j < n1; ++j) {
			int st = (int)(r1[j] >> 32), en = (int)(uint32_t)r1[j];
			MMT_CHECKF(st < en, "seq %d interval %d: %d-%d", i, j, st, en);
			MMT_CHECKF(j == 0 || st > (int)(uint32_t)r1[j-1], "seq %d interval %d starts at %d, not after the previous end", i, j, st);
			fprintf(fp, "%d\t%d\t%d\n", i, st, en);
		}
		free(r1);
	}
	sdust_buf_destroy(buf);
	km_destroy(km);
	fclose(fp);
	mmt_free_seqs(&s);
	return 0;
}

static const mmt_case_t cases[] = {
	{ "dump", case_dump },
	{ "km_append", case_km_append },
	{ "sdust", case_sdust },
	{ 0, 0 }
};

MMT_MAIN(cases)
