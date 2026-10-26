/* Area A driver: the index API of minimap.h (mm_idx_reader_*, mm_idx_is_idx,
 * mm_idx_load/dump, mm_idx_str, mm_idx_build, mm_idx_stat, mm_idx_destroy,
 * mm_idx_index_name, mm_idx_name2id, mm_idx_getseq) and the internal index
 * functions behind it (mm_idx_init, mm_idx_gen, mm_idx_get,
 * mm_idx_getseq_rev, mm_idx_getseq2, mm_idxopt_init) and splitidx.c.
 *
 * Internal API: mm_idx_get(), mm_idx_getseq2(), mm_idxopt_init() and the
 * split functions are declared in mmpriv.h; mm_idx_init(), mm_idx_gen() and
 * mm_idx_getseq_rev() are declared in no header, so their prototypes are
 * repeated below. Checks use only the documented results: the sequences,
 * names, lengths and offsets of mm_idx_t, the positions mm_idx_get() returns
 * for a minimizer from mm_sketch(), and the base codes of mm_idx_getseq().
 *
 *   reader <target> <w> <k> <hpc> <bucket_bits> <batch_size> <out> [fn_out]
 *   lookup <target> <seqs.fa> <has_dup>
 *   dump_load <target> <out1.mmi> <out2.mmi> <not_an_index>...
 *   idx_str
 *   build <seqs.fa> <w> <k> <hpc> <missing>
 *   gen <seqs.fa> <w> <k> <bucket_bits> <batch_size> <expected_parts> <mini_batch_size>
 *   idxopt
 *   stat <target>
 *   split <seqs.fa> <batch_size> <prefix> <expected_parts>
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "mmpriv.h"
#include "bseq.h"
#include "minimap2_test.h"

mm_idx_t *mm_idx_init(int w, int k, int b, int flag);
mm_idx_t *mm_idx_gen(mm_bseq_file_t *fp, int w, int k, int b, int flag, int mini_batch_size, int n_threads, uint64_t batch_size);
int mm_idx_getseq_rev(const mm_idx_t *mi, uint32_t rid, uint32_t st, uint32_t en, uint8_t *seq);

static int nt4(int c)
{
	return c == 'A' || c == 'a'? 0 : c == 'C' || c == 'c'? 1 : c == 'G' || c == 'g'? 2 : c == 'T' || c == 't' || c == 'U' || c == 'u'? 3 : 4;
}

static long file_size(const char *fn)
{
	FILE *fp = fopen(fn, "rb");
	long sz;
	if (fp == 0) return -1;
	fseek(fp, 0, SEEK_END);
	sz = ftell(fp);
	fclose(fp);
	return sz;
}

static int cmp_u64(const void *a, const void *b)
{
	uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b;
	return x < y? -1 : x > y? 1 : 0;
}

/* every minimizer of every sequence is found by mm_idx_get() at its own position, each position list
 * is sorted, and the lists hold no other position: their total size is the number of minimizers */
static void check_lookup(const mm_idx_t *mi, const mmt_seqs_t *s, const char *what)
{
	mm128_v a = {0,0,0};
	uint64_t *h = 0;
	size_t n_h = 0, m_h = 0, total = 0, sum = 0, j;
	int32_t i, n;
	for (i = 0; i < s->n && i < (int32_t)mi->n_seq; ++i) {
		if (s->len[i] == 0) continue;
		a.n = 0;
		mm_sketch(0, s->seq[i], s->len[i], mi->w, mi->k, i, mi->flag & MM_I_HPC, &a);
		total += a.n;
		for (j = 0; j < a.n; ++j) {
			const uint64_t *p = mm_idx_get(mi, a.a[j].x >> 8, &n);
			int t, hit = 0, sorted = 1;
			for (t = 0; t < n; ++t) {
				if (p[t] == a.a[j].y) ++hit;
				if (t > 0 && p[t] <= p[t-1]) sorted = 0;
			}
			MMT_CHECKF(p != 0 && hit == 1, "%s: seq %d minimizer at %u: found %d times among %d positions", what, i, (uint32_t)a.a[j].y >> 1, hit, n);
			MMT_CHECKF(sorted, "%s: positions of the minimizer at %u are not sorted", what, (uint32_t)a.a[j].y >> 1);
			if (n_h == m_h) m_h = m_h? m_h * 2 : 1024, h = (uint64_t*)realloc(h, m_h * 8);
			h[n_h++] = a.a[j].x >> 8;
		}
	}
	if (n_h > 0) qsort(h, n_h, 8, cmp_u64);
	for (j = 0; j < n_h; ++j) {
		if (j > 0 && h[j] == h[j-1]) continue;
		mm_idx_get(mi, h[j], &n);
		sum += n;
	}
	MMT_CHECKF(sum == total, "%s: position lists hold %ld entries, the sequences have %ld minimizers", what, (long)sum, (long)total);
	for (j = 0; j < (1ULL<<20); ++j) { // the smallest hash with no minimizer
		size_t lo = 0, hi = n_h;
		while (lo < hi) { size_t mid = (lo + hi) / 2; if (h[mid] < j) lo = mid + 1; else hi = mid; }
		if (lo == n_h || h[lo] != j) break;
	}
	MMT_CHECKF(mm_idx_get(mi, j, &n) == 0 && n == 0, "%s: absent hash %ld returns %d positions", what, (long)j, n);
	free(h);
	free(a.a);
}

static void check_seqs(const mm_idx_t *mi, const mmt_seqs_t *s, const char *what)
{
	int32_t i, j, n;
	uint8_t *buf;
	uint64_t off = 0;
	int32_t max_len = 0;
	for (i = 0; i < s->n; ++i) max_len = max_len > s->len[i]? max_len : s->len[i];
	buf = (uint8_t*)malloc(max_len + 128); // room for the unclipped end a faulty mm_idx_getseq() would write
	MMT_CHECK_INT(mi->n_seq, ==, s->n);
	for (i = 0; i < s->n && i < (int32_t)mi->n_seq; ++i) {
		int32_t L = s->len[i], st, en;
		MMT_CHECKF(mi->seq[i].len == (uint32_t)L, "%s: seq %d length %u, want %d", what, i, mi->seq[i].len, L);
		MMT_CHECKF(mi->seq[i].offset == off, "%s: seq %d offset %llu, want %llu", what, i, (unsigned long long)mi->seq[i].offset, (unsigned long long)off);
		MMT_CHECKF(mi->seq[i].is_alt == 0, "%s: seq %d is_alt %u", what, i, mi->seq[i].is_alt);
		off += L;
		if (mi->S == 0 || L == 0) continue;
		n = mm_idx_getseq(mi, i, 0, L, buf);
		MMT_CHECK_INT(n, ==, L);
		for (j = 0; j < n; ++j)
			if (buf[j] != nt4(s->seq[i][j])) break;
		MMT_CHECKF(j == n, "%s: seq %d base %d is %d, want %d", what, i, j, j < n? buf[j] : -1, j < n? nt4(s->seq[i][j]) : -1);
		st = L / 3, en = 2 * L / 3; // an inner interval
		n = mm_idx_getseq(mi, i, st, en, buf);
		MMT_CHECK_INT(n, ==, en - st);
		for (j = 0; j < n; ++j)
			if (buf[j] != nt4(s->seq[i][st + j])) break;
		MMT_CHECKF(j == n, "%s: seq %d [%d,%d) base %d differs", what, i, st, en, j);
		n = mm_idx_getseq2(mi, 0, i, st, en, buf);
		MMT_CHECK_INT(n, ==, en - st);
		n = mm_idx_getseq(mi, i, L > 5? L - 5 : 0, L + 100, buf); // the end is clipped to the length
		MMT_CHECK_INT(n, ==, L > 5? 5 : L);
		MMT_CHECK_INT(mm_idx_getseq(mi, i, L, L + 10, buf), ==, -1);
		n = mm_idx_getseq_rev(mi, i, st, en, buf); // [st,en) of the reverse complement
		MMT_CHECK_INT(n, ==, en - st);
		for (j = 0; j < n; ++j) {
			int c = nt4(s->seq[i][L - 1 - (st + j)]);
			if (buf[j] != (c < 4? 3 - c : 4)) break;
		}
		MMT_CHECKF(j == n, "%s: seq %d reverse [%d,%d) base %d differs", what, i, st, en, j);
		n = mm_idx_getseq2(mi, 1, i, st, en, buf);
		MMT_CHECK_INT(n, ==, en - st);
		for (j = 0; j < n; ++j) {
			int c = nt4(s->seq[i][L - 1 - (st + j)]);
			if (buf[j] != (c < 4? 3 - c : 4)) break;
		}
		MMT_CHECKF(j == n, "%s: getseq2 reverse seq %d base %d differs", what, i, j);
		MMT_CHECK_INT(mm_idx_getseq_rev(mi, i, L, L + 1, buf), ==, -1);
	}
	MMT_CHECK_INT(mm_idx_getseq(mi, mi->n_seq, 0, 1, buf), ==, -1);
	MMT_CHECK_INT(mm_idx_getseq_rev(mi, mi->n_seq, 0, 1, buf), ==, -1);
	free(buf);
}

/* args: <target> <w> <k> <hpc> <bucket_bits> <batch_size> <out> [fn_out]; writes one line per part and per
 * sequence to out, and with fn_out the index (every part) as mm_idx_reader_open() dumps it */
static int case_reader(int argc, char *argv[])
{
	mm_idxopt_t iopt;
	mm_idx_reader_t *r;
	mm_idx_t *mi;
	int64_t is_idx;
	int32_t part = 0;
	FILE *fp;
	if (argc < 7) return MMT_FAIL("usage: reader <target> <w> <k> <hpc> <bucket_bits> <batch_size> <out>"), 1;
	mm_idxopt_init(&iopt);
	iopt.w = atoi(argv[1]), iopt.k = atoi(argv[2]), iopt.bucket_bits = atoi(argv[4]);
	if (atoi(argv[3])) iopt.flag |= MM_I_HPC;
	iopt.batch_size = strtoull(argv[5], 0, 10);
	is_idx = mm_idx_is_idx(argv[0]);
	MMT_REQUIRE(is_idx >= 0);
	r = mm_idx_reader_open(argv[0], &iopt, argc > 7? argv[7] : 0);
	MMT_REQUIRE(r != 0);
	MMT_CHECK_INT(r->is_idx != 0, ==, is_idx > 0);
	if (is_idx > 0) MMT_CHECK_INT(r->idx_size, ==, file_size(argv[0]));
	MMT_CHECK((r->fp_out != 0) == (argc > 7));
	MMT_CHECK(r->opt.k == iopt.k && r->opt.w == iopt.w && r->opt.flag == iopt.flag && r->opt.bucket_bits == iopt.bucket_bits && r->opt.batch_size == iopt.batch_size);
	MMT_CHECK(!mm_idx_reader_eof(r));
	MMT_REQUIRE((fp = fopen(argv[6], "w")) != 0);
	while ((mi = mm_idx_reader_read(r, 2)) != 0) {
		uint32_t i;
		uint64_t off = 0;
		MMT_CHECK_INT(mi->index, ==, part);
		MMT_CHECK_INT(r->n_parts, ==, part + 1);
		MMT_CHECK(mi->h == 0 && mi->n_alt == 0);
		fprintf(fp, "part\t%d\t%u\t%d\t%d\t%d\t%d\t%d\n", mi->index, mi->n_seq, mi->w, mi->k, mi->b, mi->flag, mm_idx_reader_eof(r));
		for (i = 0; i < mi->n_seq; ++i) {
			MMT_CHECK(mi->seq[i].offset == off);
			off += mi->seq[i].len;
			fprintf(fp, "seq\t%s\t%u\n", mi->seq[i].name? mi->seq[i].name : "", mi->seq[i].len);
		}
		mm_idx_destroy(mi);
		++part;
	}
	MMT_CHECK(mm_idx_reader_eof(r));
	MMT_CHECK(mm_idx_reader_read(r, 1) == 0);
	MMT_CHECK_INT(r->n_parts, ==, part);
	mm_idx_reader_close(r);
	fclose(fp);
	MMT_CHECK(mm_idx_reader_open("/nonexistent/dir/no-such-file.fa", &iopt, 0) == 0);
	return 0;
}

/* args: <target> <seqs.fa> <has_dup>: positions, bases and names of a single-part index */
static int case_lookup(int argc, char *argv[])
{
	mm_idx_reader_t *r;
	mm_idx_t *mi;
	mmt_seqs_t s;
	int32_t i, j;
	if (argc < 3) return MMT_FAIL("usage: lookup <target> <seqs.fa> <has_dup>"), 1;
	MMT_REQUIRE(mmt_read_seqs(argv[1], &s) > 0);
	r = mm_idx_reader_open(argv[0], 0, 0); // NULL options: the defaults
	MMT_REQUIRE(r != 0);
	if (!r->is_idx) MMT_CHECK(r->opt.k == 15 && r->opt.w == 10 && r->opt.bucket_bits == 14);
	mi = mm_idx_reader_read(r, 3);
	MMT_REQUIRE(mi != 0);
	MMT_CHECK(mm_idx_reader_eof(r));
	check_seqs(mi, &s, "lookup");
	check_lookup(mi, &s, "lookup");
	MMT_CHECK_INT(mm_idx_name2id(mi, s.name[0]), ==, -2); // no name hash yet
	MMT_CHECK_INT(mm_idx_index_name(mi), ==, atoi(argv[2]));
	MMT_CHECK(mi->h != 0);
	MMT_CHECK_INT(mm_idx_index_name(mi), ==, 0); // already indexed
	for (i = 0; i < s.n; ++i) {
		for (j = 0; j < i; ++j)
			if (strcmp(s.name[j], s.name[i]) == 0) break; // the first sequence with this name
		MMT_CHECKF(strcmp(mi->seq[i].name, s.name[i]) == 0, "name %d: %s, want %s", i, mi->seq[i].name, s.name[i]);
		MMT_CHECKF(mm_idx_name2id(mi, s.name[i]) == j, "name2id(%s) = %d, want %d", s.name[i], mm_idx_name2id(mi, s.name[i]), j);
	}
	MMT_CHECK_INT(mm_idx_name2id(mi, "no-such-sequence-name"), ==, -1);
	MMT_CHECK_INT(mm_idx_name2id(mi, ""), ==, -1);
	mm_idx_destroy(mi);
	mm_idx_reader_close(r);
	mmt_free_seqs(&s);
	return 0;
}

static int same_index(const mm_idx_t *a, const mm_idx_t *b, const mmt_seqs_t *s)
{
	uint32_t i;
	int ok = 1;
	mm128_v v = {0,0,0};
	if (a->w != b->w || a->k != b->k || a->b != b->b || a->flag != b->flag || a->n_seq != b->n_seq) return 0;
	for (i = 0; i < a->n_seq; ++i) {
		size_t j;
		if (strcmp(a->seq[i].name, b->seq[i].name) != 0 || a->seq[i].len != b->seq[i].len || a->seq[i].offset != b->seq[i].offset) ok = 0;
		if (s->len[i] == 0) continue;
		v.n = 0;
		mm_sketch(0, s->seq[i], s->len[i], a->w, a->k, i, a->flag & MM_I_HPC, &v);
		for (j = 0; j < v.n; ++j) {
			int na, nb;
			const uint64_t *pa = mm_idx_get(a, v.a[j].x >> 8, &na), *pb = mm_idx_get(b, v.a[j].x >> 8, &nb);
			if (na != nb || na == 0 || memcmp(pa, pb, na * 8) != 0) ok = 0;
		}
	}
	free(v.a);
	return ok;
}

/* args: <target.fa> <out1.mmi> <out2.mmi> <not_an_index>...: dump, load, dump again; the other files
 * are not indexes (e.g. a FASTA, an empty file, a 3-byte "MMI", an index of another format version) */
static int case_dump_load(int argc, char *argv[])
{
	mm_idx_reader_t *r;
	mm_idx_t *mi, *mi2;
	mmt_seqs_t s;
	FILE *fp;
	int i;
	if (argc < 4) return MMT_FAIL("usage: dump_load <target.fa> <out1.mmi> <out2.mmi> <not_an_index>..."), 1;
	MMT_REQUIRE(mmt_read_seqs(argv[0], &s) > 0);
	r = mm_idx_reader_open(argv[0], 0, 0);
	MMT_REQUIRE(r != 0);
	mi = mm_idx_reader_read(r, 1);
	MMT_REQUIRE(mi != 0);
	mm_idx_reader_close(r);
	MMT_REQUIRE((fp = fopen(argv[1], "wb")) != 0);
	mm_idx_dump(fp, mi);
	fclose(fp);
	MMT_CHECK_INT(mm_idx_is_idx(argv[1]), ==, file_size(argv[1]));
	MMT_REQUIRE((fp = fopen(argv[1], "rb")) != 0);
	mi2 = mm_idx_load(fp);
	MMT_REQUIRE(mi2 != 0);
	MMT_CHECK(fgetc(fp) == EOF); // the whole file was read
	fclose(fp);
	MMT_CHECK(same_index(mi, mi2, &s));
	check_seqs(mi2, &s, "loaded");
	check_lookup(mi2, &s, "loaded");
	MMT_CHECK_INT(mi2->index, ==, 0);
	MMT_REQUIRE((fp = fopen(argv[2], "wb")) != 0);
	mm_idx_dump(fp, mi2);
	fclose(fp);
	mm_idx_destroy(mi);
	mm_idx_destroy(mi2);
	for (i = 3; i < argc; ++i) { // not indexes
		MMT_CHECKF(mm_idx_is_idx(argv[i]) == 0, "mm_idx_is_idx(%s) = %ld", argv[i], (long)mm_idx_is_idx(argv[i]));
		MMT_REQUIRE((fp = fopen(argv[i], "rb")) != 0);
		MMT_CHECKF(mm_idx_load(fp) == 0, "mm_idx_load(%s) returned an index", argv[i]);
		fclose(fp);
	}
	MMT_CHECK_INT(mm_idx_is_idx("-"), ==, 0);
	MMT_CHECK_INT(mm_idx_is_idx("/nonexistent/dir/no-such-file.mmi"), ==, -1);
	mm_idx_destroy(0);
	mmt_free_seqs(&s);
	return 0;
}

static int case_idx_str(int argc, char *argv[])
{
	const char *seq[4] = {
		"GATTACAGCTTGACCGTAGGCATTCGATCCGTTAGCAATGCTCAGGTACCTGAGTTCAGCTAGG",
		"ATCCAGTCGATTGCAATCGTACGATCGGATTCAGCGTTTTTTTTGGGGCCCCAAACCGTTAGAAC",
		"",
		"ccgtaGGCTTAAcgATCGAUUCGGATCAGGCATTTACGCGGAATTCC"
	};
	const char *name[4] = { "s1", "s2", "empty", "s4" };
	mmt_seqs_t s;
	mm_idx_t *mi, *mi2;
	FILE *fp;
	char fn[] = "idx_str.mmi";
	uint8_t buf[128];
	int i;
	s.n = 4, s.name = (char**)name, s.seq = (char**)seq;
	s.len = (int32_t*)calloc(4, 4);
	for (i = 0; i < 4; ++i) s.len[i] = strlen(seq[i]);
	MMT_CHECK(mm_idx_str(10, 15, 0, 14, 0, seq, name) == 0); // no sequence
	mi = mm_idx_str(5, 11, 1, -1, 4, seq, 0); // HPC, default bucket bits, no names
	MMT_REQUIRE(mi != 0);
	MMT_CHECK_INT(mi->b, ==, 14);
	MMT_CHECK_INT(mi->flag, ==, MM_I_HPC | MM_I_NO_NAME);
	MMT_CHECK(mi->seq[0].name == 0 && mi->seq[3].name == 0);
	check_seqs(mi, &s, "idx_str hpc");
	check_lookup(mi, &s, "idx_str hpc");
	MMT_CHECK_INT(mm_idx_getseq(mi, 2, 0, 1, buf), ==, -1); // the empty sequence
	MMT_REQUIRE((fp = fopen(fn, "wb")) != 0);
	mm_idx_dump(fp, mi);
	fclose(fp);
	MMT_REQUIRE((fp = fopen(fn, "rb")) != 0);
	mi2 = mm_idx_load(fp);
	fclose(fp);
	MMT_REQUIRE(mi2 != 0);
	MMT_CHECK_INT(mi2->flag, ==, MM_I_HPC | MM_I_NO_NAME);
	MMT_CHECK(mi2->seq[0].name == 0 && mi2->seq[1].name == 0);
	MMT_CHECK_INT(mi2->n_seq, ==, 4);
	check_seqs(mi2, &s, "idx_str hpc reloaded");
	check_lookup(mi2, &s, "idx_str hpc reloaded");
	mm_idx_destroy(mi2);
	mm_idx_destroy(mi);
	mi = mm_idx_str(3, 5, 0, 14, 4, seq, name); // bucket bits are capped at 2k
	MMT_REQUIRE(mi != 0);
	MMT_CHECK_INT(mi->b, ==, 10);
	MMT_CHECK_INT(mi->flag, ==, 0);
	for (i = 0; i < 4; ++i) MMT_CHECK_STR(mi->seq[i].name, name[i]);
	check_seqs(mi, &s, "idx_str k5");
	check_lookup(mi, &s, "idx_str k5");
	mm_idx_destroy(mi);
	mi = mm_idx_init(0, 5, 14, MM_I_HPC); // w < 1 becomes 1; b is capped at 2k
	MMT_REQUIRE(mi != 0);
	MMT_CHECK(mi->w == 1 && mi->k == 5 && mi->b == 10 && mi->flag == MM_I_HPC && mi->n_seq == 0);
	mm_idx_destroy(mi);
	free(s.len);
	return 0;
}

/* args: <seqs.fa> <w> <k> <hpc> <missing>: the deprecated mm_idx_build() */
static int case_build(int argc, char *argv[])
{
	mm_idx_t *mi;
	mmt_seqs_t s;
	int32_t i, w, k, hpc;
	if (argc < 5) return MMT_FAIL("usage: build <seqs.fa> <w> <k> <hpc> <missing>"), 1;
	w = atoi(argv[1]), k = atoi(argv[2]), hpc = atoi(argv[3]);
	MMT_REQUIRE(mmt_read_seqs(argv[0], &s) > 0);
	mi = mm_idx_build(argv[0], w, k, hpc? MM_I_HPC : 0, 2);
	MMT_REQUIRE(mi != 0);
	MMT_CHECK(mi->w == w && mi->k == k && mi->b == 14 && mi->flag == (hpc? MM_I_HPC : 0));
	for (i = 0; i < s.n && i < (int32_t)mi->n_seq; ++i)
		MMT_CHECK_STR(mi->seq[i].name, s.name[i]);
	check_seqs(mi, &s, "build");
	check_lookup(mi, &s, "build");
	mm_idx_destroy(mi);
	MMT_CHECK(mm_idx_build(argv[4], w, k, 0, 1) == 0);
	mmt_free_seqs(&s);
	return 0;
}

/* args: <seqs.fa> <w> <k> <bucket_bits> <batch_size> <expected_parts> <mini_batch_size>: mm_idx_gen() part by part */
static int case_gen(int argc, char *argv[])
{
	mm_bseq_file_t *fp;
	mm_idx_t *mi;
	mmt_seqs_t s;
	int32_t w, k, b, mini, n_parts = 0, n_seq = 0;
	uint64_t batch;
	if (argc < 7) return MMT_FAIL("usage: gen <seqs.fa> <w> <k> <bucket_bits> <batch_size> <expected_parts> <mini_batch_size>"), 1;
	w = atoi(argv[1]), k = atoi(argv[2]), b = atoi(argv[3]), batch = strtoull(argv[4], 0, 10), mini = atoi(argv[6]);
	MMT_REQUIRE(mmt_read_seqs(argv[0], &s) > 0);
	MMT_CHECK(mm_idx_gen(0, w, k, b, 0, 1000, 1, batch) == 0);
	fp = mm_bseq_open(argv[0]);
	MMT_REQUIRE(fp != 0);
	while ((mi = mm_idx_gen(fp, w, k, b, 0, mini, 3, batch)) != 0) {
		mmt_seqs_t t;
		uint32_t i;
		t.n = mi->n_seq, t.name = s.name + n_seq, t.seq = s.seq + n_seq, t.len = s.len + n_seq; // this part's sequences
		MMT_CHECK(mi->w == w && mi->k == k && mi->b == b);
		for (i = 0; i < mi->n_seq; ++i)
			MMT_CHECK_STR(mi->seq[i].name, s.name[n_seq + i]);
		check_seqs(mi, &t, "gen");
		check_lookup(mi, &t, "gen");
		n_seq += mi->n_seq;
		++n_parts;
		mm_idx_destroy(mi);
	}
	MMT_CHECK(mm_bseq_eof(fp));
	MMT_CHECK(mm_idx_gen(fp, w, k, b, 0, 1<<18, 2, batch) == 0); // at the end of the file
	mm_bseq_close(fp);
	MMT_CHECK_INT(n_seq, ==, s.n);
	MMT_CHECK_INT(n_parts, ==, atoi(argv[5]));
	mmt_free_seqs(&s);
	return 0;
}

static int case_idxopt(int argc, char *argv[])
{
	mm_idxopt_t o;
	memset(&o, 0xff, sizeof(o));
	mm_idxopt_init(&o);
	MMT_CHECK_INT(o.k, ==, 15); // minimap2.1: -k [15]
	MMT_CHECK_INT(o.w, ==, 10); // minimap2.1: -w [10]
	MMT_CHECK_INT(o.flag, ==, 0);
	MMT_CHECK_INT(o.bucket_bits, ==, 14);
	MMT_CHECK_INT(o.mini_batch_size, ==, 50000000);
	MMT_CHECK(o.batch_size == 8000000000ULL); // minimap2.1: -I [8G]
	return 0;
}

/* args: <target>: mm_idx_stat() writes its summary to stderr; pytest checks the numbers */
static int case_stat(int argc, char *argv[])
{
	mm_idx_reader_t *r;
	mm_idx_t *mi;
	if (argc < 1) return MMT_FAIL("usage: stat <target>"), 1;
	r = mm_idx_reader_open(argv[0], 0, 0);
	MMT_REQUIRE(r != 0);
	mi = mm_idx_reader_read(r, 1);
	MMT_REQUIRE(mi != 0);
	mm_idx_stat(mi);
	mm_idx_destroy(mi);
	mm_idx_reader_close(r);
	return 0;
}

/* args: <seqs.fa> <batch_size> <prefix> <expected_parts>: mm_split_init() for each part of a multi-part
 * index, mm_split_merge_prep() reads back the names and lengths, mm_split_rm_tmp() removes the files
 * (pytest checks that none is left; their names are not checked here) */
static int case_split(int argc, char *argv[])
{
	mm_idxopt_t iopt;
	mm_idx_reader_t *r;
	mm_idx_t *mi;
	mmt_seqs_t s;
	FILE *fp[16];
	uint32_t n_seq_part[16];
	int32_t n_parts = 0, i, j, n_exp;
	const char *prefix;
	if (argc < 4) return MMT_FAIL("usage: split <seqs.fa> <batch_size> <prefix> <expected_parts>"), 1;
	prefix = argv[2], n_exp = atoi(argv[3]);
	MMT_REQUIRE(mmt_read_seqs(argv[0], &s) > 0);
	mm_idxopt_init(&iopt);
	iopt.batch_size = strtoull(argv[1], 0, 10);
	r = mm_idx_reader_open(argv[0], &iopt, 0);
	MMT_REQUIRE(r != 0);
	while ((mi = mm_idx_reader_read(r, 1)) != 0) {
		FILE *f;
		MMT_REQUIRE(n_parts < 16);
		f = mm_split_init(prefix, mi);
		MMT_REQUIRE(f != 0);
		fclose(f);
		mm_idx_destroy(mi);
		++n_parts;
	}
	mm_idx_reader_close(r);
	MMT_CHECK_INT(n_parts, ==, n_exp);
	MMT_CHECK(mm_split_merge_prep(prefix, 0, fp, n_seq_part) == 0);
	MMT_CHECK(mm_split_merge_prep(prefix, n_parts + 1, fp, n_seq_part) == 0); // one file too many
	mi = mm_split_merge_prep(prefix, n_parts, fp, n_seq_part);
	MMT_REQUIRE(mi != 0);
	MMT_CHECK_INT(mi->n_seq, ==, s.n);
	MMT_CHECK_INT(mi->k, ==, 15);
	for (i = j = 0; i < n_parts; ++i) j += n_seq_part[i];
	MMT_CHECK_INT(j, ==, s.n);
	MMT_CHECK(n_parts < 2 || n_seq_part[0] < (uint32_t)s.n);
	for (i = 0; i < s.n && i < (int32_t)mi->n_seq; ++i) {
		MMT_CHECK_STR(mi->seq[i].name, s.name[i]);
		MMT_CHECK_INT(mi->seq[i].len, ==, s.len[i]);
	}
	for (i = 0; i < n_parts; ++i) fclose(fp[i]);
	mm_idx_destroy(mi);
	mm_split_rm_tmp(prefix, n_parts);
	mmt_free_seqs(&s);
	return 0;
}

static const mmt_case_t cases[] = {
	{ "reader", case_reader },
	{ "lookup", case_lookup },
	{ "dump_load", case_dump_load },
	{ "idx_str", case_idx_str },
	{ "build", case_build },
	{ "gen", case_gen },
	{ "idxopt", case_idxopt },
	{ "stat", case_stat },
	{ "split", case_split },
	{ 0, 0 }
};

MMT_MAIN(cases)
