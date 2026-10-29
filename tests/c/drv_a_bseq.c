/* C driver: the sequence reader of bseq.c (internal API, bseq.h).
 *
 * Each case writes the records it reads, one line per record with the batch
 * it came in, so pytest can compare them with what it wrote:
 * batch<TAB>name<TAB>l_seq<TAB>seq<TAB>qual<TAB>comment ('*' for NULL).
 *
 *   read <file|-> <chunk_size> <with_qual> <with_comment> <frag_mode> <fn> <out>
 *        fn 3: mm_bseq_read3; 2: mm_bseq_read2 (no comments); 1: mm_bseq_read (no comments, no frag_mode)
 *   frag <file1> <file2> <chunk_size> <with_qual> <with_comment> <fn> <out>
 *        fn 2: mm_bseq_read_frag2; 1: mm_bseq_read_frag (no comments)
 *   open_missing <path>
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "bseq.h"
#include "minimap2_test.h"

static void put_rec(FILE *fp, int batch, const mm_bseq1_t *s, int with_qual, int with_comment)
{
	MMT_CHECKF(s->name != 0 && s->seq != 0, "batch %d: NULL name or sequence", batch);
	if (s->name == 0 || s->seq == 0) return;
	MMT_CHECKF((int)strlen(s->seq) == s->l_seq, "%s: l_seq %d, strlen %ld", s->name, s->l_seq, (long)strlen(s->seq));
	MMT_CHECKF(with_qual || s->qual == 0, "%s: quality without with_qual", s->name);
	MMT_CHECKF(with_comment || s->comment == 0, "%s: comment without with_comment", s->name);
	MMT_CHECKF(s->qual == 0 || (int)strlen(s->qual) == s->l_seq, "%s: quality length %ld, l_seq %d", s->name, s->qual? (long)strlen(s->qual) : -1L, s->l_seq);
	fprintf(fp, "%d\t%s\t%d\t%s\t%s\t%s\n", batch, s->name, s->l_seq, s->seq, s->qual? s->qual : "*", s->comment? s->comment : "*");
}

static void free_recs(mm_bseq1_t *a, int n)
{
	int i;
	for (i = 0; i < n; ++i)
		free(a[i].name), free(a[i].seq), free(a[i].qual), free(a[i].comment);
	free(a);
}

static int case_read(int argc, char *argv[])
{
	mm_bseq_file_t *fp;
	mm_bseq1_t *a;
	FILE *out;
	int64_t chunk;
	int n, batch = 0, with_qual, with_comment, frag, fn;
	if (argc < 7) return MMT_FAIL("usage: read <file|-> <chunk_size> <with_qual> <with_comment> <frag_mode> <fn> <out>"), 1;
	chunk = atoll(argv[1]), with_qual = atoi(argv[2]), with_comment = atoi(argv[3]), frag = atoi(argv[4]), fn = atoi(argv[5]);
	fp = mm_bseq_open(argv[0]);
	MMT_REQUIRE(fp != 0);
	MMT_REQUIRE((out = fopen(argv[6], "w")) != 0);
	while (1) {
		int i;
		if (fn == 3) a = mm_bseq_read3(fp, chunk, with_qual, with_comment, frag, &n);
		else if (fn == 2) a = mm_bseq_read2(fp, chunk, with_qual, frag, &n);
		else a = mm_bseq_read(fp, chunk, with_qual, &n);
		if (n == 0) {
			free(a);
			break;
		}
		for (i = 0; i < n; ++i)
			put_rec(out, batch, &a[i], with_qual, fn == 3 && with_comment);
		free_recs(a, n);
		++batch;
		MMT_REQUIRE(batch < 100000);
	}
	MMT_CHECK(mm_bseq_eof(fp));
	a = mm_bseq_read3(fp, chunk, with_qual, with_comment, frag, &n); // reading past the end
	MMT_CHECK_INT(n, ==, 0);
	free(a);
	mm_bseq_close(fp);
	fclose(out);
	return 0;
}

static int case_frag(int argc, char *argv[])
{
	mm_bseq_file_t *fp[2];
	mm_bseq1_t *a;
	FILE *out;
	int64_t chunk;
	int n, batch = 0, with_qual, with_comment, fn;
	if (argc < 7) return MMT_FAIL("usage: frag <file1> <file2> <chunk_size> <with_qual> <with_comment> <fn> <out>"), 1;
	chunk = atoll(argv[2]), with_qual = atoi(argv[3]), with_comment = atoi(argv[4]), fn = atoi(argv[5]);
	MMT_CHECK(mm_bseq_read_frag2(0, fp, chunk, with_qual, with_comment, &n) == 0 && n == 0);
	fp[0] = mm_bseq_open(argv[0]);
	fp[1] = mm_bseq_open(argv[1]);
	MMT_REQUIRE(fp[0] != 0 && fp[1] != 0);
	MMT_REQUIRE((out = fopen(argv[6], "w")) != 0);
	while (1) {
		int i;
		if (fn == 2) a = mm_bseq_read_frag2(2, fp, chunk, with_qual, with_comment, &n);
		else a = mm_bseq_read_frag(2, fp, chunk, with_qual, &n);
		if (n == 0) {
			free(a);
			break;
		}
		MMT_CHECKF(n % 2 == 0, "batch %d holds %d records, not whole pairs", batch, n);
		for (i = 0; i < n; ++i)
			put_rec(out, batch, &a[i], with_qual, fn == 2 && with_comment);
		free_recs(a, n);
		++batch;
		MMT_REQUIRE(batch < 100000);
	}
	mm_bseq_close(fp[0]);
	mm_bseq_close(fp[1]);
	fclose(out);
	return 0;
}

static int case_open_missing(int argc, char *argv[])
{
	if (argc < 1) return MMT_FAIL("usage: open_missing <path>"), 1;
	MMT_CHECK(mm_bseq_open(argv[0]) == 0);
	return 0;
}

static const mmt_case_t cases[] = {
	{ "read", case_read },
	{ "frag", case_frag },
	{ "open_missing", case_open_missing },
	{ 0, 0 }
};

MMT_MAIN(cases)
