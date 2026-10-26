#ifndef MINIMAP2_TEST_H
#define MINIMAP2_TEST_H

/* A tiny check-and-report framework for C test drivers.
 *
 * A driver defines cases as functions `static int name(int argc, char *argv[])`
 * (argv holds the arguments after the case name), lists them in an
 * mmt_case_t array ending with {0, 0}, and ends with MMT_MAIN(that_array).
 *
 *   driver list              print the case names, one per line
 *   driver run CASE [ARGS]   run one case
 *
 * Each failed check prints a TAP-style "not ok N - ..." line with the file,
 * line and expression. Passing checks are silent unless MMT_VERBOSE is set
 * in the environment. A run ends with "# mmt: CASE checks=N failures=M" and
 * exits 0 only if no check failed and the case returned 0. pytest collects
 * this through minimap2_helpers.run_c_case() / assert_c_case().
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdint.h>

typedef int (*mmt_func_t)(int argc, char *argv[]);

typedef struct {
	const char *name;
	mmt_func_t func;
} mmt_case_t;

static int mmt_n_check = 0, mmt_n_fail = 0, mmt_verbose = 0;

static inline int mmt_report(int ok, const char *file, int line, const char *fmt, ...)
{
	va_list ap;
	++mmt_n_check;
	if (!ok) ++mmt_n_fail;
	if (ok && !mmt_verbose) return ok;
	printf("%s %d - %s:%d: ", ok? "ok" : "not ok", mmt_n_check, file, line);
	va_start(ap, fmt);
	vprintf(fmt, ap);
	va_end(ap);
	putchar('\n');
	fflush(stdout);
	return ok;
}

/* each macro evaluates to 1 if the check passed and 0 otherwise */
#define MMT_CHECK(cond) mmt_report(!!(cond), __FILE__, __LINE__, "%s", #cond)
#define MMT_CHECKF(cond, ...) mmt_report(!!(cond), __FILE__, __LINE__, __VA_ARGS__)
#define MMT_CHECK_INT(a, op, b) do { \
		long long mmt_a_ = (long long)(a), mmt_b_ = (long long)(b); \
		mmt_report(mmt_a_ op mmt_b_, __FILE__, __LINE__, "%s %s %s (%lld vs %lld)", #a, #op, #b, mmt_a_, mmt_b_); \
	} while (0)
#define MMT_CHECK_STR(a, b) do { \
		const char *mmt_a_ = (a), *mmt_b_ = (b); \
		int mmt_ok_ = (mmt_a_ && mmt_b_ && strcmp(mmt_a_, mmt_b_) == 0) || (mmt_a_ == 0 && mmt_b_ == 0); \
		mmt_report(mmt_ok_, __FILE__, __LINE__, "%s == %s (\"%s\" vs \"%s\")", #a, #b, mmt_a_? mmt_a_ : "(null)", mmt_b_? mmt_b_ : "(null)"); \
	} while (0)
#define MMT_FAIL(...) mmt_report(0, __FILE__, __LINE__, __VA_ARGS__)
/* abort the case (after reporting) if a precondition does not hold */
#define MMT_REQUIRE(cond) do { if (!MMT_CHECK(cond)) return 1; } while (0)

static inline void mmt_note(const char *fmt, ...)
{
	va_list ap;
	fputs("# ", stdout);
	va_start(ap, fmt);
	vprintf(fmt, ap);
	va_end(ap);
	putchar('\n');
}

static inline int mmt_main(int argc, char *argv[], const mmt_case_t *cases)
{
	const mmt_case_t *c;
	const char *v;
	int ret;
	if (argc >= 2 && strcmp(argv[1], "list") == 0) {
		for (c = cases; c->name; ++c) puts(c->name);
		return 0;
	}
	if (argc < 3 || strcmp(argv[1], "run") != 0) {
		fprintf(stderr, "Usage: %s list | run <case> [args...]\n", argv[0]);
		return 2;
	}
	for (c = cases; c->name; ++c)
		if (strcmp(c->name, argv[2]) == 0) break;
	if (c->name == 0) {
		fprintf(stderr, "ERROR: unknown case '%s'\n", argv[2]);
		return 2;
	}
	v = getenv("MMT_VERBOSE");
	mmt_verbose = v && *v && strcmp(v, "0") != 0;
	ret = c->func(argc - 3, argv + 3);
	if (ret != 0) mmt_note("case returned %d", ret);
	printf("# mmt: %s checks=%d failures=%d\n", c->name, mmt_n_check, mmt_n_fail);
	return mmt_n_fail > 0 || ret != 0? 1 : 0;
}

#define MMT_MAIN(cases) int main(int argc, char *argv[]) { return mmt_main(argc, argv, (cases)); }

/* read the sequences of a FASTA/FASTQ file without kseq (plain text only); for small test inputs only */
typedef struct {
	int32_t n;
	char **name, **seq;
	int32_t *len;
} mmt_seqs_t;

static inline int mmt_read_seqs(const char *fn, mmt_seqs_t *s)
{
	FILE *fp;
	char *line = 0;
	size_t cap = 0;
	int32_t m = 0, mode = 0; // mode: 1 after a FASTA header, 2 after a FASTQ header, 3 at FASTQ '+', 4 skip quality
	long l;
	memset(s, 0, sizeof(*s));
	if ((fp = fopen(fn, "r")) == 0) return -1;
	while (1) {
		int c;
		size_t n = 0;
		while ((c = fgetc(fp)) != EOF && c != '\n') {
			if (n + 2 > cap) cap = cap? cap * 2 : 256, line = (char*)realloc(line, cap);
			line[n++] = (char)c;
		}
		if (c == EOF && n == 0) break;
		if (line == 0) cap = 256, line = (char*)malloc(cap);
		line[n] = 0;
		l = (long)n;
		if (mode == 4) { mode = 0; continue; }
		if (mode == 3) { mode = 4; continue; }
		if ((line[0] == '>' || line[0] == '@') && (mode == 0 || mode == 1)) {
			char *p;
			if (s->n == m) {
				m = m? m * 2 : 16;
				s->name = (char**)realloc(s->name, m * sizeof(char*));
				s->seq = (char**)realloc(s->seq, m * sizeof(char*));
				s->len = (int32_t*)realloc(s->len, m * sizeof(int32_t));
			}
			for (p = line + 1; *p && *p != ' ' && *p != '\t'; ++p);
			*p = 0;
			s->name[s->n] = strdup(line + 1);
			s->seq[s->n] = strdup("");
			s->len[s->n++] = 0;
			mode = line[0] == '>'? 1 : 2;
		} else if (mode == 1 || mode == 2) {
			int32_t k = s->n - 1;
			s->seq[k] = (char*)realloc(s->seq[k], s->len[k] + l + 1);
			memcpy(s->seq[k] + s->len[k], line, l + 1);
			s->len[k] += (int32_t)l;
			if (mode == 2) mode = 3;
		}
		if (c == EOF) break;
	}
	free(line);
	fclose(fp);
	return s->n;
}

static inline void mmt_free_seqs(mmt_seqs_t *s)
{
	int32_t i;
	for (i = 0; i < s->n; ++i) free(s->name[i]), free(s->seq[i]);
	free(s->name); free(s->seq); free(s->len);
	memset(s, 0, sizeof(*s));
}

#endif
