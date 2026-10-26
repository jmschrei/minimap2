/* Area D driver: the libminimap2 API (mm_set_opt, mm_mapopt_init, mm_check_opt, the thread
 * buffers, mm_map, mm_map_file, mm_gen_cs/ds/MD), the writers of format.c, kalloc, kthread and
 * misc.c.
 *
 * Internal names used besides minimap.h: mm_bseq1_t, mm_bseq_open/read3/close (bseq.h);
 * mm_write_sam_hdr, mm_write_sam/sam2/sam3, mm_write_paf/paf3/paf4, kstring_t,
 * MM_DBG_NO_KALLOC, realtime, cputime, peakrss, mm_err_* (mmpriv.h); kalloc.h; kthread.h. The
 * long_cigar case builds an mm_reg1_t and its mm_extra_t by hand (both public structs of
 * minimap.h, fixed for the speed search). Checks are on documented results (output text, return
 * values, memory contents), not on internal state.
 *
 * Cases printing alignments write them to stdout, one record per line, before the final
 * "# mmt:" line; the Python side compares them with the minimap2 command line.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include "minimap.h"
#include "mmpriv.h"
#include "bseq.h"
#include "kalloc.h"
#include "kthread.h"
#include "minimap2_test.h"

static int64_t parse_flag(const char *s)
{
	return (int64_t)strtoull(s, 0, 0);
}

/* the options of the command line `minimap2 [-x preset]`, as main.c sets them */
static int set_cli_opt(const char *preset, mm_idxopt_t *io, mm_mapopt_t *mo)
{
	if (mm_set_opt(0, io, mo) != 0) return -1;
	if (preset && strcmp(preset, "-") != 0 && mm_set_opt(preset, io, mo) != 0) return -1;
	return 0;
}

static mm_idx_t *load_index(const char *fn, const mm_idxopt_t *io)
{
	mm_idx_reader_t *r;
	mm_idx_t *mi;
	if ((r = mm_idx_reader_open(fn, io, 0)) == 0) return 0;
	mi = mm_idx_reader_read(r, 1);
	if (mi && !mm_idx_reader_eof(r)) {
		mm_idx_destroy(mi);
		mi = 0;
	}
	mm_idx_reader_close(r);
	return mi;
}

static mm_bseq1_t *read_all(const char *fn, int *n)
{
	mm_bseq_file_t *fp;
	mm_bseq1_t *s;
	*n = 0;
	if ((fp = mm_bseq_open(fn)) == 0) return 0;
	s = mm_bseq_read3(fp, 1LL<<40, 1, 1, 0, n);
	mm_bseq_close(fp);
	return s;
}

static void free_seqs(int n, mm_bseq1_t *s)
{
	int i;
	for (i = 0; i < n; ++i) {
		free(s[i].name); free(s[i].seq); free(s[i].qual); free(s[i].comment);
	}
	free(s);
}

static void free_regs(int n_reg, mm_reg1_t *reg)
{
	int j;
	for (j = 0; j < n_reg; ++j) free(reg[j].p);
	free(reg);
}

/**********************************************
 * Macros: values the CLI tests map to options *
 **********************************************/

/* args: none. Prints "name value" for the MM_F_* flags of area D; checks the CIGAR codes */
static int case_macros(int argc, char *argv[])
{
	const char *bam_ops = "MIDNSHP=X"; /* SAM spec 4.2: BAM op codes 0..8 */
	int i;
	printf("MM_F_OUT_SAM 0x%llx\n", (unsigned long long)MM_F_OUT_SAM);
	printf("MM_F_NO_QUAL 0x%llx\n", (unsigned long long)MM_F_NO_QUAL);
	printf("MM_F_OUT_CG 0x%llx\n", (unsigned long long)MM_F_OUT_CG);
	printf("MM_F_CIGAR 0x%llx\n", (unsigned long long)MM_F_CIGAR);
	printf("MM_F_OUT_CS 0x%llx\n", (unsigned long long)MM_F_OUT_CS);
	printf("MM_F_OUT_CS_LONG 0x%llx\n", (unsigned long long)MM_F_OUT_CS_LONG);
	printf("MM_F_NO_PRINT_2ND 0x%llx\n", (unsigned long long)MM_F_NO_PRINT_2ND);
	printf("MM_F_2_IO_THREADS 0x%llx\n", (unsigned long long)MM_F_2_IO_THREADS);
	printf("MM_F_LONG_CIGAR 0x%llx\n", (unsigned long long)MM_F_LONG_CIGAR);
	printf("MM_F_SOFTCLIP 0x%llx\n", (unsigned long long)MM_F_SOFTCLIP);
	printf("MM_F_OUT_MD 0x%llx\n", (unsigned long long)MM_F_OUT_MD);
	printf("MM_F_COPY_COMMENT 0x%llx\n", (unsigned long long)MM_F_COPY_COMMENT);
	printf("MM_F_EQX 0x%llx\n", (unsigned long long)MM_F_EQX);
	printf("MM_F_PAF_NO_HIT 0x%llx\n", (unsigned long long)MM_F_PAF_NO_HIT);
	printf("MM_F_SAM_HIT_ONLY 0x%llx\n", (unsigned long long)MM_F_SAM_HIT_ONLY);
	printf("MM_F_SECONDARY_SEQ 0x%llx\n", (unsigned long long)MM_F_SECONDARY_SEQ);
	printf("MM_F_OUT_DS 0x%llx\n", (unsigned long long)MM_F_OUT_DS);
	printf("MM_F_OUT_JUNC 0x%llx\n", (unsigned long long)MM_F_OUT_JUNC);
	printf("MM_F_SPLICE 0x%llx\n", (unsigned long long)MM_F_SPLICE);
	printf("MM_VERSION %s\n", MM_VERSION);
	MMT_CHECK_INT(MM_CIGAR_MATCH, ==, 0);
	MMT_CHECK_INT(MM_CIGAR_INS, ==, 1);
	MMT_CHECK_INT(MM_CIGAR_DEL, ==, 2);
	MMT_CHECK_INT(MM_CIGAR_N_SKIP, ==, 3);
	MMT_CHECK_INT(MM_CIGAR_SOFTCLIP, ==, 4);
	MMT_CHECK_INT(MM_CIGAR_HARDCLIP, ==, 5);
	MMT_CHECK_INT(MM_CIGAR_PADDING, ==, 6);
	MMT_CHECK_INT(MM_CIGAR_EQ_MATCH, ==, 7);
	MMT_CHECK_INT(MM_CIGAR_X_MISMATCH, ==, 8);
	for (i = 0; i < 9; ++i)
		MMT_CHECKF(MM_CIGAR_STR[i] == bam_ops[i], "MM_CIGAR_STR[%d] is '%c', BAM code %d is '%c'", i, MM_CIGAR_STR[i], i, bam_ops[i]);
	return 0;
}

/************************
 * Options and checking *
 ************************/

/* args: none. Defaults of mm_set_opt(0) against the usage text and minimap2.1, mm_mapopt_init,
 * the documented preset names, and presets touching only their own parameters */
static int case_set_opt(int argc, char *argv[])
{
	/* the presets of minimap2.1 and the usage text, plus lr:hqae (options.c) */
	static const char *presets[] = { "map-ont", "lr:hq", "lr:hqae", "map-hifi", "map-pb", "map-iclr", "asm5", "asm10", "asm20",
		"splice", "splice:hq", "splice:sr", "sr", "ava-pb", "ava-ont", 0 };
	static const char *bad[] = { "foo", "asm7", "asm", "", "map-ont ", "MAP-ONT", "map_ont", "splice:lr", "srr", 0 };
	mm_idxopt_t io, io2;
	mm_mapopt_t mo, mo2;
	int i;
	MMT_CHECK_INT(mm_set_opt(0, &io, &mo), ==, 0);
	MMT_CHECK_INT(io.k, ==, 15); /* -k [15] */
	MMT_CHECK_INT(io.w, ==, 10); /* -w [10] */
	MMT_CHECK(!(io.flag & MM_I_HPC));
	MMT_CHECK(io.batch_size == 8000000000ULL); /* -I [8G] */
	MMT_CHECK(mo.flag == 0); /* PAF without CIGAR, secondary alignments shown */
	MMT_CHECK_INT(mo.seed, ==, 11); /* --seed [11] */
	MMT_CHECK_INT(mo.a, ==, 2); MMT_CHECK_INT(mo.b, ==, 4); /* -A [2] -B [4] */
	MMT_CHECK_INT(mo.q, ==, 4); MMT_CHECK_INT(mo.q2, ==, 24); /* -O [4,24] */
	MMT_CHECK_INT(mo.e, ==, 2); MMT_CHECK_INT(mo.e2, ==, 1); /* -E [2,1] */
	MMT_CHECK_INT(mo.zdrop, ==, 400); MMT_CHECK_INT(mo.zdrop_inv, ==, 200); /* -z [400,200] */
	MMT_CHECK_INT(mo.bw, ==, 500); MMT_CHECK_INT(mo.bw_long, ==, 20000); /* -r [500,20k] */
	MMT_CHECK_INT(mo.min_cnt, ==, 3); /* -n [3] */
	MMT_CHECK_INT(mo.min_chain_score, ==, 40); /* -m [40] (usage) */
	MMT_CHECK_INT(mo.best_n, ==, 5); /* -N [5] (usage) */
	MMT_CHECK(mo.pri_ratio == 0.8f); /* -p [0.8] */
	MMT_CHECK(mo.mid_occ_frac == 2e-4f); /* -f [0.0002] (usage) */
	MMT_CHECK_INT(mo.min_mid_occ, ==, 10); MMT_CHECK_INT(mo.max_mid_occ, ==, 1000000); /* -U [10,1000000] */
	MMT_CHECK(mo.q_occ_frac == 0.01f); /* --q-occ-frac [0.01] */
	MMT_CHECK_INT(mo.occ_dist, ==, 500); /* -e [500] */
	MMT_CHECK_INT(mo.max_chain_skip, ==, 25); /* --max-chain-skip [25] */
	MMT_CHECK_INT(mo.sc_ambi, ==, 1); /* --score-N [1] */
	MMT_CHECK(mo.mini_batch_size == 500000000); /* -K [500M] */
	MMT_CHECK(mo.cap_kalloc == 500000000); /* --cap-kalloc [500m] */
	MMT_CHECK(mo.max_sw_mat == 100000000); /* --cap-sw-mem [100m] */
	MMT_CHECK(mo.split_prefix == 0);
	/* mm_mapopt_init() is the deprecated form of mm_set_opt(0, ...) for the mapping options */
	memset(&mo2, 0xab, sizeof(mo2));
	mm_mapopt_init(&mo2);
	MMT_CHECK(memcmp(&mo, &mo2, sizeof(mo)) == 0);
	/* presets are accepted; unknown names return -1 */
	for (i = 0; presets[i]; ++i) {
		mm_set_opt(0, &io, &mo);
		MMT_CHECKF(mm_set_opt(presets[i], &io, &mo) == 0, "preset '%s' rejected", presets[i]);
		MMT_CHECKF(mm_check_opt(&io, &mo) == 0, "preset '%s' fails mm_check_opt", presets[i]);
	}
	for (i = 0; bad[i]; ++i) {
		mm_set_opt(0, &io, &mo);
		MMT_CHECKF(mm_set_opt(bad[i], &io, &mo) == -1, "unknown preset '%s' accepted", bad[i]);
	}
	/* minimap.h: a preset changes only the parameters it affects */
	for (i = 0; presets[i]; ++i) {
		mm_set_opt(0, &io, &mo);
		mo.seed = 99, mo.cap_kalloc = 12345, mo.max_qlen = 777, mo.flag |= MM_F_OUT_SAM | MM_F_OUT_CS | MM_F_COPY_COMMENT;
		MMT_REQUIRE(mm_set_opt(presets[i], &io, &mo) == 0);
		MMT_CHECKF(mo.seed == 99 && mo.cap_kalloc == 12345 && mo.max_qlen == 777, "preset '%s' changed seed/cap_kalloc/max_qlen", presets[i]);
		MMT_CHECKF((mo.flag & (MM_F_OUT_SAM | MM_F_OUT_CS | MM_F_COPY_COMMENT)) == (MM_F_OUT_SAM | MM_F_OUT_CS | MM_F_COPY_COMMENT),
			"preset '%s' cleared output flags", presets[i]);
	}
	/* minimap2.1: sr and splice:sr include --secondary=no and -2 with -K50m and -K100m; no other preset does */
	for (i = 0; presets[i]; ++i) {
		int sr = strcmp(presets[i], "sr") == 0, ssr = strcmp(presets[i], "splice:sr") == 0;
		mm_set_opt(0, &io, &mo);
		mm_set_opt(presets[i], &io, &mo);
		MMT_CHECKF(!!(mo.flag & MM_F_NO_PRINT_2ND) == (sr || ssr), "preset '%s': MM_F_NO_PRINT_2ND", presets[i]);
		MMT_CHECKF(!!(mo.flag & MM_F_2_IO_THREADS) == (sr || ssr), "preset '%s': MM_F_2_IO_THREADS", presets[i]);
		if (sr) MMT_CHECK(mo.mini_batch_size == 50000000);
		else if (ssr) MMT_CHECK(mo.mini_batch_size == 100000000);
		else MMT_CHECKF(mo.mini_batch_size == 500000000, "preset '%s' changed -K", presets[i]);
	}
	/* mm_set_opt(0) resets everything set before */
	mm_set_opt("sr", &io, &mo);
	mo.seed = 5;
	MMT_CHECK_INT(mm_set_opt(0, &io2, &mo2), ==, 0);
	MMT_CHECK_INT(mm_set_opt(0, &io, &mo), ==, 0);
	MMT_CHECK(memcmp(&mo, &mo2, sizeof(mo)) == 0 && memcmp(&io, &io2, sizeof(io)) == 0);
	return 0;
}

static int n_bad_checks;

#define OPT_EXPECT(want_ok, ...) do { \
		mm_idxopt_t io; mm_mapopt_t mo; int ret; \
		mm_set_opt(0, &io, &mo); \
		__VA_ARGS__; \
		ret = mm_check_opt(&io, &mo); \
		if (want_ok) MMT_CHECKF(ret == 0, "rejected (%d): %s", ret, #__VA_ARGS__); \
		else MMT_CHECKF(ret < 0, "accepted: %s", #__VA_ARGS__), ++n_bad_checks; \
	} while (0)

/* args: verbosity (mm_verbose). mm_check_opt on each rule of options.c and its boundary */
static int case_check_opt(int argc, char *argv[])
{
	if (argc < 1) return MMT_FAIL("usage: check_opt <mm_verbose>"), 1;
	mm_verbose = atoi(argv[0]);
	OPT_EXPECT(1, (void)0);
	OPT_EXPECT(0, mo.bw = 1000, mo.bw_long = 999);
	OPT_EXPECT(1, mo.bw = 1000, mo.bw_long = 1000);
	OPT_EXPECT(0, mo.flag |= MM_F_RMQ | MM_F_SR);
	OPT_EXPECT(0, mo.flag |= MM_F_RMQ | MM_F_SPLICE);
	OPT_EXPECT(1, mo.flag |= MM_F_RMQ);
	OPT_EXPECT(0, mo.split_prefix = "pfx", mo.flag |= MM_F_OUT_CS);
	OPT_EXPECT(0, mo.split_prefix = "pfx", mo.flag |= MM_F_OUT_MD);
	OPT_EXPECT(1, mo.split_prefix = "pfx", mo.flag |= MM_F_OUT_SAM | MM_F_CIGAR);
	OPT_EXPECT(0, io.k = 0);
	OPT_EXPECT(0, io.w = 0);
	OPT_EXPECT(1, io.k = 1, io.w = 1);
	OPT_EXPECT(0, mo.best_n = -1);
	OPT_EXPECT(1, mo.best_n = 0);
	OPT_EXPECT(0, mo.pri_ratio = 1.001f);
	OPT_EXPECT(0, mo.pri_ratio = -0.001f);
	OPT_EXPECT(1, mo.pri_ratio = 0.0f);
	OPT_EXPECT(1, mo.pri_ratio = 1.0f);
	OPT_EXPECT(0, mo.flag |= MM_F_FOR_ONLY | MM_F_REV_ONLY);
	OPT_EXPECT(0, mo.q = mo.q2 = 0);
	OPT_EXPECT(0, mo.e = mo.e2 = 0);
	OPT_EXPECT(0, mo.e = 1, mo.e2 = 1);
	OPT_EXPECT(0, mo.q = 4, mo.q2 = 4, mo.e = 2, mo.e2 = 1);
	OPT_EXPECT(0, mo.q = 4, mo.q2 = 5, mo.e = 2, mo.e2 = 1);
	OPT_EXPECT(1, mo.q = 4, mo.q2 = 6, mo.e = 2, mo.e2 = 1);
	OPT_EXPECT(1, mo.q = mo.q2 = 60, mo.e = mo.e2 = 3);
	OPT_EXPECT(0, mo.q = mo.q2 = 60, mo.e = mo.e2 = 4);
	OPT_EXPECT(1, mo.q = 4, mo.e = 2, mo.q2 = 120, mo.e2 = 1);
	OPT_EXPECT(0, mo.q = 4, mo.e = 2, mo.q2 = 121, mo.e2 = 1);
	OPT_EXPECT(0, mo.sc_ambi = -1);
	OPT_EXPECT(0, mo.sc_ambi = mo.b);
	OPT_EXPECT(1, mo.sc_ambi = mo.b - 1);
	OPT_EXPECT(1, mo.sc_ambi = 0);
	OPT_EXPECT(0, mo.zdrop = 199, mo.zdrop_inv = 200);
	OPT_EXPECT(1, mo.zdrop = 200, mo.zdrop_inv = 200);
	OPT_EXPECT(0, mo.flag |= MM_F_NO_PRINT_2ND | MM_F_ALL_CHAINS);
	OPT_EXPECT(1, mo.flag |= MM_F_ALL_CHAINS);
	OPT_EXPECT(0, mo.flag |= MM_F_QSTRAND | MM_F_OUT_SAM);
	OPT_EXPECT(0, mo.flag |= MM_F_QSTRAND | MM_F_SPLICE);
	OPT_EXPECT(0, mo.flag |= MM_F_QSTRAND | MM_F_FRAG_MODE);
	OPT_EXPECT(0, mo.flag |= MM_F_QSTRAND, io.flag |= MM_I_HPC);
	OPT_EXPECT(1, mo.flag |= MM_F_QSTRAND | MM_F_CIGAR | MM_F_OUT_CG);
	printf("rejected %d\n", n_bad_checks);
	return 0;
}

/******************
 * Thread buffers *
 ******************/

/* args: none */
static int case_tbuf(int argc, char *argv[])
{
	mm_tbuf_t *b;
	void *km, *p;
	b = mm_tbuf_init();
	MMT_REQUIRE(b != 0);
	km = mm_tbuf_get_km(b);
	MMT_CHECK(km != 0 && km == b->km);
	p = kmalloc(km, 1000);
	MMT_CHECK(p != 0);
	memset(p, 7, 1000);
	kfree(km, p);
	mm_tbuf_destroy(b);
	mm_tbuf_destroy(0); /* minimap.h: destroying NULL is allowed by mm_tbuf_destroy's check */
	mm_dbg_flag |= MM_DBG_NO_KALLOC; /* --no-kalloc: the libc allocator */
	b = mm_tbuf_init();
	MMT_REQUIRE(b != 0);
	MMT_CHECK(mm_tbuf_get_km(b) == 0);
	mm_tbuf_destroy(b);
	mm_dbg_flag &= ~MM_DBG_NO_KALLOC;
	return 0;
}

/*****************
 * mm_map() hits *
 *****************/

static void print_tag(const char *name, int len, const char *buf)
{
	MMT_CHECKF(len == (int)strlen(buf), "%s: returned length %d, string length %d", name, len, (int)strlen(buf));
	printf("T\t%s\t%s\n", name, buf);
}

/* args: <target> <preset|-> <reads> <extra flags> <q_occ_frac|-> <no_kalloc 0/1>
 * Maps every read with mm_map() and prints, per hit, "H qname qlen qs qe strand tname tlen rs re
 * mlen blen mapq P|S inv sam_pri cigar rep_len"; for a read without hits "U qname qlen rep_len";
 * with an alignment also the cs (short and long), ds and MD strings from mm_gen_cs/ds/MD */
static int case_map(int argc, char *argv[])
{
	mm_idxopt_t io;
	mm_mapopt_t mo;
	mm_idx_t *mi;
	mm_tbuf_t *b;
	mm_bseq1_t *seqs;
	char *buf = 0;
	int max_len = 0, i, n;
	if (argc < 6) return MMT_FAIL("usage: map <target> <preset|-> <reads> <flags> <q_occ_frac|-> <no_kalloc>"), 1;
	mm_verbose = 1;
	MMT_REQUIRE(set_cli_opt(argv[1], &io, &mo) == 0);
	mo.flag |= MM_F_CIGAR | parse_flag(argv[3]);
	if (strcmp(argv[4], "-") != 0) mo.q_occ_frac = atof(argv[4]);
	if (atoi(argv[5])) mm_dbg_flag |= MM_DBG_NO_KALLOC;
	MMT_REQUIRE((mi = load_index(argv[0], &io)) != 0);
	mm_mapopt_update(&mo, mi);
	MMT_REQUIRE((seqs = read_all(argv[2], &n)) != 0);
	b = mm_tbuf_init();
	for (i = 0; i < n; ++i) {
		mm_bseq1_t *t = &seqs[i];
		mm_reg1_t *reg;
		int j, n_reg, n_pri = 0;
		reg = mm_map(mi, t->l_seq, t->seq, &n_reg, b, &mo, t->name);
		if (n_reg == 0) {
			printf("U\t%s\t%d\t%d\n", t->name, t->l_seq, b->rep_len);
			free(reg);
			continue;
		}
		for (j = 0; j < n_reg; ++j) {
			mm_reg1_t *r = &reg[j];
			uint32_t k;
			int len;
			void *km = j & 1? mm_tbuf_get_km(b) : 0;
			printf("H\t%s\t%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%c\t%d\t%d\t", t->name, t->l_seq, r->qs, r->qe, "+-"[r->rev],
				mi->seq[r->rid].name, mi->seq[r->rid].len, r->rs, r->re, r->mlen, r->blen, r->mapq, r->id == r->parent? 'P' : 'S', r->inv, r->sam_pri);
			MMT_REQUIRE(r->p != 0);
			for (k = 0; k < r->p->n_cigar; ++k)
				printf("%d%c", r->p->cigar[k]>>4, MM_CIGAR_STR[r->p->cigar[k]&0xf]);
			printf("\t%d\n", b->rep_len);
			n_pri += r->sam_pri;
			len = mm_gen_cs(km, &buf, &max_len, mi, r, t->seq, 1);
			print_tag("cs", len, buf);
			len = mm_gen_cs(km, &buf, &max_len, mi, r, t->seq, 0);
			print_tag("csl", len, buf);
			len = mm_gen_ds(km, &buf, &max_len, mi, r, t->seq, 1);
			print_tag("ds", len, buf);
			len = mm_gen_MD(km, &buf, &max_len, mi, r, t->seq);
			print_tag("MD", len, buf);
			MMT_CHECK(max_len > len);
		}
		MMT_CHECKF(n_pri == 1, "read %s: %d hits with sam_pri", t->name, n_pri);
		free_regs(n_reg, reg);
	}
	{ /* an empty query has no hit */
		mm_reg1_t *reg;
		int n_reg = -1;
		reg = mm_map(mi, 0, "", &n_reg, b, &mo, "empty");
		MMT_CHECK_INT(n_reg, ==, 0);
		free(reg);
	}
	free(buf);
	mm_tbuf_destroy(b);
	free_seqs(n, seqs);
	mm_idx_destroy(mi);
	return 0;
}

/* args: <target> <reads> <preset|-> <flags> <n_threads> <mini_batch|-> <rg|-> <seed|->
 * mm_map_file() as the CLI calls it, with the header from mm_write_sam_hdr() in SAM mode */
static int case_map_file(int argc, char *argv[])
{
	mm_idxopt_t io;
	mm_mapopt_t mo;
	mm_idx_t *mi;
	const char *rg;
	if (argc < 8) return MMT_FAIL("usage: map_file <target> <reads> <preset|-> <flags> <threads> <K|-> <rg|-> <seed|->"), 1;
	mm_verbose = 1;
	MMT_REQUIRE(set_cli_opt(argv[2], &io, &mo) == 0);
	mo.flag |= parse_flag(argv[3]);
	if (strcmp(argv[5], "-") != 0) mo.mini_batch_size = atoll(argv[5]);
	rg = strcmp(argv[6], "-") == 0? 0 : argv[6];
	if (strcmp(argv[7], "-") != 0) mo.seed = atoi(argv[7]);
	if (!(mo.flag & MM_F_CIGAR)) io.flag |= MM_I_NO_SEQ; /* as main.c */
	MMT_REQUIRE(mm_check_opt(&io, &mo) == 0);
	MMT_REQUIRE((mi = load_index(argv[0], &io)) != 0);
	mm_mapopt_update(&mo, mi);
	if (mo.flag & MM_F_OUT_SAM)
		MMT_REQUIRE(mm_write_sam_hdr(mi, rg, MM_VERSION, 0, 0) == 0);
	MMT_CHECK_INT(mm_map_file(mi, argv[1], &mo, atoi(argv[4])), ==, 0);
	fflush(stdout);
	MMT_CHECK_INT(mm_map_file(mi, "/nonexistent/no/such/reads.fq", &mo, atoi(argv[4])), ==, -1);
	mm_idx_destroy(mi);
	return 0;
}

/* args: <target> <reads> <flags> <q_occ_frac|->
 * The writers of format.c on mm_map() hits: mm_write_paf (P), mm_write_paf3 (Q), mm_write_paf4 (F,
 * one segment), mm_write_sam (S), mm_write_sam2 (T) and mm_write_sam3 (R) */
static int case_write(int argc, char *argv[])
{
	mm_idxopt_t io;
	mm_mapopt_t mo;
	mm_idx_t *mi;
	mm_tbuf_t *b;
	mm_bseq1_t *seqs;
	kstring_t s = {0,0,0};
	int64_t flag;
	int i, n;
	if (argc < 4) return MMT_FAIL("usage: write <target> <reads> <flags> <q_occ_frac|->"), 1;
	mm_verbose = 1;
	MMT_REQUIRE(set_cli_opt("-", &io, &mo) == 0);
	flag = parse_flag(argv[2]);
	mo.flag |= MM_F_CIGAR | flag;
	if (strcmp(argv[3], "-") != 0) mo.q_occ_frac = atof(argv[3]);
	MMT_REQUIRE((mi = load_index(argv[0], &io)) != 0);
	mm_mapopt_update(&mo, mi);
	MMT_REQUIRE((seqs = read_all(argv[1], &n)) != 0);
	b = mm_tbuf_init();
	for (i = 0; i < n; ++i) {
		mm_bseq1_t *t = &seqs[i];
		mm_reg1_t *reg;
		const mm_reg1_t *creg;
		int j, n_reg;
		void *km = km_init();
		reg = mm_map(mi, t->l_seq, t->seq, &n_reg, b, &mo, t->name);
		creg = reg;
		if (n_reg == 0) {
			mm_write_paf3(&s, mi, t, 0, km, mo.flag, b->rep_len); printf("Q\t%s\n", s.s);
			mm_write_sam(&s, mi, t, 0, 0, 0); printf("S\t%s\n", s.s);
			mm_write_sam3(&s, mi, t, 0, -1, 1, &n_reg, &creg, km, mo.flag, b->rep_len); printf("R\t%s\n", s.s);
		}
		for (j = 0; j < n_reg; ++j) {
			mm_reg1_t *r = &reg[j];
			mm_write_paf(&s, mi, t, r, km, mo.flag); printf("P\t%s\n", s.s);
			mm_write_paf3(&s, mi, t, r, km, mo.flag, b->rep_len); printf("Q\t%s\n", s.s);
			mm_write_paf4(&s, mi, t, r, km, mo.flag, b->rep_len, 1, 0); printf("F\t%s\n", s.s);
			mm_write_sam(&s, mi, t, r, n_reg, reg); printf("S\t%s\n", s.s);
			mm_write_sam2(&s, mi, t, 0, j, 1, &n_reg, &creg, km, mo.flag); printf("T\t%s\n", s.s);
			mm_write_sam3(&s, mi, t, 0, j, 1, &n_reg, &creg, km, mo.flag, b->rep_len); printf("R\t%s\n", s.s);
		}
		km_destroy(km);
		free_regs(n_reg, reg);
	}
	free(s.s);
	mm_tbuf_destroy(b);
	free_seqs(n, seqs);
	mm_idx_destroy(mi);
	return 0;
}

/* args: <target> <mode> [rg]. mm_write_sam_hdr() with or without index, version, argv and
 * read group; prints "ret N" after the header */
static int case_hdr(int argc, char *argv[])
{
	mm_idxopt_t io;
	mm_mapopt_t mo;
	mm_idx_t *mi = 0;
	char *args[] = { "prog", "-a", "x y", "z" };
	const char *mode, *rg;
	int ret;
	if (argc < 2) return MMT_FAIL("usage: hdr <target> <full|noidx|nover|noargs> [rg]"), 1;
	mode = argv[1];
	rg = argc > 2? argv[2] : 0;
	mm_verbose = 1;
	set_cli_opt("-", &io, &mo);
	if (strcmp(mode, "noidx") != 0) MMT_REQUIRE((mi = load_index(argv[0], &io)) != 0);
	if (strcmp(mode, "nover") == 0) ret = mm_write_sam_hdr(mi, rg, 0, 4, args);
	else if (strcmp(mode, "noargs") == 0) ret = mm_write_sam_hdr(mi, rg, "v1", 1, args);
	else ret = mm_write_sam_hdr(mi, rg, "v1", 4, args);
	printf("ret %d\n", ret);
	if (mi) mm_idx_destroy(mi);
	return 0;
}

/* args: <n_ops> <clip5> <clip3> <kind: pri|sup|sec> <flags>
 * One hand-made alignment of n_ops CIGAR operations (alternating 1=/1X) on a generated reference,
 * written with mm_write_sam3(); for -L (MM_F_LONG_CIGAR) when the CIGAR has more than 65535
 * operations with the clips */
static int case_long_cigar(int argc, char *argv[])
{
	int n_ops, clip5, clip3, qlen, rlen, i, n_reg = 1;
	int64_t flag;
	char *ref, *qry, *name = "lc";
	const char *seq[1], *names[1];
	mm_idx_t *mi;
	mm_reg1_t *r;
	const mm_reg1_t *creg;
	mm_bseq1_t t;
	kstring_t s = {0,0,0};
	uint64_t x = 11;
	if (argc < 5) return MMT_FAIL("usage: long_cigar <n_ops> <clip5> <clip3> <pri|sup|sec> <flags>"), 1;
	n_ops = atoi(argv[0]), clip5 = atoi(argv[1]), clip3 = atoi(argv[2]);
	flag = parse_flag(argv[4]);
	rlen = n_ops + 2000;
	ref = (char*)malloc(rlen + 1);
	for (i = 0; i < rlen; ++i) {
		x = x * 6364136223846793005ULL + 1442695040888963407ULL;
		ref[i] = "ACGT"[x >> 62];
	}
	ref[rlen] = 0;
	qlen = clip5 + n_ops + clip3;
	qry = (char*)malloc(qlen + 1);
	for (i = 0; i < clip5; ++i) qry[i] = "AC"[i & 1];
	for (i = 0; i < n_ops; ++i) { /* reference from 1000, a mismatch at every odd column */
		char c = ref[1000 + i];
		qry[clip5 + i] = (i & 1)? (c == 'A'? 'C' : 'A') : c;
	}
	for (i = 0; i < clip3; ++i) qry[clip5 + n_ops + i] = "GT"[i & 1];
	qry[qlen] = 0;
	seq[0] = ref, names[0] = "ref";
	mi = mm_idx_str(10, 15, 0, 14, 1, seq, names);
	MMT_REQUIRE(mi != 0);
	r = (mm_reg1_t*)calloc(1, sizeof(mm_reg1_t));
	r->id = 0, r->parent = 0, r->sam_pri = 1, r->rid = 0, r->rev = 0, r->mapq = 60;
	if (strcmp(argv[3], "sup") == 0) r->sam_pri = 0;
	else if (strcmp(argv[3], "sec") == 0) r->parent = 1, r->sam_pri = 0;
	r->qs = clip5, r->qe = clip5 + n_ops, r->rs = 1000, r->re = 1000 + n_ops;
	r->mlen = (n_ops + 1) / 2, r->blen = n_ops;
	r->p = (mm_extra_t*)calloc(sizeof(mm_extra_t) + n_ops * sizeof(uint32_t), 1);
	r->p->capacity = n_ops, r->p->n_cigar = n_ops;
	for (i = 0; i < n_ops; ++i)
		r->p->cigar[i] = 1u << 4 | ((i & 1)? MM_CIGAR_X_MISMATCH : MM_CIGAR_EQ_MATCH);
	r->p->dp_score = r->p->dp_max = r->p->dp_max0 = 1;
	memset(&t, 0, sizeof(t));
	t.name = name, t.seq = qry, t.l_seq = qlen;
	creg = r;
	mm_write_sam3(&s, mi, &t, 0, 0, 1, &n_reg, &creg, 0, flag, -1);
	puts(s.s);
	free(s.s);
	free(r->p); free(r);
	mm_idx_destroy(mi);
	free(ref); free(qry);
	return 0;
}

/**********
 * kalloc *
 **********/

static int fill_ok(const uint8_t *p, size_t n, uint8_t v)
{
	size_t i;
	for (i = 0; i < n; ++i)
		if (p[i] != v) return 0;
	return 1;
}

/* args: none. Allocation, contents and alignment under kmalloc/kcalloc/krealloc/krelocate/kfree,
 * child allocators, km_stat, and the NULL allocator (libc) */
static int case_kalloc(int argc, char *argv[])
{
	enum { N = 600 };
	void *km, *child;
	uint8_t *p[N], *q;
	size_t sz[N], live = 0;
	km_stat_t st;
	int i, bad_align = 0, bad_fill = 0;
	km = km_init();
	MMT_REQUIRE(km != 0);
	for (i = 0; i < N; ++i) {
		sz[i] = 1 + (size_t)i * 37 % 3001 + (i % 50 == 0? 200000 : 0);
		p[i] = (uint8_t*)kmalloc(km, sz[i]);
		if (p[i] == 0) { MMT_FAIL("kmalloc(%zu) returned NULL", sz[i]); return 1; }
		if ((uintptr_t)p[i] % sizeof(size_t) != 0) ++bad_align;
		memset(p[i], i & 0xff, sz[i]);
		live += sz[i];
	}
	MMT_CHECK_INT(bad_align, ==, 0);
	for (i = 0; i < N; ++i)
		if (!fill_ok(p[i], sz[i], i & 0xff)) ++bad_fill;
	MMT_CHECKF(bad_fill == 0, "%d blocks overwritten", bad_fill);
	km_stat(km, &st);
	MMT_CHECKF(st.capacity >= live, "capacity %zu below %zu live bytes", st.capacity, live);
	MMT_CHECK(st.available <= st.capacity && st.largest <= st.capacity);
	for (i = 0; i < N; i += 2) { /* free every other block, dirtying it first */
		memset(p[i], 0xab, sz[i]);
		kfree(km, p[i]);
		p[i] = 0;
	}
	for (i = 0; i < N; i += 2) { /* reuse: kcalloc must zero recycled memory */
		p[i] = (uint8_t*)kcalloc(km, sz[i], 1);
		if (!fill_ok(p[i], sz[i], 0)) ++bad_fill;
		memset(p[i], i & 0xff, sz[i]);
	}
	MMT_CHECKF(bad_fill == 0, "kcalloc returned non-zero bytes in %d blocks", bad_fill);
	for (i = 1; i < N; i += 2)
		if (!fill_ok(p[i], sz[i], i & 0xff)) ++bad_fill;
	MMT_CHECKF(bad_fill == 0, "%d surviving blocks overwritten", bad_fill);
	/* krealloc keeps the contents when growing and when shrinking */
	q = (uint8_t*)krealloc(km, p[3], sz[3] + 50000);
	MMT_CHECK(q != 0 && fill_ok(q, sz[3], 3));
	memset(q + sz[3], 0x5a, 50000);
	p[3] = q, sz[3] += 50000;
	q = (uint8_t*)krealloc(km, p[5], 3);
	MMT_CHECK(q != 0 && fill_ok(q, 3, 5));
	p[5] = q, sz[5] = 3;
	q = (uint8_t*)krealloc(km, 0, 100); /* NULL pointer: a new block */
	MMT_CHECK(q != 0);
	memset(q, 1, 100);
	MMT_CHECK(krealloc(km, q, 0) == 0); /* size 0: freed */
	/* krelocate moves a block, contents included */
	q = (uint8_t*)krelocate(km, p[7], sz[7]);
	MMT_CHECK(q != 0 && fill_ok(q, sz[7], 7));
	p[7] = q;
	MMT_CHECK(fill_ok(p[3], sz[3] - 50000, 3) && fill_ok(p[3] + sz[3] - 50000, 50000, 0x5a));
	for (i = 0; i < N; ++i)
		if (i != 5 && !fill_ok(p[i], sz[i], i & 0xff) && i != 3) ++bad_fill;
	MMT_CHECKF(bad_fill == 0, "%d blocks overwritten after realloc", bad_fill);
	/* a child allocator takes its memory from the parent */
	child = km_init2(km, 0);
	MMT_REQUIRE(child != 0);
	q = (uint8_t*)kmalloc(child, 70000);
	memset(q, 0x77, 70000);
	MMT_CHECK(fill_ok(p[1], sz[1], 1) && fill_ok(q, 70000, 0x77));
	kfree(child, q);
	km_destroy(child);
	for (i = 0; i < N; ++i) kfree(km, p[i]);
	km_stat(km, &st);
	MMT_CHECK(st.capacity > 0 && st.available <= st.capacity);
	km_stat_print(km);
	printf("capacity %zu\n", st.capacity);
	km_destroy(km);
	km_destroy(0);
	/* km == NULL: the libc allocator */
	q = (uint8_t*)kmalloc(0, 1000);
	MMT_CHECK(q != 0);
	memset(q, 3, 1000);
	q = (uint8_t*)krealloc(0, q, 200000);
	MMT_CHECK(q != 0 && fill_ok(q, 1000, 3));
	MMT_CHECK(krelocate(0, q, 10) == q);
	kfree(0, q);
	q = (uint8_t*)kcalloc(0, 100, 10);
	MMT_CHECK(q != 0 && fill_ok(q, 1000, 0));
	kfree(0, q);
	km_stat(0, &st);
	MMT_CHECK(st.capacity == 0 && st.available == 0 && st.n_blocks == 0 && st.n_cores == 0);
	return 0;
}

/***********
 * kthread *
 ***********/

typedef struct {
	int n_threads;
	int *count;
	int bad_tid;
} kf_data_t;

static void kf_func(void *data, long i, int tid)
{
	kf_data_t *d = (kf_data_t*)data;
	__sync_fetch_and_add(&d->count[i], 1);
	if (tid < 0 || tid >= d->n_threads) __sync_fetch_and_add(&d->bad_tid, 1);
}

typedef struct {
	int n_items, next, n_out, *out;
	int in0, in2, bad;
} kp_shared_t;

typedef struct {
	int id, stage;
} kp_item_t;

static void *kp_func(void *shared, int step, void *in)
{
	kp_shared_t *s = (kp_shared_t*)shared;
	kp_item_t *it = (kp_item_t*)in;
	if (step == 0) {
		if (__sync_add_and_fetch(&s->in0, 1) != 1) __sync_fetch_and_add(&s->bad, 1);
		if (s->next < s->n_items) {
			it = (kp_item_t*)calloc(1, sizeof(kp_item_t));
			it->id = s->next++, it->stage = 1;
		} else it = 0;
		__sync_fetch_and_sub(&s->in0, 1);
		return it;
	} else if (step == 1) {
		volatile long k, spin = (long)(s->n_items - it->id) * 3000; /* earlier items take longer */
		if (it->stage != 1) __sync_fetch_and_add(&s->bad, 1);
		for (k = 0; k < spin; ++k);
		it->stage = 2;
		return it;
	} else {
		if (__sync_add_and_fetch(&s->in2, 1) != 1) __sync_fetch_and_add(&s->bad, 1);
		if (it->stage != 2) __sync_fetch_and_add(&s->bad, 1);
		s->out[s->n_out++] = it->id;
		free(it);
		__sync_fetch_and_sub(&s->in2, 1);
		return 0;
	}
}

/* args: none. kt_for visits each index once with a thread id below n_threads; kt_pipeline runs
 * step 0 and the last step one item at a time, in input order, and each item through every step */
static int case_kthread(int argc, char *argv[])
{
	static const int threads[] = { 1, 2, 3, 4, 8 };
	static const long ns[] = { 0, 1, 7, 1000 };
	int a, b, i;
	for (a = 0; a < 5; ++a) {
		for (b = 0; b < 4; ++b) {
			kf_data_t d;
			int bad = 0;
			d.n_threads = threads[a], d.bad_tid = 0;
			d.count = (int*)calloc(ns[b] + 1, sizeof(int));
			kt_for(threads[a], kf_func, &d, ns[b]);
			for (i = 0; i < ns[b]; ++i)
				if (d.count[i] != 1) ++bad;
			MMT_CHECKF(bad == 0 && d.bad_tid == 0 && d.count[ns[b]] == 0, "kt_for(%d threads, n=%ld): %d indices not visited once, %d bad thread ids",
				threads[a], ns[b], bad, d.bad_tid);
			free(d.count);
		}
	}
	for (a = 0; a < 4; ++a) {
		static const int n_items[] = { 0, 1, 40 };
		for (b = 0; b < 3; ++b) {
			kp_shared_t s;
			int ordered = 1;
			memset(&s, 0, sizeof(s));
			s.n_items = n_items[b];
			s.out = (int*)calloc(n_items[b] + 1, sizeof(int));
			kt_pipeline(threads[a], kp_func, &s, 3);
			for (i = 0; i < s.n_out; ++i)
				if (s.out[i] != i) ordered = 0;
			MMT_CHECKF(s.n_out == n_items[b] && ordered && s.bad == 0, "kt_pipeline(%d threads, %d items): %d out, ordered %d, %d violations",
				threads[a], n_items[b], s.n_out, ordered, s.bad);
			free(s.out);
		}
	}
	return 0;
}

/**********
 * misc.c *
 **********/

/* args: <temporary file>. realtime(), cputime() and peakrss(); mm_err_puts/fwrite/fread */
static int case_misc(int argc, char *argv[])
{
	double t0, t1, c0, c1;
	long r0, r1;
	int i, a[1000], b[1000];
	char *big;
	FILE *fp;
	if (argc < 1) return MMT_FAIL("usage: misc <tmpfile>"), 1;
	t0 = realtime(), c0 = cputime();
	MMT_CHECKF(t0 > (double)time(0) - 2.0 && t0 < (double)time(0) + 2.0, "realtime() %f is not the time in seconds", t0);
	while ((c1 = cputime()) - c0 < 0.05 && realtime() - t0 < 30.0);
	t1 = realtime();
	MMT_CHECKF(c1 - c0 >= 0.05, "cputime() advanced %f s in a busy loop", c1 - c0);
	MMT_CHECKF(t1 - t0 >= (c1 - c0) - 0.02, "wall time %f s less than CPU time %f s", t1 - t0, c1 - c0);
	r0 = peakrss(); /* units and availability differ by platform (0 in some sandboxes): no value checks */
	MMT_CHECK(r0 >= 0);
	big = (char*)malloc(64 << 20);
	memset(big, 1, 64 << 20);
	r1 = peakrss();
	MMT_CHECK(r1 >= r0);
	free(big);
	mm_err_puts("mm_err_puts\tline");
	for (i = 0; i < 1000; ++i) a[i] = i * 7 - 3;
	MMT_REQUIRE((fp = fopen(argv[0], "wb")) != 0);
	mm_err_fwrite(a, sizeof(int), 1000, fp);
	fclose(fp);
	MMT_REQUIRE((fp = fopen(argv[0], "rb")) != 0);
	memset(b, 0, sizeof(b));
	mm_err_fread(b, sizeof(int), 1000, fp);
	fclose(fp);
	MMT_CHECK(memcmp(a, b, sizeof(a)) == 0);
	return 0;
}

static const mmt_case_t cases[] = {
	{ "macros", case_macros },
	{ "set_opt", case_set_opt },
	{ "check_opt", case_check_opt },
	{ "tbuf", case_tbuf },
	{ "map", case_map },
	{ "map_file", case_map_file },
	{ "write", case_write },
	{ "hdr", case_hdr },
	{ "long_cigar", case_long_cigar },
	{ "kalloc", case_kalloc },
	{ "kthread", case_kthread },
	{ "misc", case_misc },
	{ 0, 0 }
};

MMT_MAIN(cases)
