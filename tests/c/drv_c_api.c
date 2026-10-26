/* Area C driver: the splice and short-read parts of the libminimap2 API.
 *
 * preset_fields   option values the splice and short-read presets set,
 *                 against minimap2.1; MM_MAX_SEG; mm_max_spsc_bonus()
 * max_intron_len  mm_mapopt_max_intron_len() in and out of splice mode
 * bed_junc        mm_idx_bed_read() and mm_idx_bed_junc() on a BED file
 * spsc            mm_idx_spsc_read(), mm_idx_spsc_read2() and mm_idx_spsc_get()
 * jump_get        mm_idx_jjump_read() and mm_idx_jump_get() (mmpriv.h; reads
 *                 the off/off2/flag fields of mm_idx_jjump1_t, which a speed
 *                 change may rearrange but not rename)
 * map_splice      mm_map() in splice mode: is_spliced, trans_strand, introns
 * map_frag        mm_map_frag() on read pairs: proper_frag, seg_split, seg_id,
 *                 pe_thru; one segment equals mm_map(); MM_MAX_SEG + 1 segments
 * map_file_frag   mm_map_file_frag() writing SAM to stdout
 *
 * Expected values come from the arguments (computed in Python from the
 * simulated gene models and reads) or from minimap2.1. */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "minimap.h"
#include "mmpriv.h"
#include "ksw2.h"
#include "minimap2_test.h"

static mm_idx_t *load_index(const char *fn, const mm_idxopt_t *iopt)
{
	mm_idx_reader_t *r;
	mm_idx_t *mi;
	r = mm_idx_reader_open(fn, iopt, 0);
	if (r == 0) return 0;
	mi = mm_idx_reader_read(r, 1);
	mm_idx_reader_close(r);
	if (mi) mm_idx_index_name(mi);
	return mi;
}

static int has_flag(int64_t flag, int64_t x) { return (flag & x) == x; }

/* args: none */
static int case_preset_fields(int argc, char *argv[])
{
	mm_idxopt_t io, io2;
	mm_mapopt_t mo, mo2;

	MMT_CHECK_INT(MM_MAX_SEG, ==, 255);

	mm_set_opt(0, &io, &mo); // defaults (minimap2.1)
	MMT_CHECK_INT(mo.noncan, ==, 0); // -C [0]
	MMT_CHECK_INT(mo.jump_min_match, ==, 3); // --jump-min-match [3]
	MMT_CHECK_INT(mo.pe_ori, ==, 0);
	MMT_CHECK(!(mo.flag & (MM_F_SPLICE|MM_F_SPLICE_FOR|MM_F_SPLICE_REV|MM_F_SR|MM_F_FRAG_MODE|MM_F_SR_RNA|MM_F_WEAK_PAIRING|MM_F_INDEPEND_SEG|MM_F_SPLICE_OLD)));
	MMT_CHECK_INT(mm_max_spsc_bonus(&mo), ==, 20); // pins the current value: max((24+1)/2-1, 24-4)

	// sr: -k21 -w11 --sr --frag=yes -A2 -B8 -O12 -E2,1 -r100 -p.5 -N20 -f1000,5000 -n2 -m25 -s40 -g100 -2K50m --heap-sort=yes --secondary=no; -F [800]
	mm_set_opt(0, &io, &mo);
	MMT_REQUIRE(mm_set_opt("sr", &io, &mo) == 0);
	MMT_CHECK(io.k == 21 && io.w == 11 && !(io.flag & MM_I_HPC));
	MMT_CHECK(has_flag(mo.flag, MM_F_SR | MM_F_FRAG_MODE | MM_F_NO_PRINT_2ND | MM_F_2_IO_THREADS | MM_F_HEAP_SORT));
	MMT_CHECK(!(mo.flag & (MM_F_SPLICE|MM_F_SR_RNA|MM_F_WEAK_PAIRING|MM_F_INDEPEND_SEG)));
	MMT_CHECK(mo.a == 2 && mo.b == 8 && mo.q == 12 && mo.e == 2 && mo.e2 == 1);
	MMT_CHECK_INT(mo.bw, ==, 100);
	MMT_CHECK_INT(mo.max_frag_len, ==, 800);
	MMT_CHECK_INT(mo.max_gap, ==, 100);
	MMT_CHECK(mo.mid_occ == 1000 && mo.max_occ == 5000);
	MMT_CHECK(mo.min_cnt == 2 && mo.min_chain_score == 25 && mo.min_dp_max == 40 && mo.best_n == 20);
	MMT_CHECK(mo.pri_ratio > 0.49f && mo.pri_ratio < 0.51f);
	MMT_CHECK_INT(mo.mini_batch_size, ==, 50000000);
	MMT_CHECK_INT(mo.pe_ori, ==, 0<<1|1); // forward-reverse pairs
	mm_set_opt(0, &io2, &mo2);
	MMT_REQUIRE(mm_set_opt("short", &io2, &mo2) == 0); // undocumented alias
	MMT_CHECK(memcmp(&io, &io2, sizeof(io)) == 0 && memcmp(&mo, &mo2, sizeof(mo)) == 0);

	// splice: -k15 -w5 --splice -g2k -G200k -A1 -B2 -O2,32 -E1,0 -C9 -z200 -ub --junc-bonus=9 --cap-sw-mem=0 --splice-flank=yes; --spsc0 [5]
	mm_set_opt(0, &io, &mo);
	MMT_REQUIRE(mm_set_opt("splice", &io, &mo) == 0);
	MMT_CHECK(io.k == 15 && io.w == 5);
	MMT_CHECK(has_flag(mo.flag, MM_F_SPLICE | MM_F_SPLICE_FOR | MM_F_SPLICE_REV | MM_F_SPLICE_FLANK));
	MMT_CHECK(!(mo.flag & (MM_F_SR|MM_F_SR_RNA|MM_F_FRAG_MODE|MM_F_WEAK_PAIRING|MM_F_SPLICE_OLD)));
	MMT_CHECK(mo.a == 1 && mo.b == 2 && mo.q == 2 && mo.q2 == 32 && mo.e == 1 && mo.e2 == 0);
	MMT_CHECK_INT(mo.max_gap, ==, 2000);
	MMT_CHECK(mo.max_gap_ref == 200000 && mo.bw == 200000 && mo.bw_long == 200000);
	MMT_CHECK_INT(mo.noncan, ==, 9);
	MMT_CHECK_INT(mo.junc_bonus, ==, 9);
	MMT_CHECK_INT(mo.junc_pen, ==, 5);
	MMT_CHECK_INT(mo.zdrop, ==, 200);
	MMT_CHECK_INT(mo.max_sw_mat, ==, 0);
	MMT_CHECK_INT(mo.pe_ori, ==, 0);
	MMT_CHECK_INT(mm_max_spsc_bonus(&mo), ==, 30); // pins the current value: max((32+1)/2-1, 32-2)
	mm_set_opt(0, &io2, &mo2);
	MMT_REQUIRE(mm_set_opt("cdna", &io2, &mo2) == 0); // undocumented alias
	MMT_CHECK(memcmp(&io, &io2, sizeof(io)) == 0 && memcmp(&mo, &mo2, sizeof(mo)) == 0);

	// splice:hq: -xsplice -C5 -O6,24 -B4
	mm_set_opt(0, &io2, &mo2);
	MMT_REQUIRE(mm_set_opt("splice:hq", &io2, &mo2) == 0);
	MMT_CHECK(mo2.noncan == 5 && mo2.q == 6 && mo2.q2 == 24 && mo2.b == 4);
	MMT_CHECK_INT(mm_max_spsc_bonus(&mo2), ==, 18); // pins the current value: max((24+1)/2-1, 24-6)
	mo2.noncan = mo.noncan, mo2.q = mo.q, mo2.q2 = mo.q2, mo2.b = mo.b; // the rest equals splice
	MMT_CHECK(memcmp(&io, &io2, sizeof(io)) == 0 && memcmp(&mo, &mo2, sizeof(mo)) == 0);

	// splice:sr: -xsplice:hq --frag=yes -m25 -s40 -2K100m --heap-sort=yes --pairing=weak --sr=rna --min-dp-len=20 --secondary=no
	mm_set_opt(0, &io, &mo);
	MMT_REQUIRE(mm_set_opt("splice:sr", &io, &mo) == 0);
	MMT_CHECK(has_flag(mo.flag, MM_F_SPLICE | MM_F_SPLICE_FOR | MM_F_SPLICE_REV | MM_F_SPLICE_FLANK | MM_F_FRAG_MODE | MM_F_WEAK_PAIRING
		| MM_F_SR_RNA | MM_F_NO_PRINT_2ND | MM_F_2_IO_THREADS | MM_F_HEAP_SORT));
	MMT_CHECK(!(mo.flag & (MM_F_SR | MM_F_INDEPEND_SEG | MM_F_SPLICE_OLD)));
	MMT_CHECK(mo.noncan == 5 && mo.q == 6 && mo.q2 == 24 && mo.b == 4);
	MMT_CHECK(mo.min_chain_score == 25 && mo.min_dp_max == 40 && mo.min_ksw_len == 20);
	MMT_CHECK_INT(mo.mini_batch_size, ==, 100000000);
	MMT_CHECK_INT(mo.pe_ori, ==, 0<<1|1);
	MMT_CHECK_INT(mo.max_gap_ref, ==, 200000);

	MMT_CHECK_INT(mm_set_opt("splice:nosuch", &io, &mo), <, 0); // NEWS 2.29: wrong splice* presets are errors
	return 0;
}

/* args: none */
static int case_max_intron_len(int argc, char *argv[])
{
	mm_idxopt_t io;
	mm_mapopt_t mo;
	int bw, bw_long, gap;

	mm_set_opt(0, &io, &mo); // not in splice mode: -G has no effect
	bw = mo.bw, bw_long = mo.bw_long, gap = mo.max_gap_ref;
	mm_mapopt_max_intron_len(&mo, 5000);
	MMT_CHECK(mo.bw == bw && mo.bw_long == bw_long && mo.max_gap_ref == gap);

	mo.flag |= MM_F_SPLICE; // --splice
	mm_mapopt_max_intron_len(&mo, 5000);
	MMT_CHECK(mo.max_gap_ref == 5000 && mo.bw == 5000 && mo.bw_long == 5000);

	mm_set_opt("splice", &io, &mo);
	mm_mapopt_max_intron_len(&mo, 12345);
	MMT_CHECK(mo.max_gap_ref == 12345 && mo.bw == 12345 && mo.bw_long == 12345);
	mm_mapopt_max_intron_len(&mo, 0); // ignored
	MMT_CHECK(mo.max_gap_ref == 12345 && mo.bw == 12345 && mo.bw_long == 12345);
	mm_mapopt_max_intron_len(&mo, -5);
	MMT_CHECK(mo.max_gap_ref == 12345 && mo.bw == 12345 && mo.bw_long == 12345);
	MMT_CHECK_INT(mm_check_opt(&io, &mo), ==, 0);
	return 0;
}

/* parse "pos:val,pos:val,..." into a dense array of length n */
static int parse_expect(const char *s, int32_t off, int32_t n, uint8_t *out, uint8_t fill)
{
	const char *p = s;
	memset(out, fill, n);
	while (*p) {
		char *q;
		long pos, val;
		pos = strtol(p, &q, 10);
		if (*q != ':') return -1;
		val = strtol(q + 1, &q, 10);
		if (pos - off < 0 || pos - off >= n) return -1;
		out[pos - off] = (uint8_t)val;
		p = *q == ','? q + 1 : q;
	}
	return 0;
}

/* args: <target.fa> <bed> <read_junc> <ctg> <st> <en> <expected "pos:bits,...">
 * the expected list gives every non-zero position of the window */
static int case_bed_junc(int argc, char *argv[])
{
	mm_idxopt_t io;
	mm_mapopt_t mo;
	mm_idx_t *mi;
	int32_t st, en, i, cid, ret;
	uint8_t *s, *want;
	if (argc < 7) return MMT_FAIL("usage: bed_junc <fa> <bed> <read_junc> <ctg> <st> <en> <expected>"), 1;
	mm_set_opt(0, &io, &mo);
	MMT_REQUIRE((mi = load_index(argv[0], &io)) != 0);
	st = atoi(argv[4]), en = atoi(argv[5]);
	MMT_REQUIRE(st < en);
	cid = mm_idx_name2id(mi, argv[3]);
	MMT_REQUIRE(cid >= 0);
	s = (uint8_t*)malloc(en - st);
	want = (uint8_t*)malloc(en - st);
	memset(s, 0xaa, en - st);
	MMT_CHECK(mi->I == 0);
	ret = mm_idx_bed_junc(mi, cid, st, en, s); // no annotation yet
	MMT_CHECK_INT(ret, ==, -1);
	for (i = 0; i < en - st && s[i] == 0; ++i) {}
	MMT_CHECK_INT(i, ==, en - st);
	MMT_CHECK_INT(mm_idx_bed_read(mi, argv[1], atoi(argv[2])), ==, 0);
	MMT_REQUIRE(mi->I != 0);
	memset(s, 0xaa, en - st);
	ret = mm_idx_bed_junc(mi, cid, st, en, s);
	MMT_CHECK_INT(ret, >=, 0);
	MMT_REQUIRE(parse_expect(argv[6], st, en - st, want, 0) == 0);
	for (i = 0; i < en - st; ++i)
		MMT_CHECKF(s[i] == want[i], "%s:%d: bits %d, expected %d", argv[3], st + i, s[i], want[i]);
	MMT_CHECK_INT(mm_idx_bed_junc(mi, mi->n_seq, st, en, s), ==, -1); // contig out of range
	free(s); free(want);
	mm_idx_destroy(mi);
	return 0;
}

/* args: <target.fa> <spsc file> <max_sc> <scale|"read"> [ctg:strand:pos:type:score ...]
 * "read" calls mm_idx_spsc_read() (no scaling). Every contig and strand is
 * read back over its whole length: a listed position holds its score and type
 * (D 0, A 1) as align.c passes them to ksw_exts2_sse(), (score +
 * KSW_SPSC_OFFSET) << 1 | type; every other position holds 0xff. */
static int case_spsc(int argc, char *argv[])
{
	mm_idxopt_t io;
	mm_mapopt_t mo;
	mm_idx_t *mi;
	int32_t i, max_sc, n_expect = 0;
	uint32_t c;
	uint8_t *sc;
	struct { int32_t cid, rev, pos, byte; } *e;
	if (argc < 4) return MMT_FAIL("usage: spsc <fa> <file> <max_sc> <scale|read> [ctg:strand:pos:type:score ...]"), 1;
	mm_set_opt(0, &io, &mo);
	MMT_REQUIRE((mi = load_index(argv[0], &io)) != 0);
	max_sc = atoi(argv[2]);
	sc = (uint8_t*)malloc(16);
	MMT_CHECK_INT(mm_idx_spsc_get(mi, 0, 0, 10, 0, sc), ==, -1); // nothing loaded
	free(sc);
	if (strcmp(argv[3], "read") == 0) MMT_CHECK_INT(mm_idx_spsc_read(mi, argv[1], max_sc), ==, 0);
	else MMT_CHECK_INT(mm_idx_spsc_read2(mi, argv[1], max_sc, atof(argv[3])), ==, 0);
	MMT_REQUIRE(mi->spsc != 0);
	e = calloc(argc, sizeof(*e));
	for (i = 4; i < argc; ++i) {
		char name[256];
		char strand, type;
		int pos, score;
		MMT_REQUIRE(sscanf(argv[i], "%255[^:]:%c:%d:%c:%d", name, &strand, &pos, &type, &score) == 5);
		e[n_expect].cid = mm_idx_name2id(mi, name);
		MMT_REQUIRE(e[n_expect].cid >= 0);
		e[n_expect].rev = strand == '-';
		e[n_expect].pos = pos, e[n_expect].byte = (score + KSW_SPSC_OFFSET) << 1 | (type == 'A');
		++n_expect;
	}
	for (c = 0; c < mi->n_seq; ++c) {
		int32_t rev, len = mi->seq[c].len;
		sc = (uint8_t*)malloc(len + 1);
		for (rev = 0; rev < 2; ++rev) {
			uint8_t *want = (uint8_t*)malloc(len);
			int64_t ret;
			memset(want, 0xff, len);
			for (i = 0; i < n_expect; ++i)
				if (e[i].cid == (int32_t)c && e[i].rev == rev)
					want[e[i].pos] = (uint8_t)e[i].byte;
			ret = mm_idx_spsc_get(mi, c, 0, len, rev, sc);
			MMT_CHECK_INT(ret, ==, len);
			for (i = 0; i < len; ++i)
				MMT_CHECKF(sc[i] == want[i], "%s:%c:%d: 0x%02x, expected 0x%02x", mi->seq[c].name, "+-"[rev], i, sc[i], want[i]);
			ret = mm_idx_spsc_get(mi, c, 0, len + 100, rev, sc); // the end is clamped to the contig
			MMT_CHECK_INT(ret, ==, len);
			free(want);
		}
		free(sc);
	}
	sc = (uint8_t*)malloc(16);
	MMT_CHECK_INT(mm_idx_spsc_get(mi, mi->n_seq, 0, 10, 0, sc), ==, -1);
	MMT_CHECK_INT(mm_idx_spsc_get(mi, -1, 0, 10, 0, sc), ==, -1);
	free(sc); free(e);
	mm_idx_destroy(mi);
	return 0;
}

/* args: <target.fa> <anno.bed|-> <pass1.bed|-> <ctg> <st> <en> <expected "off:off2:flag,...">
 * -j reads anno.bed with MM_JUNC_ANNO, --pass1 reads pass1.bed with MM_JUNC_MISC
 * and a minimum score of 5, as main.c does. */
static int case_jump_get(int argc, char *argv[])
{
	mm_idxopt_t io;
	mm_mapopt_t mo;
	mm_idx_t *mi;
	const mm_idx_jjump1_t *a;
	int32_t i, j, n, n_want = 0, cid, *got;
	struct { int32_t off, off2, flag; } w[256];
	const char *p;
	if (argc < 7) return MMT_FAIL("usage: jump_get <fa> <anno|-> <pass1|-> <ctg> <st> <en> <expected>"), 1;
	mm_set_opt(0, &io, &mo);
	MMT_REQUIRE((mi = load_index(argv[0], &io)) != 0);
	cid = mm_idx_name2id(mi, argv[3]);
	MMT_REQUIRE(cid >= 0);
	a = mm_idx_jump_get(mi, cid, 0, mi->seq[cid].len, &n);
	MMT_CHECK(n == 0 && mi->J == 0);
	if (strcmp(argv[1], "-") != 0) mm_idx_jjump_read(mi, argv[1], MM_JUNC_ANNO, -1);
	if (strcmp(argv[2], "-") != 0) mm_idx_jjump_read(mi, argv[2], MM_JUNC_MISC, 5);
	MMT_REQUIRE(mi->J != 0);
	for (p = argv[6]; *p && n_want < 256; ) {
		char *q;
		w[n_want].off = strtol(p, &q, 10);
		w[n_want].off2 = strtol(q + 1, &q, 10);
		w[n_want].flag = strtol(q + 1, &q, 10);
		++n_want;
		p = *q == ','? q + 1 : q;
	}
	a = mm_idx_jump_get(mi, cid, atoi(argv[4]), atoi(argv[5]), &n);
	MMT_CHECK_INT(n, ==, n_want);
	got = (int32_t*)calloc(n_want + 1, sizeof(int32_t));
	for (i = 0; i < n; ++i) {
		int found = 0;
		MMT_CHECKF(i == 0 || a[i-1].off <= a[i].off, "entries are not sorted by off");
		for (j = 0; j < n_want; ++j)
			if (!got[j] && w[j].off == a[i].off && w[j].off2 == a[i].off2 && w[j].flag == a[i].flag)
				got[j] = 1, found = 1;
		MMT_CHECKF(found, "unexpected jump %d -> %d flag %d", a[i].off, a[i].off2, a[i].flag);
	}
	for (j = 0; j < n_want; ++j)
		MMT_CHECKF(got[j], "missing jump %d -> %d flag %d", w[j].off, w[j].off2, w[j].flag);
	free(got);
	mm_idx_destroy(mi);
	return 0;
}

static int n_cigar_op(const mm_extra_t *p, int op)
{
	uint32_t k;
	int n = 0;
	for (k = 0; k < p->n_cigar; ++k)
		if ((int)(p->cigar[k] & 0xf) == op) ++n;
	return n;
}

/* args: <target.fa> <preset> <reads.fa> <truth.tsv> <min_exact>
 * truth.tsv: name, contig, genome strand (+/-), transcript strand of the read (+/-), introns "s-e,s-e" or "." */
static int case_map_splice(int argc, char *argv[])
{
	mm_idxopt_t io;
	mm_mapopt_t mo;
	mm_idx_t *mi;
	mm_tbuf_t *tbuf;
	mmt_seqs_t s;
	FILE *fp;
	int32_t i, n_exact = 0, n_spliced = 0;
	if (argc < 5) return MMT_FAIL("usage: map_splice <fa> <preset> <reads.fa> <truth.tsv> <min_exact>"), 1;
	mm_verbose = 2;
	mm_set_opt(0, &io, &mo);
	MMT_REQUIRE(mm_set_opt(argv[1], &io, &mo) == 0);
	mo.flag |= MM_F_CIGAR;
	MMT_REQUIRE((mi = load_index(argv[0], &io)) != 0);
	mm_mapopt_update(&mo, mi);
	MMT_REQUIRE(mmt_read_seqs(argv[2], &s) > 0);
	MMT_REQUIRE((fp = fopen(argv[3], "r")) != 0);
	tbuf = mm_tbuf_init();
	for (i = 0; i < s.n; ++i) {
		char name[256], ctg[256], gstrand, tstrand, introns[4096], want[4096];
		mm_reg1_t *reg, *h = 0;
		int32_t j, n_reg;
		uint32_t k, t;
		int l = 0;
		MMT_REQUIRE(fscanf(fp, "%255s %255s %c %c %4095s", name, ctg, &gstrand, &tstrand, introns) == 5);
		MMT_CHECK_STR(name, s.name[i]);
		reg = mm_map(mi, s.len[i], s.seq[i], &n_reg, tbuf, &mo, s.name[i]);
		for (j = 0; j < n_reg; ++j)
			if (reg[j].id == reg[j].parent && reg[j].sam_pri) h = &reg[j];
		MMT_CHECKF(h != 0, "read %s: no primary hit", name);
		if (h) {
			MMT_CHECK_STR(mi->seq[h->rid].name, ctg);
			MMT_CHECKF("+-"[h->rev] == gstrand, "read %s: strand %c, truth %c", name, "+-"[h->rev], gstrand);
			MMT_REQUIRE(h->p != 0);
			MMT_CHECKF(h->is_spliced == (n_cigar_op(h->p, MM_CIGAR_N_SKIP) > 0), "read %s: is_spliced %d", name, h->is_spliced);
			want[0] = 0;
			for (k = 0, t = h->rs; k < h->p->n_cigar; ++k) {
				int op = h->p->cigar[k] & 0xf, len = h->p->cigar[k] >> 4;
				if (op == MM_CIGAR_N_SKIP) l += snprintf(want + l, sizeof(want) - l, "%s%u-%u", l? "," : "", t, t + len);
				if (op == MM_CIGAR_MATCH || op == MM_CIGAR_DEL || op == MM_CIGAR_N_SKIP || op == MM_CIGAR_EQ_MATCH || op == MM_CIGAR_X_MISMATCH) t += len;
			}
			if (l == 0) strcpy(want, ".");
			if (strcmp(want, introns) == 0) {
				++n_exact;
				if (h->is_spliced) { // ts:A of cookbook.md: 1 when the read is on the transcript strand, 2 when opposite
					++n_spliced;
					MMT_CHECKF(h->p->trans_strand == (tstrand == '+'? 1 : 2), "read %s: trans_strand %d, truth %c", name, h->p->trans_strand, tstrand);
				}
			}
		}
		for (j = 0; j < n_reg; ++j) free(reg[j].p);
		free(reg);
	}
	MMT_CHECKF(n_exact >= atoi(argv[4]), "%d reads with the true introns, expected at least %s", n_exact, argv[4]);
	MMT_CHECK(n_spliced > 0);
	fclose(fp);
	mm_tbuf_destroy(tbuf);
	mmt_free_seqs(&s);
	mm_idx_destroy(mi);
	return 0;
}

static void revcomp_str(char *s, int l)
{
	int i;
	for (i = 0; i < l>>1; ++i) {
		char t = s[i];
		s[i] = s[l - 1 - i], s[l - 1 - i] = t;
	}
	for (i = 0; i < l; ++i) {
		char c = s[i];
		s[i] = c == 'A'? 'T' : c == 'C'? 'G' : c == 'G'? 'C' : c == 'T'? 'A' : 'N';
	}
}

static int same_hit(const mm_reg1_t *a, const mm_reg1_t *b)
{
	if (a->rid != b->rid || a->rs != b->rs || a->re != b->re || a->qs != b->qs || a->qe != b->qe || a->rev != b->rev || a->mapq != b->mapq) return 0;
	if (!a->p || !b->p || a->p->n_cigar != b->p->n_cigar || a->p->dp_score != b->p->dp_score) return 0;
	return memcmp(a->p->cigar, b->p->cigar, a->p->n_cigar * 4) == 0;
}

/* args: <target.fa> <r1.fa> <r2.fa> <truth.tsv> <max_dist>
 * truth.tsv: name, contig, r1 start, r1 strand, r2 start, r2 strand, fragment length (0-based starts of the
 * reference span of each mate). As worker_for() in map.c does for forward-reverse pairs (pe_ori), mate 2 is
 * reverse-complemented before mm_map_frag() and its hits are flipped back after. */
static int case_map_frag(int argc, char *argv[])
{
	mm_idxopt_t io;
	mm_mapopt_t mo;
	mm_idx_t *mi;
	mm_tbuf_t *tbuf;
	mmt_seqs_t s1, s2;
	FILE *fp;
	int32_t i, max_dist, n_thru = 0;
	if (argc < 5) return MMT_FAIL("usage: map_frag <fa> <r1.fa> <r2.fa> <truth.tsv> <max_dist>"), 1;
	mm_verbose = 2;
	max_dist = atoi(argv[4]);
	mm_set_opt(0, &io, &mo);
	MMT_REQUIRE(mm_set_opt("sr", &io, &mo) == 0);
	mo.flag |= MM_F_CIGAR;
	MMT_REQUIRE((mi = load_index(argv[0], &io)) != 0);
	mm_mapopt_update(&mo, mi);
	MMT_REQUIRE(mmt_read_seqs(argv[1], &s1) > 0);
	MMT_REQUIRE(mmt_read_seqs(argv[2], &s2) == s1.n);
	MMT_REQUIRE((fp = fopen(argv[3], "r")) != 0);
	tbuf = mm_tbuf_init();
	for (i = 0; i < s1.n; ++i) {
		char name[256], ctg[256], st1, st2;
		long pos[2], frag;
		int qlens[2], n_regs[2], j, k, n_pri[2] = {0, 0};
		const char *seqs[2];
		mm_reg1_t *regs[2], *pri[2] = {0, 0}, *single;
		int n_single;
		MMT_REQUIRE(fscanf(fp, "%255s %255s %ld %c %ld %c %ld", name, ctg, &pos[0], &st1, &pos[1], &st2, &frag) == 7);
		MMT_CHECK_STR(name, s1.name[i]);
		qlens[0] = s1.len[i], qlens[1] = s2.len[i];
		revcomp_str(s2.seq[i], s2.len[i]);
		seqs[0] = s1.seq[i], seqs[1] = s2.seq[i];
		mm_map_frag(mi, 2, qlens, seqs, n_regs, regs, tbuf, &mo, name);
		revcomp_str(s2.seq[i], s2.len[i]);
		for (k = 0; k < n_regs[1]; ++k) { // flip mate 2 back to its own strand
			mm_reg1_t *r = &regs[1][k];
			int t = r->qs;
			r->qs = qlens[1] - r->qe, r->qe = qlens[1] - t, r->rev = !r->rev;
		}
		for (j = 0; j < 2; ++j) {
			for (k = 0; k < n_regs[j]; ++k) {
				mm_reg1_t *r = &regs[j][k];
				MMT_CHECKF(r->seg_split == 1 && r->seg_id == j, "pair %s mate %d hit %d: seg_split %d seg_id %d", name, j + 1, k, r->seg_split, r->seg_id);
				if (r->id == r->parent && r->sam_pri) pri[j] = r, ++n_pri[j];
			}
			MMT_CHECKF(n_pri[j] == 1, "pair %s mate %d: %d primary hits", name, j + 1, n_pri[j]);
		}
		if (pri[0] && pri[1]) {
			long start5[2];
			MMT_CHECK_STR(mi->seq[pri[0]->rid].name, ctg);
			MMT_CHECKF(pri[0]->rid == pri[1]->rid && "+-"[pri[0]->rev] == st1 && "+-"[pri[1]->rev] == st2, "pair %s: strands %c%c", name, "+-"[pri[0]->rev], "+-"[pri[1]->rev]);
			for (j = 0; j < 2; ++j) {
				start5[j] = (long)pri[j]->rs - (pri[j]->rev? qlens[j] - pri[j]->qe : pri[j]->qs);
				MMT_CHECKF(labs(start5[j] - pos[j]) <= max_dist, "pair %s mate %d: start %ld, truth %ld", name, j + 1, start5[j], pos[j]);
			}
			MMT_CHECKF(pri[0]->proper_frag && pri[1]->proper_frag, "pair %s: proper_frag %d %d", name, pri[0]->proper_frag, pri[1]->proper_frag);
			// pe_thru: both mates cover the same reference interval (read-through of a fragment no longer than the reads)
			MMT_CHECKF(pri[0]->pe_thru == (frag == qlens[0] && qlens[0] == qlens[1]) && pri[1]->pe_thru == pri[0]->pe_thru,
				"pair %s (fragment %ld): pe_thru %d %d", name, frag, pri[0]->pe_thru, pri[1]->pe_thru);
			n_thru += pri[0]->pe_thru;
		}
		// mate 1 alone through mm_map_frag(n_segs = 1) and mm_map(): the same hits
		single = mm_map(mi, qlens[0], seqs[0], &n_single, tbuf, &mo, name);
		{
			int n1;
			mm_reg1_t *r1;
			mm_map_frag(mi, 1, &qlens[0], &seqs[0], &n1, &r1, tbuf, &mo, name);
			MMT_CHECKF(n1 == n_single, "pair %s: mm_map %d hits, mm_map_frag(1) %d", name, n_single, n1);
			for (k = 0; k < n1 && k < n_single; ++k) {
				MMT_CHECKF(same_hit(&r1[k], &single[k]), "pair %s: hit %d differs between mm_map and mm_map_frag", name, k);
				MMT_CHECK(r1[k].seg_split == 0 && r1[k].pe_thru == 0 && r1[k].proper_frag == 0);
			}
			for (k = 0; k < n1; ++k) free(r1[k].p);
			free(r1);
		}
		for (k = 0; k < n_single; ++k) free(single[k].p);
		free(single);
		for (j = 0; j < 2; ++j) {
			for (k = 0; k < n_regs[j]; ++k) free(regs[j][k].p);
			free(regs[j]);
		}
	}
	MMT_CHECK(n_thru > 0);
	{ // more segments than MM_MAX_SEG: no hits
		int qlens[MM_MAX_SEG + 1], n_regs[MM_MAX_SEG + 1];
		const char *seqs[MM_MAX_SEG + 1];
		mm_reg1_t *regs[MM_MAX_SEG + 1];
		for (i = 0; i <= MM_MAX_SEG; ++i)
			qlens[i] = s1.len[0], seqs[i] = s1.seq[0], n_regs[i] = -1, regs[i] = (mm_reg1_t*)1;
		mm_map_frag(mi, MM_MAX_SEG + 1, qlens, seqs, n_regs, regs, tbuf, &mo, "too_many");
		for (i = 0; i <= MM_MAX_SEG; ++i)
			MMT_CHECKF(n_regs[i] == 0 && regs[i] == 0, "segment %d: %d hits", i, n_regs[i]);
	}
	fclose(fp);
	mm_tbuf_destroy(tbuf);
	mmt_free_seqs(&s1);
	mmt_free_seqs(&s2);
	mm_idx_destroy(mi);
	return 0;
}

/* args: <target.fa> <preset> <query files...>; SAM records on stdout, as `minimap2 -a -x preset` writes them.
 * A query file that does not exist makes mm_map_file_frag() return -1 and write nothing. */
static int case_map_file_frag(int argc, char *argv[])
{
	mm_idxopt_t io;
	mm_mapopt_t mo;
	mm_idx_t *mi;
	const char *missing[2];
	if (argc < 3) return MMT_FAIL("usage: map_file_frag <fa> <preset> <files...>"), 1;
	mm_verbose = 1;
	mm_set_opt(0, &io, &mo);
	MMT_REQUIRE(mm_set_opt(argv[1], &io, &mo) == 0);
	mo.flag |= MM_F_OUT_SAM | MM_F_CIGAR;
	MMT_REQUIRE((mi = load_index(argv[0], &io)) != 0);
	mm_mapopt_update(&mo, mi);
	missing[0] = argv[2], missing[1] = "no-such-file.fq";
	fflush(stdout);
	MMT_CHECK_INT(mm_map_file_frag(mi, 2, missing, &mo, 2), ==, -1);
	MMT_CHECK_INT(mm_map_file_frag(mi, 0, missing, &mo, 2), ==, -1);
	fflush(stdout);
	MMT_CHECK_INT(mm_map_file_frag(mi, argc - 2, (const char**)&argv[2], &mo, 2), ==, 0);
	fflush(stdout);
	mm_idx_destroy(mi);
	return 0;
}

static const mmt_case_t cases[] = {
	{ "preset_fields", case_preset_fields },
	{ "max_intron_len", case_max_intron_len },
	{ "bed_junc", case_bed_junc },
	{ "spsc", case_spsc },
	{ "jump_get", case_jump_get },
	{ "map_splice", case_map_splice },
	{ "map_frag", case_map_frag },
	{ "map_file_frag", case_map_file_frag },
	{ 0, 0 }
};

MMT_MAIN(cases)
