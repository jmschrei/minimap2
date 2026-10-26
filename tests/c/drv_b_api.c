/* Area B driver: mapping through the libminimap2 API (minimap.h) with the
 * mm_mapopt_t fields and MM_F_* flags set directly, printing every mm_reg1_t
 * and mm_extra_t field so that Python can compare them with the command
 * line's PAF for the corresponding options. Uses only minimap.h.
 *
 * map <target.fa> <reads.fa> <preset|-> [setting...]
 *   setting: field=value (an mm_mapopt_t field, or k= / w= for the index),
 *            +FLAG / -FLAG (set / clear MM_F_FLAG)
 *   output, one line per hit (TAB-separated):
 *     qname qlen qs qe strand tname tlen rs re mlen blen mapq id parent
 *     sam_pri inv split split_inv is_alt cnt score score0 subsc n_sub div
 *     rep_len has_p dp_score dp_max dp_max0 dp_max2 n_ambi trans_strand
 *     n_cigar capacity cigar
 * opt_flags <preset> FLAG... : the flags mm_set_opt(preset) sets (checks that
 *   exactly the named ones among the mapping-mode flags are set) */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "minimap.h"
#include "minimap2_test.h"

static const struct { const char *name; int64_t bit; } flags[] = {
	{ "NO_DIAG", MM_F_NO_DIAG }, { "NO_DUAL", MM_F_NO_DUAL }, { "CIGAR", MM_F_CIGAR }, { "NO_LJOIN", MM_F_NO_LJOIN },
	{ "SR", MM_F_SR }, { "FRAG_MODE", MM_F_FRAG_MODE }, { "NO_PRINT_2ND", MM_F_NO_PRINT_2ND }, { "SPLICE", MM_F_SPLICE },
	{ "FOR_ONLY", MM_F_FOR_ONLY }, { "REV_ONLY", MM_F_REV_ONLY }, { "HEAP_SORT", MM_F_HEAP_SORT }, { "ALL_CHAINS", MM_F_ALL_CHAINS },
	{ "NO_END_FLT", MM_F_NO_END_FLT }, { "HARD_MLEVEL", MM_F_HARD_MLEVEL }, { "RMQ", MM_F_RMQ }, { "QSTRAND", MM_F_QSTRAND },
	{ "NO_INV", MM_F_NO_INV }, { "NO_HASH_NAME", MM_F_NO_HASH_NAME }, { "EQX", MM_F_EQX }, { 0, 0 }
};

static int64_t flag_bit(const char *name)
{
	int i;
	for (i = 0; flags[i].name; ++i)
		if (strcmp(flags[i].name, name) == 0) return flags[i].bit;
	return 0;
}

#define INT_FIELDS(X) X(seed) X(max_qlen) X(bw) X(bw_long) X(max_gap) X(max_gap_ref) X(max_chain_skip) X(max_chain_iter) X(min_cnt) \
	X(min_chain_score) X(rmq_size_cap) X(rmq_inner_dist) X(rmq_rescue_size) X(mask_len) X(best_n) X(a) X(b) X(q) X(e) X(q2) X(e2) \
	X(transition) X(sc_ambi) X(zdrop) X(zdrop_inv) X(end_bonus) X(min_dp_max) X(min_ksw_len) X(anchor_ext_len) X(anchor_ext_shift) \
	X(rank_min_len) X(min_mid_occ) X(max_mid_occ) X(mid_occ) X(max_occ) X(max_max_occ) X(occ_dist)
#define FLOAT_FIELDS(X) X(chain_gap_scale) X(chain_skip_scale) X(rmq_rescue_ratio) X(mask_level) X(pri_ratio) X(alt_drop) \
	X(max_clip_ratio) X(rank_frac) X(mid_occ_frac) X(q_occ_frac)

/* apply one setting; returns 0 if it is not recognized */
static int apply(const char *s, mm_idxopt_t *io, mm_mapopt_t *mo)
{
	const char *eq = strchr(s, '=');
	if (s[0] == '+' || s[0] == '-') {
		int64_t bit = flag_bit(s + 1);
		if (bit == 0) return 0;
		if (s[0] == '+') mo->flag |= bit;
		else mo->flag &= ~bit;
		return 1;
	}
	if (eq == 0) return 0;
#define SET_INT(f) if ((size_t)(eq - s) == strlen(#f) && strncmp(s, #f, eq - s) == 0) { mo->f = atoi(eq + 1); return 1; }
#define SET_FLOAT(f) if ((size_t)(eq - s) == strlen(#f) && strncmp(s, #f, eq - s) == 0) { mo->f = atof(eq + 1); return 1; }
	INT_FIELDS(SET_INT)
	FLOAT_FIELDS(SET_FLOAT)
	if (strncmp(s, "max_sw_mat=", 11) == 0) { mo->max_sw_mat = atoll(eq + 1); return 1; }
	if (strncmp(s, "k=", 2) == 0) { io->k = atoi(eq + 1); return 1; }
	if (strncmp(s, "w=", 2) == 0) { io->w = atoi(eq + 1); return 1; }
	if (strncmp(s, "H=", 2) == 0) { if (atoi(eq + 1)) io->flag |= 1; else io->flag &= ~1; return 1; } /* MM_I_HPC */
	return 0;
}

static int case_map(int argc, char *argv[])
{
	mm_idxopt_t iopt;
	mm_mapopt_t mopt;
	mm_idx_reader_t *r;
	mm_idx_t *mi;
	mm_tbuf_t *tbuf;
	mmt_seqs_t s;
	int i, j;
	if (argc < 3) return MMT_FAIL("usage: map <target> <reads> <preset|-> [settings]"), 1;
	mm_verbose = 1;
	MMT_REQUIRE(mm_set_opt(0, &iopt, &mopt) == 0);
	if (strcmp(argv[2], "-") != 0) MMT_REQUIRE(mm_set_opt(argv[2], &iopt, &mopt) == 0);
	for (i = 3; i < argc; ++i)
		MMT_CHECKF(apply(argv[i], &iopt, &mopt), "unknown setting %s", argv[i]);
	MMT_REQUIRE(mmt_read_seqs(argv[1], &s) > 0);
	r = mm_idx_reader_open(argv[0], &iopt, 0);
	MMT_REQUIRE(r != 0);
	mi = mm_idx_reader_read(r, 1);
	MMT_REQUIRE(mi != 0);
	mm_idx_reader_close(r);
	mm_mapopt_update(&mopt, mi);
	tbuf = mm_tbuf_init();
	for (i = 0; i < s.n; ++i) {
		mm_reg1_t *reg;
		int n_reg, n_pri = 0;
		reg = mm_map(mi, s.len[i], s.seq[i], &n_reg, tbuf, &mopt, s.name[i]);
		for (j = 0; j < n_reg; ++j) {
			mm_reg1_t *h = &reg[j];
			uint32_t k;
			MMT_CHECKF(h->id == j, "%s hit %d: id %d", s.name[i], j, h->id);
			MMT_CHECKF(h->parent == -1 || (h->parent >= 0 && h->parent < n_reg), "%s hit %d: parent %d", s.name[i], j, h->parent);
			if (h->sam_pri) ++n_pri;
			printf("%s\t%d\t%d\t%d\t%c\t%s\t%u\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%d",
				   s.name[i], s.len[i], h->qs, h->qe, "+-"[h->rev], mi->seq[h->rid].name, mi->seq[h->rid].len, h->rs, h->re, h->mlen, h->blen,
				   h->mapq, h->id, h->parent, h->sam_pri, h->inv, h->split, h->split_inv, h->is_alt, h->cnt, h->score, h->score0, h->subsc,
				   h->n_sub, h->div, tbuf->rep_len);
			if (h->p) {
				printf("\t1\t%d\t%d\t%d\t%d\t%u\t%u\t%u\t%u\t", h->p->dp_score, h->p->dp_max, h->p->dp_max0, h->p->dp_max2, h->p->n_ambi,
					   h->p->trans_strand, h->p->n_cigar, h->p->capacity);
				for (k = 0; k < h->p->n_cigar; ++k)
					printf("%d%c", h->p->cigar[k]>>4, MM_CIGAR_STR[h->p->cigar[k]&0xf]);
				MMT_CHECKF(h->p->capacity >= h->p->n_cigar + sizeof(mm_extra_t) / 4, "%s hit %d: capacity %u for %u operators", s.name[i], j, h->p->capacity, h->p->n_cigar);
			} else printf("\t0\t0\t0\t0\t0\t0\t0\t0\t0\t*");
			putchar('\n');
			free(h->p);
		}
		MMT_CHECKF(n_reg == 0 || n_pri <= 1, "%s: %d hits flagged sam_pri", s.name[i], n_pri);
		free(reg);
	}
	mm_tbuf_destroy(tbuf);
	mm_idx_destroy(mi);
	mmt_free_seqs(&s);
	return 0;
}

/* args: preset FLAG... ; among the flags listed in this driver, mm_set_opt(preset) on top of the
 * defaults sets exactly the named ones */
static int case_opt_flags(int argc, char *argv[])
{
	mm_idxopt_t iopt;
	mm_mapopt_t mopt;
	int i, j;
	if (argc < 1) return MMT_FAIL("usage: opt_flags <preset> [FLAG...]"), 1;
	mm_set_opt(0, &iopt, &mopt);
	MMT_REQUIRE(mm_set_opt(argv[0], &iopt, &mopt) == 0);
	for (i = 0; flags[i].name; ++i) {
		int want = 0;
		for (j = 1; j < argc; ++j)
			if (strcmp(argv[j], flags[i].name) == 0) want = 1;
		MMT_CHECKF(!!(mopt.flag & flags[i].bit) == want, "%s: MM_F_%s is %s", argv[0], flags[i].name, want? "not set" : "set");
	}
	return 0;
}

static const mmt_case_t cases[] = {
	{ "map", case_map },
	{ "opt_flags", case_opt_flags },
	{ 0, 0 }
};

MMT_MAIN(cases)
