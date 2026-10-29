/* C driver for test_pilot.py: the libminimap2 API used the way example.c (minimap2-lite)
 * uses it. map_truth maps simulated reads with mm_map() and checks each
 * read's primary hit against the simulator's truth, and name lookup on the
 * loaded index; idx_str builds an index from strings with mm_idx_str() and
 * reads names, lengths and bases back. */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "minimap.h"
#include "minimap2_test.h"

typedef struct {
	char name[256], ctg[256], strand;
	long pos, end;
} truth_t;

static int read_truth(const char *fn, int n_max, truth_t *t)
{
	FILE *fp;
	int n = 0;
	if ((fp = fopen(fn, "r")) == 0) return -1;
	while (n < n_max && fscanf(fp, "%255s %255s %ld %ld %c", t[n].name, t[n].ctg, &t[n].pos, &t[n].end, &t[n].strand) == 5)
		++n;
	fclose(fp);
	return n;
}

/* args: <target.fa|target.mmi> <preset or "-"> <reads.fa> <truth.tsv> <max_dist>
 * truth.tsv: name, contig, 0-based start, end, strand per read */
static int case_map_truth(int argc, char *argv[])
{
	mm_idxopt_t iopt;
	mm_mapopt_t mopt;
	mm_idx_reader_t *r;
	mm_idx_t *mi;
	mm_tbuf_t *tbuf;
	mmt_seqs_t s;
	truth_t *t;
	int32_t i, n_t, max_dist, n_parts = 0;
	if (argc < 5) return MMT_FAIL("usage: map_truth <target> <preset|-> <reads.fa> <truth.tsv> <max_dist>"), 1;
	max_dist = atoi(argv[4]);
	mm_verbose = 2;
	MMT_REQUIRE(mm_set_opt(0, &iopt, &mopt) == 0);
	if (strcmp(argv[1], "-") != 0)
		MMT_REQUIRE(mm_set_opt(argv[1], &iopt, &mopt) == 0);
	mopt.flag |= MM_F_CIGAR;
	MMT_REQUIRE(mmt_read_seqs(argv[2], &s) > 0);
	t = (truth_t*)calloc(s.n, sizeof(truth_t));
	n_t = read_truth(argv[3], s.n, t);
	MMT_CHECK_INT(n_t, ==, s.n);
	r = mm_idx_reader_open(argv[0], &iopt, 0);
	MMT_REQUIRE(r != 0);
	mi = mm_idx_reader_read(r, 1);
	MMT_REQUIRE(mi != 0);
	++n_parts;
	MMT_CHECK(mm_idx_reader_eof(r));
	mm_mapopt_update(&mopt, mi);
	MMT_CHECK(mopt.mid_occ > 0);
	MMT_CHECK_INT(mm_idx_name2id(mi, t[0].ctg), ==, -2); // no name hash yet
	MMT_CHECK_INT(mm_idx_index_name(mi), ==, 0); // 0: no duplicate names
	for (i = 0; i < (int32_t)mi->n_seq; ++i)
		MMT_CHECK_INT(mm_idx_name2id(mi, mi->seq[i].name), ==, i);
	MMT_CHECK_INT(mm_idx_name2id(mi, "no-such-contig"), ==, -1);
	tbuf = mm_tbuf_init();
	for (i = 0; i < s.n && i < n_t; ++i) {
		mm_reg1_t *reg;
		int32_t j, n_reg, n_pri = 0;
		MMT_CHECK_STR(s.name[i], t[i].name);
		reg = mm_map(mi, s.len[i], s.seq[i], &n_reg, tbuf, &mopt, s.name[i]);
		MMT_CHECKF(n_reg > 0, "read %s has no hit", s.name[i]);
		for (j = 0; j < n_reg; ++j) {
			mm_reg1_t *h = &reg[j];
			uint32_t k, qspan = 0, tspan = 0;
			const char *name;
			long start5;
			MMT_CHECKF(h->rid >= 0 && (uint32_t)h->rid < mi->n_seq, "read %s: hit %d has rid %d", s.name[i], j, h->rid);
			if (h->rid < 0 || (uint32_t)h->rid >= mi->n_seq) continue;
			name = mi->seq[h->rid].name;
			MMT_CHECKF(0 <= h->rs && h->rs < h->re && (uint32_t)h->re <= mi->seq[h->rid].len, "read %s: hit %d spans %d-%d", s.name[i], j, h->rs, h->re);
			MMT_CHECKF(0 <= h->qs && h->qs < h->qe && h->qe <= s.len[i], "read %s: hit %d query span %d-%d", s.name[i], j, h->qs, h->qe);
			MMT_CHECKF(h->p != 0 && h->p->n_cigar > 0, "read %s: hit %d has no base alignment", s.name[i], j);
			MMT_CHECKF(h->mlen <= h->blen && h->mapq <= 60, "read %s: hit %d mlen %d blen %d mapq %d", s.name[i], j, h->mlen, h->blen, h->mapq);
			if (h->p) {
				for (k = 0; k < h->p->n_cigar; ++k) {
					int op = h->p->cigar[k] & 0xf, len = h->p->cigar[k] >> 4;
					if (op == MM_CIGAR_MATCH || op == MM_CIGAR_EQ_MATCH || op == MM_CIGAR_X_MISMATCH) qspan += len, tspan += len;
					else if (op == MM_CIGAR_INS) qspan += len;
					else if (op == MM_CIGAR_DEL || op == MM_CIGAR_N_SKIP) tspan += len;
				}
				MMT_CHECKF(qspan == (uint32_t)(h->qe - h->qs) && tspan == (uint32_t)(h->re - h->rs), "read %s: hit %d CIGAR spans %u/%u, coordinates %d/%d",
					s.name[i], j, qspan, tspan, h->qe - h->qs, h->re - h->rs);
			}
			if (h->id != h->parent || !h->sam_pri) continue;
			++n_pri;
			start5 = (long)h->rs - (h->rev? s.len[i] - h->qe : h->qs); // leftmost reference position of the unclipped read
			MMT_CHECK_STR(name, t[i].ctg);
			MMT_CHECKF("+-"[h->rev] == t[i].strand, "read %s: strand %c, truth %c", s.name[i], "+-"[h->rev], t[i].strand);
			MMT_CHECKF(labs(start5 - t[i].pos) <= max_dist, "read %s: unclipped start %ld, truth %ld", s.name[i], start5, t[i].pos);
		}
		MMT_CHECKF(n_reg == 0 || n_pri == 1, "read %s: %d hits flagged primary", s.name[i], n_pri);
		for (j = 0; j < n_reg; ++j) free(reg[j].p);
		free(reg);
	}
	mm_tbuf_destroy(tbuf);
	mm_idx_destroy(mi);
	MMT_CHECK(mm_idx_reader_read(r, 1) == 0); // a single-part index or a FASTA under 8 Gbp
	mm_idx_reader_close(r);
	free(t);
	mmt_free_seqs(&s);
	return 0;
}

/* args: none. mm_idx_str() on two sequences; names, lengths and bases back */
static int case_idx_str(int argc, char *argv[])
{
	const char *seq[2] = { "ACGTACGTTTGACCAGTAGGCATCGATCGGGATCCATGCAAAGGGCCTTTAACGT", "TTTTGGGGCCCCAAAANNACGT" };
	const char *name[2] = { "first", "second" };
	mm_idx_t *mi;
	uint8_t buf[64];
	int i, n;
	mi = mm_idx_str(5, 11, 0, 14, 2, seq, name);
	MMT_REQUIRE(mi != 0);
	MMT_CHECK_INT(mi->n_seq, ==, 2);
	MMT_CHECK_INT(mi->k, ==, 11);
	MMT_CHECK_INT(mi->w, ==, 5);
	MMT_CHECK_STR(mi->seq[0].name, "first");
	MMT_CHECK_STR(mi->seq[1].name, "second");
	MMT_CHECK_INT(mi->seq[0].len, ==, strlen(seq[0]));
	MMT_CHECK_INT(mi->seq[1].len, ==, strlen(seq[1]));
	// only an absent name: mm_idx_str() puts names in the hash without setting their ids
	MMT_CHECK_INT(mm_idx_name2id(mi, "absent"), ==, -1);
	n = mm_idx_getseq(mi, 1, 2, 20, buf);
	MMT_CHECK_INT(n, ==, 18);
	for (i = 0; i < n; ++i) {
		int c = seq[1][2 + i], code = c == 'A'? 0 : c == 'C'? 1 : c == 'G'? 2 : c == 'T'? 3 : 4;
		MMT_CHECKF(buf[i] == code, "base %d: %d, want %d", 2 + i, buf[i], code);
	}
	mm_idx_destroy(mi);
	return 0;
}

static const mmt_case_t cases[] = {
	{ "map_truth", case_map_truth },
	{ "idx_str", case_idx_str },
	{ 0, 0 }
};

MMT_MAIN(cases)
