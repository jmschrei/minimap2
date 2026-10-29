"""Helpers for the tests of output formats, the command line and the libminimap2 API (test_d_*.py).

Contents:
	mixed_reads           reads of every output kind (unique, low-complexity, repeat copies,
	                      chimeric, ambiguous bases, unmapped) on the standard genome
	write_mixed           the mixed reads as FASTQ with or without comments
	junction_genome       a transcript genome whose introns carry GT-AG, GC-AG, AT-AC and other
	                      splice signals, with one gene duplicated
	lines / records_by    small text utilities
	sa_entries            parse SA:Z
	cigar_to_bam          BAM-encoded CIGAR operations (CG:B:I) of a CIGAR string
	eqx_cigar / merge_eqx expected --eqx CIGAR and its inverse
	expected_ds           the ds:Z value of an alignment, from the definition of indel ambiguity
	junction_score        the --write-junc score of an intron (minimap2.1)
	mm_version            MM_VERSION from minimap.h of the source tree under test
"""

import os
import re

import minimap2_build
import minimap2_helpers as H

_CACHE = {}


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------


class Mixed(object):
	"""Reads on the standard synthetic genome, one list per kind.

	records: [(name, seq, qual, comment)] in file order; kind: {name: kind}; seqs and quals:
	{name: as sequenced}; truth: {name: SimRead} for simulated reads; spans: {name: [(contig,
	start, end, strand)]} of every reference interval a read was copied from (one per piece).
	Kinds: 'unique' (ONT-like from unique sequence), 'ambig' (HiFi-like with N bases and a
	lowercase stretch), 'lowcomp' (covering a low-complexity tract, so rl:i > 0), 'repeat' (an
	exact copy of one repeat unit, so secondary hits), 'chimeric' (two distant pieces, so a
	supplementary hit; the first one has three N bases in each piece) and 'random' (no hit).
	"""

	def __init__(self):
		self.records, self.kind, self.seqs, self.quals, self.truth, self.spans, self.comments = [], {}, {}, {}, {}, {}, {}

	def add(self, name, seq, qual, kind, comment=None, truth=None, spans=()):
		self.records.append((name, seq, qual, comment))
		self.kind[name], self.seqs[name], self.quals[name] = kind, seq, qual
		self.comments[name] = comment
		if truth is not None:
			self.truth[name] = truth
		self.spans[name] = list(spans)

	def names(self, kind=None):
		return [r[0] for r in self.records if kind is None or self.kind[r[0]] == kind]


def _qual(rng, n):
	return "".join(chr(35 + rng.below(38)) for _ in range(n))


def mixed_reads(genome):
	"""The Mixed read set for `genome` (the standard synthetic genome of conftest.make_synth())."""
	key = ("mixed", genome.digest())
	if key in _CACHE:
		return _CACHE[key]
	rng = H.Rng(500)
	m = Mixed()
	for i, r in enumerate(H.simulate_long_reads(genome, 10, profile="ont", seed=501, unique=True, length=(1500, 4000))):
		comment = None
		if i % 3 == 0:
			comment = "BC:Z:ACGT%d" % i
		elif i % 3 == 1:
			comment = "BC:Z:TT%d\tXI:i:%d" % (i, i)
		m.add("u%d" % i, r.seq, r.qual, "unique", comment, r, [(r.contig, r.pos, r.end, r.strand)])
	# N bases and a lowercase stretch inside an otherwise ordinary read
	r = H.simulate_long_reads(genome, 1, profile="hifi", seed=502, unique=True, length=3000)[0]
	s = list(r.seq)
	for k in range(1500, 1506):
		s[k] = "N"
	for k in range(400, 520):
		s[k] = s[k].lower()
	m.add("ambig0", "".join(s), r.qual, "ambig", "BC:Z:AMB", r, [(r.contig, r.pos, r.end, r.strand)])
	# reads covering low-complexity tracts: their repetitive minimizers give rl:i > 0
	for j, f in enumerate([f for f in genome.features if f.kind == "lowcomp"][:2]):
		st, en = max(0, f.start - 900), min(genome.length(f.contig), f.end + 900)
		strand = "+" if j == 0 else "-"
		seq = H.mutate(rng, genome.fetch(f.contig, st, en, strand), 0.01)
		m.add("lc%d" % j, seq, _qual(rng, len(seq)), "lowcomp", None, None, [(f.contig, st, en, strand)])
	# one exact copy of each repeat family (plus one reverse complement): secondary hits
	fams = {}
	for f in genome.features:
		if f.kind == "repeat":
			fams.setdefault(f.info[0], []).append(f)
	j = 0
	for fam in sorted(fams):
		for f, strand in ((fams[fam][0], "+"), (fams[fam][1], "-"))[:2 if fam == 0 else 1]:
			seq = genome.fetch(f.contig, f.start + 5, f.end - 5, strand)
			copies = [(g.contig, g.start + 5, g.end - 5) for g in fams[fam]]
			m.add("rep%d" % j, seq, _qual(rng, len(seq)), "repeat", None, None, [(c, s, e, "?") for c, s, e in copies])
			j += 1
	for k, r in enumerate(H.simulate_chimeric_reads(genome, 4, profile="hifi", seed=503, length=(3000, 6000), min_part=1000)):
		seq = r.seq
		if k == 0:  # N bases in each piece: they count in NM, also in the NM of SA:Z
			s = list(seq)
			for g in r.segments:
				for j in range(g.qstart + 300, g.qstart + 303):
					s[j] = "N"
			seq = "".join(s)
		m.add("chim%s" % r.name[4:], seq, r.qual, "chimeric", "BC:Z:CH" if r.name.endswith("1") else None, r,
			[(g.contig, g.pos, g.end, g.strand) for g in r.segments])
	for i in range(2):
		seq = rng.dna(1500)
		m.add("rnd%d" % i, seq, _qual(rng, len(seq)), "random", "BC:Z:RND" if i == 0 else None)
	_CACHE[key] = m
	return m


def write_mixed(path, m, comments=True, fmt="fastq"):
	"""Write the Mixed reads to path (FASTQ, or FASTA with fmt='fasta'); comments=False drops comments."""
	recs = [(n, s, q, c if comments else None) for n, s, q, c in m.records]
	return H.write_fastx(str(path), recs, fmt=fmt)


def recode_intron(seqs, gene, idx, donor, acceptor):
	"""Write the splice signal (donor, acceptor), given on the transcript strand, into intron idx of gene."""
	s, e = gene.introns[idx]
	c = list(seqs[gene.contig])
	if gene.strand == "+":
		c[s:s + 2], c[e - 2:e] = list(donor), list(acceptor)
	else:
		c[s:s + 2], c[e - 2:e] = list(H.revcomp(acceptor)), list(H.revcomp(donor))
	seqs[gene.contig] = "".join(c)


def junction_genome():
	"""(Genome, signals) for the --write-junc tests.

	A transcript genome (6 genes of 3-6 exons). The first intron of every gene keeps GT-AG, so
	every read has a canonical junction to orient it; later introns carry GC-AG, AT-AC, GT-AC,
	AT-AG and CA-TT (a pair minimap2.1 scores 0). A third contig holds an exact copy of the first
	gene with 1 kbp of flanks, so reads from that gene map with MAPQ 0. signals maps
	(contig, intron start, intron end) to the (donor, acceptor) pair written there.
	"""
	key = ("junc",)
	if key in _CACHE:
		return _CACHE[key]
	g = H.make_genome(seed=31, contig_lengths=(50000, 40000), genes=6, exons=(3, 6), exon_len=(150, 400),
		intron_len=(200, 1500))
	seqs = dict(g.seqs)
	variants = [("GC", "AG"), ("AT", "AC"), ("GT", "AC"), ("AT", "AG"), ("CA", "TT")]
	signals, k = {}, 0
	for gene in g.genes:
		for i, (s, e) in enumerate(gene.introns):
			if i == 0:
				d, a = "GT", "AG"
			else:
				d, a = variants[k % len(variants)]
				k += 1
				recode_intron(seqs, gene, i, d, a)
			signals[(gene.contig, s, e)] = (d, a)
	dup = g.genes[0]
	st, en = dup.start - 1000, dup.end + 1000
	contigs = [(n, seqs[n]) for n in g.names] + [("dup", seqs[dup.contig][st:en])]
	out = H.Genome(contigs, g.features, g.genes)
	_CACHE[key] = (out, signals, dup.name)
	return _CACHE[key]


# ---------------------------------------------------------------------------
# Text utilities
# ---------------------------------------------------------------------------


def lines(text):
	"""Non-empty lines of text."""
	return [l for l in text.split("\n") if l]


def body_lines(sam_text):
	return [l for l in lines(sam_text) if not l.startswith("@")]


def drop_tag(line, name):
	"""A SAM or PAF line without the optional field `name`."""
	return "\t".join(f for f in line.split("\t") if not f.startswith(name + ":"))


def mm_version():
	"""MM_VERSION as #defined in minimap.h of the source tree under test."""
	with open(os.path.join(minimap2_build.src_dir(), "minimap.h")) as f:
		m = re.search(r'#define\s+MM_VERSION\s+"([^"]+)"', f.read())
	return m.group(1)


# ---------------------------------------------------------------------------
# SAM fields
# ---------------------------------------------------------------------------


def sa_entries(rec):
	"""SA:Z of a SAM record as [(rname, pos, strand, cigar, mapq, nm)]."""
	v = rec.tags.get("SA")
	if v is None:
		return []
	assert v.endswith(";"), v
	out = []
	for e in v[:-1].split(";"):
		rname, pos, strand, cigar, mapq, nm = e.split(",")
		out.append((rname, int(pos), strand, cigar, int(mapq), int(nm)))
	return out


def hard_clips(cigar):
	"""(leading, trailing) hard-clip lengths of a CIGAR."""
	ops = H.parse_cigar(cigar)
	a = ops[0][0] if ops and ops[0][1] == "H" else 0
	b = ops[-1][0] if len(ops) > 1 and ops[-1][1] == "H" else 0
	return a, b


def record_seq(rec, read):
	"""What SEQ of a mapped record must be for `read` (as sequenced): the read in alignment
	orientation without its hard-clipped bases."""
	s = H.revcomp(read) if rec.is_reverse else read
	h5, h3 = hard_clips(rec.cigar)
	return s[h5:len(s) - h3]


_BAM_OPS = "MIDNSHP=X"


def cigar_to_bam(cigar):
	"""BAM encoding of a CIGAR string: len<<4 | op with op = index in 'MIDNSHP=X' (SAM spec 4.2)."""
	return [n << 4 | _BAM_OPS.index(op) for n, op in H.parse_cigar(cigar)]


def bam_to_cigar(vals):
	return "".join("%d%s" % (v >> 4, _BAM_OPS[v & 0xf]) for v in vals)


def walk_columns(ref, pos0, seq, cigar):
	"""Yield (op, ref_index, read_index, length) for each CIGAR operation; seq is the read in
	reference orientation with soft clips; hard clips are skipped."""
	r, q = pos0, 0
	for n, op in H.parse_cigar(cigar):
		if op == "H":
			continue
		yield op, r, q, n
		if op in "M=X":
			r += n
			q += n
		elif op in "IS":
			q += n
		elif op in "DN":
			r += n


def eqx_cigar(ref, pos0, seq, cigar):
	"""The CIGAR with each M run split into = (equal minimap2 base codes) and X (different) runs."""
	out = []
	for op, r, q, n in walk_columns(ref, pos0, seq, cigar):
		if op != "M":
			out.append((n, op))
			continue
		for k in range(n):
			o = "=" if H.nt4(ref[r + k]) == H.nt4(seq[q + k]) else "X"
			if out and out[-1][1] == o:
				out[-1] = (out[-1][0] + 1, o)
			else:
				out.append((1, o))
	return "".join("%d%s" % x for x in out)


def merge_eqx(cigar):
	"""=/X operations turned back into M runs."""
	out = []
	for n, op in H.parse_cigar(cigar):
		op = "M" if op in "=X" else op
		if out and out[-1][1] == op:
			out[-1] = (out[-1][0] + n, op)
		else:
			out.append((n, op))
	return "".join("%d%s" % x for x in out)


def _codes(s):
	return [H.nt4(c) for c in s]


def _indel_ds(seq_codes, start, length, lo, hi):
	"""ds text of the indel seq_codes[start:start+length] inside the aligned interval [lo, hi) of
	the same sequence (the read for an insertion, the reference for a deletion).

	right = how far the indel can slide right without changing the sequence on either side
	(bases after it repeat its first bases), left = how far it can slide left (bases before it
	repeat its last bases), each at most its length and never past the aligned interval. The
	indel is written in square brackets if right + left covers it, else its first `right` and
	last `left` bases are bracketed (the extension of cs:Z written by --ds).
	"""
	x = seq_codes
	right = 0
	while right < length and start + length + right < hi and x[start + length + right] == x[start + right]:
		right += 1
	left = 0
	while left < length and start - left - 1 >= lo and x[start + length - left - 1] == x[start - left - 1]:
		left += 1
	b = "".join("acgtn"[c] for c in x[start:start + length])
	if right + left >= length:
		return "[" + b + "]"
	mid = b[right:length - left]
	return ("[" + b[:right] + "]" if right else "") + mid + ("[" + b[length - left:] + "]" if left else "")


def expected_ds(ref, pos0, seq, cigar):
	"""The ds:Z value (without the prefix) of an alignment: the short cs:Z with each indel's
	ambiguous bases in square brackets (see _indel_ds). seq is the read in reference orientation,
	soft clips included; hard clips are ignored."""
	ops = [x for x in walk_columns(ref, pos0, seq, cigar) if x[0] != "S"]
	if not ops:
		return ""
	r_lo = ops[0][1]
	q_lo = ops[0][2]
	r_hi = r_lo + sum(n for op, _, _, n in ops if op in "MDN=X")
	q_hi = q_lo + sum(n for op, _, _, n in ops if op in "MI=X")
	rc, qc = _codes(ref[r_lo:r_hi]), _codes(seq[q_lo:q_hi])
	out, run = [], 0
	for op, r, q, n in ops:
		if op in "M=X":
			for k in range(n):
				a, b = rc[r - r_lo + k], qc[q - q_lo + k]
				if a == b:
					run += 1
				else:
					if run:
						out.append(":%d" % run)
					run = 0
					out.append("*%s%s" % ("acgtn"[a], "acgtn"[b]))
			continue
		if run:
			out.append(":%d" % run)
		run = 0
		if op == "I":
			out.append("+" + _indel_ds(qc, q - q_lo, n, 0, q_hi - q_lo))
		elif op == "D":
			out.append("-" + _indel_ds(rc, r - r_lo, n, 0, r_hi - r_lo))
		else:
			s = "".join("acgtn"[c] for c in rc[r - r_lo:r - r_lo + n])
			out.append("~%s%d%s" % (s[:2], n, s[-2:]))
	if run:
		out.append(":%d" % run)
	return "".join(out)


def strip_ds(ds):
	"""ds:Z without its square brackets: the short cs:Z of the same alignment."""
	return ds.replace("[", "").replace("]", "")


# ---------------------------------------------------------------------------
# --write-junc
# ---------------------------------------------------------------------------

DONOR_SCORE = {"GT": 3, "GC": 2, "AT": 1}
ACCEPTOR_SCORE = {"AG": 3, "AC": 1}


def junction_score(ref, start, end, strand):
	"""minimap2.1, --write-junc: donor GT 3, GC 2, AT 1; acceptor AG 3, AC 1; read on the
	transcript strand of the intron [start, end)."""
	s = ref[start:end].upper()
	if strand == "-":
		s = H.revcomp(s)
	return DONOR_SCORE.get(s[:2], 0) + ACCEPTOR_SCORE.get(s[-2:], 0)


def expected_junctions(sam, genome):
	"""The --write-junc lines the SAM records of the same run imply: for each primary or
	supplementary record with MAPQ >= 10 and a transcript strand (ts:A + or -), one line per N
	operation: contig, start, end, read name, score, strand of the intron on the genome."""
	seqs = genome.seqs if isinstance(genome, H.Genome) else genome
	out = []
	for rec in sam.records:
		if rec.is_unmapped or rec.is_secondary or rec.mapq < 10 or rec.tags.get("ts") not in ("+", "-"):
			continue
		strand = "+" if (rec.tags["ts"] == "+") == (not rec.is_reverse) else "-"
		for s, e in rec.introns:
			out.append("%s\t%d\t%d\t%s\t%d\t%s" % (rec.rname, s, e, rec.qname, junction_score(seqs[rec.rname], s, e, strand), strand))
	return out
