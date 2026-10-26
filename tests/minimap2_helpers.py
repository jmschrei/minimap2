"""Shared helpers for the minimap2 test suite.

Standard library only. Everything random here is driven by `Rng`, a
splitmix64 generator written with integer arithmetic, so a given seed gives
the same genome and the same reads on every Python version and platform.

Contents:
	Rng                         deterministic random numbers
	make_genome / Genome / Gene synthetic genomes with repeats, inverted
	                            repeats, N runs, low-complexity sequence, short
	                            contigs and genes with GT-AG introns
	simulate_reads / _pairs /   reads with recorded truth (contig, position,
	_long_reads / _chimeric_reads / strand, CIGAR, edits, exons)
	_spliced_reads / _spliced_pairs / _contigs
	write_fastx / write_pairs / FASTA/FASTQ writers, plain or gzip
	write_interleaved / genes_bed12 / genes_gtf / junctions_bed
	read_fastx / read_fasta     FASTA/FASTQ readers, plain or gzip
	index_genome                `minimap2 -x PRESET -d`, cached per genome and preset
	kmer_hash / SketchRef /     reference (w,k)-minimizers (minimap2.1's rule and
	minimizers / expected_index sketch.c's hash, with HPC)
	read_mmi / MmiPart          a reader for the .mmi format
	run_mm2 / map_reads         run the binary with a timeout and no core dumps
	parse_sam / parse_paf       strict SAM and PAF parsers with typed tags
	compute_nm / compute_md /   recompute alignment tags from the genome
	compute_cs / compute_mlen_blen / compute_de / cigar_*
	sam_invariant_errors /      generic SAM and PAF consistency checks
	paf_invariant_errors / pair_invariant_errors
	check_pinned / strip_pg     compare output against a stored pin
	run_c_case / assert_c_case  run a case of a C driver in tests/c/
	run_k8 / run_paftools       run misc/*.js with k8
"""

import difflib
import gzip
import hashlib
import io
import os
import re
import shutil
import subprocess
import sys
import tempfile

import minimap2_build

TESTS_DIR = minimap2_build.TESTS_DIR
REPO_ROOT = minimap2_build.REPO_ROOT
DATA_DIR = os.path.join(REPO_ROOT, "test")
PINS_DIR = os.path.join(TESTS_DIR, "pins")

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------


def test_threads():
	"""Worker threads (-t) for tests that do not pin it: $MINIMAP2_TEST_THREADS, default 1."""
	try:
		n = int(os.environ.get("MINIMAP2_TEST_THREADS", "1"))
	except ValueError:
		raise ValueError("MINIMAP2_TEST_THREADS must be an integer")
	if n < 1:
		raise ValueError("MINIMAP2_TEST_THREADS must be at least 1")
	return n


def run_timeout():
	"""Timeout in seconds for one minimap2 (or driver, k8, python) process: $MINIMAP2_TEST_TIMEOUT, default 120."""
	return float(os.environ.get("MINIMAP2_TEST_TIMEOUT", "120"))


# ---------------------------------------------------------------------------
# Deterministic random numbers
# ---------------------------------------------------------------------------

_M64 = (1 << 64) - 1


class Rng(object):
	"""splitmix64. Integer arithmetic only, so streams match across Python versions."""

	def __init__(self, seed):
		self.state = int(seed) & _M64

	def u64(self):
		self.state = (self.state + 0x9E3779B97F4A7C15) & _M64
		z = self.state
		z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & _M64
		z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & _M64
		return z ^ (z >> 31)

	def below(self, n):
		"""Uniform integer in [0, n)."""
		n = int(n)
		if n <= 0:
			raise ValueError("below() needs n > 0")
		lim = (1 << 64) - ((1 << 64) % n)
		while True:
			x = self.u64()
			if x < lim:
				return x % n

	def randint(self, lo, hi):
		"""Uniform integer in [lo, hi], both inclusive."""
		return lo + self.below(hi - lo + 1)

	def random(self):
		"""Uniform float in [0, 1) with 53 random bits; exact, hence portable."""
		return (self.u64() >> 11) * (1.0 / (1 << 53))

	def chance(self, p):
		return p > 0 and self.random() < p

	def choice(self, seq):
		return seq[self.below(len(seq))]

	def weighted_index(self, weights):
		"""Index i with probability weights[i] / sum(weights); integer weights."""
		tot = sum(weights)
		x = self.below(tot)
		for i, w in enumerate(weights):
			if x < w:
				return i
			x -= w
		raise AssertionError("unreachable")

	def normal_int(self, mean, sd):
		"""Approximately normal integer (Irwin-Hall sum of 12 uniforms)."""
		s = 0
		for _ in range(12):
			s += self.below(1 << 16)
		s -= 6 << 16
		return int(mean + int(s * sd / float(1 << 16)))

	def dna(self, n, alphabet="ACGT"):
		"""n letters drawn uniformly from alphabet (a string of any length >= 1). A four-letter
		alphabet takes 2 bits per letter, 32 letters per u64; any other size draws each letter with
		below(), so it gives a different stream from a four-letter one."""
		if not alphabet:
			raise ValueError("dna() needs a non-empty alphabet")
		if len(alphabet) != 4:
			return "".join(alphabet[self.below(len(alphabet))] for _ in range(n))
		out = []
		while len(out) < n:
			x = self.u64()
			for _ in range(32):
				out.append(alphabet[x & 3])
				x >>= 2
		return "".join(out[:n])

	def fork(self, tag):
		"""An independent generator derived from this one's state and a tag."""
		h = hashlib.sha256(("%d:%s" % (self.state, tag)).encode()).digest()
		return Rng(int.from_bytes(h[:8], "little"))


# ---------------------------------------------------------------------------
# Sequences
# ---------------------------------------------------------------------------

# the complement minimap2 uses (bseq.c seq_comp_table) for the bases a test may write
_COMP = {}
for _a, _b in zip("ACGTUNRYKMSWBDHVacgtunrykmswbdhv", "TGCAANYRMKSWVHDBtgcaanyrmkswvhdb"):
	_COMP[_a] = _b
_COMP_TABLE = str.maketrans(_COMP)


def revcomp(seq):
	return seq.translate(_COMP_TABLE)[::-1]


_NT4 = {"A": 0, "C": 1, "G": 2, "T": 3, "U": 3, "a": 0, "c": 1, "g": 2, "t": 3, "u": 3}


def nt4(c):
	"""minimap2's base code (seq_nt4_table): A/C/G/T/U in any case -> 0..3, anything else -> 4."""
	return _NT4.get(c, 4)


def nt5_char(c, lower=False):
	"""The base minimap2 prints for a sequence character: A/C/G/T, or N for anything else (U -> T)."""
	s = "ACGTN"[nt4(c)]
	return s.lower() if lower else s


def mutate(rng, seq, sub_rate):
	"""Substitutions only, at rate sub_rate; returns the new sequence."""
	out = list(seq)
	for i, c in enumerate(out):
		if c.upper() in "ACGT" and rng.chance(sub_rate):
			out[i] = _other_base(rng, c)
	return "".join(out)


def _other_base(rng, c):
	alts = [b for b in "ACGT" if b != c.upper()]
	return alts[rng.below(3)]


# ---------------------------------------------------------------------------
# Synthetic genomes
# ---------------------------------------------------------------------------


class Feature(object):
	"""An annotated interval [start, end) on a contig.

	kind is 'repeat' (info: family, copy, strand), 'inversion' (one arm of an
	inverted repeat; info: pair, arm, strand), 'N' or 'lowcomp' (info: unit).
	"""
	__slots__ = ("contig", "start", "end", "kind", "info")

	def __init__(self, contig, start, end, kind, info=None):
		self.contig, self.start, self.end, self.kind, self.info = contig, start, end, kind, info

	def __repr__(self):
		return "Feature(%r, %d, %d, %r, %r)" % (self.contig, self.start, self.end, self.kind, self.info)


class Gene(object):
	"""A gene model: exons [(start, end), ...] in ascending genome order on `contig`.

	strand '+' or '-' is the transcript strand on the genome. Every intron
	starts with GT and ends with AG on the transcript strand, so on the
	genome's forward strand a '-' gene's introns read CT...AC.
	"""
	__slots__ = ("name", "contig", "strand", "exons")

	def __init__(self, name, contig, strand, exons):
		self.name, self.contig, self.strand, self.exons = name, contig, strand, list(exons)

	@property
	def start(self):
		return self.exons[0][0]

	@property
	def end(self):
		return self.exons[-1][1]

	@property
	def introns(self):
		return [(self.exons[i][1], self.exons[i + 1][0]) for i in range(len(self.exons) - 1)]

	@property
	def length(self):
		"""Transcript length (sum of exon lengths)."""
		return sum(e - s for s, e in self.exons)

	def transcript(self, genome):
		"""The mature transcript, 5' to 3'."""
		s = "".join(genome[self.contig][a:b] for a, b in self.exons)
		return s if self.strand == "+" else revcomp(s)

	def __repr__(self):
		return "Gene(%s, %s:%d-%d%s, %d exons)" % (self.name, self.contig, self.start, self.end, self.strand, len(self.exons))


class Genome(object):
	"""Contigs in order, annotated features and genes."""

	def __init__(self, contigs, features=(), genes=()):
		self.names = [n for n, _ in contigs]
		if len(set(self.names)) != len(self.names):
			raise ValueError("duplicate contig names")
		self.seqs = dict(contigs)
		self.features = list(features)
		self.genes = list(genes)

	def __getitem__(self, name):
		return self.seqs[name]

	def __contains__(self, name):
		return name in self.seqs

	def length(self, name):
		return len(self.seqs[name])

	@property
	def total_length(self):
		return sum(len(s) for s in self.seqs.values())

	def fetch(self, contig, start, end, strand="+"):
		s = self.seqs[contig][start:end]
		return s if strand == "+" else revcomp(s)

	def fasta_text(self, line_width=60):
		out = []
		for n in self.names:
			out.append(">%s\n" % n)
			s = self.seqs[n]
			if line_width and line_width > 0:
				for i in range(0, len(s), line_width):
					out.append(s[i:i + line_width] + "\n")
			else:
				out.append(s + "\n")
		return "".join(out)

	def digest(self):
		return hashlib.sha256(self.fasta_text().encode()).hexdigest()

	def write_fasta(self, path, line_width=60, gz=None):
		_write_text(path, self.fasta_text(line_width), gz)
		return path

	def features_overlapping(self, contig, start, end, kinds=None):
		return [f for f in self.features if f.contig == contig and f.start < end and start < f.end and
			(kinds is None or f.kind in kinds)]

	def is_unique(self, contig, start, end, pad=0):
		"""True if [start-pad, end+pad) overlaps no repeat, inversion arm, N run or low-complexity feature."""
		return not self.features_overlapping(contig, start - pad, end + pad)


def make_genome(seed=1, contig_lengths=(60000, 40000), names=None, repeat_families=0, repeat_len=(300, 1000),
		repeat_copies=3, repeat_divergence=0.0, soft_mask_repeats=False, inversions=0, inversion_len=(300, 1000),
		inversion_gap=(50, 500), n_runs=0, n_run_len=(10, 500), low_complexity=0, low_complexity_len=(50, 300),
		genes=0, exons=(2, 6), exon_len=(80, 400), intron_len=(80, 2000), short_contigs=0, short_contig_len=(20, 200),
		margin=200):
	"""A deterministic random genome.

	Contigs are uniform random DNA. Then, in this order and without
	overlapping each other or coming within `margin` bp of a contig end:
	`repeat_families` families of `repeat_copies` copies each (a random unit of
	length drawn from `repeat_len`, each copy on a random strand with
	substitutions at `repeat_divergence`, lowercase if `soft_mask_repeats`);
	`inversions` inverted repeats (a unit and, `inversion_gap` bp downstream,
	its reverse complement); `genes` gene models (exon count, exon and intron
	lengths drawn from the ranges, strand at random, GT-AG introns written into
	the sequence); `low_complexity` tracts (homopolymer, di- or trinucleotide
	repeats); `n_runs` runs of N. Finally `short_contigs` extra contigs with
	lengths drawn from `short_contig_len` are appended. Every placed
	interval is recorded in Genome.features, and genes in Genome.genes.
	"""
	rng = Rng(seed)
	if names is None:
		names = ["chr%d" % (i + 1) for i in range(len(contig_lengths))]
	if len(names) != len(contig_lengths):
		raise ValueError("names and contig_lengths differ in length")
	seqs = [list(rng.dna(n)) for n in contig_lengths]
	feats, gene_list, occupied = [], [], []

	def place(length):
		ok = [i for i, s in enumerate(seqs) if len(s) >= length + 2 * margin]
		if not ok:
			raise ValueError("no contig long enough for a feature of %d bp" % length)
		for _ in range(1000):
			ci = ok[rng.weighted_index([len(seqs[i]) - length - 2 * margin + 1 for i in ok])]
			st = margin + rng.below(len(seqs[ci]) - length - 2 * margin + 1)
			if all(not (c == ci and s < st + length + margin and st - margin < e) for c, s, e in occupied):
				occupied.append((ci, st, st + length))
				return ci, st
		raise ValueError("cannot place a %d bp feature without overlaps; genome too small" % length)

	for fam in range(repeat_families):
		unit = rng.dna(rng.randint(repeat_len[0], repeat_len[1]))
		for k in range(repeat_copies):
			ci, st = place(len(unit))
			strand = "+" if rng.below(2) == 0 else "-"
			copy = mutate(rng, unit, repeat_divergence)
			if strand == "-":
				copy = revcomp(copy)
			if soft_mask_repeats:
				copy = copy.lower()
			seqs[ci][st:st + len(copy)] = list(copy)
			feats.append(Feature(names[ci], st, st + len(copy), "repeat", (fam, k, strand)))
	for inv in range(inversions):
		unit = rng.dna(rng.randint(inversion_len[0], inversion_len[1]))
		gap = rng.randint(inversion_gap[0], inversion_gap[1])
		ci, st = place(2 * len(unit) + gap)
		seqs[ci][st:st + len(unit)] = list(unit)
		st2 = st + len(unit) + gap
		seqs[ci][st2:st2 + len(unit)] = list(revcomp(unit))
		feats.append(Feature(names[ci], st, st + len(unit), "inversion", (inv, 0, "+")))
		feats.append(Feature(names[ci], st2, st2 + len(unit), "inversion", (inv, 1, "-")))
	for g in range(genes):
		ne = rng.randint(exons[0], exons[1])
		el = [rng.randint(exon_len[0], exon_len[1]) for _ in range(ne)]
		il = [rng.randint(intron_len[0], intron_len[1]) for _ in range(ne - 1)]
		span = sum(el) + sum(il)
		ci, st = place(span)
		strand = "+" if rng.below(2) == 0 else "-"
		ex, p = [], st
		for i in range(ne):
			ex.append((p, p + el[i]))
			p += el[i]
			if i < ne - 1:
				s, e = p, p + il[i]
				donor, acceptor = ("GT", "AG") if strand == "+" else ("CT", "AC")
				seqs[ci][s:s + 2] = list(donor)
				seqs[ci][e - 2:e] = list(acceptor)
				p = e
		gene_list.append(Gene("gene%d" % (g + 1), names[ci], strand, ex))
	for _ in range(low_complexity):
		length = rng.randint(low_complexity_len[0], low_complexity_len[1])
		period = rng.randint(1, 3)
		unit = rng.dna(period)
		if period > 1 and len(set(unit)) == 1:
			unit = unit[0] + _other_base(rng, unit[0]) + unit[2:]
		tract = (unit * (length // period + 1))[:length]
		ci, st = place(length)
		seqs[ci][st:st + length] = list(tract)
		feats.append(Feature(names[ci], st, st + length, "lowcomp", unit))
	for _ in range(n_runs):
		length = rng.randint(n_run_len[0], n_run_len[1])
		ci, st = place(length)
		seqs[ci][st:st + length] = ["N"] * length
		feats.append(Feature(names[ci], st, st + length, "N"))
	contigs = [(names[i], "".join(s)) for i, s in enumerate(seqs)]
	for k in range(short_contigs):
		contigs.append(("short%d" % (k + 1), rng.dna(rng.randint(short_contig_len[0], short_contig_len[1]))))
	feats.sort(key=lambda f: (names.index(f.contig), f.start))
	gene_list.sort(key=lambda g: (names.index(g.contig), g.start))
	for i, g in enumerate(gene_list):
		g.name = "gene%d" % (i + 1)
	return Genome(contigs, feats, gene_list)


# ---------------------------------------------------------------------------
# Read simulation
# ---------------------------------------------------------------------------


class Segment(object):
	"""Where the read bases [qstart, qend) (read orientation) came from.

	contig, pos (0-based leftmost reference base), end (exclusive) and
	strand ('+' if these read bases equal the reference forward strand) give
	the true placement; cigar is in reference orientation without clipping
	and may hold N for introns; edits is a list of (ref_pos, op, ref_bases,
	read_bases) in reference orientation with op 'X' (substitution), 'I'
	(insertion before ref_pos) or 'D'.
	"""
	__slots__ = ("contig", "pos", "end", "strand", "qstart", "qend", "cigar", "edits")

	def __init__(self, contig, pos, end, strand, qstart, qend, cigar, edits):
		self.contig, self.pos, self.end, self.strand = contig, pos, end, strand
		self.qstart, self.qend, self.cigar, self.edits = qstart, qend, cigar, edits

	@property
	def blocks(self):
		"""Reference intervals between introns (N operations), ascending."""
		out, s, t = [], self.pos, self.pos
		for n, op in parse_cigar(self.cigar):
			if op == "N":
				out.append((s, t))
				t += n
				s = t
			elif op in "MDX=":
				t += n
		out.append((s, t))
		return out

	@property
	def introns(self):
		b = self.blocks
		return [(b[i][1], b[i + 1][0]) for i in range(len(b) - 1)]

	def __repr__(self):
		return "Segment(%s:%d-%d%s q=%d-%d %s)" % (self.contig, self.pos, self.end, self.strand, self.qstart, self.qend,
			self.cigar if len(self.cigar) < 60 else self.cigar[:57] + "...")


class SimRead(object):
	"""A simulated read and its truth.

	seq/qual are as sequenced (5'->3'). segments lists where the read came
	from in read order: one Segment for an ordinary read, several for a
	chimeric read or an assembly contig with inversions. contig, pos, end,
	strand, cigar and edits are those of the first segment. Spliced reads
	also carry gene (the Gene) and ts_strand ('+' if the read is in the
	transcript's sense orientation, as minimap2's ts:A tag reports it).
	"""

	def __init__(self, name, seq, qual, segments, comment=None, gene=None, ts_strand=None):
		self.name, self.seq, self.qual, self.segments, self.comment = name, seq, qual, segments, comment
		self.gene, self.ts_strand = gene, ts_strand

	@property
	def contig(self):
		return self.segments[0].contig

	@property
	def pos(self):
		return self.segments[0].pos

	@property
	def end(self):
		return self.segments[0].end

	@property
	def strand(self):
		return self.segments[0].strand

	@property
	def cigar(self):
		return self.segments[0].cigar

	@property
	def edits(self):
		return self.segments[0].edits

	@property
	def exons(self):
		"""Reference blocks of a single-segment read between its introns, ascending."""
		return self.segments[0].blocks

	@property
	def introns(self):
		return self.segments[0].introns

	@property
	def is_chimeric(self):
		return len(self.segments) > 1

	@property
	def fwd_seq(self):
		"""The read in the first segment's reference orientation (SAM SEQ of a single-segment read)."""
		return self.seq if self.strand == "+" else revcomp(self.seq)

	@property
	def nm(self):
		"""Number of simulated edits: substitutions plus indel lengths.

		Bases copied from a reference N run are not edits, so a SAM NM (which
		counts every N as a mismatch) can exceed this for such reads.
		"""
		n = 0
		for seg in self.segments:
			for _, op, rb, qb in seg.edits:
				n += 1 if op == "X" else len(rb) + len(qb)
		return n

	def __len__(self):
		return len(self.seq)

	def __repr__(self):
		return "SimRead(%s, %d bp, %r)" % (self.name, len(self.seq), self.segments)


class SimPair(object):
	"""A read pair. For FR pairs, contig/frag_start/frag_end/strand give the fragment."""

	def __init__(self, name, r1, r2, contig=None, frag_start=None, frag_end=None, strand=None):
		self.name, self.r1, self.r2 = name, r1, r2
		self.contig, self.frag_start, self.frag_end, self.strand = contig, frag_start, frag_end, strand

	@property
	def insert(self):
		return None if self.frag_start is None else self.frag_end - self.frag_start

	def __iter__(self):
		return iter((self.r1, self.r2))

	def __repr__(self):
		return "SimPair(%s, %r, %r)" % (self.name, self.r1, self.r2)


# Error profiles for simulate_long_reads and friends. Rates are per base;
# indel lengths are uniform in 1..max_indel.
PROFILES = {
	"short": dict(length=150, sub_rate=0.005, indel_rate=0.0005, max_indel=3),
	"ont": dict(length=(2000, 8000), sub_rate=0.03, indel_rate=0.02, max_indel=3),  # about 7% error
	"hifi": dict(length=(3000, 10000), sub_rate=0.002, indel_rate=0.001, max_indel=2),  # about 0.4% error
	"clr": dict(length=(2000, 8000), sub_rate=0.01, indel_rate=0.04, max_indel=2),  # about 10%, indel-rich
	"asm": dict(length=(20000, 50000), sub_rate=0.001, indel_rate=0.0002, max_indel=10),
}


def _cigar_str(ops):
	merged = []
	for n, op in ops:
		if n <= 0:
			continue
		if merged and merged[-1][1] == op:
			merged[-1][0] += n
		else:
			merged.append([n, op])
	return "".join("%d%s" % (n, op) for n, op in merged)


def _sim_fwd(rng, ref, start, length, sub_rate, indel_rate, max_indel, end_margin, n_rate):
	"""Walk the reference from start emitting `length` read bases with edits.

	Returns (fwd_seq, end, cigar, edits) or None if the walk ran off the contig.
	Indels stay at least end_margin bases from either read end and never
	touch another indel.
	"""
	out, ops, edits = [], [], []
	i = start
	last_indel = -10
	while len(out) < length:
		if i >= len(ref):
			return None
		k = len(out)
		if indel_rate > 0 and end_margin <= k < length - end_margin - max_indel and k - last_indel > 1 and rng.chance(indel_rate):
			n = rng.randint(1, max_indel)
			if rng.below(2) == 0:  # insertion before ref[i]
				ins = rng.dna(n)
				out.extend(ins)
				ops.append((n, "I"))
				edits.append((i, "I", "", ins))
			else:
				if i + n >= len(ref):
					return None
				ops.append((n, "D"))
				edits.append((i, "D", ref[i:i + n].upper(), ""))
				i += n
			last_indel = len(out)
			continue
		c = ref[i].upper()
		if c in "ACGT" and rng.chance(sub_rate):
			b = _other_base(rng, c)
			edits.append((i, "X", c, b))
			out.append(b)
		elif c in "ACGT" and rng.chance(n_rate):
			edits.append((i, "X", c, "N"))
			out.append("N")
		else:
			out.append(c if c in "ACGT" else "N")
		ops.append((1, "M"))
		i += 1
	return "".join(out), i, _cigar_str(ops), edits


def _sim_blocks(rng, ref, blocks, sub_rate, indel_rate, max_indel, end_margin, junction_margin, n_rate=0.0):
	"""Copy the reference blocks [(s, e), ...] (ascending) into one read with edits and N ops between blocks.

	Indels stay end_margin bases from the read ends and junction_margin bases
	from every block boundary. Returns (fwd_seq, cigar, edits).
	"""
	total = sum(e - s for s, e in blocks)
	out, ops, edits = [], [], []
	done = 0  # reference bases of the blocks consumed so far
	last_indel = -10
	for bi, (s, e) in enumerate(blocks):
		if bi > 0:
			ops.append((s - blocks[bi - 1][1], "N"))
		i = s
		while i < e:
			k_ref = done + (i - s)
			ok = (indel_rate > 0 and end_margin <= k_ref < total - end_margin - max_indel and
				i - s >= junction_margin and e - i > junction_margin + max_indel and len(out) - last_indel > 1)
			if ok and rng.chance(indel_rate):
				n = rng.randint(1, max_indel)
				if rng.below(2) == 0:
					ins = rng.dna(n)
					out.extend(ins)
					ops.append((n, "I"))
					edits.append((i, "I", "", ins))
				else:
					ops.append((n, "D"))
					edits.append((i, "D", ref[i:i + n].upper(), ""))
					i += n
				last_indel = len(out)
				continue
			c = ref[i].upper()
			if c in "ACGT" and rng.chance(sub_rate):
				b = _other_base(rng, c)
				edits.append((i, "X", c, b))
				out.append(b)
			elif c in "ACGT" and rng.chance(n_rate):
				edits.append((i, "X", c, "N"))
				out.append("N")
			else:
				out.append(c if c in "ACGT" else "N")
			ops.append((1, "M"))
			i += 1
		done += e - s
	return "".join(out), _cigar_str(ops), edits


def _quals(rng, n, qual):
	if qual is None:
		return None
	if qual == "random":
		return "".join(chr(33 + 2 + rng.below(39)) for _ in range(n))
	if len(qual) == 1:
		return qual * n
	raise ValueError("qual must be None, 'random' or one character")


def _pick_contig(rng, genome, need, contigs):
	names = [n for n in (contigs or genome.names) if genome.length(n) >= need]
	if not names:
		raise ValueError("no contig of at least %d bp" % need)
	return names[rng.weighted_index([genome.length(n) - need + 1 for n in names])]


def _draw_len(rng, length):
	if isinstance(length, (tuple, list)):
		return rng.randint(length[0], length[1])
	return int(length)


def _sim_one(rng, genome, length, sub_rate, indel_rate, max_indel, end_margin, n_rate, unique, contigs,
		strand=None, contig=None, start=None, unique_pad=0):
	"""One single-segment read: returns (fwd_seq, Segment) in read orientation."""
	slack = length // 10 + 2 * max_indel + 10
	for _ in range(10000):
		c = contig or _pick_contig(rng, genome, length + slack, contigs)
		ref = genome[c]
		st = start if start is not None else rng.below(len(ref) - length - slack + 1)
		r = _sim_fwd(rng, ref, st, length, sub_rate, indel_rate, max_indel, end_margin, n_rate)
		if r is None:
			if start is not None:
				raise ValueError("read at %s:%d runs off the contig" % (c, st))
			continue
		seq, en, cigar, edits = r
		if unique and not genome.is_unique(c, st, en, unique_pad):
			if start is not None:
				raise ValueError("read at %s:%d is not in unique sequence" % (c, st))
			continue
		s = strand if strand is not None else ("+" if rng.below(2) == 0 else "-")
		return seq, Segment(c, st, en, s, 0, length, cigar, edits)
	raise ValueError("could not place a read; genome too small or too repetitive")


def simulate_reads(genome, n, length=150, seed=1, sub_rate=0.01, indel_rate=0.0005, max_indel=3, end_margin=10,
		n_rate=0.0, unique=False, unique_pad=0, contigs=None, strand=None, qual="random", name_prefix="r"):
	"""n single-end reads with recorded truth (see SimRead).

	length is an int or a (min, max) range. sub_rate and indel_rate are per
	base; indel lengths are uniform in 1..max_indel. With unique=True, reads
	come only from sequence overlapping no genome feature (within unique_pad).
	strand None picks each read's strand at random. qual is 'random', a single
	character, or None for FASTA.
	"""
	rng = Rng(seed)
	out = []
	for i in range(n):
		L = _draw_len(rng, length)
		fwd, seg = _sim_one(rng, genome, L, sub_rate, indel_rate, max_indel, end_margin, n_rate, unique, contigs,
			strand=strand, unique_pad=unique_pad)
		seq = fwd if seg.strand == "+" else revcomp(fwd)
		out.append(SimRead("%s%d" % (name_prefix, i), seq, _quals(rng, len(seq), qual), [seg]))
	return out


def simulate_long_reads(genome, n, profile="ont", seed=1, **kw):
	"""Long reads with an error profile from PROFILES ('ont' about 7% error, 'hifi' about 0.4%,
	'clr', 'asm'); any simulate_reads argument overrides the profile's. Default end_margin 50."""
	args = dict(PROFILES[profile])
	args.update(end_margin=50, name_prefix=profile)
	args.update(kw)
	return simulate_reads(genome, n, seed=seed, **args)


def simulate_pairs(genome, n, read_length=150, insert_mean=300, insert_sd=30, seed=1, sub_rate=0.01,
		indel_rate=0.0005, max_indel=3, end_margin=10, n_rate=0.0, unique=False, unique_pad=0, contigs=None,
		qual="random", name_prefix="p", min_insert=None):
	"""n FR read pairs. Inserts are ~normal(insert_mean, insert_sd), at least min_insert (default read_length,
	or its lower bound for a (min, max) range of read lengths) and at least the longer mate.

	For a '+' fragment, r1 is forward from frag_start and r2 is reverse
	ending near frag_end; for a '-' fragment the roles swap. frag_end is the
	end of the reverse read's true reference span, so frag_end - frag_start is
	the template length a correct aligner reports for an unclipped pair.
	Both reads carry the pair's name (no /1 /2; the writers can add them).
	"""
	rng = Rng(seed)
	if min_insert is not None:
		lo = min_insert
	else:
		lo = read_length[0] if isinstance(read_length, (tuple, list)) else read_length
	out = []
	for i in range(n):
		L1 = _draw_len(rng, read_length)
		L2 = _draw_len(rng, read_length)
		for _ in range(10000):
			ins = max(lo, max(L1, L2), rng.normal_int(insert_mean, insert_sd))
			slack = ins // 10 + 2 * max_indel + 10
			c = _pick_contig(rng, genome, ins + slack, contigs)
			fs = rng.below(genome.length(c) - ins - slack + 1)
			fstrand = "+" if rng.below(2) == 0 else "-"
			ref = genome[c]
			a = _sim_fwd(rng, ref, fs, L1 if fstrand == "+" else L2, sub_rate, indel_rate, max_indel, end_margin, n_rate)
			b_start = fs + ins - (L2 if fstrand == "+" else L1)
			b = _sim_fwd(rng, ref, b_start, L2 if fstrand == "+" else L1, sub_rate, indel_rate, max_indel, end_margin, n_rate)
			if a is None or b is None:
				continue
			if unique and not (genome.is_unique(c, fs, a[1], unique_pad) and genome.is_unique(c, b_start, b[1], unique_pad)):
				continue
			left = Segment(c, fs, a[1], "+", 0, len(a[0]), a[2], a[3])
			right = Segment(c, b_start, b[1], "-", 0, len(b[0]), b[2], b[3])
			lseq, rseq = a[0], revcomp(b[0])
			name = "%s%d" % (name_prefix, i)
			if fstrand == "+":
				r1 = SimRead(name, lseq, _quals(rng, len(lseq), qual), [left])
				r2 = SimRead(name, rseq, _quals(rng, len(rseq), qual), [right])
			else:
				r1 = SimRead(name, rseq, _quals(rng, len(rseq), qual), [right])
				r2 = SimRead(name, lseq, _quals(rng, len(lseq), qual), [left])
			out.append(SimPair(name, r1, r2, c, fs, b[1], fstrand))
			break
		else:
			raise ValueError("could not place a pair; genome too small")
	return out


def _join_segments(rng, genome, pieces, err, end_margin, unique, unique_pad, contigs):
	"""Simulate each piece (length, strand, contig or None, start or None) and join them in read order.

	Returns (seq, [Segment]) with qstart/qend set, or None if a piece failed.
	"""
	seq, segs = "", []
	for L, strand, c, st in pieces:
		try:
			fwd, seg = _sim_one(rng, genome, L, err["sub_rate"], err["indel_rate"], err["max_indel"], end_margin, 0.0,
				unique, contigs, strand=strand, contig=c, start=st, unique_pad=unique_pad)
		except ValueError:
			return None
		part = fwd if seg.strand == "+" else revcomp(fwd)
		seg.qstart, seg.qend = len(seq), len(seq) + len(part)
		seq += part
		segs.append(seg)
	return seq, segs


def simulate_chimeric_reads(genome, n, profile="hifi", length=(2000, 6000), parts=2, min_part=500, seed=1,
		unique=True, unique_pad=0, contigs=None, qual="random", name_prefix="chim", min_distance=10000, **err):
	"""Chimeric reads: `parts` pieces from independent loci (different contigs, or the same contig at
	least min_distance apart), each on a random strand, concatenated. SimRead.segments records each
	piece with its read interval [qstart, qend). Error rates come from PROFILES[profile] unless given."""
	e = dict(PROFILES[profile])
	e.update(err)
	rng = Rng(seed)
	out = []
	for i in range(n):
		for _ in range(1000):
			L = _draw_len(rng, length)
			cuts = sorted(rng.randint(min_part, L - min_part) for _ in range(parts - 1))
			lens = [b - a for a, b in zip([0] + cuts, cuts + [L])]
			if min(lens) < min_part:
				continue
			r = _join_segments(rng, genome, [(x, None, None, None) for x in lens], e, 30, unique, unique_pad, contigs)
			if r is None:
				continue
			seq, segs = r
			far = all(a.contig != b.contig or abs(a.pos - b.pos) >= min_distance
				for j, a in enumerate(segs) for b in segs[j + 1:])
			if far:
				out.append(SimRead("%s%d" % (name_prefix, i), seq, _quals(rng, len(seq), qual), segs))
				break
		else:
			raise ValueError("could not place a chimeric read; genome too small")
	return out


def simulate_contigs(genome, n, length=(20000, 50000), seed=1, inversions=0, inversion_len=(1000, 3000),
		qual=None, name_prefix="ctg", contigs=None, unique=False, **err):
	"""Assembly-like contigs: long, low-divergence copies of the genome (PROFILES['asm'] rates unless
	given). With inversions=k, each contig carries k inverted pieces of `inversion_len` bp: the
	segments are adjacent on one reference contig and alternate strands (+, -, +, ...)."""
	e = dict(PROFILES["asm"])
	e.update(err)
	rng = Rng(seed)
	out = []
	for i in range(n):
		for _ in range(1000):
			L = _draw_len(rng, length)
			slack = L // 10 + 100
			c = _pick_contig(rng, genome, L + slack, contigs)
			st = rng.below(genome.length(c) - L - slack + 1)
			strand = "+" if rng.below(2) == 0 else "-"
			if inversions == 0:
				pieces_ref = [(st, L)]
			else:
				inv = [_draw_len(rng, inversion_len) for _ in range(inversions)]
				flank = (L - sum(inv)) // (inversions + 1)
				if flank < 1000:
					raise ValueError("contig too short for %d inversions" % inversions)
				pieces_ref, p = [], st
				for k in range(inversions):
					pieces_ref.append((p, flank))
					p += flank
					pieces_ref.append((p, inv[k]))
					p += inv[k]
				pieces_ref.append((p, L - sum(x[1] for x in pieces_ref)))
			# exact reference blocks per piece, so pieces abut on the reference
			segs, parts = [], []
			ok = True
			for k, (p, plen) in enumerate(pieces_ref):
				s = "+" if k % 2 == 0 else "-"
				if strand == "-":
					s = "-" if s == "+" else "+"
				ref = genome[c]
				if p + plen > len(ref) or (unique and not genome.is_unique(c, p, p + plen)):
					ok = False
					break
				fwd, cig, edits = _sim_blocks(rng, ref, [(p, p + plen)], e["sub_rate"], e["indel_rate"], e["max_indel"],
					30, 30)
				parts.append((fwd, Segment(c, p, p + plen, s, 0, 0, cig, edits)))
			if not ok:
				continue
			if strand == "-":
				parts = parts[::-1]
			seq = ""
			for fwd, seg in parts:
				part = fwd if seg.strand == "+" else revcomp(fwd)
				seg.qstart, seg.qend = len(seq), len(seq) + len(part)
				seq += part
				segs.append(seg)
			out.append(SimRead("%s%d" % (name_prefix, i), seq, _quals(rng, len(seq), qual), segs))
			break
		else:
			raise ValueError("could not place a contig; genome too small")
	return out


def transcript_blocks(gene, a, b, min_end_exon=0):
	"""Genome blocks (ascending) of the transcript interval [a, b), transcript coordinates 5'->3'.

	First or last pieces shorter than min_end_exon are dropped (an aligner
	cannot place a few bases across an intron), as long as one block remains.
	"""
	tl = gene.length
	if gene.strand == "-":
		a, b = tl - b, tl - a  # the same interval in genome-forward exon order
	blocks, off = [], 0
	for s, en in gene.exons:
		lo, hi = max(s, s + a - off), min(en, s + b - off)
		if lo < hi:
			blocks.append((lo, hi))
		off += en - s
	while len(blocks) > 1 and blocks[0][1] - blocks[0][0] < min_end_exon:
		blocks = blocks[1:]
	while len(blocks) > 1 and blocks[-1][1] - blocks[-1][0] < min_end_exon:
		blocks = blocks[:-1]
	return blocks


def _spliced_read(rng, genome, g, blocks, sense, e, end_margin, junction_margin, qual, name):
	fwd, cig, edits = _sim_blocks(rng, genome[g.contig], blocks, e["sub_rate"], e["indel_rate"], e["max_indel"],
		end_margin, junction_margin)
	read_strand = g.strand if sense else ("-" if g.strand == "+" else "+")
	seq = fwd if read_strand == "+" else revcomp(fwd)
	seg = Segment(g.contig, blocks[0][0], blocks[-1][1], read_strand, 0, len(seq), cig, edits)
	return SimRead(name, seq, _quals(rng, len(seq), qual), [seg], gene=g, ts_strand="+" if sense else "-")


def simulate_spliced_reads(genome, n, seed=1, profile="hifi", full_length_frac=0.5, min_len=300, antisense_frac=0.5,
		junction_margin=15, end_margin=20, min_end_exon=30, length=None, genes=None, qual="random", name_prefix="tx",
		**err):
	"""cDNA reads from genome.genes (or the given genes) with recorded exons.

	By default a fraction full_length_frac of reads covers the whole
	transcript and the others a random sub-interval of at least min_len
	transcript bases. With length (an int or a (min, max) range), every read
	is a random sub-interval of that many transcript bases instead (genes
	shorter than that are not used), e.g. length=100 for splice:sr. First or
	last exon pieces shorter than min_end_exon bases are dropped. A fraction
	antisense_frac is the reverse complement of the transcript (ts_strand
	'-'), the rest sense (ts_strand '+'). Errors follow PROFILES[profile]
	unless given; indels stay junction_margin bases from every exon boundary.
	Truth: SimRead.exons (genome blocks, ascending), .introns, .gene,
	.ts_strand, and .strand (read strand on the genome).
	"""
	genes = genes if genes is not None else genome.genes
	if not genes:
		raise ValueError("the genome has no genes; use make_genome(genes=...)")
	e = dict(PROFILES[profile])
	e.update(err)
	rng = Rng(seed)
	out = []
	for i in range(n):
		if length is None:
			g = genes[rng.below(len(genes))]
			tl = g.length
			if rng.chance(full_length_frac) or tl <= min_len:
				a, b = 0, tl
			else:
				L = rng.randint(min_len, tl)
				a = rng.below(tl - L + 1)
				b = a + L
		else:
			L = _draw_len(rng, length)
			ok = [x for x in genes if x.length >= L]
			if not ok:
				raise ValueError("no gene has a %d bp transcript" % L)
			g = ok[rng.below(len(ok))]
			a = rng.below(g.length - L + 1)
			b = a + L
		blocks = transcript_blocks(g, a, b, min_end_exon)
		fwd, cig, edits = _sim_blocks(rng, genome[g.contig], blocks, e["sub_rate"], e["indel_rate"], e["max_indel"],
			end_margin, junction_margin)
		sense = not rng.chance(antisense_frac)
		read_strand = g.strand if sense else ("-" if g.strand == "+" else "+")
		seq = fwd if read_strand == "+" else revcomp(fwd)
		seg = Segment(g.contig, blocks[0][0], blocks[-1][1], read_strand, 0, len(seq), cig, edits)
		out.append(SimRead("%s%d" % (name_prefix, i), seq, _quals(rng, len(seq), qual), [seg], gene=g,
			ts_strand="+" if sense else "-"))
	return out


def simulate_spliced_pairs(genome, n, read_length=100, insert_mean=250, insert_sd=30, seed=1, profile="short",
		stranded=False, junction_margin=10, end_margin=10, min_end_exon=10, genes=None, qual="random", name_prefix="sp",
		**err):
	"""Paired RNA-seq reads (FR on the transcript) with recorded exons, for splice:sr.

	Fragments are ~normal(insert_mean, insert_sd) transcript bases (at least
	read_length, at most the transcript). Read 1 is the fragment's 5' end in
	transcript sense, read 2 its 3' end in antisense; with stranded=False
	the two swap at random. Each mate is a SimRead with .exons, .introns,
	.gene, .ts_strand and .strand like simulate_spliced_reads; mates lose
	end pieces shorter than min_end_exon, so they can be shorter than
	read_length. SimPair.contig is the gene's contig; frag_start/frag_end
	are None (a spliced fragment has no single reference span).
	"""
	genes = genes if genes is not None else genome.genes
	ok = [g for g in genes if g.length >= read_length]
	if not ok:
		raise ValueError("no gene has a transcript of %d bp" % read_length)
	e = dict(PROFILES[profile])
	e.update(err)
	rng = Rng(seed)
	out = []
	for i in range(n):
		g = ok[rng.below(len(ok))]
		tl = g.length
		F = min(tl, max(read_length, rng.normal_int(insert_mean, insert_sd)))
		a = rng.below(tl - F + 1)
		name = "%s%d" % (name_prefix, i)
		left = _spliced_read(rng, genome, g, transcript_blocks(g, a, a + read_length, min_end_exon), True, e,
			end_margin, junction_margin, qual, name)
		right = _spliced_read(rng, genome, g, transcript_blocks(g, a + F - read_length, a + F, min_end_exon), False, e,
			end_margin, junction_margin, qual, name)
		if not stranded and rng.below(2) == 1:
			left, right = right, left
		out.append(SimPair(name, left, right, g.contig))
	return out


# ---------------------------------------------------------------------------
# FASTA/FASTQ and annotation I/O
# ---------------------------------------------------------------------------


def _open_out(path, gz):
	if gz is None:
		gz = path.endswith(".gz")
	if gz:
		raw = open(path, "wb")
		return gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0), raw
	return open(path, "wb"), None


def _write_text(path, text, gz=None):
	f, raw = _open_out(path, gz)
	try:
		f.write(text.encode())
	finally:
		f.close()
		if raw is not None:
			raw.close()


def write_text(path, text, gz=None):
	"""Write text to path (gzip with a fixed header if gz or the name ends in .gz); returns path."""
	_write_text(str(path), text, gz)
	return str(path)


def _as_record(r):
	"""(name, seq, qual, comment) from a SimRead, a tuple or a dict."""
	if isinstance(r, SimRead):
		return r.name, r.seq, r.qual, r.comment
	if isinstance(r, dict):
		return r["name"], r["seq"], r.get("qual"), r.get("comment")
	r = tuple(r)
	return (r + (None, None))[:4]


def fastx_text(records, fmt="fastq", suffix=None, line_width=None):
	"""FASTA or FASTQ text. suffix is appended to every name (e.g. '/1')."""
	out = []
	for r in records:
		name, seq, qual, comment = _as_record(r)
		if suffix:
			name += suffix
		head = name if not comment else "%s %s" % (name, comment)
		if fmt == "fasta":
			out.append(">%s\n" % head)
			if line_width:
				for i in range(0, len(seq), line_width):
					out.append(seq[i:i + line_width] + "\n")
				if not seq:
					out.append("\n")
			else:
				out.append(seq + "\n")
		elif fmt == "fastq":
			if qual is None:
				qual = "I" * len(seq)
			if len(qual) != len(seq):
				raise ValueError("quality length differs from sequence length for %s" % name)
			out.append("@%s\n%s\n+\n%s\n" % (head, seq, qual))
		else:
			raise ValueError("fmt must be 'fasta' or 'fastq'")
	return "".join(out)


def _fmt_from_path(path, fmt):
	if fmt:
		return fmt
	p = path[:-3] if path.endswith(".gz") else path
	return "fasta" if p.endswith((".fa", ".fasta", ".fna")) else "fastq"


def write_fastx(path, records, fmt=None, gz=None, suffix=None, line_width=None):
	"""Write records; the format and compression follow the extension unless given."""
	path = str(path)
	_write_text(path, fastx_text(records, _fmt_from_path(path, fmt), suffix, line_width), gz)
	return path


def write_fasta(path, records, gz=None, line_width=None):
	return write_fastx(path, records, "fasta", gz, line_width=line_width)


def write_fastq(path, records, gz=None):
	return write_fastx(path, records, "fastq", gz)


def write_pairs(path1, path2, pairs, fmt=None, gz=None, suffix=False):
	"""Read 1 and read 2 into two files; suffix=True appends /1 and /2 to names."""
	write_fastx(path1, [p.r1 for p in pairs], fmt, gz, "/1" if suffix else None)
	write_fastx(path2, [p.r2 for p in pairs], fmt, gz, "/2" if suffix else None)
	return str(path1), str(path2)


def write_interleaved(path, pairs, fmt=None, gz=None, suffix=False):
	"""Both mates of each pair, adjacent, in one file."""
	path = str(path)
	fmt = _fmt_from_path(path, fmt)
	parts = []
	for p in pairs:
		parts.append(fastx_text([p.r1], fmt, "/1" if suffix else None))
		parts.append(fastx_text([p.r2], fmt, "/2" if suffix else None))
	_write_text(path, "".join(parts), gz)
	return path


def genes_bed12(genes):
	"""BED12 text of gene models (thickStart/End = start/end, name = gene name, score 0)."""
	out = []
	for g in genes:
		sizes = ",".join(str(e - s) for s, e in g.exons) + ","
		starts = ",".join(str(s - g.start) for s, _ in g.exons) + ","
		out.append("%s\t%d\t%d\t%s\t0\t%s\t%d\t%d\t0\t%d\t%s\t%s\n" % (g.contig, g.start, g.end, g.name, g.strand,
			g.start, g.end, len(g.exons), sizes, starts))
	return "".join(out)


def genes_gtf(genes, source="sim"):
	"""GTF text (transcript and exon lines; gene_id = transcript_id = gene name), 1-based closed."""
	out = []
	for g in genes:
		attr = 'gene_id "%s"; transcript_id "%s";' % (g.name, g.name)
		out.append("%s\t%s\ttranscript\t%d\t%d\t.\t%s\t.\t%s\n" % (g.contig, source, g.start + 1, g.end, g.strand, attr))
		for s, e in g.exons:
			out.append("%s\t%s\texon\t%d\t%d\t.\t%s\t.\t%s\n" % (g.contig, source, s + 1, e, g.strand, attr))
	return "".join(out)


def junctions_bed(genes):
	"""One BED6 line per intron: contig, start, end (0-based half-open intron), gene name, 0, strand."""
	out = []
	for g in genes:
		for s, e in g.introns:
			out.append("%s\t%d\t%d\t%s\t0\t%s\n" % (g.contig, s, e, g.name, g.strand))
	return "".join(out)


def _open_in(path):
	with open(path, "rb") as f:
		magic = f.read(2)
	if magic == b"\x1f\x8b":
		return io.TextIOWrapper(gzip.open(path, "rb"), encoding="ascii", errors="surrogateescape")
	return open(path, "r", encoding="ascii", errors="surrogateescape")


def read_text(path):
	"""Contents of a plain or gzip text file."""
	with _open_in(str(path)) as f:
		return f.read()


def read_fastx(path):
	"""[(name, seq, qual, comment)] from FASTA or FASTQ, plain or gzip; qual None for FASTA."""
	out = []
	with _open_in(str(path)) as f:
		lines = f.read().split("\n")
	i = 0
	while i < len(lines):
		line = lines[i]
		if not line:
			i += 1
			continue
		if line[0] == ">":
			head = line[1:].split(None, 1)
			seq = []
			i += 1
			while i < len(lines) and not lines[i].startswith(">"):
				seq.append(lines[i].strip())
				i += 1
			out.append((head[0] if head else "", "".join(seq), None, head[1] if len(head) > 1 else None))
		elif line[0] == "@":
			head = line[1:].split(None, 1)
			seq = lines[i + 1]
			qual = lines[i + 3]
			out.append((head[0] if head else "", seq, qual, head[1] if len(head) > 1 else None))
			i += 4
		else:
			raise ValueError("%s: unexpected line %r" % (path, line[:50]))
	return out


def read_fasta(path):
	"""{name: sequence} from a FASTA file, plain or gzip, in file order (dict keeps order)."""
	return dict((r[0], r[1]) for r in read_fastx(path))


# ---------------------------------------------------------------------------
# Running processes
# ---------------------------------------------------------------------------


def _no_core():
	"""preexec_fn: no core dump from the child. The limit is 1 byte, not 0: when the kernel pipes
	cores to a crash handler (Linux core_pattern starting with '|', e.g. apport), it still runs the
	handler at a limit of 0, about a second per abort, but skips it at a limit of 1."""
	try:
		import resource
	except ImportError:
		return
	if hasattr(resource, "RLIMIT_CORE"):
		try:
			resource.setrlimit(resource.RLIMIT_CORE, (1, 1))
		except (ValueError, OSError):  # a hard limit of 0 cannot be raised to 1
			try:
				resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
			except (ValueError, OSError):
				pass


class Result(object):
	"""Outcome of one process: argv, returncode, stdout/stderr (str) and the raw bytes."""

	def __init__(self, argv, returncode, stdout_bytes, stderr_bytes):
		self.argv, self.returncode = argv, returncode
		self.stdout_bytes, self.stderr_bytes = stdout_bytes, stderr_bytes
		self.stdout = stdout_bytes.decode("utf-8", "surrogateescape")
		self.stderr = stderr_bytes.decode("utf-8", "surrogateescape")

	def check(self):
		if self.returncode != 0:
			raise AssertionError("%s exited with %d\nstderr:\n%s" % (" ".join(self.argv), self.returncode, self.stderr[-4000:]))
		return self

	def sam(self, strict=True):
		return parse_sam(self.stdout, strict)

	def paf(self):
		return parse_paf(self.stdout)

	@property
	def messages(self):
		"""stderr lines other than minimap2's '[M::...]' progress and timing lines."""
		return [l for l in self.stderr.splitlines() if l and not l.startswith("[M::")]

	def __repr__(self):
		return "Result(%r, returncode=%d)" % (" ".join(self.argv), self.returncode)


def run_process(argv, timeout=None, stdin=None, env=None, cwd=None, check=False):
	"""Run argv with LC_ALL=C, no core dumps and a timeout, in `cwd` or a fresh temporary directory."""
	argv = [str(a) for a in argv]
	if isinstance(stdin, str):
		stdin = stdin.encode()
	full_env = dict(os.environ)
	full_env["LC_ALL"] = "C"
	if env:
		full_env.update(env)
	kw = {}
	if os.name == "posix":
		kw["preexec_fn"] = _no_core
	tmp = None
	if cwd is None:
		tmp = cwd = tempfile.mkdtemp(prefix="mm2-run-")
	try:
		p = subprocess.run(argv, input=stdin, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
			timeout=timeout or run_timeout(), env=full_env, cwd=str(cwd), **kw)
	except subprocess.TimeoutExpired as e:
		raise AssertionError("%s timed out after %s s" % (" ".join(argv), e.timeout))
	finally:
		if tmp is not None:
			shutil.rmtree(tmp, ignore_errors=True)
	r = Result(argv, p.returncode, p.stdout, p.stderr)
	return r.check() if check else r


def run_mm2(args, binary=None, timeout=None, stdin=None, env=None, cwd=None, check=False):
	"""Run `minimap2 <args>`; stdin may be bytes, str or None. Core dumps are disabled.

	The process runs in `cwd`, or by default in a new temporary directory that is
	removed afterwards, so a relative output path (or a binary that writes a stray
	file) can never land in the source tree. Give paths in `args` as absolute paths."""
	return run_process([binary or minimap2_build.minimap2_bin()] + list(args), timeout, stdin, env, cwd, check)


def _has_thread_opt(args):
	for a in args:
		a = str(a)
		if a == "-t" or (a.startswith("-t") and a[2:].isdigit()):
			return True
	return False


def _abs(path):
	path = str(path)
	return path if path == "-" else os.path.abspath(path)


def map_reads(target, reads=(), args=(), preset=None, threads=None, binary=None, timeout=None, stdin=None, env=None,
		check=True, cwd=None):
	"""Run `minimap2 [-x preset] [args] -t N <target> <reads...>`.

	target is a FASTA or .mmi path; reads is a path or a list of paths ('-'
	reads stdin). The thread count is `threads` if given, else whatever
	`args` sets with -t N/-tN, else $MINIMAP2_TEST_THREADS. check=True raises
	on a non-zero exit. The process runs in `cwd`, by default a temporary
	directory; target and read paths are made absolute first.
	"""
	args = [str(a) for a in args]
	if isinstance(reads, (str, bytes)) or hasattr(reads, "__fspath__"):
		reads = [reads]
	reads = [_abs(r) for r in reads]
	if threads is not None:
		if _has_thread_opt(args):
			raise ValueError("give the thread count either in args or as threads=, not both")
		args = args + ["-t", str(threads)]
	elif not _has_thread_opt(args):
		args = args + ["-t", str(test_threads())]
	pre = ["-x", preset] if preset else []
	return run_mm2(pre + args + [_abs(target)] + reads, binary=binary, timeout=timeout, stdin=stdin, env=env, cwd=cwd,
		check=check)


_INDEX_CACHE = {}


def index_genome(fa, tmpdir, preset=None, extra_args=(), threads=1, binary=None):
	"""Build `minimap2 [-x preset] [extra_args] -d <mmi> <fa>`; return the .mmi path.

	fa is a FASTA path or a Genome. Indexes are cached per (genome bytes,
	preset, extra args, binary) for the life of the Python process, so tmpdir
	should outlive the test session (use the session-scoped `index_dir`
	fixture).
	"""
	binary = binary or minimap2_build.minimap2_bin()
	tmpdir = str(tmpdir)
	if isinstance(fa, Genome):
		data = fa.fasta_text().encode()
		gh = hashlib.sha256(data).hexdigest()[:16]
		path = os.path.join(tmpdir, "genome-%s.fa" % gh)
		if not os.path.exists(path):
			with open(path, "wb") as f:
				f.write(data)
		fa = path
	else:
		fa = os.path.abspath(str(fa))
		with open(fa, "rb") as f:
			gh = hashlib.sha256(f.read()).hexdigest()[:16]
	key = (gh, preset, binary, tuple(str(a) for a in extra_args))
	if key in _INDEX_CACHE and os.path.exists(_INDEX_CACHE[key]):
		return _INDEX_CACHE[key]
	kh = hashlib.sha256(repr(key).encode()).hexdigest()[:12]
	mmi = os.path.join(tmpdir, "idx-%s-%s.mmi" % ((preset or "default").replace(":", "_"), kh))
	argv = (["-x", preset] if preset else []) + [str(a) for a in extra_args] + ["-t", str(threads), "-d", mmi, fa]
	run_mm2(argv, binary=binary, check=True)
	if not os.path.exists(mmi):
		raise AssertionError("minimap2 -d did not write %s" % mmi)
	_INDEX_CACHE[key] = mmi
	return mmi


# ---------------------------------------------------------------------------
# Reference minimizers and the .mmi format
# ---------------------------------------------------------------------------

# sketch.c's seq_nt4_table: A/C/G/T/U in either case -> 0..3, anything else -> 4
_SKETCH_NT4 = {"A": 0, "C": 1, "G": 2, "T": 3, "U": 3, "a": 0, "c": 1, "g": 2, "t": 3, "u": 3}


def kmer_hash(key, k):
	"""The invertible integer hash minimap2 applies to a 2k-bit k-mer (Li 2016, 2018; sketch.c),
	restricted to 2k bits: a bijection on [0, 4^k)."""
	mask = (1 << (2 * k)) - 1
	key = (~key + (key << 21)) & mask
	key = key ^ (key >> 24)
	key = ((key + (key << 3)) + (key << 8)) & mask
	key = key ^ (key >> 14)
	key = ((key + (key << 2)) + (key << 4)) & mask
	key = key ^ (key >> 28)
	key = (key + (key << 31)) & mask
	return key


class Kmer(tuple):
	"""(x, pos, strand): x = kmer_hash << 8 | span, pos the 0-based position of the k-mer's last base,
	strand 0 if the forward k-mer has the smaller 2-bit encoding, 1 otherwise."""
	__slots__ = ()

	def __new__(cls, x, pos, strand):
		return tuple.__new__(cls, (x, pos, strand))

	x = property(lambda self: self[0])
	pos = property(lambda self: self[1])
	strand = property(lambda self: self[2])

	@property
	def hash(self):
		return self[0] >> 8


KMER_SKIP = (1 << 64) - 1  # x of a k-mer that cannot be a minimizer (an HPC k-mer spanning 256 bases or more)


def kmer_runs(seq, k, hpc=False):
	"""The canonical k-mers of seq, one list of Kmer per maximal run of A/C/G/T/U (any other character
	ends a run). span is k, or with hpc (-H) the length in seq of the k homopolymer-compressed bases,
	and pos is then the last base of the k-mer's last homopolymer run. An HPC k-mer spanning 256 bases
	or more keeps its place in the windows but is never a minimizer (x = KMER_SKIP). Only odd k is
	modelled (no k-mer is its own reverse complement); even k raises ValueError."""
	if k % 2 == 0 or not 0 < k <= 28:
		raise ValueError("the reference models odd k in 1..28 only")
	mask = (1 << 2 * k) - 1
	shift = 2 * (k - 1)
	runs, cur = [], None
	n, i = len(seq), 0
	fwd = rev = l = span = 0
	spans = []
	while i < n:
		c = _SKETCH_NT4.get(seq[i], 4)
		if c > 3:
			cur = None
			i += 1
			continue
		if cur is None:
			cur = []
			runs.append(cur)
			fwd = rev = l = span = 0
			spans = []
		if hpc:
			run = 1
			while i + run < n and _SKETCH_NT4.get(seq[i + run], 4) == c:
				run += 1
			i += run - 1
			spans.append(run)
			span += run
			if len(spans) > k:
				span -= spans.pop(0)
		else:
			span = k
		fwd = ((fwd << 2) | c) & mask
		rev = (rev >> 2) | ((3 - c) << shift)
		l += 1
		if l >= k:
			z = 0 if fwd < rev else 1
			cur.append(Kmer(kmer_hash(fwd if z == 0 else rev, k) << 8 | span if span < 256 else KMER_SKIP, i, z))
		i += 1
	return runs


def _window_minima(kms, w):
	"""[frozenset of the Kmers with the smallest x] for each window of w consecutive k-mers of kms;
	a run with 1..w-1 k-mers has one window holding all of them."""
	if not kms:
		return []
	if len(kms) < w:
		m = min(q.x for q in kms)
		return [frozenset(q for q in kms if q.x == m and m != KMER_SKIP)]
	out = []
	dq = []  # indices of increasing x, a monotone queue
	head = 0
	where = {}
	for j, q in enumerate(kms):
		where.setdefault(q.x, []).append(j)
	for j, q in enumerate(kms):
		while len(dq) > head and kms[dq[-1]].x >= q.x:
			dq.pop()
		dq.append(j)
		if dq[head] <= j - w:
			head += 1
		if j >= w - 1:
			m = kms[dq[head]].x
			if m == KMER_SKIP:
				out.append(frozenset())
				continue
			out.append(frozenset(kms[t] for t in where[m] if j - w < t <= j))
	return out


class SketchRef(object):
	"""Reference (w,k)-minimizers of one sequence, by the rule of minimap2.1 ("a minimizer is the
	smallest k-mer in a window of w consecutive k-mers"; every k-mer with the smallest hash on a tie).

	runs[r] are the Kmers of the r-th A/C/G/T/U run and windows[r] the minimum-x sets of its windows
	of w consecutive k-mers; a run with fewer than w k-mers has one short window."""

	def __init__(self, seq, w, k, hpc=False):
		self.seq, self.w, self.k, self.hpc = seq, w, k, hpc
		self.runs = kmer_runs(seq, k, hpc)
		self.windows = [_window_minima(r, w) for r in self.runs]

	@property
	def has_ties(self):
		"""True if some window's minimum occurs more than once in it."""
		return any(len(s) > 1 for ws in self.windows for s in ws)

	def minimizers(self):
		"""Every window minimum of every run (Kmers), sorted by position."""
		u = set()
		for ws in self.windows:
			for s in ws:
				u |= s
		return sorted(u, key=lambda q: q.pos)

	def kmer_at(self):
		"""{pos: Kmer} for every k-mer of the sequence that can be a minimizer."""
		return dict((q.pos, q) for r in self.runs for q in r if q.x != KMER_SKIP)


def minimizers(seq, w, k, hpc=False, ties=True):
	"""The (w,k)-minimizers of seq as a sorted list of (end position, hash, strand) (SketchRef).
	ties=False raises ValueError when a window's minimum occurs twice in it (a tandem repeat), for
	tests that model only windows with a single minimum."""
	ref = SketchRef(seq, w, k, hpc)
	if not ties and ref.has_ties:
		raise ValueError("a window of the sequence has its minimum twice; ties are not modelled here")
	return [(q.pos, q.hash, q.strand) for q in ref.minimizers()]


def expected_index(seqs, w, k, hpc=False):
	"""{hash: sorted tuple of rid << 32 | pos << 1 | strand} that an index of seqs (a list of strings,
	in index order) must hold by the documented rule. Only for sequences without ambiguous bases and
	without tied windows (raises ValueError otherwise), where the rule has no edge cases."""
	out = {}
	for rid, s in enumerate(seqs):
		ref = SketchRef(s, w, k, hpc)
		if len(ref.runs) > 1 or ref.has_ties:
			raise ValueError("sequence %d has an ambiguous base or a tied window; the exact reference does not apply" % rid)
		for q in ref.minimizers():
			out.setdefault(q.hash, []).append(rid << 32 | q.pos << 1 | q.strand)
	return dict((h, tuple(sorted(v))) for h, v in out.items())


MMI_MAGIC = b"MMI\x02"  # minimap.h MM_IDX_MAGIC
MM_I_HPC, MM_I_NO_SEQ, MM_I_NO_NAME = 0x1, 0x2, 0x4  # minimap.h MM_I_*


class MmiPart(object):
	"""One index of a .mmi file (a multi-part index holds several back to back).

	w, k, b, flag: the header; names (None with MM_I_NO_NAME) and lengths of the n_seq sequences;
	minimizers: {hash: sorted tuple of positions rid << 32 | pos << 1 | strand}; seq4: the 4-bit base
	codes (0-3 ACGT, 4 other) per sequence as bytes, or None with --idx-no-seq; n_singletons: hashes
	stored once; offset, end: the part's byte range in the file."""

	def __init__(self):
		self.w = self.k = self.b = self.flag = 0
		self.names, self.lengths, self.minimizers, self.seq4 = [], [], {}, None
		self.n_singletons = 0
		self.offset = self.end = 0

	@property
	def n_seq(self):
		return len(self.names)

	def seq_codes(self, rid):
		return self.seq4[rid] if self.seq4 is not None else None


def _parse_mmi_part(data, off, header_only=False):
	import struct
	p = MmiPart()
	p.offset = off
	if data[off:off + 4] != MMI_MAGIC:
		raise ValueError("no MMI magic at byte %d" % off)
	off += 4
	p.w, p.k, p.b, n_seq, p.flag = struct.unpack_from("<5I", data, off)
	off += 20
	for _ in range(n_seq):
		l = data[off]
		off += 1
		p.names.append(data[off:off + l].decode("ascii", "surrogateescape") if l else None)
		off += l
		p.lengths.append(struct.unpack_from("<I", data, off)[0])
		off += 4
	if header_only:
		return p, off
	mins = {}
	for bi in range(1 << p.b):
		n, = struct.unpack_from("<i", data, off)
		off += 4
		pos = struct.unpack_from("<%dQ" % n, data, off) if n else ()
		off += 8 * n
		size, = struct.unpack_from("<I", data, off)
		off += 4
		used = 0
		for _ in range(size):
			key, val = struct.unpack_from("<2Q", data, off)
			off += 16
			h = (key >> 1) << p.b | bi
			if h in mins:
				raise ValueError("hash %d stored twice" % h)
			if key & 1:
				mins[h] = (val,)
				p.n_singletons += 1
			else:
				st, cnt = val >> 32, val & 0xffffffff
				if cnt < 2 or st + cnt > n:
					raise ValueError("bucket %d: bad position range %d+%d of %d" % (bi, st, cnt, n))
				mins[h] = tuple(pos[st:st + cnt])
				used += cnt
		if used != n:
			raise ValueError("bucket %d: %d positions stored, %d referenced" % (bi, n, used))
	p.minimizers = mins
	if not p.flag & MM_I_NO_SEQ:
		total = sum(p.lengths)
		nw = (total + 7) // 8
		words = struct.unpack_from("<%dI" % nw, data, off)
		off += 4 * nw
		codes = bytearray(total)
		for i in range(total):
			codes[i] = (words[i >> 3] >> ((i & 7) << 2)) & 0xf
		p.seq4, o = [], 0
		for l in p.lengths:
			p.seq4.append(bytes(codes[o:o + l]))
			o += l
	p.end = off
	return p, off


def read_mmi(path, header_only=False):
	"""All parts of a .mmi file in order (the layout of index.c's mm_idx_dump: MM_IDX_MAGIC, w, k, b,
	n_seq, flag, the names and lengths, per bucket the position array and the hash table, then the
	packed sequence unless MM_I_NO_SEQ). Raises ValueError on a malformed file or trailing bytes.
	header_only=True returns the first part with its header, names and lengths only."""
	with open(str(path), "rb") as f:
		data = f.read()
	if header_only:
		return _parse_mmi_part(data, 0, True)[0]
	parts, off = [], 0
	while off < len(data):
		p, off = _parse_mmi_part(data, off)
		parts.append(p)
	if not parts:
		raise ValueError("empty .mmi file")
	return parts


# ---------------------------------------------------------------------------
# SAM parsing
# ---------------------------------------------------------------------------


class SamFormatError(ValueError):
	pass


_TAG_RE = re.compile(r"^([A-Za-z][A-Za-z0-9]):([AifZHB]):(.*)$", re.S)
_INT_RE = re.compile(r"^[-+]?[0-9]+$")
_FLOAT_RE = re.compile(r"^[-+]?([0-9]*\.?[0-9]+|[0-9]+\.)([eE][-+]?[0-9]+)?$")
_CIGAR_RE = re.compile(r"^(\*|([0-9]+[MIDNSHPX=])+)$")
_B_SUB = {"c": (-128, 127), "C": (0, 255), "s": (-32768, 32767), "S": (0, 65535),
	"i": (-2 ** 31, 2 ** 31 - 1), "I": (0, 2 ** 32 - 1), "f": None}


def parse_tag(text):
	"""'XX:T:value' -> (name, type, typed value). B arrays give (subtype, [values])."""
	m = _TAG_RE.match(text)
	if not m:
		raise SamFormatError("malformed tag %r" % text)
	name, typ, val = m.groups()
	if typ == "i":
		if not _INT_RE.match(val):
			raise SamFormatError("bad integer in tag %r" % text)
		v = int(val)
	elif typ == "f":
		if not _FLOAT_RE.match(val) and val not in ("inf", "-inf", "nan", "NaN"):
			raise SamFormatError("bad float in tag %r" % text)
		v = float(val)
	elif typ == "A":
		if len(val) != 1 or not ("!" <= val <= "~"):
			raise SamFormatError("bad character in tag %r" % text)
		v = val
	elif typ == "Z":
		if any(not (" " <= c <= "~") for c in val):
			raise SamFormatError("non-printable character in tag %r" % text)
		v = val
	elif typ == "H":
		if len(val) % 2 or not re.match(r"^[0-9A-F]*$", val):
			raise SamFormatError("bad hex in tag %r" % text)
		v = val
	else:
		parts = val.split(",")
		sub = parts[0]
		if sub not in _B_SUB:
			raise SamFormatError("bad B subtype in tag %r" % text)
		vals = []
		for x in parts[1:]:
			if sub == "f":
				vals.append(float(x))
			else:
				if not _INT_RE.match(x):
					raise SamFormatError("bad B array value in tag %r" % text)
				iv = int(x)
				lo, hi = _B_SUB[sub]
				if not (lo <= iv <= hi):
					raise SamFormatError("B array value out of range in tag %r" % text)
				vals.append(iv)
		v = (sub, vals)
	return name, typ, v


def _parse_tags(fields, line, strict):
	tags, types, order, extra = {}, {}, [], []
	for t in fields:
		try:
			name, typ, v = parse_tag(t)
		except SamFormatError:
			if strict:
				raise
			extra.append(t)
			continue
		if name in tags:
			raise SamFormatError("duplicate tag %s in %r" % (name, line[:200]))
		tags[name], types[name] = v, typ
		order.append(name)
	return tags, types, order, extra


class SamRecord(object):
	"""One alignment line. pos is 1-based as in the file; pos0 is 0-based.

	extra holds optional fields that are not well-formed tags (only with
	strict=False, e.g. FASTA comments copied by -y).
	"""
	__slots__ = ("qname", "flag", "rname", "pos", "mapq", "cigar", "rnext", "pnext", "tlen", "seq", "qual",
		"tags", "tag_types", "tag_order", "extra", "line")

	def __init__(self, line, strict=True):
		self.line = line
		f = line.split("\t")
		if len(f) < 11:
			raise SamFormatError("fewer than 11 fields: %r" % line[:200])
		self.qname, flag, self.rname, pos, mapq, self.cigar, self.rnext, pnext, tlen, self.seq, self.qual = f[:11]
		try:
			self.flag, self.pos, self.mapq, self.pnext, self.tlen = int(flag), int(pos), int(mapq), int(pnext), int(tlen)
		except ValueError:
			raise SamFormatError("non-integer FLAG/POS/MAPQ/PNEXT/TLEN: %r" % line[:200])
		self.tags, self.tag_types, self.tag_order, self.extra = _parse_tags(f[11:], line, strict)
		if strict:
			self._validate()

	def _validate(self):
		if not re.match(r"^[!-?A-~]{1,254}$", self.qname) and self.qname != "*":
			raise SamFormatError("bad QNAME %r" % self.qname)
		if not (0 <= self.flag < 1 << 16):
			raise SamFormatError("FLAG out of range in %r" % self.line[:200])
		if not (0 <= self.pos < 1 << 31) or not (0 <= self.pnext < 1 << 31):
			raise SamFormatError("POS/PNEXT out of range in %r" % self.line[:200])
		if not (0 <= self.mapq <= 255):
			raise SamFormatError("MAPQ out of range in %r" % self.line[:200])
		if not _CIGAR_RE.match(self.cigar):
			raise SamFormatError("bad CIGAR %r" % self.cigar)
		if self.seq != "*" and not re.match(r"^[A-Za-z=.]+$", self.seq):
			raise SamFormatError("bad SEQ in %r" % self.line[:200])
		if self.qual != "*":
			if not re.match(r"^[!-~]+$", self.qual):
				raise SamFormatError("bad QUAL in %r" % self.line[:200])
			if self.seq != "*" and len(self.qual) != len(self.seq):
				raise SamFormatError("QUAL and SEQ lengths differ in %r" % self.line[:200])
		if self.cigar != "*" and self.seq != "*" and cigar_query_span(self.cigar) != len(self.seq):
			raise SamFormatError("CIGAR query length %d != SEQ length %d in %r" % (cigar_query_span(self.cigar), len(self.seq), self.line[:200]))

	def tag(self, name, default=None):
		return self.tags.get(name, default)

	@property
	def is_paired(self):
		return bool(self.flag & 0x1)

	@property
	def is_proper_pair(self):
		return bool(self.flag & 0x2)

	@property
	def is_unmapped(self):
		return bool(self.flag & 0x4)

	@property
	def mate_unmapped(self):
		return bool(self.flag & 0x8)

	@property
	def is_reverse(self):
		return bool(self.flag & 0x10)

	@property
	def mate_reverse(self):
		return bool(self.flag & 0x20)

	@property
	def is_read1(self):
		return bool(self.flag & 0x40)

	@property
	def is_read2(self):
		return bool(self.flag & 0x80)

	@property
	def is_secondary(self):
		return bool(self.flag & 0x100)

	@property
	def is_supplementary(self):
		return bool(self.flag & 0x800)

	@property
	def is_primary(self):
		return not (self.flag & 0x900)

	@property
	def segment(self):
		"""1 for read 1, 2 for read 2, 0 for an unpaired read."""
		return 1 if self.flag & 0x40 else 2 if self.flag & 0x80 else 0

	@property
	def strand(self):
		return "-" if self.flag & 0x10 else "+"

	@property
	def pos0(self):
		return self.pos - 1

	@property
	def cigar_ops(self):
		return parse_cigar(self.cigar)

	@property
	def ref_span(self):
		return cigar_ref_span(self.cigar)

	@property
	def ref_end0(self):
		"""0-based exclusive end on the reference."""
		return self.pos0 + self.ref_span

	@property
	def left_clip(self):
		"""Soft plus hard clipping at the leftmost (reference-orientation) end."""
		n = 0
		for k, op in self.cigar_ops:
			if op not in "SH":
				break
			n += k
		return n

	@property
	def right_clip(self):
		n = 0
		for k, op in reversed(self.cigar_ops):
			if op not in "SH":
				break
			n += k
		return n

	@property
	def unclipped_start0(self):
		return self.pos0 - self.left_clip

	@property
	def unclipped_end0(self):
		return self.ref_end0 + self.right_clip

	@property
	def query_length(self):
		"""Read length implied by the CIGAR, hard clips included."""
		return cigar_query_span(self.cigar, include_hard=True)

	@property
	def introns(self):
		"""Reference intervals of the CIGAR's N operations, 0-based half-open, ascending."""
		out, t = [], self.pos0
		for n, op in self.cigar_ops:
			if op == "N":
				out.append((t, t + n))
			if op in "MDN=X":
				t += n
		return out

	def __repr__(self):
		return "SamRecord(%s)" % self.line[:200]


class SamFile(object):
	"""Header lines and records. hd, sq, rg and pg are parsed views of the header."""

	def __init__(self, header, records):
		self.header, self.records = header, records
		self.hd, self.sq, self.rg, self.pg, self.co = {}, [], [], [], []
		for h in header:
			code = h[1:3]
			if code == "CO":
				self.co.append(h[4:])
				continue
			fields = {}
			for x in h.split("\t")[1:]:
				if len(x) < 3 or x[2] != ":":
					raise SamFormatError("malformed header field %r in %r" % (x, h))
				if x[:2] in fields:
					raise SamFormatError("duplicate header field %s in %r" % (x[:2], h))
				fields[x[:2]] = x[3:]
			if code == "HD":
				self.hd = fields
			elif code == "SQ":
				self.sq.append((fields["SN"], int(fields["LN"])))
			elif code == "RG":
				self.rg.append(fields)
			elif code == "PG":
				self.pg.append(fields)

	def __iter__(self):
		return iter(self.records)

	def __len__(self):
		return len(self.records)

	@property
	def contig_lengths(self):
		return dict(self.sq)

	def by_qname(self):
		"""{qname: [records in file order]}."""
		d = {}
		for r in self.records:
			d.setdefault(r.qname, []).append(r)
		return d

	def primaries(self):
		"""{(qname, segment): the primary record}; raises if a segment has two primaries."""
		d = {}
		for r in self.records:
			if r.is_primary:
				k = (r.qname, r.segment)
				if k in d:
					raise SamFormatError("two primary records for %s segment %d" % k)
				d[k] = r
		return d


def parse_sam(text, strict=True):
	"""Parse SAM text. strict=True also validates field syntax per the SAM spec."""
	header, records = [], []
	for line in text.split("\n"):
		if not line:
			continue
		if line[0] == "@":
			if records:
				raise SamFormatError("header line after records: %r" % line)
			if not re.match(r"^@[A-Za-z][A-Za-z](\t|$)", line):
				raise SamFormatError("malformed header line %r" % line)
			header.append(line)
		else:
			records.append(SamRecord(line, strict))
	if text and not text.endswith("\n"):
		raise SamFormatError("output does not end with a newline")
	return SamFile(header, records)


def strip_pg(sam_text):
	"""SAM text without @PG lines (they carry the command line and the version)."""
	return "".join(l for l in sam_text.splitlines(True) if not l.startswith("@PG\t"))


def sam_body(sam_text):
	"""SAM text without any header line."""
	return "".join(l for l in sam_text.splitlines(True) if not l.startswith("@"))


# ---------------------------------------------------------------------------
# PAF parsing
# ---------------------------------------------------------------------------


class PafRecord(object):
	"""One PAF line. Unmapped queries (--paf-no-hit) have strand and tname None.

	minimap2 writes an unmapped query as `name qlen 0 0 * * 0 0 0 0 0 0`.
	"""
	__slots__ = ("qname", "qlen", "qstart", "qend", "strand", "tname", "tlen", "tstart", "tend", "nmatch", "alen",
		"mapq", "tags", "tag_types", "tag_order", "extra", "line")

	def __init__(self, line, strict=True):
		self.line = line
		f = line.split("\t")
		if len(f) < 12:
			raise SamFormatError("PAF line with fewer than 12 fields: %r" % line[:200])
		self.qname = f[0]

		def num(x):
			if x == "*":
				return None
			if not _INT_RE.match(x):
				raise SamFormatError("non-integer PAF field %r in %r" % (x, line[:200]))
			return int(x)
		self.qlen, self.qstart, self.qend = num(f[1]), num(f[2]), num(f[3])
		if f[4] not in ("+", "-", "*"):
			raise SamFormatError("bad PAF strand %r in %r" % (f[4], line[:200]))
		self.strand = None if f[4] == "*" else f[4]
		self.tname = None if f[5] == "*" else f[5]
		self.tlen, self.tstart, self.tend = num(f[6]), num(f[7]), num(f[8])
		self.nmatch, self.alen, self.mapq = num(f[9]), num(f[10]), num(f[11])
		self.tags, self.tag_types, self.tag_order, self.extra = _parse_tags(f[12:], line, strict)

	@property
	def is_unmapped(self):
		return self.tname is None

	@property
	def is_primary(self):
		return self.tags.get("tp") == "P"

	@property
	def cigar(self):
		return self.tags.get("cg")

	def tag(self, name, default=None):
		return self.tags.get(name, default)

	def __repr__(self):
		return "PafRecord(%s)" % self.line[:200]


def parse_paf(text, strict=True):
	if text and not text.endswith("\n"):
		raise SamFormatError("output does not end with a newline")
	out = []
	for line in text.split("\n"):
		if line:
			out.append(PafRecord(line, strict))
	return out


# ---------------------------------------------------------------------------
# CIGAR and alignment-tag recomputation
# ---------------------------------------------------------------------------

_CIGAR_OP_RE = re.compile(r"([0-9]+)([MIDNSHPX=])")


def parse_cigar(cigar):
	"""'10M2I' -> [(10, 'M'), (2, 'I')]; '*' -> []."""
	if cigar in ("*", ""):
		return []
	ops = [(int(n), op) for n, op in _CIGAR_OP_RE.findall(cigar)]
	if "".join("%d%s" % x for x in ops) != cigar:
		raise SamFormatError("bad CIGAR %r" % cigar)
	return ops


def cigar_ref_span(cigar):
	"""Reference bases consumed: M, D, N, = and X."""
	return sum(n for n, op in parse_cigar(cigar) if op in "MDN=X")


def cigar_query_span(cigar, include_hard=False):
	"""Read bases in SEQ: M, I, S, = and X (plus H with include_hard)."""
	ops = "MIS=XH" if include_hard else "MIS=X"
	return sum(n for n, op in parse_cigar(cigar) if op in ops)


def cigar_clips(cigar):
	"""(leading clip, trailing clip) lengths, soft and hard together."""
	ops = parse_cigar(cigar)
	a = b = 0
	for n, op in ops:
		if op not in "SH":
			break
		a += n
	for n, op in reversed(ops):
		if op not in "SH":
			break
		b += n
	return a, b


def strip_clips(cigar):
	"""The CIGAR without S/H operations: what minimap2 writes as cg:Z in PAF."""
	return "".join("%d%s" % (n, op) for n, op in parse_cigar(cigar) if op not in "SH")


def _walk(ref, pos0, seq, cigar):
	"""Yield ('M', ref_base, read_base), ('I', read_bases), ('D', ref_bases), ('N', ref_bases).

	seq is the read in reference orientation; soft clips are skipped, hard clips ignored."""
	r, q = pos0, 0
	for n, op in parse_cigar(cigar):
		if op in "M=X":
			if r + n > len(ref) or q + n > len(seq):
				raise ValueError("CIGAR %s runs past the reference or the read" % cigar)
			for k in range(n):
				yield "M", ref[r + k], seq[q + k]
			r += n
			q += n
		elif op == "I":
			yield "I", seq[q:q + n]
			q += n
		elif op == "S":
			q += n
		elif op in "DN":
			if r + n > len(ref):
				raise ValueError("CIGAR %s runs past the reference" % cigar)
			yield op, ref[r:r + n]
			r += n


def compute_nm(ref, pos0, seq, cigar):
	"""Edit distance per the SAM spec: mismatches (any non-ACGT base on either side is a mismatch)
	plus I and D lengths; introns (N) do not count."""
	nm = 0
	for x in _walk(ref, pos0, seq, cigar):
		if x[0] == "M":
			a, b = x[1].upper(), x[2].upper()
			if a not in "ACGT" or b not in "ACGT" or a != b:
				nm += 1
		elif x[0] in "ID":
			nm += len(x[1])
	return nm


def compute_md(ref, pos0, seq, cigar, ambiguous_equal=True):
	"""MD:Z value in the SAM-spec form [0-9]+(([A-Z]|\\^[A-Z]+)[0-9]+)*.

	Bases compare by minimap2's 0-4 code, so with ambiguous_equal=True (the
	default, minimap2's convention) any two non-ACGT bases are equal;
	ambiguous_equal=False is the htslib calmd convention (N never matches).
	Reference bases are written as minimap2 stores them: uppercase, non-ACGT
	as N. Introns (N operations) are skipped. Use md_as_minimap2() to compare
	with minimap2's own MD, which omits a final 0.
	"""
	out, run = [], 0
	for x in _walk(ref, pos0, seq, cigar):
		if x[0] == "M":
			a, b = nt4(x[1]), nt4(x[2])
			if a == b and (a < 4 or ambiguous_equal):
				run += 1
			else:
				out.append("%d%s" % (run, nt5_char(x[1])))
				run = 0
		elif x[0] == "D":
			out.append("%d^%s" % (run, "".join(nt5_char(c) for c in x[1])))
			run = 0
	out.append("%d" % run)
	return "".join(out)


def md_as_minimap2(md):
	"""The MD value minimap2 writes for a spec-form MD: no final 0 after a mismatch or a deletion.

	The SAMtags regex [0-9]+(([A-Z]|\\^[A-Z]+)[0-9]+)* requires the final
	number; base minimap2 (format.c, write_MD_core) omits it when it is 0.
	Use only to compare against minimap2's MD, never to excuse a different
	discrepancy.
	"""
	return md[:-1] if md.endswith("0") and len(md) > 1 and not md[-2].isdigit() else md


def compute_cs(ref, pos0, seq, cigar, long_form=False):
	"""cs:Z value (without the 'cs:Z:' prefix), as minimap2's man page defines it.

	Short form: ':N' identical run, '*xy' substitution (reference x, read y),
	'+seq' insertion, '-seq' deletion, '~xxNyy' intron with its first and
	last two reference bases; bases lowercase. Bases compare by minimap2's
	0-4 code (N equals N) and print as acgtn. long_form writes '=SEQ' (read
	bases, uppercase) instead of ':N'.
	"""
	out, run = [], []

	def flush():
		if run:
			out.append("=" + "".join(run) if long_form else ":%d" % len(run))
			del run[:]
	for x in _walk(ref, pos0, seq, cigar):
		if x[0] == "M":
			if nt4(x[1]) == nt4(x[2]):
				run.append(nt5_char(x[2]))
			else:
				flush()
				out.append("*%s%s" % (nt5_char(x[1], True), nt5_char(x[2], True)))
		elif x[0] == "I":
			flush()
			out.append("+" + "".join(nt5_char(c, True) for c in x[1]))
		elif x[0] == "D":
			flush()
			out.append("-" + "".join(nt5_char(c, True) for c in x[1]))
		else:
			flush()
			s = "".join(nt5_char(c, True) for c in x[1])
			out.append("~%s%d%s" % (s[:2], len(s), s[-2:]))
	flush()
	return "".join(out)


def compute_mlen_blen(ref, pos0, seq, cigar):
	"""(matches, block length, ambiguous bases) as minimap2 reports them for a base-level alignment:
	PAF columns 10 and 11 and the nn:i tag. A column with a non-ACGT base on either side counts in
	neither column 10 nor 11 but in nn; insertions and deletions count in column 11 (their non-ACGT
	bases in nn instead); introns count nowhere."""
	mlen = blen = ambi = 0
	for x in _walk(ref, pos0, seq, cigar):
		if x[0] == "M":
			a, b = nt4(x[1]), nt4(x[2])
			if a > 3 or b > 3:
				ambi += 1
			else:
				blen += 1
				mlen += a == b
		elif x[0] in "ID":
			k = sum(1 for c in x[1] if nt4(c) > 3)
			ambi += k
			blen += len(x[1]) - k
	return mlen, blen, ambi


def compute_de(ref, pos0, seq, cigar):
	"""Gap-compressed per-base divergence (the de:f tag): 1 - matches / (aligned columns + gap opens),
	where an aligned column is an M/=/X position (ambiguous ones included) and every I or D run is one
	gap open. minimap2 prints it with %.4f, or 0 when it is exactly 0."""
	mlen, _, _ = compute_mlen_blen(ref, pos0, seq, cigar)
	cols = sum(n for n, op in parse_cigar(cigar) if op in "M=X")
	gapo = sum(1 for n, op in parse_cigar(cigar) if op in "ID")
	return 1.0 - float(mlen) / (cols + gapo)


_TRANSITIONS = frozenset([(0, 2), (2, 0), (1, 3), (3, 1)])  # A<->G, C<->T as 0-3 codes


def cigar_score(ref, pos0, seq, cigar, a=2, b=4, q=4, e=2, q2=24, e2=1, sc_ambi=1, transition=0):
	"""Alignment score of a CIGAR under minimap2's two-piece affine gap model: +a per match (-A), -b
	per mismatch (-B), -sc_ambi for a column with a non-ACGT base (--score-N), and min(q+k*e,
	q2+k*e2) per k-long gap (the cost the usage text gives for -O/-E). transition is minimap2.1's
	-b (the mismatch penalty for transitions): when non-zero, an A<->G or C<->T mismatch costs
	-transition instead of -b. Introns are not scored; seq is in reference orientation with soft
	clips included."""
	s = 0
	for x in _walk(ref, pos0, seq, cigar):
		if x[0] == "M":
			c1, c2 = nt4(x[1]), nt4(x[2])
			if c1 > 3 or c2 > 3:
				s -= sc_ambi
			elif c1 == c2:
				s += a
			elif transition and (c1, c2) in _TRANSITIONS:
				s -= transition
			else:
				s -= b
		elif x[0] in "ID":
			k = len(x[1])
			s -= min(q + k * e, q2 + k * e2)
	return s


# ---------------------------------------------------------------------------
# Generic SAM/PAF invariants
# ---------------------------------------------------------------------------


def cs_is_long(cs):
	"""True if a cs:Z value is in the long form (--cs=long): identical runs are '=SEQ', never ':N'.
	A long-form value may start with a substitution, an indel or an intron, so the first character
	does not tell; one without any identical run reads the same in both forms."""
	return ":" not in cs


def _seqs_of(genome):
	return genome.seqs if isinstance(genome, Genome) else genome


def sam_invariant_errors(sam, genome, reads=None, check_tags=("NM", "MD", "cs")):
	"""List of violated invariants in a parsed SAM file; empty if all hold.

	genome is a Genome or a {name: sequence} dict; reads (optional) is
	{qname: [read sequences in segment order]} as written to the input.
	Checks: @SQ matches the genome; one primary per read segment; every
	mapped record lies within its contig; the CIGAR (hard clips included)
	spans the whole read; SEQ is the read in alignment orientation (only the
	aligned part for hard-clipped records); and each tag in check_tags that
	is present (NM, MD, cs short or long form) equals its recomputation from
	the genome. Secondary records with SEQ '*' are checked against the
	primary's sequence.
	"""
	seqs = _seqs_of(genome)
	errs = []
	if sam.sq != [(n, len(s)) for n, s in seqs.items()]:
		errs.append("@SQ lines %r differ from the genome" % (sam.sq[:5],))
	try:
		prim = sam.primaries()
	except SamFormatError as e:
		errs.append(str(e))
		prim = {}
	seen = set()
	for r in sam.records:
		seen.add((r.qname, r.segment))
	for k in seen:
		if k not in prim:
			errs.append("no primary record for %s segment %d" % k)
	fwd = {}
	for k, p in prim.items():
		if p.seq != "*" and "H" not in p.cigar:
			fwd[k] = revcomp(p.seq) if p.is_reverse else p.seq
	if reads is not None:
		for name, rs in reads.items():
			for i, s in enumerate(rs):
				seg = 0 if len(rs) == 1 else i + 1
				k = (name, seg)
				if k not in prim:
					errs.append("read %s segment %d missing from the output" % k)
				elif k in fwd and fwd[k].upper() != s.upper():
					errs.append("SEQ of %s segment %d differs from the input read" % k)
				else:
					fwd[k] = s
	for r in sam.records:
		if r.is_unmapped:
			continue
		if r.rname not in seqs:
			errs.append("%s: unknown RNAME %s" % (r.qname, r.rname))
			continue
		clen = len(seqs[r.rname])
		if r.pos < 1 or r.ref_end0 > clen:
			errs.append("%s: alignment %s:%d-%d outside contig of length %d" % (r.qname, r.rname, r.pos, r.ref_end0, clen))
			continue
		read = fwd.get((r.qname, r.segment))
		if read is None:
			continue
		if r.query_length != len(read):
			errs.append("%s: CIGAR %s implies %d bases, read has %d" % (r.qname, r.cigar[:60], r.query_length, len(read)))
			continue
		s = read if not r.is_reverse else revcomp(read)
		ops = r.cigar_ops
		h5 = ops[0][0] if ops and ops[0][1] == "H" else 0
		h3 = ops[-1][0] if len(ops) > 1 and ops[-1][1] == "H" else 0
		s = s[h5:len(s) - h3]
		if r.seq != "*" and r.seq.upper() != s.upper():
			errs.append("%s: SEQ is not the read in alignment orientation" % r.qname)
			continue
		ref = seqs[r.rname]
		if "NM" in check_tags and "NM" in r.tags:
			nm = compute_nm(ref, r.pos0, s, r.cigar)
			if nm != r.tags["NM"]:
				errs.append("%s: NM:i:%d but the alignment has %d differences" % (r.qname, r.tags["NM"], nm))
		if "MD" in check_tags and "MD" in r.tags:
			md = md_as_minimap2(compute_md(ref, r.pos0, s, r.cigar))
			if md != r.tags["MD"]:
				errs.append("%s: MD:Z:%s, recomputed %s" % (r.qname, r.tags["MD"][:80], md[:80]))
		if "cs" in check_tags and "cs" in r.tags:
			cs = compute_cs(ref, r.pos0, s, r.cigar, long_form=cs_is_long(r.tags["cs"]))
			if cs != r.tags["cs"]:
				errs.append("%s: cs:Z:%s, recomputed %s" % (r.qname, r.tags["cs"][:80], cs[:80]))
	return errs


def paf_query_part(rec, read):
	"""The aligned part of `read` (as sequenced) in the target's forward orientation."""
	part = read[rec.qstart:rec.qend]
	return part if rec.strand == "+" else revcomp(part)


def paf_invariant_errors(paf, genome, reads=None, check_tags=("NM", "MD", "cs", "de")):
	"""List of violated invariants in parsed PAF records; empty if all hold.

	reads (optional) is {qname: sequence as sequenced}. Checks: target name
	and length match the genome; 0 <= qstart < qend <= qlen and
	0 <= tstart < tend <= tlen; nmatch <= alen; mapq in 0..255; tp:A is one
	of P/S/I/i; and, when cg:Z is present, its query and target spans match
	the coordinates, columns 10 and 11 equal compute_mlen_blen, and each tag
	of check_tags present (NM, MD, cs, de) equals its recomputation. Not for
	--qstrand output, whose target coordinates follow the query strand.
	"""
	seqs = _seqs_of(genome)
	errs = []
	for p in paf:
		tag = "%s %s:%s-%s" % (p.qname, p.tname, p.tstart, p.tend)
		if reads is not None and p.qname in reads and len(reads[p.qname]) != p.qlen:
			errs.append("%s: qlen %d, read has %d bases" % (tag, p.qlen, len(reads[p.qname])))
		if p.is_unmapped:
			continue
		if p.tname not in seqs:
			errs.append("%s: unknown target" % tag)
			continue
		if p.tlen != len(seqs[p.tname]):
			errs.append("%s: target length %d, genome has %d" % (tag, p.tlen, len(seqs[p.tname])))
		if not (0 <= p.qstart < p.qend <= p.qlen) or not (0 <= p.tstart < p.tend <= p.tlen):
			errs.append("%s: coordinates out of range" % tag)
			continue
		if not (0 <= p.nmatch <= p.alen) or not (0 <= p.mapq <= 255):
			errs.append("%s: nmatch %d, alen %d, mapq %d" % (tag, p.nmatch, p.alen, p.mapq))
		if p.tag("tp") not in ("P", "S", "I", "i"):
			errs.append("%s: tp:A:%s" % (tag, p.tag("tp")))
		cg = p.cigar
		if cg is None:
			continue
		if cigar_query_span(cg) != p.qend - p.qstart or cigar_ref_span(cg) != p.tend - p.tstart:
			errs.append("%s: cg:Z spans %d/%d, coordinates %d/%d" % (tag, cigar_query_span(cg), cigar_ref_span(cg),
				p.qend - p.qstart, p.tend - p.tstart))
			continue
		if reads is None or p.qname not in reads:
			continue
		s = paf_query_part(p, reads[p.qname])
		ref = seqs[p.tname]
		mlen, blen, ambi = compute_mlen_blen(ref, p.tstart, s, cg)
		if (p.nmatch, p.alen) != (mlen, blen):
			errs.append("%s: columns 10/11 are %d/%d, recomputed %d/%d" % (tag, p.nmatch, p.alen, mlen, blen))
		if "nn" in p.tags and p.tags["nn"] != ambi:
			errs.append("%s: nn:i:%d, recomputed %d" % (tag, p.tags["nn"], ambi))
		if "NM" in check_tags and "NM" in p.tags:
			nm = compute_nm(ref, p.tstart, s, cg)
			if nm != p.tags["NM"]:
				errs.append("%s: NM:i:%d, recomputed %d" % (tag, p.tags["NM"], nm))
		if "MD" in check_tags and "MD" in p.tags:
			md = md_as_minimap2(compute_md(ref, p.tstart, s, cg))
			if md != p.tags["MD"]:
				errs.append("%s: MD:Z:%s, recomputed %s" % (tag, p.tags["MD"][:80], md[:80]))
		if "cs" in check_tags and "cs" in p.tags:
			cs = compute_cs(ref, p.tstart, s, cg, long_form=cs_is_long(p.tags["cs"]))
			if cs != p.tags["cs"]:
				errs.append("%s: cs:Z:%s, recomputed %s" % (tag, p.tags["cs"][:80], cs[:80]))
		if "de" in check_tags and "de" in p.tags:
			de = compute_de(ref, p.tstart, s, cg)
			if abs(de - p.tags["de"]) > 0.00006:
				errs.append("%s: de:f:%g, recomputed %.6f" % (tag, p.tags["de"], de))
	return errs


def pair_invariant_errors(sam):
	"""List of violated mate-field invariants for paired records; empty if all hold.

	For every record of a two-segment read, against the mate's primary record M:
	flags 0x1, 0x40/0x80 by segment; 0x8 iff M is unmapped; 0x20 iff M is
	reverse; RNEXT/PNEXT give M's RNAME/POS ('=' when equal to this RNAME);
	'*'/0 when neither has a position. The two primaries carry opposite
	TLENs, zero unless both are mapped to one contig, and agree on flag 0x2.
	"""
	errs = []
	prim = {}
	for r in sam.records:
		if r.is_primary:
			prim[(r.qname, r.segment)] = r
	for r in sam.records:
		if r.segment == 0:
			if r.flag & 0xC9:
				errs.append("%s: unpaired record with pair flags (%d)" % (r.qname, r.flag))
			continue
		if not r.is_paired:
			errs.append("%s: segment flag without 0x1 (%d)" % (r.qname, r.flag))
		if r.is_read1 and r.is_read2:
			errs.append("%s: both 0x40 and 0x80 (%d)" % (r.qname, r.flag))
		m = prim.get((r.qname, 3 - r.segment))
		if m is None:
			errs.append("%s: no primary record for the mate of segment %d" % (r.qname, r.segment))
			continue
		if r.mate_unmapped != m.is_unmapped:
			errs.append("%s/%d: flag 0x8 is %d but the mate is %smapped" % (r.qname, r.segment, r.mate_unmapped, "un" if m.is_unmapped else ""))
		if not m.is_unmapped and r.mate_reverse != m.is_reverse:
			errs.append("%s/%d: flag 0x20 disagrees with the mate's strand" % (r.qname, r.segment))
		if m.rname == "*":
			if r.rnext != "*" or r.pnext != 0:
				errs.append("%s/%d: RNEXT/PNEXT %s/%d for a mate without position" % (r.qname, r.segment, r.rnext, r.pnext))
		else:
			want = "=" if m.rname == r.rname else m.rname
			if r.rnext != want or r.pnext != m.pos:
				errs.append("%s/%d: RNEXT/PNEXT %s/%d, mate is at %s/%d" % (r.qname, r.segment, r.rnext, r.pnext, m.rname, m.pos))
	for (q, s), r in prim.items():
		if s != 1 or (q, 2) not in prim:
			continue
		m = prim[(q, 2)]
		if r.tlen != -m.tlen:
			errs.append("%s: TLENs %d and %d are not opposite" % (q, r.tlen, m.tlen))
		if (r.is_unmapped or m.is_unmapped or r.rname != m.rname) and r.tlen != 0:
			errs.append("%s: TLEN %d without both mates on one contig" % (q, r.tlen))
		if r.is_proper_pair != m.is_proper_pair:
			errs.append("%s: only one mate has flag 0x2" % q)
	return errs


def tlen_5prime(rec, mate):
	"""TLEN as the 5'-to-5' distance ("TLEN#1" in the SAM spec), the form minimap2 writes
	(format.c): from this record's 5' mapped base to the mate's, plus one in the direction of
	the sign, 0 when the two 5' ends coincide. Both records mapped on one contig."""
	this5 = rec.ref_end0 - 1 if rec.is_reverse else rec.pos0
	mate5 = mate.ref_end0 - 1 if mate.is_reverse else mate.pos0
	d = mate5 - this5
	return d + 1 if d > 0 else d - 1 if d < 0 else 0


def tlen_outer(rec, mate):
	"""TLEN as leftmost to rightmost mapped base ("TLEN#2" in the SAM spec); positive for the
	leftmost segment (ties: the forward one). Agrees with tlen_5prime for an FR pair whose
	forward mate starts and ends no later than the reverse mate."""
	left = min(rec.pos0, mate.pos0)
	right = max(rec.ref_end0, mate.ref_end0)
	n = right - left
	first = rec.pos0 < mate.pos0 or (rec.pos0 == mate.pos0 and not rec.is_reverse)
	return n if first else -n


def placement_error(rec, read):
	"""|unclipped start - true start| for a SAM record of a single-segment read, or None on wrong contig/strand."""
	if rec.is_unmapped or rec.rname != read.contig or rec.strand != read.strand:
		return None
	return abs(rec.unclipped_start0 - read.pos)


def paf_placement_error(rec, read):
	"""For a PAF record of a single-segment read: the larger of |start - truth| and |end - truth| after
	extending the target interval by the query clipping, or None on wrong target/strand."""
	if rec.is_unmapped or rec.tname != read.contig or rec.strand != read.strand:
		return None
	lc, rc = (rec.qstart, rec.qlen - rec.qend) if rec.strand == "+" else (rec.qlen - rec.qend, rec.qstart)
	return max(abs(rec.tstart - lc - read.pos), abs(rec.tend + rc - read.end))


# ---------------------------------------------------------------------------
# Pinned outputs
# ---------------------------------------------------------------------------


def check_pinned(name, text, max_diff_lines=40):
	"""Compare text with tests/pins/<name>; MINIMAP2_TEST_UPDATE_PINS=1 rewrites the pin.

	A name ending in .gz is stored gzip-compressed. Pins are generated on
	base minimap2; a test using this must say "pins the current output" in
	its docstring and carry the `pinned` marker.
	"""
	path = os.path.join(PINS_DIR, name)
	if os.environ.get("MINIMAP2_TEST_UPDATE_PINS") == "1":
		os.makedirs(os.path.dirname(path), exist_ok=True)
		_write_text(path, text, name.endswith(".gz"))
		return
	if not os.path.exists(path):
		raise AssertionError("missing pin %s; create it on base minimap2 with MINIMAP2_TEST_UPDATE_PINS=1" % path)
	with _open_in(path) as f:
		want = f.read()
	if want != text:
		diff = list(difflib.unified_diff(want.splitlines(), text.splitlines(), "pinned/" + name, "actual", lineterm="", n=1))
		more = "" if len(diff) <= max_diff_lines else "\n... (%d more diff lines)" % (len(diff) - max_diff_lines)
		raise AssertionError("output differs from pin %s:\n%s%s" % (name, "\n".join(diff[:max_diff_lines]), more))


# ---------------------------------------------------------------------------
# C drivers
# ---------------------------------------------------------------------------


class CResult(object):
	def __init__(self, argv, returncode, stdout, stderr):
		self.argv, self.returncode, self.stdout, self.stderr = argv, returncode, stdout, stderr
		self.failures = [l for l in stdout.splitlines() if l.startswith("not ok")]
		m = re.search(r"^# mmt: \S+ checks=(\d+) failures=(\d+)$", stdout, re.M)
		self.checks = int(m.group(1)) if m else None
		self.n_failures = int(m.group(2)) if m else None

	@property
	def passed(self):
		return self.returncode == 0 and not self.failures and self.checks is not None and self.n_failures == 0

	def report(self):
		return "$ %s\nexit %d\n--- stdout ---\n%s\n--- stderr ---\n%s" % (" ".join(self.argv), self.returncode,
			self.stdout[-6000:], self.stderr[-3000:])


def c_cases(driver):
	"""Names of the cases a driver in tests/c/<driver>.c defines."""
	exe = minimap2_build.c_driver(driver)
	p = subprocess.run([exe, "list"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60)
	if p.returncode != 0:
		raise AssertionError("%s list failed: %s" % (driver, p.stderr.decode("utf-8", "replace")))
	return p.stdout.decode().split()


def run_c_case(driver, case, args=(), timeout=None, stdin=None, env=None):
	"""Run one case of a C driver (in a temporary directory, LC_ALL=C, no core dumps); return a CResult."""
	exe = minimap2_build.c_driver(driver)
	r = run_process([exe, "run", case] + [str(a) for a in args], timeout=timeout, stdin=stdin, env=env)
	return CResult(r.argv, r.returncode, r.stdout, r.stderr)


def assert_c_case(driver, case, *args, **kw):
	"""Run a C driver case and raise AssertionError with its output unless every check passed."""
	r = run_c_case(driver, case, args, **kw)
	if not r.passed:
		raise AssertionError("C case %s/%s failed\n%s" % (driver, case, r.report()))
	return r


# ---------------------------------------------------------------------------
# k8 scripts and python/minimap2.py
# ---------------------------------------------------------------------------

PAFTOOLS = os.path.join(REPO_ROOT, "misc", "paftools.js")
PAFCLUSTER = os.path.join(REPO_ROOT, "misc", "pafcluster.js")


def k8_path():
	"""The k8 binary on PATH, or None (tests that need it skip)."""
	return shutil.which("k8")


def run_k8(script, args=(), stdin=None, timeout=None, cwd=None, check=False):
	"""Run `k8 <script> <args>` (LC_ALL=C, fresh temporary directory unless cwd)."""
	k8 = k8_path()
	if k8 is None:
		raise AssertionError("k8 is not on PATH")
	src = minimap2_build.src_dir()
	script = os.path.join(src, "misc", os.path.basename(script)) if not os.path.isabs(script) else script
	return run_process([k8, script] + [str(a) for a in args], timeout=timeout, stdin=stdin, cwd=cwd, check=check)


def run_paftools(args, stdin=None, timeout=None, cwd=None, check=True):
	"""Run `k8 misc/paftools.js <args>` from $MINIMAP2_TEST_SRC_DIR."""
	return run_k8("paftools.js", args, stdin=stdin, timeout=timeout, cwd=cwd, check=check)


def run_minimap2_py(args, stdin=None, timeout=None, check=True):
	"""Run python/minimap2.py from $MINIMAP2_TEST_SRC_DIR with this interpreter and the cached mappy build."""
	lib = minimap2_build.mappy_dir()
	env = {"PYTHONPATH": lib + os.pathsep + os.environ.get("PYTHONPATH", "")}
	script = os.path.join(minimap2_build.src_dir(), "python", "minimap2.py")
	return run_process([sys.executable, script] + [str(a) for a in args], timeout=timeout, stdin=stdin, env=env,
		check=check)
