"""Helpers for the area C tests (splice and short reads).

Hand-built loci whose correct alignment is known by construction:

	tie_locus           two GT-AG junctions that align a read equally well
	                    (the intron shifted by 4 bp), with control over the
	                    flanking base of each splice site
	long_intron_locus   two 400 bp exons around a long GT-AG intron
	noncanonical_locus  two exons around an intron with no GT/GC/AT donor and
	                    no AG/AC acceptor on either strand
	alt_acceptor_locus  one donor and two acceptors whose downstream exons
	                    start with the same bases, for ambiguous jumps

and readers for what the tests need from them:

	overhang_reads      reads that end a few bases past an exon boundary
	spliced_truth_tsv   truth file for the drv_c_api map_splice case
	intron_bits         the bit pattern mm_idx_bed_junc() writes for introns
"""

import os

import minimap2_helpers as H


class Locus(object):
	"""A genome with one engineered gene: .genome (Genome), .fasta (path once written), .exons and
	.introns in genome coordinates, .strand of the gene, .read (the spliced transcript in genome
	forward orientation) and extra attributes per construction."""

	def __init__(self, genome, exons, strand, **kw):
		self.genome, self.exons, self.strand = genome, exons, strand
		self.fasta = None
		self.__dict__.update(kw)

	@property
	def introns(self):
		return [(self.exons[i][1], self.exons[i + 1][0]) for i in range(len(self.exons) - 1)]

	@property
	def read(self):
		ctg = self.genome["chr1"]
		return "".join(ctg[s:e] for s, e in self.exons)

	def write(self, path):
		self.fasta = self.genome.write_fasta(path)
		return self.fasta


def _assemble(rng, parts, strand, flank=3000, decoy=5000):
	"""chr1 = random flank + parts (joined, reverse-complemented for a '-' gene) + random flank, and a
	random chr2. Returns (Genome, offset of parts[0] in chr1, total length of parts)."""
	body = "".join(parts)
	left, right = rng.dna(flank), rng.dna(flank)
	ctg = left + (body if strand == "+" else H.revcomp(body)) + right
	g = H.Genome([("chr1", ctg), ("chr2", rng.fork("decoy").dna(decoy))])
	return g, len(left), len(body)


def _flip(iv, off, total, strand):
	"""Map an interval of the gene body (transcript orientation) to chr1 coordinates."""
	s, e = iv
	if strand == "+":
		return (off + s, off + e)
	return (off + total - e, off + total - s)


def tie_locus(seed, strand="+", flank_a=True, flank_b=True, exon_len=300, intron_len=600):
	"""Exon 1 ends with GTAG and the intron ends with AGGTAG, so the read (exon 1 + exon 2) aligns
	equally well with intron A = [s, e) or with intron B, A shifted 4 bp towards the 5' end of the
	transcript; both are GT..AG. Every other shift breaks the GT..AG signal.

	The donor of B is always followed by A (GTA..) and the acceptor of A is always preceded by T
	(..TAG), both the conserved bases of README.md (GT[A/G]..[C/T]AG). flank_a sets the base after
	A's donor (A/G if True, C/T if False); flank_b the base before B's acceptor (C/T if True,
	A/G if False). With both True, A and B tie under every splice model.

	Returns a Locus with .A and .B (genome intervals of the two introns) and .exons for A.
	"""
	rng = H.Rng(seed)
	e1 = rng.dna(exon_len - 5)
	e1 += rng.choice("ACT") + "GTAG"  # the base before GTAG is not G: no 5 bp shift
	third = rng.choice("AG") if flank_a else rng.choice("CT")
	b7 = rng.choice("CT") if flank_b else rng.choice("AG")
	intron = "GT" + third + rng.dna(intron_len - 10) + b7 + "AGGTAG"
	e2 = rng.choice("ACT") + rng.dna(exon_len - 1)  # not G: no shift towards exon 2
	assert len(intron) == intron_len
	g, off, total = _assemble(rng, [e1, intron, e2], strand)
	A = (len(e1), len(e1) + intron_len)
	B = (A[0] - 4, A[1] - 4)
	gA, gB = _flip(A, off, total, strand), _flip(B, off, total, strand)
	exons = [(off, gA[0]), (gA[1], off + total)]
	loc = Locus(g, exons, strand, A=gA, B=gB)
	ctg = g["chr1"]
	for s, e in (gA, gB):  # both are canonical on the gene strand
		sig = ctg[s:s + 2] + ctg[e - 2:e]
		assert sig == ("GTAG" if strand == "+" else "CTAC"), sig
	return loc


def long_intron_locus(seed, intron_len=6000, strand="+"):
	"""Two 400 bp exons around a GT..AG intron of intron_len bp."""
	rng = H.Rng(seed)
	e1 = rng.dna(400)
	e2 = rng.dna(400)
	intron = "GTA" + rng.dna(intron_len - 6) + "CAG"
	g, off, total = _assemble(rng, [e1, intron, e2], strand)
	iv = _flip((400, 400 + intron_len), off, total, strand)
	return Locus(g, [(off, iv[0]), (iv[1], off + total)], strand)


def noncanonical_locus(seed, intron_len=800):
	"""Two 400 bp exons around an intron that starts with CA and ends with TC: neither a donor
	(GT, GC, AT) nor an acceptor (AG, AC) on either strand. The last exon-1 base differs from the
	last intron base and the first exon-2 base from the first intron base, so the intron cannot
	slide by one base without a mismatch."""
	rng = H.Rng(seed)
	intron = "CA" + rng.dna(intron_len - 4) + "TC"
	e1 = rng.dna(399) + rng.choice([b for b in "ACGT" if b != intron[-1]])
	e2 = rng.choice([b for b in "ACGT" if b != intron[0]]) + rng.dna(399)
	g, off, total = _assemble(rng, [e1, intron, e2], "+")
	return Locus(g, [(off, off + 400), (off + 400 + intron_len, off + total)], "+")


def alt_acceptor_locus(seed, shared=20, exon1_len=300, intron1=500, exon2_len=200, intron2=500):
	"""exon 1 | intron | exon 2a | spacer | exon 2b, '+' strand: exon 2a and exon 2b start with the
	same `shared` bases, and exon 1 can be spliced to either. A read whose 3' end runs a few
	bases past exon 1 matches both junctions. Returns a Locus with .ja and .jb (the two introns,
	genome coordinates) and .exons = [exon 1, exon 2a]."""
	rng = H.Rng(seed)
	head = rng.choice("ACT") + rng.dna(shared - 1)  # exon 2 does not start with G (no shift into the intron)
	e1 = rng.dna(exon1_len - 1) + rng.choice("ACT")
	i1 = "GTA" + rng.dna(intron1 - 6) + "CAG"
	e2a = head + rng.dna(exon2_len - shared)
	i2 = "GTA" + rng.dna(intron2 - 6) + "CAG"
	e2b = head + rng.dna(exon2_len - shared)
	parts = [e1, i1, e2a, i2, e2b]
	g, off, total = _assemble(rng, parts, "+")
	s = off + exon1_len
	a2 = s + intron1
	b2 = a2 + exon2_len + intron2
	loc = Locus(g, [(off, s), (a2, a2 + exon2_len)], "+", ja=(s, a2), jb=(s, b2), e1=(off, s))
	return loc


def overhang_reads(genome, genes, overhangs=(4, 8, 12), length=100, prefix="oh"):
	"""For every exon boundary of every gene and every k in overhangs, two reads of `length` bases
	(genome forward orientation): one ending k bases into the next exon ('end') and one starting k
	bases before the boundary ('start'). Returns [(name, seq, truth)] with truth = (contig,
	[(start, end) of the two genome blocks], k)."""
	out = []
	for gi, gene in enumerate(genes):
		ex = gene.exons
		for j in range(len(ex) - 1):
			left, right = ex[j], ex[j + 1]
			for k in overhangs:
				for side in ("end", "start"):
					if side == "end":  # the read ends k bases into the right exon
						a = length - k
						if right[1] - right[0] < k or left[1] - left[0] < a:
							continue
						blocks = [(left[1] - a, left[1]), (right[0], right[0] + k)]
					else:  # the read starts k bases before the end of the left exon
						a = length - k
						if left[1] - left[0] < k or right[1] - right[0] < a:
							continue
						blocks = [(left[1] - k, left[1]), (right[0], right[0] + a)]
					seq = "".join(genome[gene.contig][s:e] for s, e in blocks)
					out.append(("%s%d_%d_%s%d" % (prefix, gi, j, side, k), seq, (gene.contig, blocks, k)))
	return out


def spliced_truth_tsv(path, reads):
	"""name, contig, genome strand, transcript strand, introns "s-e,..." ("." if none), one line per read."""
	with open(path, "w") as f:
		for r in reads:
			intr = ",".join("%d-%d" % iv for iv in r.introns) or "."
			f.write("%s\t%s\t%s\t%s\t%s\n" % (r.name, r.contig, r.strand, r.ts_strand, intr))
	return path


def intron_bits(introns, st, en):
	"""{position: bits} that mm_idx_bed_junc() writes for oriented introns [(s, e, strand), ...]
	inside the window [st, en): 1 at a '+' intron's first base, 2 at its last base, 8 and 4 for a
	'-' intron. Introns not entirely inside the window are left out."""
	bits = {}
	for s, e, strand in introns:
		if s < st or e > en:
			continue
		a, b = (1, 2) if strand == "+" else (8, 4)
		bits[s] = bits.get(s, 0) | a
		bits[e - 1] = bits.get(e - 1, 0) | b
	return bits


def encode_bits(d):
	"""'pos:val,pos:val' for a C driver argument (empty string for nothing)."""
	return ",".join("%d:%d" % (k, d[k]) for k in sorted(d))


def write_reads(path, records):
	"""FASTA or FASTQ of (name, seq) tuples or SimReads, chosen by the extension."""
	return H.write_fastx(path, records)


def tmp(tmp_path, name):
	return str(tmp_path / name) if hasattr(tmp_path, "__truediv__") else os.path.join(str(tmp_path), name)
