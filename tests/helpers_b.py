"""Helpers for the seeding, chaining, alignment and hit selection tests (test_b_*.py).

Reference implementations written from minimap2's documentation and papers,
not from its code paths, so tests can derive expected values independently:

	mm_hash64 / minimizers      the invertible integer hash and the (w,k)-minimizer
	                            rule (minimap2_helpers.kmer_hash and .minimizers, with
	                            a tied window raising ValueError)
	MinimizerIndex / anchors    reference minimizers by k-mer, and the anchors
	                            (matching minimizers) of a query against them
	chain_dp / chain_gap_cost   the chaining score of Li (2018), eq. 1-2, with the
	                            gap cost scaled as minimap2 documents for
	                            --chain-gap-scale/--chain-skip-scale
	covered_bases               the chaining score of gap-free collinear anchors
	cigar_score_full            CIGAR score with -A/-B/-b/-O/-E/--score-N (minimap2_helpers.cigar_score)
	peak_score                  the max-scoring segment of an alignment (ms:i, -s)
	affine_global / affine_ext  brute-force two-piece affine alignment (Gotoh)
	many_copy_genome            a 160 kbp genome with a 120-copy family (high-occurrence seeds)
	large_genome                a 1.5 Mbp genome for thread-invariance and pinned output
"""

import math

import minimap2_helpers as H

def mm_hash64(key, k):
	"""Thomas Wang's invertible 64-bit integer hash restricted to 2k bits, as minimap2 uses it
	(Li 2016, minimap paper; Li 2018): a bijection on [0, 4^k) (minimap2_helpers.kmer_hash)."""
	return H.kmer_hash(key, k)


def kmers(seq, k):
	"""[(end position, hash of the canonical k-mer, strand)] for every k-mer of seq without an ambiguous
	base (minimap2_helpers.kmer_runs; odd k only)."""
	return [(q.pos, q.hash, q.strand) for run in H.kmer_runs(seq, k) for q in run]


def minimizers(seq, w, k):
	"""The (w,k)-minimizers of seq as a sorted list of (end position, hash, strand)
	(minimap2_helpers.minimizers). A window holding its minimum twice (a tandem repeat) raises
	ValueError: B's brute-force chaining assumes no ties."""
	return H.minimizers(seq, w, k, ties=False)


class MinimizerIndex(object):
	"""Reference minimizers of a {name: seq} dict (or Genome): hash -> [(contig, end pos, strand)]."""

	def __init__(self, seqs, w, k):
		seqs = seqs.seqs if isinstance(seqs, H.Genome) else seqs
		self.w, self.k = w, k
		self.names = list(seqs)
		self.lens = dict((n, len(s)) for n, s in seqs.items())
		self.occ = {}
		for n, s in seqs.items():
			for pos, h, z in minimizers(s, w, k):
				self.occ.setdefault(h, []).append((n, pos, z))

	def count(self, h):
		return len(self.occ.get(h, ()))

	def mid_occ(self, frac):
		"""-f FLOAT: one more than the count at rank floor((1 - frac) * n) of the sorted per-minimizer
		counts, i.e. minimizers in the top `frac` fraction by count are ignored (minimap2.1, -f)."""
		counts = sorted(len(v) for v in self.occ.values())
		return counts[int((1.0 - frac) * len(counts))] + 1


class Anchor(object):
	"""A seed match: contig, strand ('+'/'-'), x = target end position, y = query end position on the
	query strand that matches the target's forward strand, span = k."""
	__slots__ = ("contig", "strand", "x", "y", "span", "qpos")

	def __init__(self, contig, strand, x, y, span, qpos):
		self.contig, self.strand, self.x, self.y, self.span, self.qpos = contig, strand, x, y, span, qpos

	def __repr__(self):
		return "Anchor(%s%s x=%d y=%d)" % (self.contig, self.strand, self.x, self.y)


def anchors(query, index, max_occ=None):
	"""Anchors of query against a MinimizerIndex, skipping minimizers with more than max_occ
	occurrences. For a reverse-strand anchor, y is on the reverse complement of the query."""
	qlen = len(query)
	out = []
	for pos, h, z in minimizers(query, index.w, index.k):
		occ = index.occ.get(h, ())
		if not occ or (max_occ is not None and len(occ) > max_occ):
			continue
		for contig, tpos, tz in occ:
			if tz == z:
				out.append(Anchor(contig, "+", tpos, pos, index.k, pos))
			else:
				out.append(Anchor(contig, "-", tpos, qlen - (pos + 1 - index.k) - 1, index.k, pos))
	return out


def repeat_len(query, index, max_occ):
	"""rl:i: the length of the union of the query intervals of minimizers occurring more than max_occ
	times in the reference."""
	iv = [(pos + 1 - index.k, pos + 1) for pos, h, _ in minimizers(query, index.w, index.k)
		if index.count(h) > max_occ]
	return union_length(iv)


def union_length(intervals):
	tot, end = 0, None
	for s, e in sorted(intervals):
		if end is None or s >= end:
			tot += e - s
			end = e
		elif e > end:
			tot += e - end
			end = e
	return tot


def chain_gap_cost(dd, dg, k, gap_scale=0.8, skip_scale=0.0):
	"""The concave gap cost between two anchors (Li 2018, eq. 2): for a diagonal difference dd != 0,
	gap_scale * 0.01 * k * dd + skip_scale * 0.01 * k * dg + 0.5 * log2(dd + 1), truncated to an
	integer; 0 when dd == 0 and only the skip term applies otherwise."""
	lin = gap_scale * 0.01 * k * dd + skip_scale * 0.01 * k * dg
	log = math.log(dd + 1, 2) if dd >= 1 else 0.0
	return lin + 0.5 * log


def chain_dp(anc, k, max_gap=5000, bw=500, gap_scale=0.8, skip_scale=0.0):
	"""Exact chaining DP over anchors of one contig and strand (Li 2018, eq. 1):
	f(i) = max(span, max_j f(j) + min(dq, dr, span) - gap cost) over predecessors j with
	0 < dq <= max_gap, 0 < dr <= max_gap and |dr - dq| <= bw. Returns (best score, anchors of the best
	chain, list of (f, predecessor)). The fractional gap cost is truncated as minimap2 truncates it;
	`margins` of the transitions used are recorded for callers that must avoid rounding edges."""
	a = sorted(anc, key=lambda z: (z.x, z.y))
	n = len(a)
	f = [0] * n
	p = [-1] * n
	edge = [None] * n
	for i in range(n):
		f[i] = a[i].span
		for j in range(i - 1, -1, -1):  # on a tie the nearest predecessor wins, so collinear anchors all chain
			dq = a[i].y - a[j].y
			dr = a[i].x - a[j].x
			if dq <= 0 or dr <= 0 or dq > max_gap or dr > max_gap:
				continue
			dd = abs(dr - dq)
			if dd > bw:
				continue
			dg = min(dr, dq)
			sc = min(a[j].span, dg)
			cost = 0.0
			if dd or dg > a[j].span:
				cost = chain_gap_cost(dd, dg, k, gap_scale, skip_scale)
			s = f[j] + sc - int(cost)
			if s > f[i]:
				f[i], p[i], edge[i] = s, j, cost
	if n == 0:
		return 0, [], []
	best = max(range(n), key=lambda i: (f[i], -i))
	chain, i = [], best
	costs = []
	while i >= 0:
		chain.append(a[i])
		if edge[i] is not None:
			costs.append(edge[i])
		i = p[i]
	return f[best], chain[::-1], costs


def covered_bases(chain):
	"""Chaining score of anchors on one diagonal with no gap cost: the number of query bases covered
	by their k-mers."""
	return union_length([(c.y + 1 - c.span, c.y + 1) for c in chain])


# ---------------------------------------------------------------------------
# Alignment scores
# ---------------------------------------------------------------------------

_TRANSITION = {(0, 2), (2, 0), (1, 3), (3, 1)}  # A<->G, C<->T


def column_score(t, q, a, b, transition, sc_ambi):
	ct, cq = H.nt4(t), H.nt4(q)
	if ct > 3 or cq > 3:
		return -sc_ambi
	if ct == cq:
		return a
	if transition and (ct, cq) in _TRANSITION:
		return -transition
	return -b


# -A/-B/-b/-O/-E/--score-N score of a CIGAR: the shared scorer
cigar_score_full = H.cigar_score


def peak_score(ref, pos0, seq, cigar, a=2, b=4, q=4, e=2, sc_ambi=1, transition=0):
	"""The score of the max-scoring segment of a gap-free alignment (ms:i and the -s filter for
	alignments without indels): the largest sum of column scores over any run of columns."""
	best = run = 0
	for x in H._walk(ref, pos0, seq, cigar):
		if x[0] != "M":
			raise ValueError("peak_score handles gap-free alignments only")
		run += column_score(x[1], x[2], a, b, transition, sc_ambi)
		if run < 0:
			run = 0
		best = max(best, run)
	return best


NEG = -(1 << 40)


def affine_matrices(qs, ts, a, b, q, e, q2, e2, ext):
	"""Gotoh DP of query qs against target ts with gap cost min(q+k*e, q2+k*e2). With ext=False the
	whole of both is aligned; with ext=True the alignment starts at (0, 0) and may end anywhere.
	Returns H (best score of an alignment of qs[:i] and ts[:j])."""
	n, m = len(qs), len(ts)
	Hm = [[NEG] * (m + 1) for _ in range(n + 1)]
	E1 = [[NEG] * (m + 1) for _ in range(n + 1)]
	F1 = [[NEG] * (m + 1) for _ in range(n + 1)]
	E2 = [[NEG] * (m + 1) for _ in range(n + 1)]
	F2 = [[NEG] * (m + 1) for _ in range(n + 1)]
	Hm[0][0] = 0
	for i in range(n + 1):
		for j in range(m + 1):
			if i == 0 and j == 0:
				continue
			if j > 0:  # deletion (target base without query)
				E1[i][j] = max(Hm[i][j - 1] - q - e, E1[i][j - 1] - e)
				E2[i][j] = max(Hm[i][j - 1] - q2 - e2, E2[i][j - 1] - e2)
			if i > 0:
				F1[i][j] = max(Hm[i - 1][j] - q - e, F1[i - 1][j] - e)
				F2[i][j] = max(Hm[i - 1][j] - q2 - e2, F2[i - 1][j] - e2)
			best = max(E1[i][j], E2[i][j], F1[i][j], F2[i][j])
			if i > 0 and j > 0:
				best = max(best, Hm[i - 1][j - 1] + (a if qs[i - 1] == ts[j - 1] else -b))
			Hm[i][j] = best
	return Hm


def affine_ext(qs, ts, a, b, q, e, q2=None, e2=None):
	"""Extension alignment from the start of both: (max, mqe) where max is the best score over all
	(query prefix, target prefix) pairs including the empty one (0) and mqe the best score of an
	alignment consuming the whole query."""
	q2 = q if q2 is None else q2
	e2 = e if e2 is None else e2
	Hm = affine_matrices(qs, ts, a, b, q, e, q2, e2, True)
	mx = max(max(r) for r in Hm)
	mqe = max(Hm[len(qs)])
	return mx, mqe


# ---------------------------------------------------------------------------
# Genomes for seeding tests
# ---------------------------------------------------------------------------


def copy_number_genome(seed=31, contig_lengths=(40000, 30000), families=((400, 3), (300, 6), (200, 12)), gap=300):
	"""Random contigs with repeat families of different copy numbers: for each (unit length, copies)
	a random unit is written `copies` times (alternating strands, exact copies) at non-overlapping
	places at least `gap` bases apart. Returns (Genome, [(family, contig, start, end, strand)])."""
	rng = H.Rng(seed)
	seqs = [list(rng.dna(n)) for n in contig_lengths]
	names = ["c%d" % (i + 1) for i in range(len(seqs))]
	placed, occupied = [], []
	for fam, (ulen, copies) in enumerate(families):
		unit = rng.dna(ulen)
		for c in range(copies):
			for _ in range(1000):
				ci = rng.below(len(seqs))
				st = 200 + rng.below(len(seqs[ci]) - ulen - 400)
				if all(not (oc == ci and os - gap < st + ulen and st < oe + gap) for oc, os, oe in occupied):
					break
			else:
				raise ValueError("cannot place repeat copies")
			strand = "+" if c % 2 == 0 else "-"
			seqs[ci][st:st + ulen] = list(unit if strand == "+" else H.revcomp(unit))
			occupied.append((ci, st, st + ulen))
			placed.append((fam, names[ci], st, st + ulen, strand))
	contigs = [(names[i], "".join(s)) for i, s in enumerate(seqs)]
	return H.Genome(contigs), placed


def collect_expectation(query, index, max_occ, max_max_occ=4095, dist=0):
	"""What mm_collect_matches should report for a query, from brute-force minimizers: (number of
	anchors, number of kept query minimizers, repetitive query length or None).

	Query minimizers absent from the reference are ignored. A minimizer occurring more than max_occ
	times is dropped, except that with dist > 0 and max_max_occ > max_occ, in every run of such
	minimizers between two kept ones (spanning ps..pe on the query) the round((pe - ps) / dist) with
	the fewest occurrences (at most 128) are kept (-e: sample a high-frequency minimizer every dist
	bp) unless they occur more than max_max_occ times. rep_len (the union of the dropped minimizers'
	k-mers) is returned only when no sampling happens, since which of equally frequent minimizers are
	sampled is not documented."""
	m = [(pos, index.count(h)) for pos, h, _ in minimizers(query, index.w, index.k) if index.count(h) > 0]
	keep = [c <= max_occ for _, c in m]
	sampling = dist > 0 and max_max_occ > max_occ
	if sampling:
		i, n = 0, len(m)
		last0 = -1
		while i <= n:
			if i == n or m[i][1] <= max_occ:
				if i - last0 > 1:
					ps = 0 if last0 < 0 else m[last0][0]
					pe = len(query) if i == n else m[i][0]
					k = min(int((pe - ps) / float(dist) + .499), 128)
					run = sorted(range(last0 + 1, i), key=lambda j: m[j][1])
					for j in run[:k]:
						keep[j] = m[j][1] <= max_max_occ
				last0 = i
			i += 1
	n_a = sum(c for (_, c), kp in zip(m, keep) if kp)
	n_kept = sum(keep)
	rep = None
	if not sampling:
		rep = union_length([(p + 1 - index.k, p + 1) for (p, c), kp in zip(m, keep) if not kp])
	return n_a, n_kept, rep


def max_zdrop(ref, pos0, seq, cigar, a=2, b=4, q=4, e=2, q2=24, e2=1, sc_ambi=1):
	"""The largest Z-drop along an alignment (minimap2.1, -z): over pairs of points on the path, the
	running score at the earlier point minus that at the later one, minus e times the difference of
	their diagonals; the running score adds the column scores and subtracts min(q+k*e, q2+k*e2) at
	each k-long gap."""
	s = 0
	best, bi, bj = 0, 0, 0
	i = j = 0  # query and target offsets
	drop = 0
	for x in H._walk(ref, pos0, seq, cigar):
		if x[0] == "M":
			s += column_score(x[1], x[2], a, b, 0, sc_ambi)
			i += 1
			j += 1
		elif x[0] in "ID":
			n = len(x[1])
			s -= min(q + n * e, q2 + n * e2)
			if x[0] == "I":
				i += n
			else:
				j += n
		else:
			continue
		if s > best:
			best, bi, bj = s, i, j
		else:
			drop = max(drop, best - s - e * abs((i - bi) - (j - bj)))
	return drop


def local_score(qs, ts, a=2, b=4, q=4, e=2):
	"""Best Smith-Waterman score of qs against ts, a k-long gap costing q + k*e."""
	n, m = len(qs), len(ts)
	Hp = [0] * (m + 1)
	Ep = [NEG] * (m + 1)
	best = 0
	for i in range(1, n + 1):
		Hc = [0] * (m + 1)
		Fc = NEG
		for j in range(1, m + 1):
			Ep[j] = max(Hp[j] - q - e, Ep[j] - e)  # gap in the target direction (query base unaligned)
			Fc = max(Hc[j - 1] - q - e, Fc - e)
			h = max(0, Hp[j - 1] + (a if qs[i - 1] == ts[j - 1] else -b), Ep[j], Fc)
			Hc[j] = h
			best = max(best, h)
		Hp = Hc
	return best


def many_copy_genome():
	"""A 160 kbp random genome (100 and 60 kbp contigs) with one family of 120 copies of a 300 bp unit at
	2% divergence, each on a random strand: the unit's minimizers occur up to 96 times, so the -f
	thresholds of the default fractions (92 at 0.0002, 80 at 0.002) fall inside the family's counts.
	Seed 37 is one whose brute-force minimizers have no tied window (helpers_b.minimizers raises on one).
	Returns (Genome, [(contig, start, end, strand)] of the copies)."""
	g = H.make_genome(seed=37, contig_lengths=(100000, 60000), repeat_families=1, repeat_copies=120, repeat_len=(300, 300),
		repeat_divergence=0.02)
	copies = [(f.contig, f.start, f.end, f.info[2]) for f in g.features if f.kind == "repeat"]
	return g, copies


def large_genome():
	"""1.5 Mbp: chr1 (600 kbp) and chr2 (400 kbp) with six families of three 1-4 kbp segmental duplications
	at 0.5% divergence, four inverted repeats, eight low-complexity tracts and three N runs; chr3 (500 kbp)
	with a family of 150 copies of a 280-320 bp unit at 3% divergence (an Alu-like family whose
	minimizers exceed the default occurrence thresholds, so high-occurrence seeds are filtered and
	sampled) and four low-complexity tracts. Changing it changes the pins built on it."""
	a = H.make_genome(seed=91, contig_lengths=(600000, 400000), repeat_families=6, repeat_copies=3, repeat_len=(1000, 4000),
		repeat_divergence=0.005, inversions=4, low_complexity=8, n_runs=3)
	b = H.make_genome(seed=92, contig_lengths=(500000,), names=["chr3"], repeat_families=1, repeat_copies=150,
		repeat_len=(280, 320), repeat_divergence=0.03, low_complexity=4)
	return H.Genome([(n, a[n]) for n in a.names] + [("chr3", b["chr3"])], a.features + b.features)
