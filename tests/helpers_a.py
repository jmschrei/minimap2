"""Helpers for the area A tests (input, sketching, index).

Standard library only.

	hash64, kmer_runs,          a reference (w,k)-minimizer scan written from
	SketchRef                   the rule in minimap2.1 ("a minimizer is the
	                            smallest k-mer in a window of w consecutive
	                            k-mers"), sketch.c's documented hash and the
	                            HPC definition of -H
	read_mmi / MmiPart          a reader for the .mmi format (minimap.h's
	                            MM_IDX_MAGIC, mm_idx_dump's layout) that
	                            decodes the minimizer table into
	                            {hash: positions}
	sdust_reference             the symmetric DUST rule of Morgulis et al.
	                            (2006): the union of perfect intervals of at
	                            most W bases
	random_dna, write_fa        small sequence utilities
"""

import bisect
import collections
import os
import struct

import minimap2_helpers as H

# ---------------------------------------------------------------------------
# Minimizers
# ---------------------------------------------------------------------------

# sketch.c's seq_nt4_table: A/C/G/T/U in either case -> 0..3, anything else -> 4
_NT4 = {"A": 0, "C": 1, "G": 2, "T": 3, "U": 3, "a": 0, "c": 1, "g": 2, "t": 3, "u": 3}


def hash64(key, mask):
	"""sketch.c's invertible integer hash of a 2k-bit k-mer, restricted to `mask`."""
	key = (~key + (key << 21)) & mask
	key = key ^ (key >> 24)
	key = ((key + (key << 3)) + (key << 8)) & mask
	key = key ^ (key >> 14)
	key = ((key + (key << 2)) + (key << 4)) & mask
	key = key ^ (key >> 28)
	key = (key + (key << 31)) & mask
	return key


Kmer = collections.namedtuple("Kmer", "x pos strand")  # x = hash << 8 | span; pos = last base; strand 0 fwd < rev
SKIP = (1 << 64) - 1  # x of a k-mer that cannot be a minimizer (HPC span >= 256)


def kmer_runs(seq, k, hpc=False):
	"""The canonical k-mers of seq, one list per maximal run of A/C/G/T/U.

	A k-mer is canonical on the strand whose 2-bit encoding is smaller
	(strand 0: forward); it is keyed by x = hash64(kmer) << 8 | span, where
	span is k, or with hpc the length in the original sequence of the k
	homopolymer-compressed bases. pos is the 0-based position of the k-mer's
	last base (with hpc, the last base of its last homopolymer run). Any
	other character ends a run. An HPC k-mer spanning 256 bases or more
	cannot be stored in x's 8 bits: it keeps its place in the window but is
	never a minimizer (x = SKIP). Only odd k (no k-mer is its own reverse
	complement) is modelled; even k raises.
	"""
	if k % 2 == 0 or not 0 < k <= 28:
		raise ValueError("the reference models odd k in 1..28 only")
	mask = (1 << 2 * k) - 1
	shift = 2 * (k - 1)
	runs, cur = [], None
	n, i = len(seq), 0
	fwd = rev = l = span = 0
	spans = collections.deque()
	while i < n:
		c = _NT4.get(seq[i], 4)
		if c > 3:
			cur = None
			i += 1
			continue
		if cur is None:
			cur = []
			runs.append(cur)
			fwd = rev = l = span = 0
			spans.clear()
		run = 1
		if hpc:
			while i + run < n and _NT4.get(seq[i + run], 4) == c:
				run += 1
			i += run - 1
			spans.append(run)
			span += run
			if len(spans) > k:
				span -= spans.popleft()
		else:
			span = k
		fwd = ((fwd << 2) | c) & mask
		rev = (rev >> 2) | ((3 - c) << shift)
		l += 1
		if l >= k:
			z = 0 if fwd < rev else 1
			cur.append(Kmer(hash64(fwd if z == 0 else rev, mask) << 8 | span if span < 256 else SKIP, i, z))
		i += 1
	return runs


def _window_minima(kms, w):
	"""[frozenset of the k-mers with the smallest x] for each window of w consecutive k-mers of kms;
	a run with 1..w-1 k-mers has one window holding all of them."""
	if not kms:
		return []
	if len(kms) < w:
		m = min(q.x for q in kms)
		return [frozenset(q for q in kms if q.x == m and m != SKIP)]
	where = collections.defaultdict(list)
	for j, q in enumerate(kms):
		where[q.x].append(j)
	out, dq = [], collections.deque()
	for j, q in enumerate(kms):
		while dq and kms[dq[-1]].x >= q.x:
			dq.pop()
		dq.append(j)
		if dq[0] <= j - w:
			dq.popleft()
		if j >= w - 1:
			m = kms[dq[0]].x
			if m == SKIP:
				out.append(frozenset())
				continue
			idx = where[m]
			lo = bisect.bisect_left(idx, j - w + 1)
			hi = bisect.bisect_right(idx, j)
			out.append(frozenset(kms[t] for t in idx[lo:hi]))
	return out


class SketchRef(object):
	"""Reference (w,k)-minimizers of one sequence.

	runs[r] are the k-mers of the r-th A/C/G/T/U run and windows[r] the
	minimum-x sets of its windows of w consecutive k-mers. A run with fewer
	than w k-mers has one short window; its minimum is the run's minimizer.
	"""

	def __init__(self, seq, w, k, hpc=False):
		self.seq, self.w, self.k, self.hpc = seq, w, k, hpc
		self.runs = kmer_runs(seq, k, hpc)
		self.windows = [_window_minima(r, w) for r in self.runs]

	@property
	def has_ties(self):
		"""True if some window's minimum occurs more than once in it."""
		return any(len(s) > 1 for ws in self.windows for s in ws)

	def minimizers(self):
		"""Every window minimum of every run, sorted by position: the documented rule."""
		u = set()
		for ws in self.windows:
			for s in ws:
				u |= s
		return sorted(u, key=lambda q: q.pos)

	def kmer_at(self):
		"""{pos: Kmer} for every k-mer of the sequence."""
		return dict((q.pos, q) for r in self.runs for q in r if q.x != SKIP)


def expected_index(seqs, w, k, hpc=False):
	"""{hash: sorted tuple of rid << 32 | pos << 1 | strand} that an index of seqs must hold, from the
	documented minimizer rule (for sequences without ambiguous bases and without ties in a window)."""
	out = collections.defaultdict(list)
	for rid, s in enumerate(seqs):
		ref = SketchRef(s, w, k, hpc)
		if len(ref.runs) > 1 or ref.has_ties:
			raise ValueError("sequence %d has an ambiguous base or a tied window; the exact reference does not apply" % rid)
		for q in ref.minimizers():
			out[q.x >> 8].append(rid << 32 | q.pos << 1 | q.strand)
	return dict((h, tuple(sorted(v))) for h, v in out.items())


# ---------------------------------------------------------------------------
# .mmi files
# ---------------------------------------------------------------------------

MMI_MAGIC = b"MMI\x02"  # minimap.h MM_IDX_MAGIC
I_HPC, I_NO_SEQ, I_NO_NAME = 0x1, 0x2, 0x4  # minimap.h MM_I_*


class MmiPart(object):
	"""One index of a .mmi file (a multi-part index holds several back to back).

	w, k, b, flag: the header; names and lengths of the n_seq sequences;
	minimizers: {hash: tuple of positions rid << 32 | pos << 1 | strand,
	sorted}; seq4: the 4-bit base codes (0-3 ACGT, 4 other) per sequence,
	or None with --idx-no-seq; n_singletons: hashes stored once.
	"""

	def __init__(self):
		self.w = self.k = self.b = self.flag = 0
		self.names, self.lengths, self.minimizers, self.seq4 = [], [], {}, None
		self.n_singletons = 0
		self.offset = self.end = 0  # byte range in the file

	@property
	def n_seq(self):
		return len(self.names)

	def seq_codes(self, rid):
		return self.seq4[rid] if self.seq4 is not None else None


def _parse_part(data, off):
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
	mins = {}
	for bi in range(1 << p.b):
		n, = struct.unpack_from("<i", data, off)
		off += 4
		pos = struct.unpack_from("<%dQ" % n, data, off) if n else ()
		off += 8 * n
		size, = struct.unpack_from("<I", data, off)
		off += 4
		used = 0
		for j in range(size):
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
	if not p.flag & I_NO_SEQ:
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


def read_mmi(path):
	"""All parts of a .mmi file, in order; raises ValueError on a malformed file or trailing bytes."""
	with open(str(path), "rb") as f:
		data = f.read()
	parts, off = [], 0
	while off < len(data):
		p, off = _parse_part(data, off)
		parts.append(p)
	if not parts:
		raise ValueError("empty .mmi file")
	return parts


def seq_codes(seq):
	"""4-bit codes minimap2 stores for a target sequence: A/C/G/T/U in either case 0..3, else 4."""
	return bytes(_NT4.get(c, 4) for c in seq)


# ---------------------------------------------------------------------------
# SDUST
# ---------------------------------------------------------------------------

_SD = {"A": 0, "C": 1, "G": 2, "T": 3, "a": 0, "c": 1, "g": 2, "t": 3}  # sdust.c's own table: U is not a base


def _better(a, b):
	"""True if score a = (r, l) is larger than b (r / l compared exactly)."""
	return a[0] * b[1] > b[0] * a[1]


def sdust_reference(seq, T=20, W=64):
	"""Low-complexity intervals [(start, end)] of seq by symmetric DUST.

	A string of n >= 4 bases holds n - 2 triplets with counts c_t; its score
	is sum_t c_t (c_t - 1) / 2 divided by n - 3. An interval is perfect if
	its score exceeds T / 10 and no sub-interval scores higher. The result is
	the union of perfect intervals of at most W bases inside A/C/G/T runs,
	with overlapping or adjacent intervals merged. O(len * W).
	"""
	ivs = []
	n, i = len(seq), 0
	while i < n:
		if seq[i] not in _SD:
			i += 1
			continue
		j = i
		while j < n and seq[j] in _SD:
			j += 1
		piece = seq[i:j]
		trip = [_SD[piece[p]] << 4 | _SD[piece[p + 1]] << 2 | _SD[piece[p + 2]] for p in range(len(piece) - 2)]
		m = len(trip)
		prev = {}  # b -> best score among intervals inside [a + 1, b] (triplet indices)
		for a in range(m - 1, -1, -1):
			cnt = {trip[a]: 1}
			r, cur = 0, {}
			for b in range(a + 1, min(m, a + W - 2)):
				t = trip[b]
				c = cnt.get(t, 0)
				r += c
				cnt[t] = c + 1
				me = (r, b - a)
				sub = None
				if b - a >= 2:
					c1, c2 = prev[b], cur[b - 1]
					sub = c1 if _better(c1, c2) else c2
				if r * 10 > T * (b - a) and (sub is None or not _better(sub, me)):
					ivs.append((i + a, i + b + 3))
				cur[b] = me if sub is None or _better(me, sub) else sub
			prev = cur
		i = j
	ivs.sort()
	out = []
	for s, e in ivs:
		if out and s <= out[-1][1]:
			out[-1][1] = max(out[-1][1], e)
		else:
			out.append([s, e])
	return [tuple(x) for x in out]


# ---------------------------------------------------------------------------
# Sequences
# ---------------------------------------------------------------------------


def random_dna(seed, n, alphabet="ACGT"):
	return H.Rng(seed).dna(n, alphabet)


def write_fa(path, records, line_width=None):
	"""FASTA from [(name, seq)]; returns the path."""
	return H.write_fasta(str(path), records, line_width=line_width)


def low_complexity(rng, n, units=("AAG", "ACT", "AGT")):
	"""n bases of a random concatenation of short units: low complexity by DUST, yet with many
	distinct 15-mers (unlike a perfect microsatellite)."""
	out = []
	while sum(len(u) for u in out) < n:
		out.append(rng.choice(units))
	return "".join(out)[:n]


def file_sha256(path):
	import hashlib
	with open(str(path), "rb") as f:
		return hashlib.sha256(f.read()).hexdigest()


def clean_genome(seed=7):
	"""[(name, seq)]: random contigs without ambiguous bases or tied windows for the presets' (w, k), so
	the exact minimizer reference applies to every sequence: 20 kbp, 8 kbp with a lowercase stretch, 3 kbp
	written with U for T, 30 bp (fewer k-mers than some windows), 20 bp and 12 bp (shorter than k=15/19/21/25)."""
	rng = H.Rng(seed)
	b = rng.dna(8000)
	c = rng.dna(3000)
	return [("ctgA", rng.dna(20000)), ("ctgB", b[:3000] + b[3000:5000].lower() + b[5000:]), ("rnaC", c.replace("T", "U")),
		("short30", rng.dna(30)), ("short20", rng.dna(20)), ("tiny12", rng.dna(12))]


def stored_minimizers(part):
	"""{rid: set of (hash, pos, strand)} held by a decoded index part."""
	out = collections.defaultdict(set)
	for h, ys in part.minimizers.items():
		for y in ys:
			out[y >> 32].add((h, (y & 0xffffffff) >> 1, y & 1))
	return out


def index_problems(part, records, w, k, hpc):
	"""Problems with the minimizers an index part holds for sequences [(name, seq)] that may contain
	ambiguous bases and tied windows (where only properties are checked, see test_a_sketch.py):
	every stored entry is a real k-mer at its position; in each run, after its first window, the stored
	k-mers are exactly the window minima; every window has a stored minimum except the windows whose
	minimum is that of the last window of a run followed by an ambiguous base; runs with fewer than w
	k-mers are not checked."""
	probs = []
	stored = stored_minimizers(part)
	for rid, (name, seq) in enumerate(records):
		ref = SketchRef(seq, w, k, hpc)
		at = ref.kmer_at()
		got = stored.get(rid, set())
		for h, pos, z in got:
			q = at.get(pos)
			if q is None or q.x >> 8 != h or q.strand != z:
				probs.append("%s: stored hash %d at %d strand %d is not the k-mer there" % (name, h, pos, z))
		gpos = set(pos for _, pos, _ in got)
		for r, (run, ws) in enumerate(zip(ref.runs, ref.windows)):
			if len(run) < w:
				continue
			cut, end = run[w - 1].pos, run[-1].pos
			mins = set()
			for s in ws:
				mins |= s
			want = set(q.pos for q in mins if q.pos > cut)
			have = set(p for p in gpos if cut < p <= end)
			if r < len(ref.runs) - 1:
				lost = set(q.pos for q in ws[-1])
				extra = have - want
				missing = want - have - lost
			else:
				extra, missing = have - want, want - have
			if extra or missing:
				probs.append("%s run %d: %d stored k-mers are not window minima, %d window minima are missing (e.g. %s)" %
					(name, r, len(extra), len(missing), sorted(extra or missing)[:3]))
			first = ws[0]
			if first and not (set(q.pos for q in first) & gpos) and not (r < len(ref.runs) - 1 and first == ws[-1]):
				probs.append("%s run %d: no minimizer stored for the first window" % (name, r))
	return probs
