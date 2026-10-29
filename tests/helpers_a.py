"""Helpers for the input, sketching and index tests (test_a_*.py).

Standard library only.

	hash64, kmer_runs,          the reference (w,k)-minimizer scan and the .mmi
	SketchRef, read_mmi,        reader, aliases of the shared ones in
	MmiPart, expected_index     minimap2_helpers (kmer_hash, kmer_runs, SketchRef,
	                            read_mmi, ...)
	sdust_reference             the symmetric DUST rule of Morgulis et al.
	                            (2006): the union of perfect intervals of at
	                            most W bases
	random_dna, write_fa        small sequence utilities
"""

import collections

import minimap2_helpers as H

# ---------------------------------------------------------------------------
# Minimizers and .mmi files: the shared reference in minimap2_helpers, under the names the
# test_a_*.py tests use
# ---------------------------------------------------------------------------

_NT4 = H._SKETCH_NT4  # sketch.c's seq_nt4_table: A/C/G/T/U in either case -> 0..3, anything else -> 4


def hash64(key, mask):
	"""sketch.c's invertible integer hash of a 2k-bit k-mer, restricted to `mask` = 4^k - 1."""
	return H.kmer_hash(key, mask.bit_length() // 2)


Kmer = H.Kmer  # (x, pos, strand), x = hash << 8 | span
SKIP = H.KMER_SKIP
kmer_runs = H.kmer_runs
_window_minima = H._window_minima
SketchRef = H.SketchRef
expected_index = H.expected_index
MMI_MAGIC = H.MMI_MAGIC
I_HPC, I_NO_SEQ, I_NO_NAME = H.MM_I_HPC, H.MM_I_NO_SEQ, H.MM_I_NO_NAME
MmiPart = H.MmiPart
_parse_part = H._parse_mmi_part
read_mmi = H.read_mmi


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
