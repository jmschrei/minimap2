"""Minimizer sketching (sketch.c, mm_sketch).

mm_sketch() runs through tests/c/drv_a_sketch.c and its output is compared
with helpers_a.SketchRef, a brute-force scan written from minimap2.1 ("a
minimizer is the smallest k-mer in a window of w consecutive k-mers"), the
hash in sketch.c, -H's definition (a k-mer of the homopolymer-compressed
sequence) and the layout documented above mm_sketch() (x = hash<<8 | span,
y = rid<<32 | lastPos<<1 | strand).

Two edge behaviours of minimap2 that contradict this definition are kept out
of the exact comparisons: in the first window of a
run, a k-mer tied with the running minimum can be dropped or written when it
is not a window minimum; and the minimizer of the last window before an
ambiguous base is never written.
"""


import pytest

import helpers_a as A
import minimap2_helpers as H

# (w, k, hpc) of every preset (options.c, minimap2.1) and edge values
PRESET_WKH = [(10, 15, False), (5, 15, False), (10, 19, True), (5, 19, True), (19, 19, False), (51, 25, False),
	(10, 19, False), (11, 21, False)]
EDGE_WKH = [(1, 15, False), (2, 11, False), (255, 15, False), (3, 27, False), (7, 27, True), (1, 13, True)]


def _dump(tmp_path, records, w, k, hpc, tag="d"):
	"""Run drv_a_sketch dump on [(name, seq)]; {seq index: [(x, pos, strand)]} in output order."""
	fa = A.write_fa(tmp_path / ("%s.fa" % tag), records)
	out = str(tmp_path / ("%s.tsv" % tag))
	H.assert_c_case("drv_a_sketch", "dump", fa, w, k, int(hpc), out)
	got = dict((i, []) for i in range(len(records)))
	with open(out) as f:
		for line in f:
			rid, pos, strand, h, span = (int(v) for v in line.split("\t"))
			got[rid].append(A.Kmer(h << 8 | span, pos, strand))
	return got


def _random_records(seed, w, k, hpc=False):
	"""Sequences around the length thresholds of the rule (no k-mer, one k-mer, fewer than w k-mers,
	exactly one window) and longer ones, without ties in any window."""
	rng = H.Rng(seed)
	lengths = [k - 1, k, k + 1, k + w - 2, k + w - 1, k + w, 2 * (k + w)] + [rng.randint(300, 3000) for _ in range(6)]
	recs = []
	for i, L in enumerate(lengths):
		if L < 1:
			continue
		for attempt in range(20):
			s = rng.dna(L)
			if not A.SketchRef(s, w, k, hpc).has_ties:
				break
		recs.append(("s%d" % i, s))
	return recs


@pytest.mark.parametrize("w,k,hpc", PRESET_WKH + EDGE_WKH)
def test_minimizers_equal_the_reference_scan(tmp_path, w, k, hpc):
	"""On random sequences without ties in a window, mm_sketch() returns exactly the window minima of the
	reference scan, sorted by position, with each k-mer's hash, span (k, or the HPC run lengths), last
	position and strand: none for a sequence shorter than k, the smallest k-mer for a sequence with fewer
	than w k-mers, and the union of window minima otherwise."""
	recs = _random_records(1000 + 7 * w + k + hpc, w, k, hpc)
	got = _dump(tmp_path, recs, w, k, hpc)
	n_checked = 0
	for i, (name, s) in enumerate(recs):
		ref = A.SketchRef(s, w, k, hpc)
		if ref.has_ties:
			continue
		n_checked += 1
		assert got[i] == ref.minimizers(), (name, len(s), w, k, hpc)
		n_kmers = sum(len(r) for r in ref.runs)
		if n_kmers == 0:
			assert got[i] == []
		elif n_kmers < w:
			assert len(got[i]) == 1
	assert n_checked == len(recs)


def test_hpc_minimizers_on_homopolymer_rich_sequence(tmp_path):
	"""-H: on a sequence of homopolymer runs (1-9 bases), k-mers are taken on the compressed sequence,
	each minimizer's span is the length of its k runs, its position the last base of its last run, and
	the minimizers equal those of the compressed sequence with positions mapped back."""
	rng = H.Rng(31)
	runs = []
	while len(runs) < 900:
		c = rng.choice("ACGT")
		if runs and runs[-1][0] == c:
			continue
		runs.append((c, rng.randint(1, 9)))
	seq = "".join(c * n for c, n in runs)
	comp = "".join(c for c, _ in runs)
	ends, p = [], -1
	for _, n in runs:
		p += n
		ends.append(p)
	for w, k in ((10, 19), (5, 19), (10, 15)):
		got = _dump(tmp_path, [("hp", seq)], w, k, True, tag="hp%d_%d" % (w, k))[0]
		ref = A.SketchRef(seq, w, k, True)
		assert got == ref.minimizers()
		plain = A.SketchRef(comp, w, k, False)  # the compressed sequence without -H
		assert [(q.x >> 8, ends[q.pos], q.strand) for q in plain.minimizers()] == [(q.x >> 8, q.pos, q.strand) for q in got]
		for q in got:
			j = ends.index(q.pos)
			assert q.x & 0xff == sum(n for _, n in runs[j - k + 1:j + 1])
		assert any(q.x & 0xff > k for q in got)


def test_minimizers_with_tied_kmers(tmp_path):
	"""Low-complexity sequence where the smallest k-mer occurs several times in a window: after the first
	window every occurrence of a window minimum is a minimizer and nothing else is; in the first window
	every k-mer written is a real k-mer of the sequence and the rightmost minimum is written, and when the
	window's last k-mer is not its minimum, every occurrence of the minimum is written."""
	rng = H.Rng(41)
	recs = []
	for i in range(8):
		unit = rng.dna(rng.randint(1, 4))
		s = (unit * 400)[:rng.randint(100, 600)]
		s = H.mutate(rng, s, 0.05) + rng.dna(rng.randint(0, 200))
		recs.append(("t%d" % i, s))
	n_tied = n_first = 0
	for w, k, hpc in ((10, 15, False), (5, 11, False), (19, 19, False), (10, 15, True)):
		got = _dump(tmp_path, recs, w, k, hpc, tag="tie%d_%d_%d" % (w, k, hpc))
		for i, (name, s) in enumerate(recs):
			ref = A.SketchRef(s, w, k, hpc)
			n_tied += ref.has_ties
			km = ref.runs[0]
			at = ref.kmer_at()
			if len(km) < w:
				continue
			cut = km[w - 1].pos  # the first window's last k-mer
			assert [q for q in got[i] if q.pos > cut] == [q for q in ref.minimizers() if q.pos > cut], (name, w, k, hpc)
			first = [q for q in got[i] if q.pos <= cut]
			assert all(at.get(q.pos) == q for q in first), (name, w, k, hpc)
			m = min(q.x for q in km[:w])
			rightmost = [q for q in km[:w] if q.x == m][-1]
			assert rightmost in got[i], (name, w, k, hpc)
			if km[w - 1].x > m:
				mins = [q for q in km[:w] if q.x == m]
				assert set(mins) <= set(got[i]), (name, w, k, hpc)
				n_first += len(mins) > 1
	assert n_tied >= 10 and n_first >= 3, (n_tied, n_first)


def test_minimizers_around_ambiguous_bases(tmp_path):
	"""N and IUPAC codes end a k-mer run: no minimizer spans one, every minimizer is a window minimum of
	its own run, and every window of every run has a minimizer, except the windows whose minimum is
	that of the last window of a run followed by an ambiguous base (minimap2 never writes that minimizer)."""
	rng = H.Rng(51)
	recs = []
	for i in range(10):
		parts = []
		for j in range(rng.randint(2, 5)):
			parts.append(rng.dna(rng.randint(5, 400)))
			parts.append(rng.choice(["N", "NNNNN", "R", "y", "N" * 60, "-", "X"]))
		parts.append(rng.dna(600))  # a long final run
		recs.append(("n%d" % i, "".join(parts)))
	for w, k, hpc in ((10, 15, False), (5, 15, False), (11, 21, False), (10, 19, True), (51, 25, False)):
		got = _dump(tmp_path, recs, w, k, hpc, tag="amb%d_%d_%d" % (w, k, hpc))
		for i, (name, s) in enumerate(recs):
			ref = A.SketchRef(s, w, k, hpc)
			if ref.has_ties:
				continue
			gs = set(got[i])
			allmin = set()
			for ws in ref.windows:
				for x in ws:
					allmin |= x
			assert gs <= allmin, (name, sorted(gs - allmin, key=lambda q: q.pos)[:3])
			for r, ws in enumerate(ref.windows):
				last = ws[-1] if ws else None
				for x in ws:
					if not x or x & gs:
						continue
					assert r < len(ref.windows) - 1 and x == last, (name, w, k, hpc, r)
			for q in got[i]:  # the k-mer ends at q.pos and covers only A/C/G/T/U
				assert all(c in "ACGTUacgtu" for c in s[q.pos - (q.x & 0xff) + 1:q.pos + 1])


def test_case_and_u_give_the_same_minimizers(tmp_path):
	"""Lowercase bases and U (RNA) sketch as their uppercase DNA bases (sketch.c's seq_nt4_table)."""
	rng = H.Rng(61)
	s = rng.dna(3000)
	lower = s[:1000].lower() + s[1000:]
	rna = s[:2000].replace("T", "U") + s[2000:].lower().replace("t", "u")
	got = _dump(tmp_path, [("a", s), ("b", lower), ("c", rna)], 10, 15, False)
	assert got[0] == got[1] == got[2] == A.SketchRef(s, 10, 15).minimizers()


def test_sketch_with_kalloc_pool_and_appending(tmp_path):
	"""mm_sketch() gives the same result with a kalloc pool as with malloc (km NULL), stores rid in the
	high 32 bits of y, and appends to a vector that is not empty."""
	recs = _random_records(71, 10, 15) + [("lc", "ACACACACAC" * 30 + "GATTACA" * 20)]
	for w, k, hpc in ((10, 15, 0), (5, 19, 1)):
		H.assert_c_case("drv_a_sketch", "km_append", A.write_fa(tmp_path / "km.fa", recs), w, k, hpc)


def _palindrome(rng, half):
	x = rng.dna(half)
	return x + H.revcomp(x)


@pytest.mark.parametrize("k", [16, 14])
def test_palindromic_kmers_are_never_minimizers(tmp_path, k):
	"""Even k: a k-mer equal to its reverse complement has no strand, and mm_sketch() skips it (the
	comment in sketch.c: symmetric k-mers are skipped as their strand is unknown); it takes no place in
	any window. With w 1, where every other k-mer is a minimizer, no minimizer ends at a planted
	palindromic k-mer, and the k-mers just before and after each palindrome are minimizers."""
	rng = H.Rng(71 + k)
	parts, pal_ends = [], []
	pos = 0
	for i in range(12):
		flank = rng.dna(rng.randint(40, 90))
		p = _palindrome(rng, k // 2)
		parts += [flank, p]
		pos += len(flank) + k
		pal_ends.append(pos - 1)
	parts.append(rng.dna(60))
	seq = "".join(parts)
	got = _dump(tmp_path, [("pal", seq)], 1, k, False, tag="pal%d" % k)[0]
	ends = set(q.pos for q in got)
	for e in pal_ends:
		assert seq[e - k + 1:e + 1] == H.revcomp(seq[e - k + 1:e + 1])
		assert e not in ends, (k, e)
		assert e - 1 in ends and e + 1 in ends, (k, e)
	assert all(seq[q.pos - k + 1:q.pos + 1] != H.revcomp(seq[q.pos - k + 1:q.pos + 1]) for q in got)


def test_hpc_kmers_spanning_256_bases_or_more(tmp_path):
	"""-H stores a k-mer's span in 8 bits (x = hash<<8 | span, documented above mm_sketch()), so a
	homopolymer-compressed k-mer spanning 256 bases or more is never a minimizer, while one spanning 255
	is. Runs of 10-25 bases make the spans of k = 15 compressed k-mers range around 255; with w 1 every
	k-mer spanning at most 255 bases is a minimizer, and with w 5 the windows follow the reference scan
	(such long k-mers keep their place in a window without being a minimum)."""
	rng = H.Rng(83)
	runs = []
	while len(runs) < 400:
		c = rng.choice("ACGT")
		if runs and runs[-1][0] == c:
			continue
		runs.append((c, rng.randint(10, 25)))
	seq = "".join(c * n for c, n in runs)
	for w in (1, 5):
		ref = A.SketchRef(seq, w, 15, True)
		if ref.has_ties:
			continue
		got = _dump(tmp_path, [("hp", seq)], w, 15, True, tag="hp256_%d" % w)[0]
		assert got == ref.minimizers(), w
	spans = [sum(n for _, n in runs[j - 14:j + 1]) for j in range(14, len(runs))]
	assert 255 in spans and any(s >= 256 for s in spans) and any(s < 255 for s in spans)
	got1 = _dump(tmp_path, [("hp", seq)], 1, 15, True, tag="hp256_w1")[0]
	assert 255 in [q.x & 0xff for q in got1]
	assert len(got1) == sum(1 for s in spans if s < 256)
