"""Area A: low-complexity masking with SDUST (sdust.c), the `sdust` program of
`make extra`, and -T (SDUST masking of query minimizers, map.c).

Expected intervals come from helpers_a.sdust_reference, a direct
implementation of the symmetric DUST definition (Morgulis et al. 2006): the
union of perfect intervals of at most W bases, where a perfect interval
scores above T/10 and no sub-interval scores higher. Ambiguous bases appear
only at the ends of test sequences: base minimap2 carries the window state
across an ambiguous base inside a sequence (known-bugs/a-sdust-state-across-n.md).
"""

import pytest

import helpers_a as A
import minimap2_helpers as H


def _lcr_sequences(seed, n=24):
	"""Random sequences with planted low-complexity tracts: homopolymers, imperfect di- to hexanucleotide
	repeats, two-letter random sequence, and runs of N at either end of some sequences."""
	rng = H.Rng(seed)
	out = []
	for i in range(n):
		parts = []
		for j in range(rng.randint(1, 6)):
			kind = rng.below(4)
			if kind == 0:
				parts.append(rng.dna(rng.randint(5, 150)))
			elif kind == 1:
				unit = rng.dna(rng.randint(1, 6))
				s = (unit * 200)[:rng.randint(8, 120)]
				parts.append(H.mutate(rng, s, 0.08))
			elif kind == 2:
				two = "AC" if rng.below(2) else "GT"
				parts.append("".join(rng.choice(two) for _ in range(rng.randint(10, 80))))
			else:
				parts.append(A.low_complexity(rng, rng.randint(20, 120)))
		s = "".join(parts)
		if i % 5 == 1:
			s = "N" * rng.randint(1, 30) + s
		if i % 5 == 2:
			s = s + "n" * rng.randint(1, 30)
		if i % 7 == 3:
			s = s.lower()
		out.append(("s%d" % i, s))
	out.append(("tiny", "ACG"))
	out.append(("allA", "A" * 200))
	return out


@pytest.fixture(scope="module")
def lcr(tmp_path_factory):
	d = tmp_path_factory.mktemp("sdust")
	recs = _lcr_sequences(501)
	return recs, A.write_fa(d / "lcr.fa", recs, line_width=60)


def _program_intervals(text):
	out = {}
	for line in text.splitlines():
		name, st, en = line.split("\t")
		out.setdefault(name, []).append((int(st), int(en)))
	return out


@pytest.mark.parametrize("opts,T,W", [([], 20, 64), (["-t", "15", "-w", "32"], 15, 32), (["-w", "100", "-t", "30"], 30, 100),
	(["-t", "10"], 10, 64)])
def test_sdust_program_matches_reference(sdust_bin, lcr, opts, T, W):
	"""`sdust [-w W] [-t T] in.fa` prints name, start and end of each masked interval (0-based, end
	exclusive), equal to the reference SDUST with the defaults -w 64 -t 20 or the given values."""
	recs, fa = lcr
	r = H.run_process([sdust_bin] + opts + [fa], check=True)
	got = _program_intervals(r.stdout)
	n_masked = 0
	for name, s in recs:
		want = A.sdust_reference(s, T, W)
		assert got.get(name, []) == want, (name, T, W)
		n_masked += len(want)
	assert n_masked >= len(recs)


def test_sdust_api_matches_reference(lcr, tmp_path):
	"""sdust() (length from the NUL) and sdust_core() with one buffer reused across sequences return the
	same intervals, equal to the reference."""
	recs, fa = lcr
	for T, W in ((20, 64), (25, 48)):
		out = str(tmp_path / "sd.tsv")
		H.assert_c_case("drv_a_sketch", "sdust", fa, T, W, out)
		got = {}
		for line in H.read_text(out).splitlines():
			i, st, en = (int(x) for x in line.split("\t"))
			got.setdefault(i, []).append((st, en))
		for i, (name, s) in enumerate(recs):
			assert got.get(i, []) == A.sdust_reference(s, T, W), (name, T, W)


def test_sdust_usage_stdin_and_gzip(sdust_bin, lcr, tmp_path):
	"""Without an input file, sdust prints its usage, naming -w and -t, to stderr and exits 1; '-' reads
	standard input and a gzip file is read like the plain one."""
	r = H.run_process([sdust_bin])
	assert r.returncode == 1 and r.stdout == "" and "-w" in r.stderr and "-t" in r.stderr
	recs, fa = lcr
	plain = H.run_process([sdust_bin, fa], check=True).stdout
	with open(fa, "rb") as f:
		assert H.run_process([sdust_bin, "-"], stdin=f.read(), check=True).stdout == plain
	gz = A.write_fa(tmp_path / "lcr.fa.gz", recs, line_width=60)
	assert H.run_process([sdust_bin, gz], check=True).stdout == plain


@pytest.fixture(scope="module")
def lcr_target(tmp_path_factory):
	"""A 21.2 kbp target: random sequence around a 600 bp low-complexity tract made of the units AAG, ACT
	and AGT (DUST masks it, yet its 15-mers are diverse enough to seed), and reads that start inside or
	at the start of the tract and run into unique sequence."""
	d = tmp_path_factory.mktemp("lcrt")
	rng = H.Rng(511)
	lc = A.low_complexity(rng, 600)
	seq = rng.dna(10000) + lc + rng.dna(10600)
	fa = A.write_fa(d / "t.fa", [("t", seq)])
	reads = [("lcr_then_unique%d" % i, seq[10000 + o:10000 + o + 2100 + 37 * i]) for i, o in enumerate((0, 150, 320))]
	reads.append(("lcr_only", lc[20:580]))
	return seq, fa, reads


@pytest.mark.parametrize("T", [20, 10, 30])
def test_T_masks_low_complexity_query_minimizers(lcr_target, tmp_path, T):
	"""-T T masks query minimizers of which more than half lies in a region SDUST (window 64, threshold
	T) masks on the query. For exact substrings of the target, PAF without -c then has cm:i equal to the
	number of reference minimizers of the read that are not masked, and starts at the first of them; the
	read made only of the tract has no unmasked minimizer and stays unmapped. Without -T nothing is
	masked: every minimizer is on the chain."""
	seq, fa, reads = lcr_target
	rf = H.write_fasta(str(tmp_path / "r.fa"), reads)
	plain = dict((p.qname, p) for p in H.map_reads(fa, rf).paf())
	masked = dict((p.qname, p) for p in H.map_reads(fa, rf, ["-T", str(T)]).paf())
	n_cut = 0
	for name, s in reads:
		mins = A.SketchRef(s, 10, 15).minimizers()
		dust = A.sdust_reference(s, T, 64)
		keep = []
		for q in mins:
			st = q.pos - 14
			overlap = sum(max(0, min(q.pos + 1, e) - max(st, b)) for b, e in dust)
			if overlap <= 7:  # at most half of the 15-mer
				keep.append(q)
		if name == "lcr_only":
			assert name in plain and plain[name].tags["cm"] == len(mins)
			if not keep:
				assert name not in masked
			continue
		assert plain[name].tags["cm"] == len(mins) and plain[name].qstart == mins[0].pos - 14
		p = masked[name]
		assert p.tags["cm"] == len(keep), (name, T)
		assert p.qstart == keep[0].pos - 14
		n_cut += len(keep) < len(mins)
	assert n_cut >= 2
