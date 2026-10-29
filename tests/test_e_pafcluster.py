"""misc/pafcluster.js: greedy clustering of sequences from all-vs-all PAF.

Five sequences: A-B identical over their whole length, A-C covering 95% of A
at 0.2% divergence, A-D covering 98% of A at 3% divergence, C-E a short
overlap. Each pair is given in both directions, as `minimap2 -x ava-ont
--dual=yes` writes it. pafcluster.js picks the sequence with the largest
summed s1, groups the sequences whose hit to it passes -c (min coverage of
both), -d (max divergence) and -e (max unaligned bases), prints the group
(`SD name size`, one `CL name length` per member, `//`), drops every hit
involving a member and repeats; sequences never grouped are printed as
singletons. Skipped when k8 is not on PATH.
"""

import pytest

import minimap2_helpers as H

pytestmark = pytest.mark.k8


def _hit(n1, l1, s1, e1, n2, l2, s2, e2, score, dv):
	return "\t".join(str(x) for x in (n1, l1, s1, e1, "+", n2, l2, s2, e2, e2 - s2 - 10, e2 - s2, 60, "tp:A:P",
		"s1:i:%d" % score, "dv:f:%.4f" % dv))


PAF = [
	_hit("A", 10000, 0, 10000, "B", 10000, 0, 10000, 9000, 0.001),
	_hit("B", 10000, 0, 10000, "A", 10000, 0, 10000, 9000, 0.001),
	_hit("A", 10000, 250, 9750, "C", 9000, 0, 9000, 8000, 0.002),
	_hit("C", 9000, 0, 9000, "A", 10000, 250, 9750, 8000, 0.002),
	_hit("A", 10000, 0, 9800, "D", 9800, 0, 9800, 7000, 0.03),
	_hit("D", 9800, 0, 9800, "A", 10000, 0, 9800, 7000, 0.03),
	_hit("C", 9000, 0, 4000, "E", 5000, 1000, 5000, 3000, 0.001),
	_hit("E", 5000, 1000, 5000, "C", 9000, 0, 4000, 3000, 0.001),
	"F\t3000\t0\t3000\t-\tA\t10000\t0\t3000\t2990\t3000\t60\ttp:A:P\ts1:i:2000\tdv:f:0.0010",
]
LEN = {"A": 10000, "B": 10000, "C": 9000, "D": 9800, "E": 5000, "G": 10000}


def _groups(*names):
	out = []
	for g in names:
		out.append("SD\t%s\t%d" % (g[0], len(g)))
		out += ["CL\t%s\t%d" % (n, LEN[n]) for n in g]
		out.append("//")
	return out


@pytest.fixture
def paf(tmp_path):
	p = tmp_path / "ava.paf"
	p.write_text("\n".join(PAF) + "\n")
	return str(p)


def _run(*args):
	return H.run_k8("pafcluster.js", [str(a) for a in args], check=True)


def test_default_thresholds(k8, paf):
	"""-c 0.9 -d 0.015 -e 20000: A (largest summed s1) groups B and C; D (3% divergence) and E (no hit to A)
	stay singletons; the reverse-strand line is ignored; the hit count read goes to stderr."""
	r = _run(paf)
	assert r.stdout.splitlines() == _groups(["A", "B", "C"], ["D"], ["E"])
	assert "Read 8 hits" in r.stderr


def test_min_coverage(k8, paf):
	"""-c 0.96 excludes C (95% of A covered); C then leads its own one-member group since its hit to E is
	short, and D and E remain singletons."""
	assert _run("-c", "0.96", paf).stdout.splitlines() == _groups(["A", "B"], ["C"], ["D"], ["E"])


def test_max_divergence(k8, paf):
	"""-d 0.05 admits D's 3% divergence: A, B, C and D form one group."""
	assert _run("-d", "0.05", paf).stdout.splitlines() == _groups(["A", "B", "C", "D"], ["E"])


def test_max_difference(k8, paf):
	"""-e 400: the A-C hit leaves 500 bases of A unaligned, more than 400, so C is not grouped with A."""
	assert _run("-e", "400", paf).stdout.splitlines() == _groups(["A", "B"], ["C"], ["D"], ["E"])


def test_split_hits_are_chained(k8, tmp_path):
	"""Several hits between the same two sequences are chained in order of their start on the first one and
	scored as one: G matches A in two collinear pieces, A[0,4900) and A[5000,10000) (listed out of order),
	each 0.1% divergent. Together they cover 99% of both, so G joins A's group at the default -c 0.9 (each
	piece alone covers about half); at -c 0.995 it stays a singleton. Divergence is averaged over the
	aligned length."""
	g = [_hit("A", 10000, 5000, 10000, "G", 10000, 5000, 10000, 4000, 0.001),
		_hit("A", 10000, 0, 4900, "G", 10000, 0, 4900, 3900, 0.001),
		_hit("G", 10000, 5000, 10000, "A", 10000, 5000, 10000, 4000, 0.001),
		_hit("G", 10000, 0, 4900, "A", 10000, 0, 4900, 3900, 0.001)]
	p = tmp_path / "ava.paf"
	p.write_text("\n".join(PAF + g) + "\n")
	r = _run(p)
	assert r.stdout.splitlines() == _groups(["A", "B", "C", "G"], ["D"], ["E"]) and "Read 12 hits" in r.stderr
	assert _run("-c", "0.995", p).stdout.splitlines() == _groups(["A", "B"], ["C"], ["D"], ["E"], ["G"])


def test_usage(k8):
	"""Without an input file pafcluster.js prints its usage with the option defaults."""
	r = H.run_k8("pafcluster.js", [], check=False)
	assert r.returncode == 0 and r.stdout.startswith("Usage: pafcluster.js [options] <ava.paf>")
	assert "-c FLOAT     min coverage [0.9]" in r.stdout and "-d FLOAT     max divergence [0.015]" in r.stdout
	assert "-e FLOAT     max difference [20000]" in r.stdout
