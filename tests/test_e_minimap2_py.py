"""python/minimap2.py, the example command-line wrapper around mappy (area E).

It prints `name  qlen  str(Alignment)` per alignment: the first twelve
columns are PAF's, so they are compared with `minimap2 -c --no-hash-name
--secondary=yes` (minimap2.py does not pass read names to map()), and every
option is compared with the minimap2 option it stands for. `-m` is parsed but
never reaches mappy (known-bugs/minimap2-py-m-ignored.md) and is not tested.
"""

import pytest

import minimap2_helpers as H

pytestmark = pytest.mark.mappy  # the `mappy` fixture builds the module (or skips without Cython) before the script runs


@pytest.fixture(scope="module")
def wd(tmp_path_factory):
	return tmp_path_factory.mktemp("e_minimap2_py")


def _expected(target, reads, args=(), preset=None, tags=()):
	"""minimap2.py's lines as the command line computes them: PAF columns 1-12, tp:A (P or S), ts:A, cg:Z, then
	the cs/ds/MD tags in `tags` (one command-line run per tag; it writes only one of them per run)."""
	base = H.map_reads(target, reads, ["-c", "--no-hash-name", "--secondary=yes"] + list(args), preset=preset).paf()
	lines = []
	for p in base:
		ts = p.tags.get("ts")
		lines.append(p.line.split("\t")[:12] + ["tp:A:" + ("P" if p.is_primary else "S"),
			"ts:A:" + (ts if ts in ("+", "-") else "."), "cg:Z:" + p.cigar])
	for opt, tag in (("--cs", "cs"), ("--ds", "ds"), ("--MD", "MD")):
		if tag not in tags:
			continue
		run = H.map_reads(target, reads, ["-c", "--no-hash-name", "--secondary=yes", opt] + list(args), preset=preset).paf()
		assert len(run) == len(lines)
		for fields, p in zip(lines, run):
			fields.append("%s:Z:%s" % (tag, p.tags[tag]))
	return lines


def _got(args):
	r = H.run_minimap2_py(args)
	return [l.split("\t") for l in r.stdout.splitlines()]


def test_usage_with_fewer_than_two_arguments(mappy):
	"""Fewer than two positional arguments: the usage text (listing every option) on stdout and exit 1."""
	for args in ([], ["only_one.fa"], ["-c", "only_one.fa"]):
		r = H.run_minimap2_py(args, check=False)
		assert r.returncode == 1
		assert r.stdout.startswith("Usage: minimap2.py [options] <ref.fa>|<ref.mmi> <query.fq>")
		for opt in ("-x STR", "-n INT", "-m INT", "-k INT", "-w INT", "-r INT", "-c ", "-d ", "-M "):
			assert "\n  %s" % opt in r.stdout, opt


def test_default_output_is_paf_like(mappy, synth, wd):
	"""No options: map-ont defaults, one line per alignment (secondaries included), equal to the command line,
	on a FASTQ of HiFi-like reads from both strands and repeats."""
	reads = H.simulate_long_reads(synth.genome, 12, profile="hifi", seed=601)
	fq = H.write_fastq(str(wd / "r.fq"), reads)
	got = _got([synth.fasta, fq])
	assert len(got) >= len(reads)
	assert got == _expected(synth.fasta, fq)


def test_preset_k_w_and_tags(mappy, synth, wd):
	"""-x map-hifi -k 17 -w 9 with -c, -d and -M: the alignments of `minimap2 -x map-hifi -k17 -w9` followed
	by cs:Z, ds:Z and MD:Z equal to the command line's --cs, --ds and --MD values."""
	reads = H.simulate_long_reads(synth.genome, 8, profile="hifi", seed=602, unique=True)
	fa = H.write_fasta(str(wd / "r.fa"), reads)
	got = _got(["-x", "map-hifi", "-k", "17", "-w", "9", "-c", "-d", "-M", synth.fasta, fa])
	assert len(got) >= len(reads)
	assert got == _expected(synth.fasta, fa, ["-k17", "-w9"], preset="map-hifi", tags=("cs", "ds", "MD"))
	only_md = _got(["-x", "map-hifi", "-k", "17", "-w", "9", "-M", synth.fasta, fa])
	assert only_md == [g[:15] + g[17:] for g in got]


def test_min_count(mappy, synth, wd):
	"""-n INT is the minimum number of minimizers on a chain (-n of minimap2): -n 40 drops short reads' chains
	exactly as `minimap2 -n 40` does, and fewer lines remain than without it."""
	reads = H.simulate_long_reads(synth.genome, 30, profile="ont", seed=507, unique=True, length=(300, 1500))
	fq = H.write_fastq(str(wd / "short.fq"), reads)
	got = _got(["-n", "40", synth.fasta, fq])
	assert got == _expected(synth.fasta, fq, ["-n", "40"])
	assert len(got) < len(_got([synth.fasta, fq]))


def test_k_and_w_set_the_seed_density(mappy, synth, wd):
	"""-k and -w are minimap2's: a window of w k-mers keeps about 2/(w+1) of the positions, so with a minimum
	of 40 minimizers per chain, which short reads survive depends on both. `-k 12 -w 20 -n 40` equals
	`minimap2 -k12 -w20 -n40`, and differs from changing either value alone."""
	reads = H.simulate_long_reads(synth.genome, 30, profile="ont", seed=507, unique=True, length=(300, 1500))
	fq = H.write_fastq(str(wd / "short_kw.fq"), reads)
	got = _got(["-k", "12", "-w", "20", "-n", "40", synth.fasta, fq])
	assert got == _expected(synth.fasta, fq, ["-k12", "-w20", "-n40"])
	assert 0 < len(got) < len(reads)
	for alone in (["-k12", "-n40"], ["-w20", "-n40"], ["-k20", "-n40"]):
		assert _expected(synth.fasta, fq, alone) != got, alone


def test_band_width(mappy, synth, wd):
	"""-r INT is the band width (-r of minimap2): with -x sr, 150 bp reads with a 40 bp deletion align in one
	piece with a 40D by default and in two pieces with -r 10, as `minimap2 -x sr -r 10` does."""
	c = synth.genome["chr2"]
	reads = [("del%d" % i, c[st:st + 75] + c[st + 115:st + 190]) for i, st in enumerate((3000, 9000, 15000))]
	fa = H.write_fasta(str(wd / "del.fa"), reads)
	plain = _got(["-x", "sr", synth.fasta, fa])
	assert [l[0] for l in plain] == ["del0", "del1", "del2"] and all("40D" in l[14] for l in plain)
	narrow = _got(["-x", "sr", "-r", "10", synth.fasta, fa])
	assert sorted(l[0] for l in narrow) == ["del0", "del0", "del1", "del1", "del2", "del2"]
	assert narrow == _expected(synth.fasta, fa, ["-r", "10"], preset="sr")
