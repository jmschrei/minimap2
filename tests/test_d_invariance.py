"""Output that must not depend on how the work is split: thread counts (-t), mini-batch sizes (-K),
two I/O threads (-2), the allocator (--no-kalloc, --cap-kalloc), and the hash seed (--seed) for
reads without equal hits.

minimap2.1: -t sets the worker threads, -K the bases per mini-batch ("A large NUM helps load
balancing ... at the cost of increased memory"), -2 an extra I/O thread, --no-kalloc the libc
allocator, --cap-kalloc when thread-local memory is freed; none of them is documented to change
the alignments, and tests/README.md records byte-identical output across them.
"""

import pytest

import helpers_d as D
import minimap2_helpers as H


@pytest.fixture(scope="module")
def inputs(synth, synth_tx, tmp_path_factory):
	d = tmp_path_factory.mktemp("d_inv")
	m = D.mixed_reads(synth.genome)
	fq = D.write_mixed(d / "mixed.fq", m)
	pairs = H.simulate_pairs(synth.genome, 150, seed=701, insert_mean=350, insert_sd=40)
	p1, p2 = H.write_pairs(str(d / "p1.fq"), str(d / "p2.fq"), pairs)
	tx = H.simulate_spliced_reads(synth_tx.genome, 40, seed=702, profile="hifi")
	txq = H.write_fastq(str(d / "tx.fq"), tx)
	return dict(mixed=(m, fq), pairs=(p1, p2), tx=txq)


def _body(r):
	return H.strip_pg(r.stdout)


SPLITS = [["-t", "1"], [], ["-t", "2", "-K", "2k"], ["-t", "4", "-K", "10k"], ["-t", "4", "-2"], ["-t", "8", "-K", "1"],
	["-t", "3", "-K", "3k", "-2"]]


@pytest.mark.thread_invariance
def test_long_reads_sam(synth, inputs):
	"""map-ont SAM with --cs and secondary and supplementary records: the same bytes for -t 1, 2,
	3 (the CLI default, no -t), 4 and 8, -K 1 (a batch per read), 2k, 3k, 10k and with -2."""
	fq = inputs["mixed"][1]
	outs = []
	for args in SPLITS:
		argv = ["-a", "--cs"] + args
		r = H.map_reads(synth.fasta, fq, argv, threads=None) if args else H.run_mm2(["-a", "--cs", synth.fasta, fq], check=True)
		outs.append(_body(r))
	assert len(D.body_lines(outs[0])) > 30
	for a, o in zip(SPLITS[1:], outs[1:]):
		assert o == outs[0], a


@pytest.mark.thread_invariance
def test_long_reads_paf(synth, inputs):
	"""PAF with -c --ds and --paf-no-hit: the same bytes for -t 1 and -t 4 -K 3k -2."""
	fq = inputs["mixed"][1]
	a = H.map_reads(synth.fasta, fq, ["-c", "--ds", "--paf-no-hit", "-t", "1"]).stdout
	b = H.map_reads(synth.fasta, fq, ["-c", "--ds", "--paf-no-hit", "-t", "4", "-K", "3k", "-2"]).stdout
	assert a == b and len(D.lines(a)) > 30


@pytest.mark.thread_invariance
def test_paired_reads(synth, index_of, inputs):
	"""sr on read pairs from two files (SAM with --MD) and the same pairs interleaved in one file:
	the same bytes for -t 1, -t 4 -K 5k and -t 3 -K 1k (sr turns on two I/O threads)."""
	p1, p2 = inputs["pairs"]
	idx = index_of(synth.genome, "sr")
	outs = [_body(H.map_reads(idx, [p1, p2], ["-a", "--MD"] + a, preset="sr"))
		for a in (["-t", "1"], ["-t", "4", "-K", "5k"], ["-t", "3", "-K", "1k"])]
	assert outs[0] == outs[1] == outs[2] and len(D.body_lines(outs[0])) >= 300


@pytest.mark.thread_invariance
def test_spliced_reads(synth_tx, inputs):
	"""splice:hq SAM with --cs: the same bytes for -t 1 and -t 4 -K 2k."""
	a = _body(H.map_reads(synth_tx.fasta, inputs["tx"], ["-a", "--cs", "-t", "1"], preset="splice:hq"))
	b = _body(H.map_reads(synth_tx.fasta, inputs["tx"], ["-a", "--cs", "-t", "4", "-K", "2k"], preset="splice:hq"))
	assert a == b and len(D.body_lines(a)) >= 40


@pytest.mark.parametrize("args,threads", [(["--no-kalloc"], 2), (["--no-kalloc"], 1), (["--cap-kalloc", "1k"], 2),
	(["--cap-kalloc", "0"], 2)])
def test_allocator_options(synth, inputs, args, threads):
	"""--no-kalloc (libc allocator) and --cap-kalloc (1k: thread-local memory is reset after nearly
	every read; 0: never) leave the alignments unchanged (minimap2.1 describes both as memory
	controls)."""
	fq = inputs["mixed"][1]
	base = H.map_reads(synth.fasta, fq, ["-a", "--cs"], threads=threads)
	r = H.map_reads(synth.fasta, fq, ["-a", "--cs"] + args, threads=threads)
	assert _body(r) == _body(base)


def test_seed(synth, inputs):
	"""--seed (minimap2.1: integer seed for randomizing equally best hits, hashed with the read
	name; default 11): --seed 11 gives the default output; reads with one best hit come out the
	same under every seed; reads copied from an exact repeat unit have equally good hits, and
	across seeds 1-8 their primary moves between the copies (always onto one of them)."""
	m, fq = inputs["mixed"]
	base = H.map_reads(synth.fasta, fq, ["-a"])
	assert _body(H.map_reads(synth.fasta, fq, ["-a", "--seed", "11"])) == _body(base)
	keep = set(m.names("unique") + m.names("ambig") + m.names("chimeric") + m.names("random"))
	base_by = base.sam().by_qname()
	chosen = dict((n, set()) for n in m.names("repeat"))
	for seed in range(1, 9):
		by = H.map_reads(synth.fasta, fq, ["-a", "--seed", str(seed)]).sam().by_qname()
		for n in keep:
			assert [r.line for r in by[n]] == [r.line for r in base_by[n]], (seed, n)
		for n in chosen:
			p = [r for r in by[n] if r.is_primary][0]
			assert (p.rname, p.pos0) in set((c, s) for c, s, e, _ in m.spans[n])
			chosen[n].add((p.rname, p.pos0))
	assert any(len(v) > 1 for v in chosen.values()), chosen


@pytest.fixture(scope="module")
def larger(synth, tmp_path_factory):
	"""A larger input: 2000 sr pairs (600 kbp; N bases, indels, repeats included) and 180 ONT-like and
	HiFi-like reads (about 750 kbp) from the whole standard genome, so that -K splits the input into
	tens of batches and every thread gets reads in each."""
	d = tmp_path_factory.mktemp("d_inv_large")
	pairs = H.simulate_pairs(synth.genome, 2000, seed=711, insert_mean=350, insert_sd=60, n_rate=0.001, indel_rate=0.002)
	p1, p2 = H.write_pairs(str(d / "p1.fq"), str(d / "p2.fq"), pairs)
	lr = H.simulate_long_reads(synth.genome, 120, profile="ont", seed=712, length=(1000, 6000))
	lr += H.simulate_long_reads(synth.genome, 60, profile="hifi", seed=713, length=(2000, 8000))
	return dict(pairs=(p1, p2), long=H.write_fastq(str(d / "long.fq"), lr), n_pairs=len(pairs), n_long=len(lr))


LARGER_SPLITS = [["-t", "8", "-K", "50k"], ["-t", "4", "-K", "5k", "-2"], ["-t", "8", "-K", "1"]]


@pytest.mark.thread_invariance
def test_larger_input(synth, larger):
	"""On the larger input, the same bytes for -t 1 (one batch) and for -t 8 -K 50k, -t 4 -K 5k -2 and
	-t 8 -K 1 (a batch per read or pair): paired sr SAM with --MD, map-ont SAM with --cs and PAF
	with -c --ds."""
	p1, p2 = larger["pairs"]
	runs = [(["-a", "--MD"], [p1, p2], "sr"), (["-a", "--cs"], [larger["long"]], None), (["-c", "--ds"], [larger["long"]], None)]
	for args, reads, preset in runs:
		one = H.map_reads(synth.fasta, reads, args + ["-t", "1"], preset=preset)
		ref = _body(one)
		n = len(D.body_lines(ref))
		assert n >= (2 * larger["n_pairs"] if preset else larger["n_long"]), (args, n)
		for split in LARGER_SPLITS:
			assert _body(H.map_reads(synth.fasta, reads, args + split, preset=preset)) == ref, (args, split)
