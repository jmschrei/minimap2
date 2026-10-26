"""Session fixtures shared by every test file.

mm2_bin         path of the binary under test (built out of tree unless $MINIMAP2_BIN)
sdust_bin       the `sdust` program of `make extra` (built from $MINIMAP2_TEST_SRC_DIR)
lite_bin        the `minimap2-lite` program of `make extra` (example.c)
threads         $MINIMAP2_TEST_THREADS (default 1); map_reads() uses it unless told otherwise
index_dir       a session-lived directory for indexes (pass it to index_genome)
index_of        index_of(genome_or_fasta, preset=None, extra_args=()) -> .mmi path, cached
mt              the shipped test/MT-human.fa and test/MT-orang.fa
synth           the standard synthetic genome (make_synth) and its FASTA
synth_tx        the standard transcript genome (make_synth_tx), its FASTA and annotation files
mappy           the mappy module built from the sources under test (skips without Cython)
k8              path of k8 (skips when k8 is not on PATH)
tree_before     build products present in the source tree before the session started
"""

import os

import pytest

import minimap2_build
import minimap2_helpers as H

_TREE_BEFORE = {}


def pytest_configure(config):
	global _TREE_BEFORE
	_TREE_BEFORE = minimap2_build.tree_products()


def pytest_report_header(config):
	lines = ["minimap2 source: %s" % minimap2_build.src_dir()]
	try:
		lines.append("minimap2 binary: %s" % minimap2_build.minimap2_bin())
	except minimap2_build.BuildError as e:
		lines.append("minimap2 binary: BUILD FAILED (%s)" % str(e).splitlines()[0])
	lines.append("CC=%s  MINIMAP2_TEST_MAKE_ARGS=%s  MINIMAP2_TEST_THREADS=%d" % (minimap2_build.compiler(),
		" ".join(minimap2_build.make_args()) or "-", H.test_threads()))
	why = minimap2_build.mappy_missing()
	lines.append("mappy: %s" % ("SKIPPED (%s)" % why.split(";")[0] if why else "built with " + os.path.basename(H.sys.executable)))
	lines.append("k8: %s" % (H.k8_path() or "not on PATH (paftools tests skip)"))
	return lines


@pytest.fixture(scope="session")
def mm2_bin():
	try:
		return minimap2_build.minimap2_bin()
	except minimap2_build.BuildError as e:
		pytest.fail("cannot build minimap2: %s" % e, pytrace=False)


@pytest.fixture(scope="session")
def sdust_bin():
	try:
		return minimap2_build.sdust_bin()
	except minimap2_build.BuildError as e:
		pytest.fail("cannot build sdust (make extra): %s" % e, pytrace=False)


@pytest.fixture(scope="session")
def lite_bin():
	try:
		return minimap2_build.lite_bin()
	except minimap2_build.BuildError as e:
		pytest.fail("cannot build minimap2-lite (make extra): %s" % e, pytrace=False)


@pytest.fixture(scope="session")
def mappy():
	"""The mappy module built out of tree from the sources under test with this interpreter."""
	try:
		return minimap2_build.import_mappy()
	except minimap2_build.MappyUnavailable as e:
		pytest.skip(str(e))
	except minimap2_build.BuildError as e:
		pytest.fail("cannot build mappy: %s" % e, pytrace=False)


@pytest.fixture(scope="session")
def k8():
	p = H.k8_path()
	if p is None:
		pytest.skip("k8 is not on PATH")
	return p


@pytest.fixture(scope="session")
def threads():
	return H.test_threads()


@pytest.fixture(scope="session")
def tree_before():
	return dict(_TREE_BEFORE)


@pytest.fixture(scope="session")
def index_dir(tmp_path_factory):
	return str(tmp_path_factory.mktemp("indexes"))


@pytest.fixture(scope="session")
def index_of(mm2_bin, index_dir):
	def f(genome_or_fasta, preset=None, extra_args=()):
		return H.index_genome(genome_or_fasta, index_dir, preset=preset, extra_args=extra_args, binary=mm2_bin)
	return f


class _Data(object):
	def __init__(self, **kw):
		self.__dict__.update(kw)


@pytest.fixture(scope="session")
def mt():
	"""test/MT-human.fa (MT_human, 16,569 bp) and test/MT-orang.fa (MT_orang, 16,499 bp), the README example."""
	human = os.path.join(H.DATA_DIR, "MT-human.fa")
	orang = os.path.join(H.DATA_DIR, "MT-orang.fa")
	return _Data(human=human, orang=orang, human_seqs=H.read_fasta(human), orang_seqs=H.read_fasta(orang))


def make_synth():
	"""The standard synthetic genome: two contigs (80 and 50 kbp) with two exact-repeat families of
	three 300-1000 bp copies, one inverted repeat, three low-complexity tracts, two N runs and two short
	contigs (20-200 bp). Changing it changes every pin built on it."""
	return H.make_genome(seed=11, contig_lengths=(80000, 50000), repeat_families=2, repeat_copies=3, inversions=1,
		low_complexity=3, n_runs=2, short_contigs=2)


def make_synth_tx():
	"""The standard transcript genome: two contigs (60 and 40 kbp) with 8 genes of 2-6 exons
	(80-400 bp) and GT-AG introns of 80-2000 bp, plus one exact-repeat family. Changing it changes
	every pin built on it."""
	return H.make_genome(seed=12, contig_lengths=(60000, 40000), genes=8, exons=(2, 6), exon_len=(80, 400),
		intron_len=(80, 2000), repeat_families=1, repeat_copies=2)


@pytest.fixture(scope="session")
def synth(index_dir):
	g = make_synth()
	fasta = g.write_fasta(os.path.join(index_dir, "synth.fa"))
	return _Data(genome=g, fasta=fasta)


@pytest.fixture(scope="session")
def synth_tx(index_dir):
	g = make_synth_tx()
	fasta = g.write_fasta(os.path.join(index_dir, "synth_tx.fa"))
	bed12 = H.write_text(os.path.join(index_dir, "synth_tx.bed"), H.genes_bed12(g.genes))
	gtf = H.write_text(os.path.join(index_dir, "synth_tx.gtf"), H.genes_gtf(g.genes))
	return _Data(genome=g, fasta=fasta, bed12=bed12, gtf=gtf)
