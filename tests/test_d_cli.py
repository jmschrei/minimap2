"""The minimap2 command line (main.c, options.c): --version, --help and the usage text, -x and its
errors, unknown options, missing files, -o, -v, option validation (mm_check_opt) with its
boundaries, deprecated options, debugging options, output write errors and the lifted address
space limit.

Expected values come from the usage text, minimap2.1 and the error messages themselves (which
state the rule they enforce). Messages are checked for presence and for the option they name,
never for their wording.
"""

import os
import subprocess
import sys

import pytest

import helpers_d as D
import minimap2_build
import minimap2_helpers as H


@pytest.fixture(scope="module")
def small(synth, tmp_path_factory):
	"""(target FASTA, reads FASTQ, read names): eight unique ONT-like reads."""
	reads = H.simulate_long_reads(synth.genome, 8, profile="ont", seed=601, unique=True, length=(1000, 2500))
	fq = H.write_fastq(str(tmp_path_factory.mktemp("d_cli") / "r.fq"), reads)
	return synth.fasta, fq, [r.name for r in reads]


def test_version(mm2_bin):
	"""-V and --version print MM_VERSION of minimap.h and a newline to stdout and exit 0."""
	for opt in ("-V", "--version"):
		r = H.run_mm2([opt])
		assert r.returncode == 0 and r.stdout == D.mm_version() + "\n" and r.stderr == ""


def test_help_and_usage(mm2_bin):
	"""-h and --help print the usage to stdout and exit 0; without arguments the same text goes to
	stderr and the exit status is 1 (tests/README.md)."""
	h = H.run_mm2(["-h"])
	assert h.returncode == 0 and h.stdout.startswith("Usage: minimap2 [options] <target.fa>|<target.idx>")
	assert h.stderr == ""
	assert H.run_mm2(["--help"]).stdout == h.stdout
	n = H.run_mm2([])
	assert n.returncode == 1 and n.stdout == "" and n.stderr == h.stdout
	for opt in ("-x STR", "-a ", "-o FILE", "-t INT", "-K NUM", "--cs[=STR]", "--MD", "--eqx", "--version"):
		assert opt in h.stdout


def _usage_default(text, opt):
	"""The [default] the usage text shows for an option line such as '    -k INT ...'."""
	for l in text.splitlines():
		if l.strip().startswith(opt + " "):
			return l[l.rindex("[") + 1:l.rindex("]")]
	raise AssertionError("no usage line for %s" % opt)


def test_usage_shows_current_values_and_preset_first(mm2_bin):
	"""The usage text shows each option's current value in brackets: -k 15 and -w 10 by default,
	21 and 11 under -x sr (minimap2.1); -x is applied before other options whatever its position,
	so `-k 17 -x sr` and `-x sr -k 17` both show 17."""
	d = H.run_mm2(["-h"]).stdout
	assert (_usage_default(d, "-k"), _usage_default(d, "-w"), _usage_default(d, "-t")) == ("15", "10", "3")
	sr = H.run_mm2(["-x", "sr", "-h"]).stdout
	assert (_usage_default(sr, "-k"), _usage_default(sr, "-w")) == ("21", "11")
	for args in (["-k", "17", "-x", "sr", "-h"], ["-x", "sr", "-k", "17", "-h"]):
		assert _usage_default(H.run_mm2(args).stdout, "-k") == "17"


def test_preset_applied_first(small):
	"""-x applies its options before the others wherever it appears (minimap2.1): `-k 17 -x map-hifi`
	maps like `-x map-hifi -k 17`, and both differ from plain map-hifi (k 19)."""
	fa, fq, _ = small
	a = H.map_reads(fa, fq, ["-k", "17", "-x", "map-hifi"])
	b = H.map_reads(fa, fq, ["-x", "map-hifi", "-k", "17"])
	c = H.map_reads(fa, fq, ["-x", "map-hifi"])
	assert a.stdout == b.stdout and a.stdout != c.stdout and len(a.paf()) >= 8


@pytest.mark.parametrize("preset", ["foo", "asm7", "", "map_ont", "MAP-ONT"])
def test_unknown_preset(small, preset):
	"""An unknown preset (e.g. -x foo or -x asm7) is an error: exit 1, a message naming it and no output."""
	fa, fq, _ = small
	r = H.run_mm2(["-x", preset, fa, fq])
	assert r.returncode == 1 and r.stdout == ""
	assert r.messages and (preset == "" or any(preset in l for l in r.messages))


@pytest.mark.parametrize("args,named", [(["--no-such-option"], "--no-such-option"), (["-Z"], "-Z"),
	(["--cs-long"], "--cs-long"), (["-k"], None)])
def test_unknown_option_and_missing_argument(small, args, named):
	"""An unknown option, or an option missing its argument, is an error: exit 1 and a message
	before anything is read, before or after the file arguments. When the unknown option comes
	first the message names it (after the files it names another argument:
	known-bugs/d-unknown-option-message.md)."""
	fa, fq, _ = small
	for argv in ((args + [fa, fq], [fa, fq] + args) if named else ([fa, fq] + args,)):
		r = H.run_mm2(argv)
		assert r.returncode == 1 and r.stdout == "" and r.messages
		if named and argv[0] == args[0]:
			assert any(named in l for l in r.messages)


def test_missing_files(small, tmp_path):
	"""A missing target is an error naming the file (exit 1, nothing written); a missing query
	file is an error after the SAM header is written (tests/README.md); a target without a query
	file and without -d is an error."""
	fa, fq, _ = small
	nofile = str(tmp_path / "absent.fa")
	r = H.run_mm2(["-a", nofile, fq])
	assert r.returncode == 1 and r.stdout == "" and any(nofile in l for l in r.messages)
	r = H.run_mm2(["-a", fa, str(tmp_path / "absent.fq")])
	assert r.returncode == 1 and r.messages
	assert r.stdout.startswith("@HD") and D.body_lines(r.stdout) == []
	r = H.run_mm2([fa])
	assert r.returncode == 1 and r.stdout == "" and r.messages


def test_sr_takes_at_most_two_query_files(small):
	"""In the sr mode more than two query files is an error (exit 1, a message, no output)."""
	fa, fq, _ = small
	r = H.run_mm2(["-x", "sr", fa, fq, fq, fq])
	assert r.returncode == 1 and r.stdout == "" and r.messages


def test_output_file(small, tmp_path):
	"""-o FILE writes what stdout would get (the @PG line differs only by the -o argument) and
	leaves stdout empty; -o - writes to stdout; an output path in a missing directory is an
	error naming the path (exit 1)."""
	fa, fq, _ = small
	base = H.map_reads(fa, fq, ["-a"])
	out = str(tmp_path / "o.sam")
	r = H.map_reads(fa, fq, ["-a", "-o", out])
	assert r.stdout == ""
	text = H.read_text(out)
	assert H.strip_pg(text) == H.strip_pg(base.stdout)
	assert H.parse_sam(text).pg[0]["CL"] == "minimap2 " + " ".join(r.argv[1:])
	dash = H.map_reads(fa, fq, ["-a", "-o", "-"])
	assert H.strip_pg(dash.stdout) == H.strip_pg(base.stdout)
	paf = str(tmp_path / "o.paf")
	H.map_reads(fa, fq, ["-o", paf])
	assert H.read_text(paf) == H.map_reads(fa, fq).stdout
	bad = str(tmp_path / "no" / "such" / "dir" / "o.sam")
	r = H.map_reads(fa, fq, ["-a", "-o", bad], check=False)
	assert r.returncode == 1 and any(bad in l for l in r.messages) and r.stdout == ""


def test_verbosity(small):
	"""-v sets the verbosity (mm_verbose: 0 none, 1 errors, 2 warnings, 3 messages, the CLI's
	default 3): [M::] progress lines only at 3, the --cs warning only from 2, mm_check_opt errors
	only from 1; the alignments do not change."""
	fa, fq, _ = small
	outs = []
	for v in ("0", "1", "2", "3"):
		r = H.map_reads(fa, fq, ["-v", v, "--cs=bogus"])
		outs.append(r.stdout)
		has_m = any(l.startswith("[M::") for l in r.stderr.splitlines())
		assert has_m == (v == "3")
		assert any("--cs" in l for l in r.messages) == (v in "23")
		e = H.run_mm2(["-v", v, "-p", "2", "-h"])
		assert e.returncode == 1 and (e.stderr != "") == (v != "0") and e.stdout == ""
	assert len(set(outs)) == 1
	assert H.map_reads(fa, fq, ["--cs"]).stdout == outs[0]


# (arguments, option the message names or None): each violates a rule mm_check_opt states in its message
CHECK_ERRORS = [
	(["-r", "1000,500"], "-r"),
	(["-r", "600,599"], "-r"),
	(["-r", "20001"], "-r"),
	(["--rmq", "-x", "sr"], "--rmq"),
	(["--rmq=yes", "--splice"], "--rmq"),
	(["-x", "lr:hqae", "-x", "splice"], "--rmq"),
	(["--split-prefix", "pfx", "--cs"], "--split-prefix"),
	(["--split-prefix", "pfx", "--MD"], "--split-prefix"),
	(["-k", "0"], "-k"),
	(["-w", "0"], "-w"),
	(["-N", "-1"], "-N"),
	(["-p", "1.01"], "-p"),
	(["-p", "-0.01"], "-p"),
	(["--for-only", "--rev-only"], "--for-only"),
	(["-O", "0"], "-O"),
	(["-E", "0"], "-E"),
	(["-O", "4,24", "-E", "1,2"], None),
	(["-O", "4,24", "-E", "2,2"], None),
	(["-O", "4,3", "-E", "2,1"], None),
	(["-O", "4,121", "-E", "2,1"], "-O"),
	(["-O", "60", "-E", "4"], "-O"),
	(["--score-N", "4"], "--score-N"),
	(["-B", "6", "--score-N", "6"], "--score-N"),
	(["--score-N", "-1"], "--score-N"),
	(["-z", "100,200"], None),
	(["-z", "199,200"], None),
	(["-X", "--secondary=no"], "-X"),
	(["-P", "--secondary=no"], "-P"),
	(["-x", "ava-ont", "--secondary=no"], "-X"),
	(["--qstrand", "-a"], "--qstrand"),
	(["--qstrand", "-H"], "--qstrand"),
	(["--qstrand", "--frag=yes"], "--qstrand"),
	(["--qstrand", "--splice"], "--qstrand"),
	(["--qstrand", "-x", "map-pb"], "--qstrand"),
]

# the same rules at their boundaries, which hold
CHECK_OK = [
	["-r", "500,500"], ["-r", "20000"], ["-p", "0"], ["-p", "1"], ["-N", "0"], ["-k", "1"], ["-w", "1"],
	["-O", "4,120", "-E", "2,1"], ["-O", "59", "-E", "4"], ["-O", "5,24", "-E", "2,1"], ["-O", "4,23", "-E", "2,1"],
	["-O", "4,24", "-E", "2,1"], ["-O", "6", "-E", "1"], ["--score-N", "3"], ["--score-N", "0"],
	["-z", "200,200"], ["-z", "200"], ["--for-only"], ["--rev-only"], ["--rmq"], ["--qstrand"], ["-X"],
	["--split-prefix", "pfx"], ["--cs", "--MD"],
]


@pytest.mark.parametrize("args,named", CHECK_ERRORS, ids=[" ".join(a) for a, _ in CHECK_ERRORS])
def test_option_check_errors(args, named):
	"""Each option combination mm_check_opt rejects ends minimap2 with exit 1 and a message (naming
	the option where the message carries one) before any input is read; -h does not bypass it."""
	r = H.run_mm2(args + ["-h"])
	assert r.returncode == 1 and r.stdout == "" and r.messages, r.stderr
	if named:
		assert any(named in l for l in r.messages), r.messages


@pytest.mark.parametrize("args", CHECK_OK, ids=[" ".join(a) for a in CHECK_OK])
def test_option_check_boundaries(args):
	"""The same rules hold at their boundaries (-r NUM1 equal to NUM2; -p 0 and 1 inclusive; -N 0
	allowed with a warning; -k/-w 1; O1+E1 < O2+E2 with E1 > E2; ({-O}+{-E})+({-O2}+{-E2}) equal to
	127 or below; --score-N within [0,{-B}); Z-drop equal to the inversion Z-drop): minimap2 prints
	the usage and exits 0."""
	r = H.run_mm2(args + ["-h"])
	assert r.returncode == 0 and r.stdout.startswith("Usage:"), r.stderr
	assert not any("ERROR" in l for l in r.messages)


def test_deprecated_lj_min_ratio(small):
	"""--lj-min-ratio is deprecated: a warning naming it, and the output is unchanged."""
	fa, fq, _ = small
	r = H.map_reads(fa, fq, ["--lj-min-ratio", "0.9"])
	assert any("--lj-min-ratio" in l for l in r.messages)
	assert r.stdout == H.map_reads(fa, fq).stdout


@pytest.mark.parametrize("opt", ["--print-qname", "--print-seeds", "--print-aln-seq", "--print-chains", "--dbg-seed-occ"])
def test_debug_options(small, opt):
	"""Debugging options write to stderr only: the SAM output is unchanged and stderr gets lines
	besides the [M::] messages; --print-qname (minimap2.1: print query names to stderr) names
	every read there."""
	fa, fq, names = small
	base = H.map_reads(fa, fq, ["-a"])
	r = H.map_reads(fa, fq, ["-a", opt])
	assert H.strip_pg(r.stdout) == H.strip_pg(base.stdout)
	assert len(r.messages) > len(base.messages)
	if opt == "--print-qname":
		for n in names:
			assert any(n in l.split("\t") for l in r.messages), n


def _have_dev_full():
	return os.path.exists("/dev/full") and sys.platform.startswith("linux")


@pytest.mark.skipif(not _have_dev_full(), reason="needs /dev/full (Linux)")
@pytest.mark.parametrize("args", [["-a"], [], ["--paf-no-hit"]])
def test_write_error_exits_nonzero(small, tmp_path, args):
	"""A failed write of the alignments (here to /dev/full, which is always full) is an error:
	non-zero exit and a message (NEWS: exit with non-zero code on file writing errors). Both an
	output larger than the stdio buffer (the failure shows in mm_err_puts) and a tiny one (it
	shows only in the final fflush) fail."""
	fa, fq, _ = small
	tiny = H.write_fastq(str(tmp_path / "tiny.fq"), [H.read_fastx(fq)[0][:3]])
	bin_ = minimap2_build.minimap2_bin()
	for reads in (fq, tiny):
		with open("/dev/full", "wb") as out:
			p = subprocess.run([bin_] + args + ["-t", "2", fa, reads], stdout=out, stderr=subprocess.PIPE, timeout=120,
				env=dict(os.environ, LC_ALL="C"))
		assert p.returncode != 0
		assert any(l and not l.startswith(b"[M::") for l in p.stderr.split(b"\n"))


@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="liftrlimit() acts on Linux only")
def test_address_space_limit_is_lifted(small):
	"""On Linux minimap2 raises its soft RLIMIT_AS to the hard limit before working (main.c
	liftrlimit()), so a low soft limit with an unlimited hard limit does not stop it: with a
	40 MB soft limit, mapping with four threads (each thread stack alone reserves 8 MB) succeeds
	and gives the usual output."""
	import resource
	fa, fq, _ = small
	soft, hard = resource.getrlimit(resource.RLIMIT_AS)
	if hard != resource.RLIM_INFINITY:
		pytest.skip("the hard RLIMIT_AS is limited here")

	def limit():
		resource.setrlimit(resource.RLIMIT_AS, (40 << 20, resource.RLIM_INFINITY))
	bin_ = minimap2_build.minimap2_bin()
	p = subprocess.run([bin_, "-a", "-t", "4", fa, fq], stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=limit,
		timeout=120, env=dict(os.environ, LC_ALL="C"))
	assert p.returncode == 0, p.stderr[-2000:]
	base = H.map_reads(fa, fq, ["-a", "-t", "4"])
	assert H.strip_pg(p.stdout.decode()) == H.strip_pg(base.stdout)
