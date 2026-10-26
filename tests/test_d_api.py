"""The libminimap2 API through the C driver tests/c/drv_d_api.c (mm_set_opt, mm_mapopt_init,
mm_check_opt, the thread buffers, mm_map, mm_map_file, mm_gen_cs/ds/MD, the writers of
format.c, mm_write_sam_hdr, -L's CG tag), kalloc, kthread and misc.c, and minimap2-lite
(example.c).

The API must give what the command line gives: each case is compared with the minimap2 run
that uses the same options (main.c is a thin layer over these calls). Option defaults come
from the usage text and minimap2.1, the CG tag and BAM operation codes from the SAM spec.
"""

import re

import pytest

import helpers_d as D
import minimap2_helpers as H

DRV = "drv_d_api"


def _driver_out(r):
	"""Driver stdout without the framework's final '# mmt:' line."""
	return "".join(l for l in r.stdout.splitlines(True) if not l.startswith("# mmt:"))


@pytest.fixture(scope="module")
def flags():
	"""{MM_F_* name: value} and MM_VERSION as the driver sees them in minimap.h."""
	r = H.assert_c_case(DRV, "macros")
	out = {}
	for l in r.stdout.splitlines():
		if l.startswith("MM_F_"):
			k, v = l.split()
			out[k] = int(v, 16)
		elif l.startswith("MM_VERSION "):
			out["MM_VERSION"] = l.split(" ", 1)[1]
	return out


@pytest.fixture(scope="module")
def mixed(synth, tmp_path_factory):
	m = D.mixed_reads(synth.genome)
	d = tmp_path_factory.mktemp("d_api")
	return m, D.write_mixed(d / "mixed.fq", m), d


def test_macros_and_version(flags):
	"""MM_CIGAR_* are the BAM operation codes 0-8 of 'MIDNSHP=X' and MM_CIGAR_STR spells them (SAM
	spec 4.2); MM_VERSION is what --version prints."""
	assert flags["MM_VERSION"] == D.mm_version() == H.run_mm2(["--version"]).stdout.strip()


def test_set_opt_and_mapopt_init():
	"""mm_set_opt(0) gives the defaults of the usage text and minimap2.1 (-k 15 -w 10 -I 8G -A 2 -B 4
	-O 4,24 -E 2,1 -z 400,200 -r 500,20k -n 3 -m 40 -N 5 -p 0.8 -f 0.0002 -U 10,1000000
	--q-occ-frac 0.01 -e 500 --max-chain-skip 25 --score-N 1 -K 500M --cap-kalloc 500m
	--cap-sw-mem 100m --seed 11, PAF without CIGAR); mm_mapopt_init() fills the same mapping
	options; every documented preset is accepted and passes mm_check_opt, unknown names return -1;
	a preset leaves the parameters it does not set alone (minimap.h); sr and splice:sr include
	--secondary=no, -2 and their -K (minimap2.1)."""
	H.assert_c_case(DRV, "set_opt")


def test_check_opt():
	"""mm_check_opt() returns a negative value for each rule options.c enforces and 0 at each rule's
	boundary (the same rules as test_d_cli.py through the API); with mm_verbose 0 it prints
	nothing, with 1 one error line per rejection."""
	r = H.assert_c_case(DRV, "check_opt", 0)
	assert r.stderr == ""
	n = int(re.search(r"^rejected (\d+)$", r.stdout, re.M).group(1))
	r = H.assert_c_case(DRV, "check_opt", 1)
	assert len([l for l in r.stderr.splitlines() if l.strip()]) == n


def test_thread_buffer():
	"""mm_tbuf_init() gives a buffer whose kalloc pool (mm_tbuf_get_km) serves allocations;
	mm_tbuf_destroy(NULL) is harmless; with --no-kalloc's debug flag the buffer has no pool, so the
	libc allocator is used (minimap2.1 --no-kalloc)."""
	H.assert_c_case(DRV, "tbuf")


def _parse_map(text):
	hits, tags, unmapped = [], [], {}
	for l in text.splitlines():
		f = l.split("\t")
		if f[0] == "H":
			hits.append(f[1:])
			tags.append({})
		elif f[0] == "T":
			tags[-1][f[1]] = f[2]
		elif f[0] == "U":
			unmapped[f[1]] = (int(f[2]), int(f[3]))
	return hits, tags, unmapped


@pytest.mark.parametrize("preset,no_kalloc", [("-", 0), ("map-hifi", 0), ("-", 1)])
def test_mm_map_matches_cli(synth, mixed, preset, no_kalloc):
	"""mm_map() on every read gives the hits `minimap2 -c` prints, in the same order: query and
	target intervals, strand, columns 10-12, primary or secondary, CIGAR and the repetitive length
	(mm_tbuf_t.rep_len is rl:i); mm_gen_cs (short and long), mm_gen_ds and mm_gen_MD give the tags
	of --cs, --cs=long, --ds and --MD and return their lengths; exactly one hit per read has
	sam_pri, and it is the read's SAM primary record (the others are supplementary, or secondary
	when not their own parent); reads without hits match --paf-no-hit. With --q-occ-frac 0, reads
	on low-complexity tracts have rl:i > 0. no_kalloc=1 repeats it with the libc allocator."""
	m, fq, _ = mixed
	r = H.assert_c_case(DRV, "map", synth.fasta, preset, fq, "0", "0", no_kalloc)
	hits, tags, unmapped = _parse_map(r.stdout)
	pre = [] if preset == "-" else ["-x", preset]
	base = pre + ["--q-occ-frac", "0"]
	paf = H.map_reads(synth.fasta, fq, base + ["-c", "--paf-no-hit"]).paf()
	mapped = [p for p in paf if not p.is_unmapped]
	assert len(hits) == len(mapped) > 0
	for h, p in zip(hits, mapped):
		assert h[:12] == p.line.split("\t")[:12]
		assert h[12] == ("P" if p.tags["tp"] in "PI" else "S") and int(h[13]) == (p.tags["tp"] in "Ii")
		assert h[15] == p.cigar and int(h[16]) == p.tags["rl"]
	assert any(int(h[16]) > 0 for h in hits)
	assert unmapped == dict((p.qname, (p.qlen, p.tags["rl"])) for p in paf if p.is_unmapped)
	for opt, key in (("--cs", "cs"), ("--cs=long", "csl"), ("--ds", "ds"), ("--MD", "MD")):
		tag = "cs" if key == "csl" else key
		got = H.map_reads(synth.fasta, fq, base + ["-c", opt]).paf()
		assert [t[key] for t in tags] == [p.tags[tag] for p in got], opt
	sam = H.map_reads(synth.fasta, fq, base + ["-a"]).sam()
	recs = [x for x in sam.records if not x.is_unmapped]
	assert len(recs) == len(hits)
	for h, x in zip(hits, recs):
		assert (h[5], int(h[7]), h[4]) == (x.rname, x.pos0, x.strand)
		assert int(h[14]) == (x.is_primary)
		assert x.is_secondary == (h[12] == "S")


@pytest.fixture(scope="module")
def tx_reads(synth_tx, tmp_path_factory):
	reads = H.simulate_spliced_reads(synth_tx.genome, 25, seed=801, profile="hifi")
	return H.write_fastq(str(tmp_path_factory.mktemp("d_api_tx") / "tx.fq"), reads)


# (CLI options, MM_F_* flags main.c sets for them, threads, -K, read group, seed, preset)
MAP_FILE_CASES = [
	([], [], 1, None, None, None, None),
	(["-a"], ["OUT_SAM", "CIGAR"], 1, None, None, None, None),
	(["-c", "--ds", "--paf-no-hit"], ["OUT_CG", "CIGAR", "OUT_DS", "PAF_NO_HIT"], 2, None, None, None, None),
	(["-a", "--cs", "-Y", "--secondary-seq", "-y"], ["OUT_SAM", "CIGAR", "OUT_CS", "SOFTCLIP", "SECONDARY_SEQ", "COPY_COMMENT"],
		3, None, None, None, None),
	(["-a", "--MD", "-Q", "-L", "--eqx", "--sam-hit-only"], ["OUT_SAM", "CIGAR", "OUT_MD", "NO_QUAL", "LONG_CIGAR", "EQX",
		"SAM_HIT_ONLY"], 2, None, None, None, None),
	(["-a", "--cs=long", "--secondary=no", "-2"], ["OUT_SAM", "CIGAR", "OUT_CS", "OUT_CS_LONG", "NO_PRINT_2ND", "2_IO_THREADS"],
		4, None, None, None, None),
	(["-c", "-y", "-K", "3000"], ["OUT_CG", "CIGAR", "COPY_COMMENT"], 3, 3000, None, None, None),
	(["-a", "-R", "@RG\\tID:grp7\\tSM:s1"], ["OUT_SAM", "CIGAR"], 1, None, "@RG\\tID:grp7\\tSM:s1", None, None),
	(["-a", "--seed", "4"], ["OUT_SAM", "CIGAR"], 2, None, None, 4, None),
	(["--write-junc"], ["OUT_JUNC", "CIGAR"], 2, None, None, None, "splice:hq"),
	(["-a", "--cs"], ["OUT_SAM", "CIGAR", "OUT_CS"], 3, None, None, None, "splice:hq"),
]


@pytest.mark.parametrize("case", MAP_FILE_CASES, ids=[" ".join(c[0]) + (" -x " + c[6] if c[6] else "") or "paf" for c in MAP_FILE_CASES])
def test_mm_map_file_matches_cli(synth, synth_tx, mixed, tx_reads, flags, case):
	"""mm_map_file() (after mm_write_sam_hdr() in SAM mode, as main.c does) with the MM_F_* flags
	main.c sets for each output option prints what `minimap2` prints for that option (without
	@PG); minimap2 needs the same flag for the same output, so each option maps to its flag
	(TARGETS.md: MM_F_OUT_SAM set by -a, and so on). It returns 0, and -1 for a missing file."""
	args, names, threads, K, rg, seed, preset = case
	m, fq, _ = mixed
	fa, reads = (synth_tx.fasta, tx_reads) if preset else (synth.fasta, fq)
	flag = 0
	for n in names:
		flag |= flags["MM_F_" + n]
	r = H.assert_c_case(DRV, "map_file", fa, reads, preset or "-", hex(flag), threads, K or "-", rg or "-",
		"-" if seed is None else seed)
	got = _driver_out(r)
	cli = H.map_reads(fa, reads, (["-x", preset] if preset else []) + args, threads=threads)
	assert H.strip_pg(got) == H.strip_pg(cli.stdout)
	assert len(D.lines(got)) > 10


def test_writers_match_cli(synth, mixed, flags):
	"""On the hits of mm_map(): mm_write_paf3 and mm_write_paf4 (one segment) write the lines of
	`minimap2 -c --cs -y --paf-no-hit`, and mm_write_paf the same without rl:i (it has no repetitive
	length); mm_write_sam3 writes the records of `minimap2 -a --cs -y`, mm_write_sam2 the same
	without rl:i, and mm_write_sam (no options) those of `minimap2 -a` without rl:i."""
	m, fq, _ = mixed
	flag = 0
	for n in ("OUT_CG", "CIGAR", "OUT_CS", "COPY_COMMENT", "PAF_NO_HIT"):
		flag |= flags["MM_F_" + n]
	r = H.assert_c_case(DRV, "write", synth.fasta, fq, hex(flag), "0")
	out = {}
	for l in r.stdout.splitlines():
		if len(l) > 2 and l[1] == "\t":
			out.setdefault(l[0], []).append(l[2:])
	base = ["--q-occ-frac", "0"]
	paf = D.lines(H.map_reads(synth.fasta, fq, base + ["-c", "--cs", "-y", "--paf-no-hit"]).stdout)
	sam_cs = D.body_lines(H.map_reads(synth.fasta, fq, base + ["-a", "--cs", "-y"]).stdout)
	sam = D.body_lines(H.map_reads(synth.fasta, fq, base + ["-a"]).stdout)
	assert out["Q"] == paf
	assert out["F"] == [l for l in paf if "\t*\t*\t" not in l]
	assert out["P"] == [D.drop_tag(l, "rl") for l in paf if "\t*\t*\t" not in l]
	assert out["R"] == sam_cs
	assert out["T"] == [D.drop_tag(l, "rl") for l in sam_cs if not (int(l.split("\t")[1]) & 4)]
	assert out["S"] == [D.drop_tag(l, "rl") for l in sam]
	assert any("SA:Z:" in l for l in out["S"]) and any("\t256\t" in l or "\t272\t" in l for l in out["S"])


def test_sam_header_writer(synth):
	"""mm_write_sam_hdr(): @HD VN:1.6 SO:unsorted GO:query, an @SQ line per index sequence, the read
	group line with its escaped TABs turned into TABs, and @PG with ID and PN minimap2, VN when a
	version is given and CL when there are arguments (the program name replaced by minimap2);
	without an index no @SQ. A read group line not starting with @RG, holding a literal TAB or
	without ID returns -1 (the header is written without it)."""
	sq = "".join("@SQ\tSN:%s\tLN:%d\n" % (n, synth.genome.length(n)) for n in synth.genome.names)
	hd = "@HD\tVN:1.6\tSO:unsorted\tGO:query\n"
	cl = "\tCL:minimap2 -a x y z"

	def run(mode, *rg):
		r = H.assert_c_case(DRV, "hdr", synth.fasta, mode, *rg)
		out = _driver_out(r)
		ret = int(re.search(r"^ret (-?\d+)$", out, re.M).group(1))
		return out[:out.index("ret ")], ret, r
	assert run("full")[:2] == (hd + sq + "@PG\tID:minimap2\tPN:minimap2\tVN:v1" + cl + "\n", 0)
	assert run("noidx")[:2] == (hd + "@PG\tID:minimap2\tPN:minimap2\tVN:v1" + cl + "\n", 0)
	assert run("nover")[:2] == (hd + sq + "@PG\tID:minimap2\tPN:minimap2" + cl + "\n", 0)
	assert run("noargs")[:2] == (hd + sq + "@PG\tID:minimap2\tPN:minimap2\tVN:v1\n", 0)
	assert run("full", "@RG\\tID:a1\\tSM:x y")[:2] == (hd + sq + "@RG\tID:a1\tSM:x y\n@PG\tID:minimap2\tPN:minimap2\tVN:v1" + cl + "\n", 0)
	for bad in ("RG\\tID:a1", "@RG\tID:a1", "@RG\\tSM:x"):
		text, ret, r = run("full", bad)
		assert ret == -1 and text == hd + sq + "@PG\tID:minimap2\tPN:minimap2\tVN:v1" + cl + "\n" and r.stderr.strip()


def _cg_expect(n_ops, clip5, clip3, clip_op):
	ops = [(1 << 4 | (7 if i % 2 == 0 else 8)) for i in range(n_ops)]
	pre = [clip5 << 4 | clip_op] if clip5 else []
	post = [clip3 << 4 | clip_op] if clip3 else []
	return pre + ops + post


# (operations, 5' clip, 3' clip, kind, flags, CIGAR moved to CG:B:I)
LONG_CASES = [
	(65533, 1, 1, "pri", ["LONG_CIGAR"], False),
	(65534, 1, 1, "pri", ["LONG_CIGAR"], True),
	(65535, 0, 0, "pri", ["LONG_CIGAR"], False),
	(65536, 0, 0, "pri", ["LONG_CIGAR"], True),
	(65534, 1, 0, "pri", ["LONG_CIGAR"], False),
	(65540, 1, 1, "pri", [], False),
	(70000, 5, 7, "sup", ["LONG_CIGAR"], True),
	(70000, 5, 7, "sup", ["LONG_CIGAR", "SOFTCLIP"], True),
	(70000, 5, 7, "sec", ["LONG_CIGAR"], True),
	(70000, 5, 7, "sec", ["LONG_CIGAR", "SECONDARY_SEQ"], True),
]


@pytest.mark.parametrize("case", LONG_CASES, ids=["%d-%d-%d-%s-%s" % (c[0], c[1], c[2], c[3], "+".join(c[4]) or "none") for c in LONG_CASES])
def test_long_cigar(flags, case):
	"""-L (minimap2.1; NEWS: the fake CIGAR is <readLen>S<refLen>N): when an alignment's CIGAR has
	more than 65535 operations counting the clips (the BAM limit), the CIGAR field becomes
	<SEQ length>S<reference span>N and CG:B:I holds the real CIGAR in BAM encoding, clips included
	(hard on supplementary records and on secondary records with --secondary-seq, unless -Y); at
	65535 or fewer, or without -L, the CIGAR stays in its field. Through mm_write_sam3() on a
	hand-made hit of alternating 1=/1X operations."""
	n_ops, c5, c3, kind, names, in_tag = case
	flag = 0
	for n in names:
		flag |= flags["MM_F_" + n]
	r = H.assert_c_case(DRV, "long_cigar", n_ops, c5, c3, kind, hex(flag))
	line = [l for l in r.stdout.splitlines() if l.startswith("lc\t")]
	assert len(line) == 1
	rec = H.parse_sam(line[0] + "\n", strict=False).records[0]
	f = line[0].split("\t")
	hard = (kind == "sup" or (kind == "sec" and "SECONDARY_SEQ" in names)) and "SOFTCLIP" not in names
	text_cigar = ("%dH" % c5 if c5 else "") if hard else ("%dS" % c5 if c5 else "")
	text_cigar += "".join("1=" if i % 2 == 0 else "1X" for i in range(n_ops))
	text_cigar += (("%dH" % c3) if hard else ("%dS" % c3)) if c3 else ""
	qlen = c5 + n_ops + c3
	if kind == "pri" or "SOFTCLIP" in names:
		slen = qlen
	elif kind == "sec" and "SECONDARY_SEQ" not in names:
		slen = 0
	else:
		slen = n_ops
	assert (len(f[9]) if f[9] != "*" else 0) == slen
	if in_tag:
		assert f[5] == "%dS%dN" % (slen, n_ops)
		assert rec.tag_types.get("CG") == "B" and rec.tags["CG"][0] == "I"
		assert rec.tags["CG"][1] == _cg_expect(n_ops, c5, c3, 5 if hard else 4)
		assert D.bam_to_cigar(rec.tags["CG"][1]) == text_cigar
	else:
		assert f[5] == text_cigar and "CG" not in rec.tags


def test_kalloc():
	"""kmalloc, kcalloc, krealloc, krelocate and kfree keep every block's contents, return memory
	aligned for 8-byte words, zero kcalloc memory even when recycled, follow the libc semantics
	for a NULL pointer or size 0, and fall back to libc with a NULL pool; a child pool (km_init2)
	draws on its parent; km_stat reports a capacity covering the live bytes. km_stat_print runs, but
	what it prints is kalloc's debug output and is not checked."""
	r = H.assert_c_case(DRV, "kalloc")
	assert re.search(r"^capacity \d+$", r.stdout, re.M)


def test_kthread():
	"""kt_for visits each index exactly once with a thread id below the thread count, for 1-8
	threads and 0 to 1000 indices; kt_pipeline passes every item through its three steps, runs the
	first and last step one item at a time, and finishes items in input order even when later
	items finish their middle step first (what keeps minimap2's output in input order)."""
	H.assert_c_case(DRV, "kthread")


def test_misc(tmp_path):
	"""realtime() is the wall-clock time in seconds; cputime() advances in a busy loop and never
	faster than the wall clock; peakrss() is non-negative and does not decrease (no values: units
	differ by platform); mm_err_puts writes the string and a newline to stdout; mm_err_fwrite and
	mm_err_fread round-trip data through a file."""
	r = H.assert_c_case(DRV, "misc", str(tmp_path / "f.bin"))
	assert "mm_err_puts\tline\n" in r.stdout


def test_minimap2_lite(synth, mixed, lite_bin, index_of):
	"""minimap2-lite (example.c), the documented way to use the library: for each read, each hit as
	12 PAF columns and cg:Z, exactly as `minimap2 -c` with the same default options prints them
	(--no-hash-name, since example.c passes no read name to mm_map), from a FASTA target or a
	prebuilt index; without two arguments it prints its usage to stderr and exits 1."""
	m, fq, _ = mixed
	cli = H.map_reads(synth.fasta, fq, ["-c", "--no-hash-name"]).paf()
	want = ["\t".join(p.line.split("\t")[:12]) + "\tcg:Z:" + p.cigar for p in cli]
	for target in (synth.fasta, index_of(synth.genome)):
		r = H.run_process([lite_bin, target, fq], check=True)
		assert D.lines(r.stdout) == want
	for args in ([], [synth.fasta]):
		r = H.run_process([lite_bin] + args)
		assert r.returncode == 1 and r.stdout == "" and "Usage" in r.stderr
