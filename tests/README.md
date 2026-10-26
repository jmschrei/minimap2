# minimap2 test suite

Unit and end-to-end tests of minimap2: reading sequences, sketching and
index building, seeding, chaining, alignment, splice and paired-end
mapping, SAM/PAF output, every command-line option and preset, the
libminimap2 API in `minimap.h`, the mappy Python binding and the `misc/`
scripts. `TARGETS.md` lists every target, grouped by area, and which test
covers it.

## Running

Requirements: Python 3.8 or later, pytest, make, a C compiler and zlib.
Optional: Cython and setuptools (mappy tests), k8 (paftools.js and
pafcluster.js tests), samtools and pysam (cross-checks). Tests that need a
missing optional tool are skipped rather than reported as errors.

```sh
python3 -m pytest tests                           # from the repository root
cd tests && python3 -m pytest                     # same
MINIMAP2_TEST_THREADS=4 python3 -m pytest tests   # run minimap2 with -t 4
python3 -m pytest tests -k pilot                  # a subset
```

The first run builds minimap2 with the Makefile's own flags in a cache
directory (about 2 s), and `make extra` and the mappy extension when a test
needs them; later runs reuse the builds until a source file or the compiler
changes. Nothing is written to the source tree: sources are copied to the
cache and built there, and every child process runs in a temporary working
directory.

### mappy (Cython)

The mappy tests build the extension with `python3 setup.py build_ext` from a
copy of the sources, with the interpreter running pytest, and import it from
the cache; nothing is installed. That interpreter needs Cython and
setuptools. Without installing anything into a shared environment, run the
suite through uv:

```sh
uv run --no-project --with pytest --with cython --with setuptools python3 -m pytest tests
uv run --no-project --python 3.8 --with pytest --with cython --with setuptools python3 -m pytest tests
```

Without Cython or setuptools the mappy tests are skipped (the pytest header
says why). If both are present and the build fails, the tests fail.

### k8 (paftools.js, pafcluster.js)

Tests of the `misc/` scripts run `k8` from `PATH` and are skipped when it
is not there. k8 1.x release binaries (MIT, github.com/attractivechaos/k8)
work; for example:

```sh
curl -L https://github.com/attractivechaos/k8/releases/download/v1.2/k8-1.2.tar.bz2 | tar -jxf -
PATH=$PWD/k8-1.2:$PATH python3 -m pytest tests      # after copying k8-x86_64-Linux to k8-1.2/k8
```

### Compilers and builds

`CC` selects the compiler for the build, the C drivers and mappy, for
example `CC=gcc-11` or `CC=clang`. `MINIMAP2_TEST_MAKE_ARGS` passes extra
arguments to make, e.g. `sse2only=1` (the SSE2-only ksw2 kernels instead of
the SSE2/SSE4.1 runtime dispatch) or `arm_neon=1 aarch64=1` on arm64; it is
part of the cache key, so each combination gets its own build. A clang
installed in a conda environment does not find its zlib by itself; point the
compiler, the linker and the loader at the environment:

```sh
E=/path/to/conda/env    # holds bin/clang, include/zlib.h, lib/libz.so
CC=$E/bin/clang CPATH=$E/include LIBRARY_PATH=$E/lib LD_LIBRARY_PATH=$E/lib \
	python3 -m pytest tests
```

| Variable | Default | Meaning |
|:---------|:--------|:--------|
| `MINIMAP2_BIN` | build from source | binary under test |
| `MINIMAP2_TEST_SRC_DIR` | repository root | source tree to build minimap2, `libminimap2.a`, `make extra`, mappy and the C drivers from, and to take `misc/` and `python/minimap2.py` from |
| `CC` | `gcc` | compiler for the build, the C drivers and mappy |
| `MINIMAP2_TEST_MAKE_ARGS` | unset | extra `make` arguments, e.g. `sse2only=1` or `asan=1` (part of the cache key) |
| `MINIMAP2_TEST_MAKE_JOBS` | 4 | `make -j` (and `build_ext -j`) |
| `MINIMAP2_TEST_CACHE_DIR` | `$XDG_CACHE_HOME/minimap2-tests` or `~/.cache/minimap2-tests` | build cache |
| `MINIMAP2_TEST_THREADS` | 1 | `-t` for tests that do not set it themselves |
| `MINIMAP2_TEST_TIMEOUT` | 120 | seconds before one child process is killed |
| `MINIMAP2_TEST_UPDATE_PINS` | unset | `1` rewrites pinned outputs (only on base minimap2) |
| `MMT_VERBOSE` | unset | C drivers also print passing checks |

`MINIMAP2_BIN` replaces only the `minimap2` binary. `sdust`,
`minimap2-lite`, the C drivers and mappy are always built from
`MINIMAP2_TEST_SRC_DIR`, so to test a modified tree point
`MINIMAP2_TEST_SRC_DIR` at it. The build cache is keyed by a hash of the
Makefile, the top-level `*.c`/`*.h` files, `sse2neon/`, the compiler's path
and version and `MINIMAP2_TEST_MAKE_ARGS` (mappy's key adds setup.py,
`python/`, the interpreter and the Cython and setuptools versions), and
builds hold a file lock, so parallel workers and concurrent runs are safe.
Child processes run with `LC_ALL=C`, in a fresh temporary working directory
unless a test passes `cwd=`, and with the core-file limit at 1 byte: no core
is written, and on Linux a crash handler named in `core_pattern` (such as
apport, about a second per abort) is skipped, which a limit of 0 would not
do.

## Layout

| File | Contents |
|:-----|:---------|
| `TARGETS.md` | every target, grouped by area, with a Tested column |
| `minimap2_build.py` | locate or build the binary, `libminimap2.a`, `make extra`, C drivers and mappy out of tree |
| `minimap2_helpers.py` | genome and read simulation, FASTA/FASTQ and annotation I/O, running minimap2, SAM/PAF parsers, tag recomputation, invariant checks, pins, C driver, k8 and minimap2.py runners |
| `conftest.py` | session fixtures |
| `pytest.ini` | test paths, `-ra --strict-markers`, the markers |
| `c/minimap2_test.h` | check-and-report framework for C drivers |
| `c/drv_*.c` | C drivers, one binary each, compiled on demand |
| `pins/<area>/` | pinned outputs |
| `test_harness.py` | tests of the harness itself |
| `test_pilot.py` | whole-program sanity checks, one API, one mappy and one paftools case |
| `test_<area>_*.py`, `helpers_<area>.py` | tests and helpers per area of `TARGETS.md` (a to e) |

Markers (`pytest.ini`): `pinned` (pins base output), `external`
(samtools/pysam), `thread_invariance`, `mappy`, `k8`.

## Fixtures (conftest.py)

| Fixture | Scope | Value |
|:--------|:------|:------|
| `mm2_bin` | session | path of the binary under test |
| `sdust_bin`, `lite_bin` | session | `sdust` and `minimap2-lite` from `make extra` |
| `mappy` | session | the mappy module built from the sources under test (skips without Cython/setuptools) |
| `k8` | session | path of k8 (skips when it is not on `PATH`) |
| `threads` | session | `MINIMAP2_TEST_THREADS` |
| `index_dir` | session | a directory that lives for the session |
| `index_of` | session | `index_of(genome_or_fasta, preset=None, extra_args=())` -> `.mmi` path, built once per genome, preset and arguments |
| `mt` | session | `.human`, `.orang` (paths) and `.human_seqs`, `.orang_seqs` ({name: seq}) for `test/MT-*.fa` |
| `synth` | session | `.genome` (a `Genome`) and `.fasta` of the standard genome (`conftest.make_synth()`): 80 and 50 kbp contigs, two repeat families of three exact copies, one inverted repeat, three low-complexity tracts, two N runs, two short contigs |
| `synth_tx` | session | `.genome`, `.fasta`, `.bed12`, `.gtf` of the transcript genome (`conftest.make_synth_tx()`): 60 and 40 kbp, 8 genes of 2-6 exons with GT-AG introns of 80-2000 bp, one repeat family |
| `tree_before` | session | build products in the source tree before the session |

## Helpers (minimap2_helpers.py)

Simulation, deterministic for a given seed on every Python version (the
generator is splitmix64 in integer arithmetic; never use `random` or
`hash()` for inputs):

* `make_genome(seed, contig_lengths, names, repeat_families, repeat_len, repeat_copies, repeat_divergence, soft_mask_repeats, inversions, inversion_len, inversion_gap, n_runs, n_run_len, low_complexity, low_complexity_len, genes, exons, exon_len, intron_len, short_contigs, short_contig_len)` -> `Genome` with `.names`, `[name]`, `.fetch(contig, s, e, strand)`, `.features` (`Feature(contig, start, end, kind, info)`, kind `repeat`/`inversion`/`N`/`lowcomp`), `.genes` (`Gene(name, contig, strand, exons)` with `.introns`, `.transcript(genome)`), `.is_unique(contig, start, end, pad)`, `.write_fasta(path)`.
* `simulate_reads(genome, n, length, seed, sub_rate, indel_rate, max_indel, end_margin, n_rate, unique, contigs, strand, qual, name_prefix)` -> `[SimRead]`. A `SimRead` has `.name .seq .qual` (as sequenced) and the truth `.segments` (`Segment(contig, pos, end, strand, qstart, qend, cigar, edits)` with `.blocks` and `.introns`), plus shortcuts for the first segment: `.contig .pos .end .strand .cigar .edits .fwd_seq .nm .exons .introns`.
* `simulate_long_reads(genome, n, profile, seed, **overrides)` with `PROFILES`: `ont` (2-8 kbp, about 7% error), `hifi` (3-10 kbp, about 0.4%), `clr` (indel-rich, about 7%), `asm`, `short`.
* `simulate_pairs(genome, n, read_length, insert_mean, insert_sd, ...)` -> `[SimPair]` (FR; `.r1 .r2 .contig .frag_start .frag_end .strand .insert`).
* `simulate_spliced_reads(genome, n, seed, profile, full_length_frac, min_len, antisense_frac, junction_margin, min_end_exon, length)`: cDNA reads from `genome.genes` (whole transcripts or sub-intervals; `length=100` for fixed-length short reads); truth `.exons`, `.introns`, `.gene`, `.ts_strand` (`+` for a read in the transcript's sense, as ts:A reports it) and `.strand` (on the genome).
* `simulate_spliced_pairs(genome, n, read_length, insert_mean, insert_sd, stranded)`: paired RNA-seq reads (FR on the transcript) for splice:sr, each mate with the same truth fields. `transcript_blocks(gene, a, b)` maps a transcript interval to genome blocks.
* `simulate_chimeric_reads(genome, n, profile, length, parts, min_part, min_distance)`: pieces from distant loci; one `Segment` per piece with its `qstart/qend`.
* `simulate_contigs(genome, n, length, seed, inversions, inversion_len)`: assembly-like contigs (about 0.1% divergence) with optional inverted pieces (segments adjacent on the reference, alternating strands).
* `Rng(seed)`: `u64 below randint random chance choice normal_int dna fork`.

I/O: `write_fastx/write_fasta/write_fastq(path, records)` (format and gzip
from the extension; gzip bytes are stable), `write_pairs(p1, p2, pairs,
suffix=False)`, `write_interleaved(path, pairs)`, `write_text`,
`genes_bed12(genes)`, `genes_gtf(genes)`, `junctions_bed(genes)`,
`fastx_text`, `read_fastx`, `read_fasta`, `read_text`. Records may be
`SimRead`s or `(name, seq[, qual[, comment]])` tuples.

Running: `run_mm2(argv, cwd=None, check=False)` and `map_reads(target,
reads, args=(), preset=None, threads=None, check=True, cwd=None)` return a
`Result` with `.returncode .stdout .stderr .stdout_bytes .messages .sam()
.paf() .check()`; `.messages` is stderr without the `[M::...]` progress
lines. `map_reads` adds `-t $MINIMAP2_TEST_THREADS` unless `threads=` or a
`-t N`/`-tN` argument is given, puts `-x preset` first, and makes the target
and read paths absolute. Both run the process in a new temporary directory
(removed afterwards) unless `cwd` is given, so give other paths in `args`
(`-o`, `-d`, `--split-prefix`, BED files) as absolute paths.
`index_genome(fa_or_genome, dir, preset, extra_args)` runs `minimap2 -x
preset -d`. `run_process(argv)` runs any program the same way;
`run_paftools(args, stdin)` and `run_k8(script, args)` run `misc/*.js`;
`run_minimap2_py(args)` runs `python/minimap2.py` with the cached mappy. A
process that aborts returns `-6` without a core dump or crash-handler delay.

Parsing: `parse_sam(text, strict=True)` -> `SamFile` (`.hd .sq .rg .pg .co
.records .by_qname() .primaries()`); `SamRecord` fields `qname flag rname
pos mapq cigar rnext pnext tlen seq qual tags tag_types tag_order extra`,
with flag properties (`is_reverse`, `is_secondary`, `is_supplementary`,
`segment`, ...) and coordinates (`pos0 ref_end0 left_clip right_clip
unclipped_start0 query_length introns`). Tags are typed (`i` int, `f`
float, `A`/`Z`/`H` str, `B` `(subtype, [values])`); strict parsing rejects
duplicate tags, CIGAR/SEQ/QUAL length mismatches, out-of-range fields and
fields that are not tags (use `strict=False` for `-y` comments; they land
in `.extra`). `parse_paf(text)` -> `[PafRecord]` (`qname qlen qstart qend
strand tname tlen tstart tend nmatch alen mapq tags`, `.cigar` = cg:Z,
`.is_unmapped`, `.is_primary`).

Recomputation from the genome (seq in reference orientation, soft clips
included): `compute_nm`, `compute_md` (spec form; `md_as_minimap2()` drops
the final 0 the way minimap2 does), `compute_cs(..., long_form)`,
`compute_mlen_blen` (PAF columns 10/11 and nn:i), `compute_de` (de:f),
`cigar_score(..., a, b, q, e, q2, e2, sc_ambi)` (the two-piece affine score
of the usage text), `cigar_ref_span`, `cigar_query_span`, `cigar_clips`,
`strip_clips` (cg:Z from a SAM CIGAR), `parse_cigar`, `revcomp`, `nt4`,
`tlen_5prime`, `tlen_outer`, `placement_error(rec, read)`,
`paf_placement_error(rec, read)`, `paf_query_part(rec, read)`. Generic
checks returning a list of problems: `sam_invariant_errors(sam, genome,
reads)` (header, one primary per segment, SEQ, CIGAR span, NM/MD/cs),
`paf_invariant_errors(paf, genome, reads)` (columns, cg spans, columns
10/11, nn, NM/MD/cs/de) and `pair_invariant_errors(sam)`.

Pins: `check_pinned(name, text)` compares with `pins/<name>` (gzip when the
name ends in `.gz`) and shows a diff; `strip_pg(sam)` drops the `@PG` line
and `sam_body(sam)` all header lines. Area writers put pins under
`pins/<area>/`.

## Adding a C driver

1. Create `c/drv_<area>_<topic>.c` that includes `minimap2_test.h` and any of
   minimap2's headers (`minimap.h`, or internal ones such as `mmpriv.h`,
   `bseq.h`, `ksw2.h`, `sdust.h`, `kalloc.h`, `kthread.h`); the include path
   holds the copied source tree.
2. Write each case as `static int name(int argc, char *argv[])` using
   `MMT_CHECK(cond)`, `MMT_CHECKF(cond, fmt, ...)`, `MMT_CHECK_INT(a, op, b)`,
   `MMT_CHECK_STR(a, b)`, `MMT_REQUIRE(cond)` (stops the case) and
   `MMT_FAIL(fmt, ...)`. `argv` holds the arguments after the case name.
3. End the file with a `{name, function}` table terminated by `{0, 0}` and
   `MMT_MAIN(table)`.
4. From Python: `H.assert_c_case("drv_<...>", "case", arg1, ...)` compiles
   the driver if needed (`$CC -std=gnu99 -g -O2 -Wall ... libminimap2.a -lz
   -lm -lpthread`), runs it in a temporary directory and fails with its
   output unless every check passed. `H.c_cases(driver)` lists the cases.

`libminimap2.a` holds every object of `make minimap2` except `main.o`,
including bseq.o, kthread.o, kalloc.o, sdust.o and the ksw2 kernels, so a
driver can call internal functions through their headers. Each failed check
prints `not ok N - file:line: message`; a run ends with `# mmt: CASE
checks=N failures=M`. Pass inputs (FASTA or `.mmi` paths, expected values)
as arguments from Python. `mmt_read_seqs()` reads a small plain FASTA/FASTQ
file. Match minimap2's C style (tabs, C99, K&R braces). `drv_pilot` uses
only `minimap.h`. A driver that builds or reads internal structs
(`mm128_t` arrays for chaining, `ksw_extz_t`, `mm_reg1_t` by hand) should
say so in its header comment, since a speed change may alter them; the
public structs in `minimap.h` are fixed for the speed search.

## Rules for expected values

These bind every test in the suite.

* Derive expected values independently wherever practical: inputs small
  enough to work by hand; a short reference implementation inside the test
  (brute-force minimizers with minimap2's documented hash and (w,k) rule,
  a brute-force affine-gap or two-piece affine alignment, NM/MD/cs/cg
  recomputed from the genome); the simulator's recorded truth; the SAM and
  PAF specifications; minimap2.1, the README, cookbook.md, FAQ.md and
  python/README.rst. A test that only records the current output says so in
  its docstring ("pins the current output") and carries the `pinned`
  marker; keep these few.
* No test of behaviour base minimap2 gets wrong (behaviour that contradicts
  its documentation or the SAM/PAF spec). Such behaviour is reported
  separately, not tested.
* No dependence on platform, compiler, Python version or thread count unless
  the test is explicitly about thread-count invariance (`thread_invariance`
  marker). minimap2 builds on Linux x86-64 (gcc, clang), Linux aarch64
  (`arm_neon=1 aarch64=1`), macOS arm64 (clang) and with `sse2only=1`: do
  not depend on SIMD width, `ru_maxrss` units, `/proc`, `sed -i` dialects or
  GNU-only tool options.
* Each test file runs fast; the whole suite should stay under about two
  minutes at `MINIMAP2_TEST_THREADS=4` on a loaded machine, since it runs
  once per iteration of a speed search. A minimap2 process on these inputs
  takes 10-50 ms; prefer one run checked many ways over many runs, and keep
  simulated read counts in the tens to hundreds (the Python simulator costs
  about 2 ms per kbp of long reads).
* Do not pin implementation details that a speed optimization could change
  while keeping minimap2's output byte-identical: internal struct layouts
  outside `minimap.h`, allocation counts, kalloc statistics, the order of
  internal calls, stderr timing lines, debug output formats (`--print-*`).
  Pinning the bytes of `minimap2` SAM/PAF output on fixed inputs is fine
  and wanted; so are `.mmi` bytes and libminimap2 API results, which the
  speed search keeps fixed.
* Never compare stderr wholesale. For a warning or an error, check that a
  message appears (`Result.messages`) and, where it carries one, that it
  names the offending option or value; do not match its wording or its ANSI
  colour codes.

## Test data

Synthetic by default. Small real subsets (a few hundred reads at most) may
be added to `tests/data/` only from these sources, each with a NOTICE file
beside it:

* ENCODE and 4DN K562 data: a NOTICE citing ENCODE or 4DN and the ENCSR/ENCFF
  or 4DN accessions.
* NIST GIAB HG002 PacBio HiFi and ONT ultra-long (ONT-UL) reads: the NIST
  notice and a note that the data were changed (subsetted), with the date.
* GIAB HG002 HiSeq WGS short reads: NIST/GIAB attribution and a "subsetted"
  notice.

Never copy HG002 ONT Open Data reads (CC BY-NC 4.0); they may only be used
for timing outside the repository. The files already in `test/` may be used
as they are. The `__cpuidex` wrapper in `ksw2_dispatch.c` carries unknown
terms: no test may depend on editing it.

## Things to know about minimap2's output

* The `@PG` line carries `VN:` (`MM_VERSION`) and `CL:minimap2` followed by
  the command-line arguments, including temporary paths. Compare or pin
  `strip_pg(text)`.
* stderr carries `[M::...]` progress lines with timings and CPU ratios, the
  index statistics, `mid_occ = N`, and a final `Real time ... Peak RSS`
  line. `[WARNING]` and `[ERROR]` lines carry ANSI colour escapes.
* Exit status: 1 for an unknown option, an unknown preset, a missing or
  unreadable target, a missing query file (`ERROR: failed to map the query
  file`; the SAM header is already written) and option conflicts rejected by
  `mm_check_opt`; `-h` prints the usage to stdout and exits 0; no arguments
  print it to stderr and exit 1; `--version` prints `2.31-r1302`.
* The CLI's default is `-t 3`; `map_reads` always passes `-t`.
* Output was byte-identical between gcc 13, gcc 11, clang 18 and the
  `sse2only=1` build for map-ont, map-hifi, map-pb, lr:hq, lr:hqae, sr,
  splice, splice:hq, splice:sr, asm5, asm20, ava-ont and `--rmq` on the
  simulated data, so pins hold across these builds on x86-64. NEON builds
  were not compared.
* Output was byte-identical across `-t 1, 2, 4, 8` and `-K 5k, 50k` and the
  default for single-end map-ont (SAM and PAF), paired-end sr and splice on
  the simulated data, and `.mmi` files are identical across `-t` at
  indexing. A multi-part index (`-I` smaller than the target) changes the
  output: `mid_occ` is computed per part, and without `--split-prefix` the
  SAM header has no `@SQ` lines.
* Read names seed the hash that breaks ties between equal hits (unless
  `--no-hash-name`), as does `--seed`. Renaming a read can move a read from a
  repeat to another copy; do not assert which copy a repetitive read lands
  on.
* `--cs` together with `--MD` writes only MD:Z, and `--cs` with `--ds` only
  ds:Z (base behaviour; see Rules). Run minimap2 twice when a test needs
  both.
* MD:Z omits the final `0` when an alignment ends with a mismatch (`65T`
  instead of `65T0`; common with `-x sr`, whose end bonus extends alignments
  to the read end). Compare with `md_as_minimap2(compute_md(...))`.
* Bases compare by minimap2's 0-4 code: in MD and cs, N against N is a match
  and IUPAC codes in the target are stored and printed as N; NM counts every
  column with a non-ACGT base as a mismatch; PAF columns 10 and 11 leave
  such columns out and nn:i counts them.
* cs:Z long form (`--cs=long`) writes the read's bases in uppercase. An
  intron is `~` plus its first two reference bases, its length and its last
  two (e.g. `~gt364ag`); MD skips introns.
* ts:A appears only on spliced alignments (cookbook.md): `+` when the read
  is in the transcript's orientation.
* Secondary alignments are written by default (up to `-N`, 5) with SEQ and
  QUAL `*` (`--secondary-seq`: the aligned part, hard-clipped);
  supplementary alignments are hard-clipped unless `-Y`; the primary of a
  split read lists the others in SA:Z with approximate CIGARs.
* Unmapped reads: a SAM line with flag 4 (not with `--sam-hit-only`); PAF
  leaves them out unless `--paf-no-hit`, which writes `name qlen 0 0 * * 0
  0 0 0 0 0 rl:i:N`.
* PAF without `-c` has approximate columns 10/11 and the tags tp, cm, s1,
  s2, dv, rl; `-c` adds NM, ms, AS, nn, de and cg (and ts on spliced reads).
* AS:i equalled the CIGAR's score under the preset's `-A/-B/-O/-E`
  (`cigar_score`) on every alignment of simulated map-ont and map-hifi reads
  and the MT example; under `-x sr` it was a few points higher on about 4%
  of reads.
* TLEN is the 5'-to-5' template length (`tlen_5prime`). An unmapped mate is
  placed at its mapped mate's position. With `-x sr` or `--frag=yes`, a
  single query file is read as interleaved pairs: adjacent reads with the
  same name (a trailing `/1` or `/2` ignored) are mates, so single-end test
  reads need distinct names. SAM drops a trailing `/1` or `/2` from paired
  reads' names; PAF appends `/1` or `/2` to them, even to a name that already
  ends in `/1` (`p0/1/1`).
* `mm_idx_name2id()` on an index from `mm_idx_str()` returns an unset value
  (base behaviour; see Rules): look up names only on indexes read from a
  file, after `mm_idx_index_name()`.
* mappy builds a single-part index, always computes CIGARs, returns
  secondary alignments too, and computes cs and MD independently of each
  other. `Alignment.__str__` writes `ts:A:.` when the strand is unknown.
  `python/minimap2.py` ignores `-m` (base behaviour).

## Coverage

Not measured yet; see `TARGETS.md` for what each area covers.

## Local setup for the 2026-10-26 speed search

Remove this section before the suite goes upstream.

* k8 1.2: `/home/jmschrei/auto-optimization/2026-10-26-minimap2-speed/toolchains/k8/bin`
  (put it on `PATH`).
* clang 18.1.8: `E=/home/jmschrei/auto-optimization/2026-10-24-chromap-speed/toolchains/clang`
  with `CC`, `CPATH`, `LIBRARY_PATH` and `LD_LIBRARY_PATH` as above.
* gcc 11: `CC=gcc-11`. Python 3.8: `uv run --no-project --python 3.8 ...`.
* The data rules above are the user's decisions on the license record in
  `/home/jmschrei/auto-optimization/2026-10-26-minimap2-speed/TERMS.md`
  (with `USAGE.md` beside it); the NOTICE texts are in
  `/home/jmschrei/auto-optimization/2026-10-24-minibwa-speed/TERMS.md`.
* Behaviour of base minimap2 that contradicts its documentation is filed in
  `/home/jmschrei/auto-optimization/2026-10-26-minimap2-speed/known-bugs/`,
  not tested.
