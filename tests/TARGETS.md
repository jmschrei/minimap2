# Test targets

Every target of the minimap2 test suite, listed from the source before any
area test was written: every command-line option of `minimap2`, every `-x`
preset, every public name in `minimap.h`, every non-static function in the
C files the suite covers, the `make extra` programs, `misc/paftools.js` and
`misc/pafcluster.js`, and the mappy binding. Writers tick a row by replacing
`[ ]` in its Tested column with the test that covers it, e.g.
`test_b_chain.py::test_max_chain_skip`. A row that cannot be tested
(platform-only code, or behaviour that contradicts the documentation) gets
a short reason instead.

Rows are grouped by area: each writer edits only its own `## Area X`
section, so parallel edits never touch the same lines.

## Areas

Five areas, each written in parallel by one writer. A writer creates and
edits only the files with its prefix, so no two writers edit the same file:
`tests/test_<a>_*.py`, `tests/helpers_<a>.py`, `tests/c/drv_<a>_*.c` and
`tests/pins/<a>/`. The shared harness (`minimap2_build.py`,
`minimap2_helpers.py`, `conftest.py`, `c/minimap2_test.h`, `pytest.ini`,
`README.md`, `test_harness.py`, `test_pilot.py`) is not edited by area
writers; a helper one area needs goes in its `helpers_<a>.py`, and one
several areas need is reported back for the harness.

| Area | Prefix | Scope |
|:-----|:-------|:------|
| A | `a` | Input, sketching and index: bseq.c (FASTA/FASTQ, gzip, several files, stdin), sketch.c (minimizers, HPC), sdust.c and the `sdust` program, index.c (building, `-d`, `.mmi` load/dump/round trip, `-k -w -H -I`, `--idx-no-seq`, `--bucket-bits`, multi-part indexes, name and sequence lookup), splitidx.c and `--split-prefix`, the index parameters of every preset, `-T`. |
| B | `b` | Seeding, chaining, alignment and hit selection: seed.c, lchain.c (DP and RMQ chaining), hit.c, align.c, esterr.c, ksw2_*.c, `mm_mapopt_update`/`mm_idx_cal_max_occ` (`-f -U --q-occ-frac -e`), the chaining, scoring and selection options, `-c`, `--alt`, `--qstrand`, the all-vs-all options (`-X -D -P --dual --no-long-join`), and the mapping parameters of the long-read and assembly presets. |
| C | `c` | Splice and short reads: jump.c, pe.c, the junction/annotation readers in index.c, `mm_map_frag`/`mm_map_file_frag`, `mm_seg_gen`, splice options (`-u -G -C -J -j --junc-bed --junc-bonus --junc-pen --spsc* --pass1 --jump-min-match --splice --splice-flank`), short-read and pairing options (`--sr --frag -F --pairing`), paired input (two files and interleaved), proper-pair flags and mate fields, and the mapping parameters of the `sr` and splice presets. |
| D | `d` | Output, CLI and API: format.c (SAM and PAF fields and tags, `--cs --MD --ds -c`'s cg:Z, `-L -Y -y -Q --eqx -R --secondary --secondary-seq --write-junc --paf-no-hit --sam-hit-only -o`), every remaining CLI option (`--seed -t -K -2 -v -x` errors, `--version/--help`, debug and kalloc options, deprecated options), option validation (`mm_check_opt`), thread and batch invariance, the libminimap2 API through C drivers (`mm_set_opt`, `mm_map`, `mm_map_file`, `mm_gen_cs/ds/MD`, `mm_tbuf_*`), `minimap2-lite`, kthread.c, kalloc.c and misc.c. |
| E | `e` | mappy and misc/: the mappy binding built from the source under test (`Aligner` and its arguments, `map`, `seq`, `seq_names`, `Alignment` fields, `revcomp`, `fastx_read`, `ThreadBuffer`, `verbose`, options passed through `extra_flags`/`scoring`), `python/minimap2.py`, and every `paftools.js` subcommand and `pafcluster.js` (skipped without k8). |

Differences from the initially proposed split, and why:

* jump.c is C's only. The first proposal listed it under both B and C; its one
  non-static function, `mm_jump_split`, runs only for spliced short reads
  with junction annotation (`-j`, `--pass1`), which C owns.
* The junction, jump and splice-score readers in index.c
  (`mm_idx_bed_read`, `mm_idx_bed_junc`, `mm_idx_jjump_read`,
  `mm_idx_jump_get`, `mm_idx_spsc_*`) are C's, not A's: they serve
  `--junc-bed`, `-j`, `--pass1` and `--spsc`.
* `mm_idx_alt_read` and `--alt/--alt-drop` are B's: minimap2.1 lists
  `--alt` under indexing options, but it only changes hit selection and
  MAPQ (`mm_mark_alt` in hit.c).
* `mm_mapopt_update` and `mm_idx_cal_max_occ` are B's: they turn `-f` into
  the seed occurrence threshold.
* `mm_map_frag_core` (map.c) is B's: every read goes through it, and B owns
  seeding to alignment. `mm_map`, `mm_map_file` and the thread buffers are
  D's (API); `mm_map_frag` and `mm_map_file_frag` are C's (fragments).
* `mm_split_merge` (map.c) is A's with the rest of `--split-prefix`.
* `-T` (SDUST masking of queries) is A's with sdust.c.
* `mm_set_opt` and `mm_check_opt` are D's as API and validation; the
  presets themselves are split by what they set (the two Presets tables).
* `--write-junc` stays in D (mm_write_junc is in
  format.c); it needs C's spliced inputs, which the shared simulator
  provides.
* `--seed` stays in D, although it acts in map.c: it seeds
  the hash that breaks ties between equal hits.
* kalloc.c, kthread.c and misc.c are D's: reached through `--no-kalloc`,
  `--cap-kalloc`, `-t`, `-2` and I/O errors. kalloc internals are not pinned.

ksw2.h also declares `ksw_extz`, `ksw_extd`, `ksw_extf2_sse`, `ksw_gg`, `ksw_gg2` and `ksw_gg2_sse`, which no file in the tree defines; they are not targets.

| Area | Rows |
|:--|--:|
| A | 131 |
| B | 208 |
| C | 84 |
| D | 134 |
| E | 102 |
| total | 659 |

## Area A

### A: Command-line options (9)

"Documented": whether the usage text (`minimap2 -h`) and minimap2.1 describe the option.

| Option | Argument | Documented | Notes | Tested |
|:--|:--|:--|:--|:--|
| `-w` | yes | usage, man |  | [ ] |
| `-k` | yes | usage, man |  | [ ] |
| `-I` | yes | usage, man |  | [ ] |
| `-d` | yes | usage, man |  | [ ] |
| `-T` | yes | usage (commented out) | SDUST threshold for masking queries; commented out of the usage text ("never used; might be buggy") | [ ] |
| `-H` | no | usage, man |  | [ ] |
| `--bucket-bits` | yes | - |  | [ ] |
| `--idx-no-seq` | no | man |  | [ ] |
| `--split-prefix` | yes | man | errors with --cs or --MD | [ ] |

### A: Presets: index parameters (k, w, HPC) (21)

| Preset | Documented | Notes | Tested |
|:--|:--|:--|:--|
| `lr` | - | alias of map-ont (the defaults); undocumented | [ ] |
| `map-ont` | usage, man |  | [ ] |
| `ava-ont` | usage, man |  | [ ] |
| `map10k` | - | alias of map-pb; undocumented | [ ] |
| `map-pb` | usage, man |  | [ ] |
| `ava-pb` | usage, man |  | [ ] |
| `lr:hq` | usage, man |  | [ ] |
| `map-hifi` | usage, man | also sets -A1 -B4 -O6,26 -E2,1 -s200 | [ ] |
| `map-ccs` | - | alias of map-hifi; undocumented | [ ] |
| `lr:hqae` | - | high-quality assembly evaluation (options.c comment); undocumented | [ ] |
| `map-iclr-prerender` | - | undocumented | [ ] |
| `map-iclr` | usage, man |  | [ ] |
| `asm5` | usage, man |  | [ ] |
| `asm10` | usage, man |  | [ ] |
| `asm20` | usage, man |  | [ ] |
| `short` | - | alias of sr; undocumented | [ ] |
| `sr` | usage, man |  | [ ] |
| `splice` | usage, man |  | [ ] |
| `splice:hq` | usage, man |  | [ ] |
| `splice:sr` | usage, man |  | [ ] |
| `cdna` | - | alias of splice; undocumented | [ ] |

### A: minimap.h functions (14)

| Function | Tested |
|:--|:--|
| `mm_idx_reader_open` | [ ] |
| `mm_idx_reader_read` | [ ] |
| `mm_idx_reader_close` | [ ] |
| `mm_idx_reader_eof` | [ ] |
| `mm_idx_is_idx` | [ ] |
| `mm_idx_load` | [ ] |
| `mm_idx_dump` | [ ] |
| `mm_idx_str` | [ ] |
| `mm_idx_stat` | [ ] |
| `mm_idx_destroy` | [ ] |
| `mm_idx_index_name` | [ ] |
| `mm_idx_name2id` | [ ] |
| `mm_idx_getseq` | [ ] |
| `mm_idx_build` | [ ] |

### A: minimap.h macros (4)

| Macro | Meaning | Tested |
|:--|:--|:--|
| `MM_I_HPC` | -H | [ ] |
| `MM_I_NO_SEQ` | --idx-no-seq; set when no CIGAR is needed and no -d | [ ] |
| `MM_I_NO_NAME` | mm_idx_str with name NULL | [ ] |
| `MM_IDX_MAGIC` | .mmi header | [ ] |

### A: minimap.h globals and opaque types (1)

| Name | Notes | Tested |
|:--|:--|:--|
| `mm_idx_reader_t` (`is_idx`, `n_parts`, `idx_size`, `opt`, `fp_out`, `fp`) |  | [ ] |

### A: minimap.h struct fields (29)

Fields of `mm_reg1_t` and `mm_extra_t` are results a C driver can check against the CLI or the simulator's truth; fields of the option structs name the option that sets them.

| Field | Set by / notes | Tested |
|:--|:--|:--|
| `mm128_t.x` |  | [ ] |
| `mm128_t.y` |  | [ ] |
| `mm128_v.n` |  | [ ] |
| `mm128_v.m` |  | [ ] |
| `mm128_v.a` |  | [ ] |
| `mm_idx_seq_t.name` |  | [ ] |
| `mm_idx_seq_t.offset` |  | [ ] |
| `mm_idx_seq_t.len` |  | [ ] |
| `mm_idx_seq_t.is_alt` |  | [ ] |
| `mm_idx_t.b` |  | [ ] |
| `mm_idx_t.w` |  | [ ] |
| `mm_idx_t.k` |  | [ ] |
| `mm_idx_t.flag` |  | [ ] |
| `mm_idx_t.n_seq` |  | [ ] |
| `mm_idx_t.index` |  | [ ] |
| `mm_idx_t.n_alt` |  | [ ] |
| `mm_idx_t.seq` |  | [ ] |
| `mm_idx_t.S` |  | [ ] |
| `mm_idx_t.B` | hidden buckets | [ ] |
| `mm_idx_t.km` | allocator | [ ] |
| `mm_idx_t.h` | name hash (mm_idx_index_name) | [ ] |
| `mm_idxopt_t.k` | -k | [ ] |
| `mm_idxopt_t.w` | -w | [ ] |
| `mm_idxopt_t.flag` | -H, --idx-no-seq | [ ] |
| `mm_idxopt_t.bucket_bits` | --bucket-bits | [ ] |
| `mm_idxopt_t.mini_batch_size` | no option (index building batch) | [ ] |
| `mm_idxopt_t.batch_size` | -I | [ ] |
| `mm_mapopt_t.sdust_thres` | -T | [ ] |
| `mm_mapopt_t.split_prefix` | --split-prefix | [ ] |

### A: Non-static functions (38)

Defined without `static` (prototypes in mmpriv.h, bseq.h, sdust.h, ksw2.h, kalloc.h, kthread.h). A row is ticked when a test exercises the function with a checked result, through the CLI, the API or a C driver.

| File | Function | Notes | Tested |
|:--|:--|:--|:--|
| bseq.c | `mm_bseq_open` |  | [ ] |
| bseq.c | `mm_bseq_close` |  | [ ] |
| bseq.c | `mm_bseq_read3` |  | [ ] |
| bseq.c | `mm_bseq_read2` |  | [ ] |
| bseq.c | `mm_bseq_read` |  | [ ] |
| bseq.c | `mm_bseq_read_frag2` |  | [ ] |
| bseq.c | `mm_bseq_read_frag` |  | [ ] |
| bseq.c | `mm_bseq_eof` |  | [ ] |
| sketch.c | `mm_sketch` |  | [ ] |
| sdust.c | `sdust_buf_init` |  | [ ] |
| sdust.c | `sdust_buf_destroy` |  | [ ] |
| sdust.c | `sdust_core` |  | [ ] |
| sdust.c | `sdust` |  | [ ] |
| sdust.c | `main` | the `sdust` program of `make extra` (-D_SDUST_MAIN) | [ ] |
| index.c | `mm_idx_init` |  | [ ] |
| index.c | `mm_idx_destroy` |  | [ ] |
| index.c | `mm_idx_get` |  | [ ] |
| index.c | `mm_idx_stat` |  | [ ] |
| index.c | `mm_idx_index_name` |  | [ ] |
| index.c | `mm_idx_name2id` |  | [ ] |
| index.c | `mm_idx_getseq` |  | [ ] |
| index.c | `mm_idx_getseq_rev` | not declared in any header | [ ] |
| index.c | `mm_idx_getseq2` |  | [ ] |
| index.c | `mm_idx_gen` |  | [ ] |
| index.c | `mm_idx_build` | deprecated API | [ ] |
| index.c | `mm_idx_str` | known issue: ids in its name hash are unset (do not test mm_idx_name2id on it) | [ ] |
| index.c | `mm_idx_dump` |  | [ ] |
| index.c | `mm_idx_load` |  | [ ] |
| index.c | `mm_idx_is_idx` |  | [ ] |
| index.c | `mm_idx_reader_open` |  | [ ] |
| index.c | `mm_idx_reader_close` |  | [ ] |
| index.c | `mm_idx_reader_read` |  | [ ] |
| index.c | `mm_idx_reader_eof` |  | [ ] |
| splitidx.c | `mm_split_init` |  | [ ] |
| splitidx.c | `mm_split_merge_prep` |  | [ ] |
| splitidx.c | `mm_split_rm_tmp` |  | [ ] |
| options.c | `mm_idxopt_init` |  | [ ] |
| map.c | `mm_split_merge` |  | [ ] |

### A: Programs and scripts (4)

| Target | Notes | Tested |
|:--|:--|:--|
| `sdust -w` | window size [64] | [ ] |
| `sdust -t` | score threshold [20] | [ ] |
| `sdust` (no argument) | usage, exit 1 | [ ] |
| `sdust -` (stdin) |  | [ ] |

### A: Behaviours not tied to one name (11)

| Behaviour | Notes | Tested |
|:--|:--|:--|
| FASTA and FASTQ input, multi-line FASTA, lowercase bases, IUPAC codes, U |  | [ ] |
| gzip input (target and query) |  | [ ] |
| several query files in one run (mapped one after another) | main.c loops over files unless --frag | [ ] |
| query from stdin (`-`) |  | [ ] |
| empty query file; empty sequences; reads shorter than k |  | [ ] |
| index built on the fly from FASTA equals the prebuilt .mmi for the same preset |  | [ ] |
| `.mmi` round trip: -d output, reload, byte-identical re-dump | .mmi bytes may be pinned | [ ] |
| multi-part index (-I small): parts, @SQ warning without --split-prefix |  | [ ] |
| target file that is neither FASTA nor .mmi; missing target file | error, exit 1 | [ ] |
| duplicate target names | warning from mm_idx_index_name | [ ] |
| minimizers equal a brute-force (w,k) minimizer scan with minimap2's hash64 | sketch.c; HPC with -H | [ ] |

## Area B

### B: Command-line options (48)

"Documented": whether the usage text (`minimap2 -h`) and minimap2.1 describe the option.

| Option | Argument | Documented | Notes | Tested |
|:--|:--|:--|:--|:--|
| `-D` / `--no-self` | no | man |  | [ ] |
| `-r` | yes | usage, man |  | [ ] |
| `-f` | yes | usage, man |  | [ ] |
| `-g` | yes | usage, man |  | [ ] |
| `-X` | no | usage, man |  | [ ] |
| `-s` / `--min-dp-score` | yes | usage, man |  | [ ] |
| `-c` | no | usage, man | base-level alignment; PAF cg:Z (B owns alignment, D owns the cg:Z format) | [ ] |
| `-p` | yes | usage, man |  | [ ] |
| `-M` / `--mask-level` | yes | man |  | [ ] |
| `-n` / `--min-count` | yes | usage, man |  | [ ] |
| `-z` | yes | usage, man |  | [ ] |
| `-A` | yes | usage, man |  | [ ] |
| `-B` | yes | usage, man |  | [ ] |
| `-b` | yes | man |  | [ ] |
| `-O` | yes | usage, man |  | [ ] |
| `-E` | yes | usage, man |  | [ ] |
| `-m` / `--min-chain-score` | yes | usage, man |  | [ ] |
| `-N` | yes | usage, man |  | [ ] |
| `-P` / `--all-chain` | no | man |  | [ ] |
| `-e` | yes | man |  | [ ] |
| `-U` | yes | man |  | [ ] |
| `--max-chain-skip` | yes | man |  | [ ] |
| `--min-dp-len` | yes | - |  | [ ] |
| `--no-long-join` | no | man |  | [ ] |
| `--end-bonus` | yes | man |  | [ ] |
| `--end-seed-pen` | yes | man |  | [ ] |
| `--for-only` | no | man |  | [ ] |
| `--rev-only` | no | man |  | [ ] |
| `--heap-sort` | yes | man |  | [ ] |
| `--dual` | yes | man |  | [ ] |
| `--max-clip-ratio` | yes | - |  | [ ] |
| `--min-occ-floor` | yes | man |  | [ ] |
| `--score-N` | yes | man |  | [ ] |
| `--no-end-flt` | no | man |  | [ ] |
| `--hard-mask-level` | no | man |  | [ ] |
| `--cap-sw-mem` | yes | man | sets max_sw_mat (main.c comment says --cap-sw-mat) | [ ] |
| `--max-qlen` | yes | man |  | [ ] |
| `--max-chain-iter` | yes | man |  | [ ] |
| `--chain-gap-scale` | yes | man |  | [ ] |
| `--alt` | yes | man | listed under indexing options in minimap2.1; changes hit selection (mm_mark_alt) | [ ] |
| `--alt-drop` | yes | man |  | [ ] |
| `--mask-len` | yes | man |  | [ ] |
| `--rmq` | optional | man | optional argument yes/no; errors with --sr or --splice | [ ] |
| `--qstrand` | no | - | errors with -a, -H, --frag or --splice | [ ] |
| `--q-occ-frac` | yes | man |  | [ ] |
| `--chain-skip-scale` | yes | - |  | [ ] |
| `--no-hash-name` | no | man |  | [ ] |
| `--rmq-inner` | yes | man |  | [ ] |

### B: Presets: mapping parameters (15)

| Preset | Documented | Notes | Tested |
|:--|:--|:--|:--|
| `lr` | - | alias of map-ont (the defaults); undocumented | [ ] |
| `map-ont` | usage, man |  | [ ] |
| `ava-ont` | usage, man |  | [ ] |
| `map10k` | - | alias of map-pb; undocumented | [ ] |
| `map-pb` | usage, man |  | [ ] |
| `ava-pb` | usage, man |  | [ ] |
| `lr:hq` | usage, man |  | [ ] |
| `map-hifi` | usage, man | also sets -A1 -B4 -O6,26 -E2,1 -s200 | [ ] |
| `map-ccs` | - | alias of map-hifi; undocumented | [ ] |
| `lr:hqae` | - | high-quality assembly evaluation (options.c comment); undocumented | [ ] |
| `map-iclr-prerender` | - | undocumented | [ ] |
| `map-iclr` | usage, man |  | [ ] |
| `asm5` | usage, man |  | [ ] |
| `asm10` | usage, man |  | [ ] |
| `asm20` | usage, man |  | [ ] |

### B: minimap.h functions (2)

| Function | Tested |
|:--|:--|
| `mm_mapopt_update` | [ ] |
| `mm_idx_alt_read` | [ ] |

### B: minimap.h macros (14)

| Macro | Meaning | Tested |
|:--|:--|:--|
| `MM_F_NO_DIAG` | set by -D, -X, ava presets | [ ] |
| `MM_F_NO_DUAL` | set by --dual=no, -X, ava presets | [ ] |
| `MM_F_CIGAR` | set by -c, -a, --cs, --write-junc; mappy always | [ ] |
| `MM_F_NO_LJOIN` | set by --no-long-join, -X, ava presets | [ ] |
| `MM_F_FOR_ONLY` | set by --for-only | [ ] |
| `MM_F_REV_ONLY` | set by --rev-only | [ ] |
| `MM_F_HEAP_SORT` | set by --heap-sort, sr and splice:sr presets | [ ] |
| `MM_F_ALL_CHAINS` | set by -P, -X, ava presets | [ ] |
| `MM_F_NO_END_FLT` | set by --no-end-flt | [ ] |
| `MM_F_HARD_MLEVEL` | set by --hard-mask-level | [ ] |
| `MM_F_RMQ` | set by --rmq, asm and lr:hqae presets | [ ] |
| `MM_F_QSTRAND` | set by --qstrand | [ ] |
| `MM_F_NO_INV` | set by --qstrand | [ ] |
| `MM_F_NO_HASH_NAME` | set by --no-hash-name | [ ] |

### B: minimap.h struct fields (79)

Fields of `mm_reg1_t` and `mm_extra_t` are results a C driver can check against the CLI or the simulator's truth; fields of the option structs name the option that sets them.

| Field | Set by / notes | Tested |
|:--|:--|:--|
| `mm_extra_t.capacity` |  | [ ] |
| `mm_extra_t.dp_score` |  | [ ] |
| `mm_extra_t.dp_max` |  | [ ] |
| `mm_extra_t.dp_max2` |  | [ ] |
| `mm_extra_t.dp_max0` |  | [ ] |
| `mm_extra_t.n_ambi` |  | [ ] |
| `mm_extra_t.n_cigar` |  | [ ] |
| `mm_extra_t.cigar` |  | [ ] |
| `mm_reg1_t.id` |  | [ ] |
| `mm_reg1_t.cnt` |  | [ ] |
| `mm_reg1_t.rid` |  | [ ] |
| `mm_reg1_t.score` |  | [ ] |
| `mm_reg1_t.qs` |  | [ ] |
| `mm_reg1_t.qe` |  | [ ] |
| `mm_reg1_t.rs` |  | [ ] |
| `mm_reg1_t.re` |  | [ ] |
| `mm_reg1_t.parent` |  | [ ] |
| `mm_reg1_t.subsc` |  | [ ] |
| `mm_reg1_t.as` |  | [ ] |
| `mm_reg1_t.mlen` |  | [ ] |
| `mm_reg1_t.blen` |  | [ ] |
| `mm_reg1_t.n_sub` |  | [ ] |
| `mm_reg1_t.score0` |  | [ ] |
| `mm_reg1_t.mapq` |  | [ ] |
| `mm_reg1_t.split` |  | [ ] |
| `mm_reg1_t.rev` |  | [ ] |
| `mm_reg1_t.inv` |  | [ ] |
| `mm_reg1_t.split_inv` |  | [ ] |
| `mm_reg1_t.is_alt` |  | [ ] |
| `mm_reg1_t.dummy` |  | [ ] |
| `mm_reg1_t.hash` |  | [ ] |
| `mm_reg1_t.div` |  | [ ] |
| `mm_reg1_t.p` |  | [ ] |
| `mm_mapopt_t.max_qlen` | --max-qlen | [ ] |
| `mm_mapopt_t.bw` | -r | [ ] |
| `mm_mapopt_t.bw_long` | -r ,NUM2 | [ ] |
| `mm_mapopt_t.max_gap` | -g | [ ] |
| `mm_mapopt_t.max_chain_skip` | --max-chain-skip | [ ] |
| `mm_mapopt_t.max_chain_iter` | --max-chain-iter | [ ] |
| `mm_mapopt_t.min_cnt` | -n | [ ] |
| `mm_mapopt_t.min_chain_score` | -m | [ ] |
| `mm_mapopt_t.chain_gap_scale` | --chain-gap-scale | [ ] |
| `mm_mapopt_t.chain_skip_scale` | --chain-skip-scale | [ ] |
| `mm_mapopt_t.rmq_size_cap` | no option | [ ] |
| `mm_mapopt_t.rmq_inner_dist` | --rmq-inner | [ ] |
| `mm_mapopt_t.rmq_rescue_size` | no option | [ ] |
| `mm_mapopt_t.rmq_rescue_ratio` | no option | [ ] |
| `mm_mapopt_t.mask_level` | -M | [ ] |
| `mm_mapopt_t.mask_len` | --mask-len | [ ] |
| `mm_mapopt_t.pri_ratio` | -p | [ ] |
| `mm_mapopt_t.best_n` | -N | [ ] |
| `mm_mapopt_t.alt_drop` | --alt-drop | [ ] |
| `mm_mapopt_t.a` | -A | [ ] |
| `mm_mapopt_t.b` | -B | [ ] |
| `mm_mapopt_t.q` | -O | [ ] |
| `mm_mapopt_t.e` | -E | [ ] |
| `mm_mapopt_t.q2` | -O ,INT2 | [ ] |
| `mm_mapopt_t.e2` | -E ,INT2 | [ ] |
| `mm_mapopt_t.transition` | -b | [ ] |
| `mm_mapopt_t.sc_ambi` | --score-N | [ ] |
| `mm_mapopt_t.zdrop` | -z | [ ] |
| `mm_mapopt_t.zdrop_inv` | -z ,INT2 | [ ] |
| `mm_mapopt_t.end_bonus` | --end-bonus | [ ] |
| `mm_mapopt_t.min_dp_max` | -s | [ ] |
| `mm_mapopt_t.min_ksw_len` | --min-dp-len | [ ] |
| `mm_mapopt_t.anchor_ext_len` | no option | [ ] |
| `mm_mapopt_t.anchor_ext_shift` | --end-seed-pen | [ ] |
| `mm_mapopt_t.max_clip_ratio` | --max-clip-ratio | [ ] |
| `mm_mapopt_t.rank_min_len` | no option | [ ] |
| `mm_mapopt_t.rank_frac` | no option | [ ] |
| `mm_mapopt_t.mid_occ_frac` | -f FLOAT | [ ] |
| `mm_mapopt_t.q_occ_frac` | --q-occ-frac | [ ] |
| `mm_mapopt_t.min_mid_occ` | -U, --min-occ-floor | [ ] |
| `mm_mapopt_t.max_mid_occ` | -U ,INT2 | [ ] |
| `mm_mapopt_t.mid_occ` | -f INT | [ ] |
| `mm_mapopt_t.max_occ` | -f ,INT2 | [ ] |
| `mm_mapopt_t.max_max_occ` | no option | [ ] |
| `mm_mapopt_t.occ_dist` | -e | [ ] |
| `mm_mapopt_t.max_sw_mat` | --cap-sw-mem | [ ] |

### B: Non-static functions (46)

Defined without `static` (prototypes in mmpriv.h, bseq.h, sdust.h, ksw2.h, kalloc.h, kthread.h). A row is ticked when a test exercises the function with a checked result, through the CLI, the API or a C driver.

| File | Function | Notes | Tested |
|:--|:--|:--|:--|
| index.c | `mm_idx_cal_max_occ` |  | [ ] |
| index.c | `mm_idx_alt_read` |  | [ ] |
| options.c | `mm_mapopt_update` |  | [ ] |
| seed.c | `mm_seed_mz_flt` |  | [ ] |
| seed.c | `mm_seed_collect_all` |  | [ ] |
| seed.c | `mm_seed_select` |  | [ ] |
| seed.c | `mm_collect_matches` |  | [ ] |
| lchain.c | `mg_chain_backtrack` |  | [ ] |
| lchain.c | `mg_lchain_dp` |  | [ ] |
| lchain.c | `mg_lchain_rmq` |  | [ ] |
| hit.c | `mm_gen_regs` |  | [ ] |
| hit.c | `mm_mark_alt` |  | [ ] |
| hit.c | `mm_split_reg` |  | [ ] |
| hit.c | `mm_set_parent` |  | [ ] |
| hit.c | `mm_hit_sort` |  | [ ] |
| hit.c | `mm_set_sam_pri` |  | [ ] |
| hit.c | `mm_sync_regs` |  | [ ] |
| hit.c | `mm_select_sub` |  | [ ] |
| hit.c | `mm_filter_strand_retained` |  | [ ] |
| hit.c | `mm_filter_regs` |  | [ ] |
| hit.c | `mm_squeeze_a` |  | [ ] |
| hit.c | `mm_set_mapq2` |  | [ ] |
| align.c | `mm_enlarge_cigar` |  | [ ] |
| align.c | `mm_event_identity` |  | [ ] |
| align.c | `mm_update_dp_max` |  | [ ] |
| align.c | `mm_align_skeleton` |  | [ ] |
| esterr.c | `mm_est_err` |  | [ ] |
| map.c | `mm_map_frag_core` | the mapping core every path goes through | [ ] |
| ksw2_extz2_sse.c | `ksw_extz2_sse41` |  | [ ] |
| ksw2_extz2_sse.c | `ksw_extz2_sse2` |  | [ ] |
| ksw2_extz2_sse.c | `ksw_extz2_sse` |  | [ ] |
| ksw2_extd2_sse.c | `ksw_extd2_sse41` |  | [ ] |
| ksw2_extd2_sse.c | `ksw_extd2_sse2` |  | [ ] |
| ksw2_extd2_sse.c | `ksw_extd2_sse` |  | [ ] |
| ksw2_exts2_sse.c | `ksw_exts2_sse41` |  | [ ] |
| ksw2_exts2_sse.c | `ksw_exts2_sse2` |  | [ ] |
| ksw2_exts2_sse.c | `ksw_exts2_sse` | splice DP; reached through C's splice tests as well | [ ] |
| ksw2_ll_sse.c | `ksw_ll_qinit` |  | [ ] |
| ksw2_ll_sse.c | `ksw_ll_i16` |  | [ ] |
| ksw2_dispatch.c | `__cpuidex` | Intel-derived wrapper with unknown terms: off limits, test only through dispatch | [ ] |
| ksw2_dispatch.c | `ksw_extz2_sse` |  | [ ] |
| ksw2_dispatch.c | `ksw_extd2_sse` |  | [ ] |
| ksw2_dispatch.c | `ksw_exts2_sse` | splice DP; reached through C's splice tests as well | [ ] |
| misc.c | `radix_sort_128x` | macro-generated (KRADIX_SORT_INIT) | [ ] |
| misc.c | `radix_sort_64` | macro-generated (KRADIX_SORT_INIT) | [ ] |
| misc.c | `ks_ksmall_uint32_t` | macro-generated (KSORT_INIT_GENERIC); used by mm_idx_cal_max_occ | [ ] |

### B: Behaviours not tied to one name (4)

| Behaviour | Notes | Tested |
|:--|:--|:--|
| unique reads map to their true locus; repeats get MAPQ 0 and secondaries |  | [ ] |
| AS equals the CIGAR's score under -A/-B/-O/-E for long-read presets | observed on simulated reads; sr differs | [ ] |
| inversion rescue: tp:A:I/i alignments |  | [ ] |
| chimeric reads: supplementary alignments, SA:Z |  | [ ] |

## Area C

### C: Command-line options (20)

"Documented": whether the usage text (`minimap2 -h`) and minimap2.1 describe the option.

| Option | Argument | Documented | Notes | Tested |
|:--|:--|:--|:--|:--|
| `-G` / `--max-intron-len` | yes | usage, man | only takes effect in splice mode (mm_mapopt_max_intron_len) | [ ] |
| `-u` | yes | usage, man |  | [ ] |
| `-F` | yes | usage, man |  | [ ] |
| `-C` / `--cost-non-gt-ag` | yes | man |  | [ ] |
| `-J` | yes | usage, man |  | [ ] |
| `-j` | yes | usage, man |  | [ ] |
| `--splice` | no | man |  | [ ] |
| `--sr` | optional | man |  | [ ] |
| `--frag` | yes | man |  | [ ] |
| `--no-pairing` | no | - | deprecated; same as --pairing=no | [ ] |
| `--splice-flank` | yes | man |  | [ ] |
| `--junc-bed` | yes | man |  | [ ] |
| `--junc-bonus` | yes | man |  | [ ] |
| `--spsc` | yes | man |  | [ ] |
| `--junc-pen` | yes | - |  | [ ] |
| `--pairing` | yes | man |  | [ ] |
| `--jump-min-match` | yes | man |  | [ ] |
| `--pass1` | yes | man |  | [ ] |
| `--spsc-scale` | yes | man |  | [ ] |
| `--spsc0` | yes | man | same field as --junc-pen | [ ] |

### C: Presets: mapping parameters (6)

| Preset | Documented | Notes | Tested |
|:--|:--|:--|:--|
| `short` | - | alias of sr; undocumented | [ ] |
| `sr` | usage, man |  | [ ] |
| `splice` | usage, man |  | [ ] |
| `splice:hq` | usage, man |  | [ ] |
| `splice:sr` | usage, man |  | [ ] |
| `cdna` | - | alias of splice; undocumented | [ ] |

### C: minimap.h functions (9)

| Function | Tested |
|:--|:--|
| `mm_mapopt_max_intron_len` | [ ] |
| `mm_map_frag` | [ ] |
| `mm_map_file_frag` | [ ] |
| `mm_idx_bed_read` | [ ] |
| `mm_idx_bed_junc` | [ ] |
| `mm_max_spsc_bonus` | [ ] |
| `mm_idx_spsc_read` | [ ] |
| `mm_idx_spsc_read2` | [ ] |
| `mm_idx_spsc_get` | [ ] |

### C: minimap.h macros (11)

| Macro | Meaning | Tested |
|:--|:--|:--|
| `MM_F_SPLICE` | set by --splice, splice presets, -u f/b via mm_mapopt_update | [ ] |
| `MM_F_SPLICE_FOR` | set by -u f/b, splice presets | [ ] |
| `MM_F_SPLICE_REV` | set by -u r/b, splice presets | [ ] |
| `MM_F_SR` | set by --sr, sr preset | [ ] |
| `MM_F_FRAG_MODE` | set by --frag, sr and splice:sr presets | [ ] |
| `MM_F_INDEPEND_SEG` | set by --pairing=no, --no-pairing | [ ] |
| `MM_F_SPLICE_FLANK` | set by --splice-flank, splice presets | [ ] |
| `MM_F_SPLICE_OLD` | set by -J 0 | [ ] |
| `MM_F_WEAK_PAIRING` | set by --pairing=weak, splice:sr preset | [ ] |
| `MM_F_SR_RNA` | set by --sr=rna, splice:sr preset | [ ] |
| `MM_MAX_SEG` | max segments per fragment | [ ] |

### C: minimap.h struct fields (18)

Fields of `mm_reg1_t` and `mm_extra_t` are results a C driver can check against the CLI or the simulator's truth; fields of the option structs name the option that sets them.

| Field | Set by / notes | Tested |
|:--|:--|:--|
| `mm_idx_t.I` | hidden intervals (--junc-bed); C | [ ] |
| `mm_idx_t.spsc` | hidden (--spsc); C | [ ] |
| `mm_idx_t.J` | hidden (-j, --pass1); C | [ ] |
| `mm_extra_t.trans_strand` |  | [ ] |
| `mm_reg1_t.proper_frag` |  | [ ] |
| `mm_reg1_t.pe_thru` |  | [ ] |
| `mm_reg1_t.seg_split` |  | [ ] |
| `mm_reg1_t.seg_id` |  | [ ] |
| `mm_reg1_t.strand_retained` |  | [ ] |
| `mm_reg1_t.is_spliced` |  | [ ] |
| `mm_mapopt_t.max_gap_ref` | -G (splice), presets | [ ] |
| `mm_mapopt_t.max_frag_len` | -F | [ ] |
| `mm_mapopt_t.noncan` | -C | [ ] |
| `mm_mapopt_t.junc_bonus` | --junc-bonus | [ ] |
| `mm_mapopt_t.junc_pen` | --junc-pen / --spsc0 | [ ] |
| `mm_mapopt_t.pe_ori` | presets (sr, splice:sr) | [ ] |
| `mm_mapopt_t.pe_bonus` | no option | [ ] |
| `mm_mapopt_t.jump_min_match` | --jump-min-match | [ ] |

### C: Non-static functions (17)

Defined without `static` (prototypes in mmpriv.h, bseq.h, sdust.h, ksw2.h, kalloc.h, kthread.h). A row is ticked when a test exercises the function with a checked result, through the CLI, the API or a C driver.

| File | Function | Notes | Tested |
|:--|:--|:--|:--|
| index.c | `mm_idx_bed_read` |  | [ ] |
| index.c | `mm_idx_bed_junc` |  | [ ] |
| index.c | `mm_idx_jjump_read` |  | [ ] |
| index.c | `mm_idx_jump_get` |  | [ ] |
| index.c | `mm_idx_spsc_read2` |  | [ ] |
| index.c | `mm_idx_spsc_read` |  | [ ] |
| index.c | `mm_idx_spsc_get` |  | [ ] |
| options.c | `mm_mapopt_max_intron_len` |  | [ ] |
| options.c | `mm_max_spsc_bonus` |  | [ ] |
| hit.c | `mm_seg_gen` |  | [ ] |
| hit.c | `mm_seg_free` |  | [ ] |
| jump.c | `mm_jump_split` |  | [ ] |
| map.c | `mm_map_frag` |  | [ ] |
| map.c | `mm_map_file_frag` |  | [ ] |
| pe.c | `mm_select_sub_multi` |  | [ ] |
| pe.c | `mm_set_pe_thru` |  | [ ] |
| pe.c | `mm_pair` |  | [ ] |

### C: Behaviours not tied to one name (3)

| Behaviour | Notes | Tested |
|:--|:--|:--|
| paired-end from two files, interleaved one file, /1 /2 suffixes |  | [ ] |
| frag mode: reads with the same name in one file |  | [ ] |
| spliced reads: introns, ts:A, junction BED bonus |  | [ ] |

## Area D

### D: Command-line options (33)

"Documented": whether the usage text (`minimap2 -h`) and minimap2.1 describe the option.

| Option | Argument | Documented | Notes | Tested |
|:--|:--|:--|:--|:--|
| `-2` | no | man |  | [ ] |
| `-a` / `--sam` | no | usage, man |  | [ ] |
| `-S` | no | - | deprecated: same as --cs=long with a warning | [ ] |
| `-K` / `--mb-size` | yes | usage, man |  | [ ] |
| `-t` | yes | usage, man |  | [ ] |
| `-V` / `--version` | no | usage, man | prints MM_VERSION | [ ] |
| `-v` | yes | usage (commented out) | verbosity; commented out of the usage text | [ ] |
| `-x` | yes | usage, man | preset machinery: applied before other options whatever its position; unknown preset is an error (presets in their own table) | [ ] |
| `-Q` | no | man |  | [ ] |
| `-R` | yes | usage, man |  | [ ] |
| `-h` / `--help` | no | - | --help: usage to stdout, exit 0; no arguments: usage to stderr, exit 1 | [ ] |
| `-L` | no | usage, man |  | [ ] |
| `-y` | no | usage, man |  | [ ] |
| `-Y` | no | usage, man |  | [ ] |
| `-o` | yes | usage, man |  | [ ] |
| `--seed` | yes | man | seeds the hash that breaks ties between equal hits | [ ] |
| `--no-kalloc` | no | man |  | [ ] |
| `--print-qname` | no | man | debug output to stderr | [ ] |
| `--print-seeds` | no | man | debug output to stderr; forces -t 1 | [ ] |
| `--print-aln-seq` | no | - | debug output to stderr; forces -t 1 | [ ] |
| `--secondary` | yes | man |  | [ ] |
| `--cs` | optional | usage, man |  | [ ] |
| `--MD` | no | usage, man |  | [ ] |
| `--lj-min-ratio` | yes | - | deprecated; prints a warning and is ignored | [ ] |
| `--eqx` | no | usage, man |  | [ ] |
| `--paf-no-hit` | no | man |  | [ ] |
| `--sam-hit-only` | no | man |  | [ ] |
| `--cap-kalloc` | yes | man |  | [ ] |
| `--print-chains` | no | - | debug output to stderr | [ ] |
| `--secondary-seq` | no | man |  | [ ] |
| `--ds` | no | usage |  | [ ] |
| `--write-junc` | no | man |  | [ ] |
| `--dbg-seed-occ` | no | - | debug output to stderr | [ ] |

### D: Presets: mapping parameters (1)

| Preset | Documented | Notes | Tested |
|:--|:--|:--|:--|
| unknown preset, e.g. `-x foo` or `-x asm7` | - | `[ERROR] unknown preset`, exit 1 | [ ] |

### D: minimap.h functions (11)

| Function | Tested |
|:--|:--|
| `mm_set_opt` | [ ] |
| `mm_check_opt` | [ ] |
| `mm_tbuf_init` | [ ] |
| `mm_tbuf_destroy` | [ ] |
| `mm_tbuf_get_km` | [ ] |
| `mm_map` | [ ] |
| `mm_map_file` | [ ] |
| `mm_gen_cs` | [ ] |
| `mm_gen_ds` | [ ] |
| `mm_gen_MD` | [ ] |
| `mm_mapopt_init` | [ ] |

### D: minimap.h macros (28)

| Macro | Meaning | Tested |
|:--|:--|:--|
| `MM_VERSION` | `--version`, SAM @PG VN | [ ] |
| `MM_F_OUT_SAM` | set by -a | [ ] |
| `MM_F_NO_QUAL` | set by -Q | [ ] |
| `MM_F_OUT_CG` | set by -c | [ ] |
| `MM_F_OUT_CS` | set by --cs, -S | [ ] |
| `MM_F_OUT_CS_LONG` | set by --cs=long, -S | [ ] |
| `MM_F_NO_PRINT_2ND` | set by --secondary=no, -N 0, sr and splice:sr presets | [ ] |
| `MM_F_2_IO_THREADS` | set by -2, sr and splice:sr presets | [ ] |
| `MM_F_LONG_CIGAR` | set by -L | [ ] |
| `MM_F_SOFTCLIP` | set by -Y | [ ] |
| `MM_F_OUT_MD` | set by --MD | [ ] |
| `MM_F_COPY_COMMENT` | set by -y | [ ] |
| `MM_F_EQX` | set by --eqx | [ ] |
| `MM_F_PAF_NO_HIT` | set by --paf-no-hit | [ ] |
| `MM_F_SAM_HIT_ONLY` | set by --sam-hit-only | [ ] |
| `MM_F_SECONDARY_SEQ` | set by --secondary-seq | [ ] |
| `MM_F_OUT_DS` | set by --ds | [ ] |
| `MM_F_OUT_JUNC` | set by --write-junc | [ ] |
| `MM_CIGAR_MATCH` | CIGAR operation code in mm_extra_t.cigar | [ ] |
| `MM_CIGAR_INS` | CIGAR operation code in mm_extra_t.cigar | [ ] |
| `MM_CIGAR_DEL` | CIGAR operation code in mm_extra_t.cigar | [ ] |
| `MM_CIGAR_N_SKIP` | CIGAR operation code in mm_extra_t.cigar | [ ] |
| `MM_CIGAR_SOFTCLIP` | CIGAR operation code in mm_extra_t.cigar | [ ] |
| `MM_CIGAR_HARDCLIP` | CIGAR operation code in mm_extra_t.cigar | [ ] |
| `MM_CIGAR_PADDING` | CIGAR operation code in mm_extra_t.cigar | [ ] |
| `MM_CIGAR_EQ_MATCH` | CIGAR operation code in mm_extra_t.cigar | [ ] |
| `MM_CIGAR_X_MISMATCH` | CIGAR operation code in mm_extra_t.cigar | [ ] |
| `MM_CIGAR_STR` | CIGAR operator letters | [ ] |

### D: minimap.h globals and opaque types (4)

| Name | Notes | Tested |
|:--|:--|:--|
| `mm_verbose` | -v; 3 by default in the CLI | [ ] |
| `mm_dbg_flag` | --print-*, --no-kalloc | [ ] |
| `mm_realtime0` | wall-clock origin for stderr timings | [ ] |
| `mm_tbuf_t` / `struct mm_tbuf_s` (`km`, `rep_len`, `frag_gap`) |  | [ ] |

### D: minimap.h struct fields (5)

Fields of `mm_reg1_t` and `mm_extra_t` are results a C driver can check against the CLI or the simulator's truth; fields of the option structs name the option that sets them.

| Field | Set by / notes | Tested |
|:--|:--|:--|
| `mm_reg1_t.sam_pri` |  | [ ] |
| `mm_mapopt_t.flag` | MM_F_* flags | [ ] |
| `mm_mapopt_t.seed` | --seed | [ ] |
| `mm_mapopt_t.mini_batch_size` | -K | [ ] |
| `mm_mapopt_t.cap_kalloc` | --cap-kalloc | [ ] |

### D: Non-static functions (43)

Defined without `static` (prototypes in mmpriv.h, bseq.h, sdust.h, ksw2.h, kalloc.h, kthread.h). A row is ticked when a test exercises the function with a checked result, through the CLI, the API or a C driver.

| File | Function | Notes | Tested |
|:--|:--|:--|:--|
| options.c | `mm_mapopt_init` | deprecated API (mm_set_opt(0, ...)) | [ ] |
| options.c | `mm_set_opt` |  | [ ] |
| options.c | `mm_check_opt` |  | [ ] |
| map.c | `mm_tbuf_init` |  | [ ] |
| map.c | `mm_tbuf_destroy` |  | [ ] |
| map.c | `mm_tbuf_get_km` |  | [ ] |
| map.c | `mm_map` |  | [ ] |
| map.c | `mm_map_file` |  | [ ] |
| format.c | `mm_write_sam_hdr` |  | [ ] |
| format.c | `mm_write_junc` |  | [ ] |
| format.c | `mm_gen_cs_ds_or_MD` |  | [ ] |
| format.c | `mm_gen_cs_or_MD` |  | [ ] |
| format.c | `mm_gen_cs` |  | [ ] |
| format.c | `mm_gen_ds` |  | [ ] |
| format.c | `mm_gen_MD` |  | [ ] |
| format.c | `mm_write_paf4` |  | [ ] |
| format.c | `mm_write_paf3` |  | [ ] |
| format.c | `mm_write_paf` |  | [ ] |
| format.c | `mm_write_sam3` |  | [ ] |
| format.c | `mm_write_sam2` |  | [ ] |
| format.c | `mm_write_sam` |  | [ ] |
| misc.c | `gettimeofday` | Windows only (_WIN32) | [ ] |
| misc.c | `cputime` | two definitions (Windows / POSIX) | [ ] |
| misc.c | `peakrss` | ru_maxrss units differ by platform: never assert values | [ ] |
| misc.c | `realtime` |  | [ ] |
| misc.c | `mm_err_puts` |  | [ ] |
| misc.c | `mm_err_fwrite` |  | [ ] |
| misc.c | `mm_err_fread` |  | [ ] |
| kthread.c | `kt_for` |  | [ ] |
| kthread.c | `kt_pipeline` |  | [ ] |
| kalloc.c | `km_init2` | kalloc internals are not pinned (speed work may change them) | [ ] |
| kalloc.c | `km_init` |  | [ ] |
| kalloc.c | `km_destroy` |  | [ ] |
| kalloc.c | `kfree` |  | [ ] |
| kalloc.c | `kmalloc` |  | [ ] |
| kalloc.c | `kcalloc` |  | [ ] |
| kalloc.c | `krealloc` |  | [ ] |
| kalloc.c | `krelocate` |  | [ ] |
| kalloc.c | `km_stat` |  | [ ] |
| kalloc.c | `km_stat_print` |  | [ ] |
| main.c | `liftrlimit` | Linux: raises RLIMIT_AS to the hard limit | [ ] |
| main.c | `main` | the minimap2 CLI | [ ] |
| example.c | `main` | `minimap2-lite` of `make extra` | [ ] |

### D: Programs and scripts (2)

| Target | Notes | Tested |
|:--|:--|:--|
| `minimap2-lite <target> <query>` | example.c; PAF-like lines with cg:Z | [ ] |
| `minimap2-lite` (no argument) | usage, exit 1 | [ ] |

### D: Behaviours not tied to one name (7)

| Behaviour | Notes | Tested |
|:--|:--|:--|
| SAM header: @HD, @SQ, @RG (-R), @PG with CL |  | [ ] |
| output identical across -t and -K | observed for SE, PE, splice and PAF on simulated data | [ ] |
| unmapped reads in SAM; --sam-hit-only; --paf-no-hit line format |  | [ ] |
| secondary (SEQ *), supplementary (hard clips; -Y soft) |  | [ ] |
| -o FILE equals stdout; -o to an unwritable path is an error |  | [ ] |
| option errors from mm_check_opt (each returns 1 with a message) |  | [ ] |
| libminimap2 used as example.c does (C drivers) |  | [ ] |

## Area E

### E: Programs and scripts (34)

| Target | Notes | Tested |
|:--|:--|:--|
| `paftools.js view` |  | [ ] |
| `paftools.js sam2paf` |  | [ ] |
| `paftools.js delta2paf` |  | [ ] |
| `paftools.js splice2bed` |  | [ ] |
| `paftools.js gff2bed` |  | [ ] |
| `paftools.js gff2junc` |  | [ ] |
| `paftools.js stat` |  | [ ] |
| `paftools.js asmstat` |  | [ ] |
| `paftools.js asmgene` |  | [ ] |
| `paftools.js misjoin` |  | [ ] |
| `paftools.js liftover` |  | [ ] |
| `paftools.js liftOver` | alias of liftover | [ ] |
| `paftools.js vcfpair` |  | [ ] |
| `paftools.js call` |  | [ ] |
| `paftools.js mapeval` |  | [ ] |
| `paftools.js pafcmp` |  | [ ] |
| `paftools.js bedcov` |  | [ ] |
| `paftools.js mason2fq` |  | [ ] |
| `paftools.js sim2bed` |  | [ ] |
| `paftools.js pbsim2fq` |  | [ ] |
| `paftools.js badread2fa` |  | [ ] |
| `paftools.js junceval` |  | [ ] |
| `paftools.js exoneval` |  | [ ] |
| `paftools.js ov-eval` |  | [ ] |
| `paftools.js vcfstat` |  | [ ] |
| `paftools.js sveval` |  | [ ] |
| `paftools.js vcfsel` |  | [ ] |
| `paftools.js longcs2seq` |  | [ ] |
| `paftools.js paf2gff` |  | [ ] |
| `paftools.js version` |  | [ ] |
| `paftools.js` (no command) | usage, exit 1 | [ ] |
| `pafcluster.js -c` |  | [ ] |
| `pafcluster.js -d` |  | [ ] |
| `pafcluster.js -e` |  | [ ] |

### E: mappy and python/minimap2.py (68)

| Name | Notes | Tested |
|:--|:--|:--|
| `mappy.__version__` |  | [ ] |
| `Aligner(fn_idx_in=)` |  | [ ] |
| `Aligner(preset=)` |  | [ ] |
| `Aligner(k=)` |  | [ ] |
| `Aligner(w=)` |  | [ ] |
| `Aligner(min_cnt=)` |  | [ ] |
| `Aligner(min_chain_score=)` |  | [ ] |
| `Aligner(min_dp_score=)` |  | [ ] |
| `Aligner(bw=)` |  | [ ] |
| `Aligner(bw_long=)` |  | [ ] |
| `Aligner(best_n=)` |  | [ ] |
| `Aligner(n_threads=)` |  | [ ] |
| `Aligner(fn_idx_out=)` |  | [ ] |
| `Aligner(max_frag_len=)` |  | [ ] |
| `Aligner(extra_flags=)` |  | [ ] |
| `Aligner(seq=)` |  | [ ] |
| `Aligner(scoring=)` |  | [ ] |
| `Aligner(sc_ambi=)` |  | [ ] |
| `Aligner(max_chain_skip=)` |  | [ ] |
| `Aligner.__bool__` | False when the index failed to load/build | [ ] |
| `Aligner.map(seq=)` |  | [ ] |
| `Aligner.map(seq2=)` |  | [ ] |
| `Aligner.map(name=)` |  | [ ] |
| `Aligner.map(buf=)` |  | [ ] |
| `Aligner.map(cs=)` |  | [ ] |
| `Aligner.map(ds=)` |  | [ ] |
| `Aligner.map(MD=)` |  | [ ] |
| `Aligner.map(max_frag_len=)` |  | [ ] |
| `Aligner.map(extra_flags=)` |  | [ ] |
| `Aligner.seq(name, start, end)` | None for an unknown name or empty range | [ ] |
| `Aligner.k` |  | [ ] |
| `Aligner.w` |  | [ ] |
| `Aligner.n_seq` |  | [ ] |
| `Aligner.seq_names` |  | [ ] |
| `Alignment.ctg` |  | [ ] |
| `Alignment.ctg_len` |  | [ ] |
| `Alignment.r_st` |  | [ ] |
| `Alignment.r_en` |  | [ ] |
| `Alignment.strand` |  | [ ] |
| `Alignment.trans_strand` |  | [ ] |
| `Alignment.blen` |  | [ ] |
| `Alignment.mlen` |  | [ ] |
| `Alignment.NM` |  | [ ] |
| `Alignment.is_primary` |  | [ ] |
| `Alignment.q_st` |  | [ ] |
| `Alignment.q_en` |  | [ ] |
| `Alignment.mapq` |  | [ ] |
| `Alignment.cigar` |  | [ ] |
| `Alignment.read_num` |  | [ ] |
| `Alignment.cs` |  | [ ] |
| `Alignment.ds` |  | [ ] |
| `Alignment.MD` |  | [ ] |
| `Alignment.cigar_str` |  | [ ] |
| `Alignment.__str__` | PAF-like line (q_st onward) with tp, ts, cg, cs/ds/MD | [ ] |
| `ThreadBuffer()` | reusable per-thread buffer passed as map(buf=) | [ ] |
| `fastx_read(fn, read_comment=False)` | plain or gzip FASTA/FASTQ; '-' is stdin | [ ] |
| `revcomp(seq)` |  | [ ] |
| `verbose(v=None)` | get/set mm_verbose | [ ] |
| `python/minimap2.py -x` |  | [ ] |
| `python/minimap2.py -n` |  | [ ] |
| `python/minimap2.py -m` | known issue: -m is parsed but not passed on | [ ] |
| `python/minimap2.py -r` |  | [ ] |
| `python/minimap2.py -k` |  | [ ] |
| `python/minimap2.py -w` |  | [ ] |
| `python/minimap2.py -c` |  | [ ] |
| `python/minimap2.py -d` |  | [ ] |
| `python/minimap2.py -M` |  | [ ] |
| `python/minimap2.py` usage | fewer than two arguments: usage, exit 1 | [ ] |
