# Test targets

Every target of the minimap2 test suite, listed from the source before any
area test was written: every command-line option of `minimap2`, every `-x`
preset, every public name in `minimap.h`, every non-static function in the
C files the suite covers, the `make extra` programs, `misc/paftools.js` and
`misc/pafcluster.js`, and the mappy binding. Writers tick a row by replacing
`[ ]` in its Tested column with the test that covers it, e.g.
`test_b_chain.py::test_max_chain_skip`. A row that cannot be tested
(platform-only code, or behaviour that contradicts the documentation) gets
a short reason instead. In the Tested column, `partial` marks a row only part
of whose behaviour is checked, `needs k8` a `misc/` script whose tests skip
without k8, and "a defect of base minimap2" behaviour that contradicts the
documentation, which is reported rather than tested.

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
| `-w` | yes | usage, man |  | test_a_index.py::test_mmi_k_w_H_options, test_a_index.py::test_mmi_presets, test_a_index.py::test_preset_is_applied_before_other_options |
| `-k` | yes | usage, man |  | test_a_index.py::test_mmi_k_w_H_options, test_a_index.py::test_target_errors (-k 0, -k 29) |
| `-I` | yes | usage, man |  | test_a_index.py::test_multi_part_index_file, test_a_index.py::test_multi_part_index_mapping, test_a_split.py::test_split_prefix_* (the exact-equality case is a known bug, not tested) |
| `-d` | yes | usage, man |  | test_a_index.py::test_mmi_default_index_holds_the_reference_minimizers, test_a_index.py::test_dump_while_mapping, test_a_index.py::test_mmi_round_trip, test_a_index.py::test_synth_index_bytes (unwritable path is a known bug, not tested) |
| `-T` | yes | usage (commented out) | SDUST threshold for masking queries; commented out of the usage text ("never used; might be buggy") | test_a_sdust.py::test_T_masks_low_complexity_query_minimizers |
| `-H` | no | usage, man |  | test_a_index.py::test_mmi_k_w_H_options, test_a_sketch.py::test_hpc_minimizers_on_homopolymer_rich_sequence |
| `--bucket-bits` | yes | - |  | test_a_index.py::test_bucket_bits |
| `--idx-no-seq` | no | man |  | test_a_index.py::test_idx_no_seq |
| `--split-prefix` | yes | man | errors with --cs or --MD | test_a_split.py::test_split_prefix_merges_a_multi_part_index, test_a_split.py::test_split_prefix_three_parts, test_a_split.py::test_split_prefix_paired_reads, test_a_split.py::test_split_prefix_with_a_single_part, test_a_split.py::test_split_prefix_rejects_cs_and_md, test_a_split.py::test_split_prefix_output_pin (several single-end query files: known bug, not tested) |

### A: Presets: index parameters (k, w, HPC) (21)

| Preset | Documented | Notes | Tested |
|:--|:--|:--|:--|
| `lr` | - | alias of map-ont (the defaults); undocumented | each by test_a_index.py::test_mmi_presets (header w/k/HPC against minimap2.1 or options.c, minimizers against the reference, aliases byte-identical to their targets); map-ont, map-pb, sr, asm5, splice also by test_a_index.py::test_synth_index_bytes; map-ont, map-hifi, map-pb, sr, asm5, splice by test_a_index.py::test_prebuilt_index_equals_on_the_fly |
| `map-ont` | usage, man |  | each by test_a_index.py::test_mmi_presets (header w/k/HPC against minimap2.1 or options.c, minimizers against the reference, aliases byte-identical to their targets); map-ont, map-pb, sr, asm5, splice also by test_a_index.py::test_synth_index_bytes; map-ont, map-hifi, map-pb, sr, asm5, splice by test_a_index.py::test_prebuilt_index_equals_on_the_fly |
| `ava-ont` | usage, man |  | each by test_a_index.py::test_mmi_presets (header w/k/HPC against minimap2.1 or options.c, minimizers against the reference, aliases byte-identical to their targets); map-ont, map-pb, sr, asm5, splice also by test_a_index.py::test_synth_index_bytes; map-ont, map-hifi, map-pb, sr, asm5, splice by test_a_index.py::test_prebuilt_index_equals_on_the_fly |
| `map10k` | - | alias of map-pb; undocumented | each by test_a_index.py::test_mmi_presets (header w/k/HPC against minimap2.1 or options.c, minimizers against the reference, aliases byte-identical to their targets); map-ont, map-pb, sr, asm5, splice also by test_a_index.py::test_synth_index_bytes; map-ont, map-hifi, map-pb, sr, asm5, splice by test_a_index.py::test_prebuilt_index_equals_on_the_fly |
| `map-pb` | usage, man |  | each by test_a_index.py::test_mmi_presets (header w/k/HPC against minimap2.1 or options.c, minimizers against the reference, aliases byte-identical to their targets); map-ont, map-pb, sr, asm5, splice also by test_a_index.py::test_synth_index_bytes; map-ont, map-hifi, map-pb, sr, asm5, splice by test_a_index.py::test_prebuilt_index_equals_on_the_fly |
| `ava-pb` | usage, man |  | each by test_a_index.py::test_mmi_presets (header w/k/HPC against minimap2.1 or options.c, minimizers against the reference, aliases byte-identical to their targets); map-ont, map-pb, sr, asm5, splice also by test_a_index.py::test_synth_index_bytes; map-ont, map-hifi, map-pb, sr, asm5, splice by test_a_index.py::test_prebuilt_index_equals_on_the_fly |
| `lr:hq` | usage, man |  | each by test_a_index.py::test_mmi_presets (header w/k/HPC against minimap2.1 or options.c, minimizers against the reference, aliases byte-identical to their targets); map-ont, map-pb, sr, asm5, splice also by test_a_index.py::test_synth_index_bytes; map-ont, map-hifi, map-pb, sr, asm5, splice by test_a_index.py::test_prebuilt_index_equals_on_the_fly |
| `map-hifi` | usage, man | also sets -A1 -B4 -O6,26 -E2,1 -s200 | each by test_a_index.py::test_mmi_presets (header w/k/HPC against minimap2.1 or options.c, minimizers against the reference, aliases byte-identical to their targets); map-ont, map-pb, sr, asm5, splice also by test_a_index.py::test_synth_index_bytes; map-ont, map-hifi, map-pb, sr, asm5, splice by test_a_index.py::test_prebuilt_index_equals_on_the_fly |
| `map-ccs` | - | alias of map-hifi; undocumented | each by test_a_index.py::test_mmi_presets (header w/k/HPC against minimap2.1 or options.c, minimizers against the reference, aliases byte-identical to their targets); map-ont, map-pb, sr, asm5, splice also by test_a_index.py::test_synth_index_bytes; map-ont, map-hifi, map-pb, sr, asm5, splice by test_a_index.py::test_prebuilt_index_equals_on_the_fly |
| `lr:hqae` | - | high-quality assembly evaluation (options.c comment); undocumented | each by test_a_index.py::test_mmi_presets (header w/k/HPC against minimap2.1 or options.c, minimizers against the reference, aliases byte-identical to their targets); map-ont, map-pb, sr, asm5, splice also by test_a_index.py::test_synth_index_bytes; map-ont, map-hifi, map-pb, sr, asm5, splice by test_a_index.py::test_prebuilt_index_equals_on_the_fly |
| `map-iclr-prerender` | - | undocumented | each by test_a_index.py::test_mmi_presets (header w/k/HPC against minimap2.1 or options.c, minimizers against the reference, aliases byte-identical to their targets); map-ont, map-pb, sr, asm5, splice also by test_a_index.py::test_synth_index_bytes; map-ont, map-hifi, map-pb, sr, asm5, splice by test_a_index.py::test_prebuilt_index_equals_on_the_fly |
| `map-iclr` | usage, man |  | each by test_a_index.py::test_mmi_presets (header w/k/HPC against minimap2.1 or options.c, minimizers against the reference, aliases byte-identical to their targets); map-ont, map-pb, sr, asm5, splice also by test_a_index.py::test_synth_index_bytes; map-ont, map-hifi, map-pb, sr, asm5, splice by test_a_index.py::test_prebuilt_index_equals_on_the_fly |
| `asm5` | usage, man |  | each by test_a_index.py::test_mmi_presets (header w/k/HPC against minimap2.1 or options.c, minimizers against the reference, aliases byte-identical to their targets); map-ont, map-pb, sr, asm5, splice also by test_a_index.py::test_synth_index_bytes; map-ont, map-hifi, map-pb, sr, asm5, splice by test_a_index.py::test_prebuilt_index_equals_on_the_fly |
| `asm10` | usage, man |  | each by test_a_index.py::test_mmi_presets (header w/k/HPC against minimap2.1 or options.c, minimizers against the reference, aliases byte-identical to their targets); map-ont, map-pb, sr, asm5, splice also by test_a_index.py::test_synth_index_bytes; map-ont, map-hifi, map-pb, sr, asm5, splice by test_a_index.py::test_prebuilt_index_equals_on_the_fly |
| `asm20` | usage, man |  | each by test_a_index.py::test_mmi_presets (header w/k/HPC against minimap2.1 or options.c, minimizers against the reference, aliases byte-identical to their targets); map-ont, map-pb, sr, asm5, splice also by test_a_index.py::test_synth_index_bytes; map-ont, map-hifi, map-pb, sr, asm5, splice by test_a_index.py::test_prebuilt_index_equals_on_the_fly |
| `short` | - | alias of sr; undocumented | each by test_a_index.py::test_mmi_presets (header w/k/HPC against minimap2.1 or options.c, minimizers against the reference, aliases byte-identical to their targets); map-ont, map-pb, sr, asm5, splice also by test_a_index.py::test_synth_index_bytes; map-ont, map-hifi, map-pb, sr, asm5, splice by test_a_index.py::test_prebuilt_index_equals_on_the_fly |
| `sr` | usage, man |  | each by test_a_index.py::test_mmi_presets (header w/k/HPC against minimap2.1 or options.c, minimizers against the reference, aliases byte-identical to their targets); map-ont, map-pb, sr, asm5, splice also by test_a_index.py::test_synth_index_bytes; map-ont, map-hifi, map-pb, sr, asm5, splice by test_a_index.py::test_prebuilt_index_equals_on_the_fly |
| `splice` | usage, man |  | each by test_a_index.py::test_mmi_presets (header w/k/HPC against minimap2.1 or options.c, minimizers against the reference, aliases byte-identical to their targets); map-ont, map-pb, sr, asm5, splice also by test_a_index.py::test_synth_index_bytes; map-ont, map-hifi, map-pb, sr, asm5, splice by test_a_index.py::test_prebuilt_index_equals_on_the_fly |
| `splice:hq` | usage, man |  | each by test_a_index.py::test_mmi_presets (header w/k/HPC against minimap2.1 or options.c, minimizers against the reference, aliases byte-identical to their targets); map-ont, map-pb, sr, asm5, splice also by test_a_index.py::test_synth_index_bytes; map-ont, map-hifi, map-pb, sr, asm5, splice by test_a_index.py::test_prebuilt_index_equals_on_the_fly |
| `splice:sr` | usage, man |  | each by test_a_index.py::test_mmi_presets (header w/k/HPC against minimap2.1 or options.c, minimizers against the reference, aliases byte-identical to their targets); map-ont, map-pb, sr, asm5, splice also by test_a_index.py::test_synth_index_bytes; map-ont, map-hifi, map-pb, sr, asm5, splice by test_a_index.py::test_prebuilt_index_equals_on_the_fly |
| `cdna` | - | alias of splice; undocumented | each by test_a_index.py::test_mmi_presets (header w/k/HPC against minimap2.1 or options.c, minimizers against the reference, aliases byte-identical to their targets); map-ont, map-pb, sr, asm5, splice also by test_a_index.py::test_synth_index_bytes; map-ont, map-hifi, map-pb, sr, asm5, splice by test_a_index.py::test_prebuilt_index_equals_on_the_fly |

### A: minimap.h functions (14)

| Function | Tested |
|:--|:--|
| `mm_idx_reader_open` | test_a_index.py::test_c_api_reader_single_part, test_a_index.py::test_multi_part_index_file (incl. fn_out, missing file) |
| `mm_idx_reader_read` | same |
| `mm_idx_reader_close` | same |
| `mm_idx_reader_eof` | test_a_index.py::test_multi_part_index_file (per part), test_a_index.py::test_c_api_reader_single_part |
| `mm_idx_is_idx` | test_a_index.py::test_mmi_round_trip (dump_load: index size, FASTA, empty, 3-byte, "-", missing) |
| `mm_idx_load` | test_a_index.py::test_mmi_round_trip |
| `mm_idx_dump` | test_a_index.py::test_mmi_round_trip (bytes equal `-d`) |
| `mm_idx_str` | test_a_index.py::test_c_api_idx_str_init_build_gen (name2id on it not tested: a defect of base minimap2) |
| `mm_idx_stat` | test_a_index.py::test_c_api_idx_stat (numbers only) |
| `mm_idx_destroy` | every C case; NULL in dump_load |
| `mm_idx_index_name` | test_a_index.py::test_c_api_lookup_names_and_sequences, test_a_index.py::test_duplicate_target_names |
| `mm_idx_name2id` | same |
| `mm_idx_getseq` | test_a_index.py::test_c_api_lookup_names_and_sequences, idx_str, build, gen |
| `mm_idx_build` | test_a_index.py::test_c_api_idx_str_init_build_gen |

### A: minimap.h macros (4)

| Macro | Meaning | Tested |
|:--|:--|:--|
| `MM_I_HPC` | -H | test_a_index.py::test_mmi_k_w_H_options, idx_str |
| `MM_I_NO_SEQ` | --idx-no-seq; set when no CIGAR is needed and no -d | partial. test_a_index.py::test_idx_no_seq covers --idx-no-seq; the automatic flag when no CIGAR is needed and no -d changes no output and is not observable |
| `MM_I_NO_NAME` | mm_idx_str with name NULL | test_a_index.py::test_c_api_idx_str_init_build_gen (also through dump and load) |
| `MM_IDX_MAGIC` | .mmi header | test_a_index.py::test_mmi_default_index_holds_the_reference_minimizers, helpers_a.read_mmi |

### A: minimap.h globals and opaque types (1)

| Name | Notes | Tested |
|:--|:--|:--|
| `mm_idx_reader_t` (`is_idx`, `n_parts`, `idx_size`, `opt`, `fp_out`, `fp`) |  | drv_a_index reader case (is_idx, idx_size, opt copy, n_parts, fp_out with and without a file); `fp` only indirectly |

### A: minimap.h struct fields (29)

Fields of `mm_reg1_t` and `mm_extra_t` are results a C driver can check against the CLI or the simulator's truth; fields of the option structs name the option that sets them.

| Field | Set by / notes | Tested |
|:--|:--|:--|
| `mm128_t.x` |  | K (drv_a_sketch dump: hash, span, rid, position, strand) |
| `mm128_t.y` |  | K (drv_a_sketch dump: hash, span, rid, position, strand) |
| `mm128_v.n` |  | test_a_sketch.py::test_sketch_with_kalloc_pool_and_appending (n, a; m only as capacity through appending) |
| `mm128_v.m` |  | test_a_sketch.py::test_sketch_with_kalloc_pool_and_appending (n, a; m only as capacity through appending) |
| `mm128_v.a` |  | test_a_sketch.py::test_sketch_with_kalloc_pool_and_appending (n, a; m only as capacity through appending) |
| `mm_idx_seq_t.name` |  | drv_a_index check_seqs / reader (is_alt only as 0; --alt is area B) |
| `mm_idx_seq_t.offset` |  | drv_a_index check_seqs / reader (is_alt only as 0; --alt is area B) |
| `mm_idx_seq_t.len` |  | drv_a_index check_seqs / reader (is_alt only as 0; --alt is area B) |
| `mm_idx_seq_t.is_alt` |  | drv_a_index check_seqs / reader (is_alt only as 0; --alt is area B) |
| `mm_idx_t.b` |  | drv_a_index reader, idx_str, build, gen; .mmi headers in I |
| `mm_idx_t.w` |  | drv_a_index reader, idx_str, build, gen; .mmi headers in I |
| `mm_idx_t.k` |  | drv_a_index reader, idx_str, build, gen; .mmi headers in I |
| `mm_idx_t.flag` |  | drv_a_index reader, idx_str, build, gen; .mmi headers in I |
| `mm_idx_t.n_seq` |  | drv_a_index reader, idx_str, build, gen; .mmi headers in I |
| `mm_idx_t.index` |  | drv_a_index reader, idx_str, build, gen; .mmi headers in I |
| `mm_idx_t.n_alt` |  | drv_a_index reader, idx_str, build, gen; .mmi headers in I |
| `mm_idx_t.seq` |  | drv_a_index reader, idx_str, build, gen; .mmi headers in I |
| `mm_idx_t.S` |  | drv_a_index reader, idx_str, build, gen; .mmi headers in I |
| `mm_idx_t.B` | hidden buckets | through mm_idx_get (hidden struct, not read directly) |
| `mm_idx_t.km` | allocator | partial. checked non-NULL after reading; its contents are kalloc internals |
| `mm_idx_t.h` | name hash (mm_idx_index_name) | NULL before and set after mm_idx_index_name (lookup case) |
| `mm_idxopt_t.k` | -k | test_a_index.py::test_c_api_idxopt_defaults, reader case |
| `mm_idxopt_t.w` | -w | test_a_index.py::test_c_api_idxopt_defaults, reader case |
| `mm_idxopt_t.flag` | -H, --idx-no-seq | test_a_index.py::test_c_api_idxopt_defaults, reader case |
| `mm_idxopt_t.bucket_bits` | --bucket-bits | test_a_index.py::test_c_api_idxopt_defaults, reader case |
| `mm_idxopt_t.mini_batch_size` | no option (index building batch) | idxopt defaults; gen case with mini-batches of 500 and 2^18 bases |
| `mm_idxopt_t.batch_size` | -I | test_a_index.py::test_c_api_idxopt_defaults, reader case |
| `mm_mapopt_t.sdust_thres` | -T | test_a_sdust.py::test_T_masks_low_complexity_query_minimizers (through -T) |
| `mm_mapopt_t.split_prefix` | --split-prefix | P (through --split-prefix) |

### A: Non-static functions (38)

Defined without `static` (prototypes in mmpriv.h, bseq.h, sdust.h, ksw2.h, kalloc.h, kthread.h). A row is ticked when a test exercises the function with a checked result, through the CLI, the API or a C driver.

| File | Function | Notes | Tested |
|:--|:--|:--|:--|
| bseq.c | `mm_bseq_open` |  | test_a_input.py::test_bseq_read_fasta_and_fastq, test_a_input.py::test_bseq_read_frag_mode_keeps_mates_together |
| bseq.c | `mm_bseq_close` |  | test_a_input.py::test_bseq_read_fasta_and_fastq, test_a_input.py::test_bseq_read_frag_mode_keeps_mates_together |
| bseq.c | `mm_bseq_read3` |  | test_a_input.py::test_bseq_read_fasta_and_fastq, test_a_input.py::test_bseq_read_frag_mode_keeps_mates_together |
| bseq.c | `mm_bseq_read2` |  | test_a_input.py::test_bseq_read_fasta_and_fastq, test_a_input.py::test_bseq_read_frag_mode_keeps_mates_together |
| bseq.c | `mm_bseq_read` |  | test_a_input.py::test_bseq_read_fasta_and_fastq, test_a_input.py::test_bseq_read_frag_mode_keeps_mates_together |
| bseq.c | `mm_bseq_read_frag2` |  | test_a_input.py::test_bseq_read_frag_two_files |
| bseq.c | `mm_bseq_read_frag` |  | test_a_input.py::test_bseq_read_frag_two_files |
| bseq.c | `mm_bseq_eof` |  | test_a_input.py::test_bseq_read_fasta_and_fastq, test_a_input.py::test_bseq_read_frag_mode_keeps_mates_together |
| sketch.c | `mm_sketch` |  | K (all tests); through every .mmi in I; test_a_sketch.py::test_palindromic_kmers_are_never_minimizers, test_a_sketch.py::test_hpc_kmers_spanning_256_bases_or_more |
| sdust.c | `sdust_buf_init` |  | test_a_sdust.py::test_sdust_api_matches_reference |
| sdust.c | `sdust_buf_destroy` |  | test_a_sdust.py::test_sdust_api_matches_reference |
| sdust.c | `sdust_core` |  | test_a_sdust.py::test_sdust_api_matches_reference |
| sdust.c | `sdust` |  | test_a_sdust.py::test_sdust_api_matches_reference |
| sdust.c | `main` | the `sdust` program of `make extra` (-D_SDUST_MAIN) | test_a_sdust.py::test_sdust_program_matches_reference, test_a_sdust.py::test_sdust_usage_stdin_and_gzip |
| index.c | `mm_idx_init` |  | idx_str case (w < 1, b cap), and every build |
| index.c | `mm_idx_destroy` |  | every C case; NULL in dump_load |
| index.c | `mm_idx_get` |  | drv_a_index check_lookup (lookup, dump_load, idx_str, build, gen) |
| index.c | `mm_idx_stat` |  | test_a_index.py::test_c_api_idx_stat (numbers only) |
| index.c | `mm_idx_index_name` |  | test_a_index.py::test_c_api_lookup_names_and_sequences, test_a_index.py::test_duplicate_target_names |
| index.c | `mm_idx_name2id` |  | same |
| index.c | `mm_idx_getseq` |  | test_a_index.py::test_c_api_lookup_names_and_sequences, idx_str, build, gen |
| index.c | `mm_idx_getseq_rev` | not declared in any header | drv_a_index check_seqs |
| index.c | `mm_idx_getseq2` |  | drv_a_index check_seqs |
| index.c | `mm_idx_gen` |  | test_a_index.py::test_c_api_idx_str_init_build_gen (gen case: parts, mini-batches, NULL file, end of file) |
| index.c | `mm_idx_build` | deprecated API | test_a_index.py::test_c_api_idx_str_init_build_gen |
| index.c | `mm_idx_str` | known issue: ids in its name hash are unset (do not test mm_idx_name2id on it) | test_a_index.py::test_c_api_idx_str_init_build_gen (name2id on it not tested: a defect of base minimap2) |
| index.c | `mm_idx_dump` |  | test_a_index.py::test_mmi_round_trip (bytes equal `-d`) |
| index.c | `mm_idx_load` |  | test_a_index.py::test_mmi_round_trip |
| index.c | `mm_idx_is_idx` |  | test_a_index.py::test_mmi_round_trip (dump_load: index size, FASTA, empty, 3-byte, "-", missing) |
| index.c | `mm_idx_reader_open` |  | test_a_index.py::test_c_api_reader_single_part, test_a_index.py::test_multi_part_index_file (incl. fn_out, missing file) |
| index.c | `mm_idx_reader_close` |  | same |
| index.c | `mm_idx_reader_read` |  | same |
| index.c | `mm_idx_reader_eof` |  | test_a_index.py::test_multi_part_index_file (per part), test_a_index.py::test_c_api_reader_single_part |
| splitidx.c | `mm_split_init` |  | test_a_split.py::test_c_api_split_files; through --split-prefix in P |
| splitidx.c | `mm_split_merge_prep` |  | test_a_split.py::test_c_api_split_files; through --split-prefix in P |
| splitidx.c | `mm_split_rm_tmp` |  | test_a_split.py::test_c_api_split_files; through --split-prefix in P |
| options.c | `mm_idxopt_init` |  | test_a_index.py::test_c_api_idxopt_defaults |
| map.c | `mm_split_merge` |  | P (merged output against the single-part index, two and three parts, paired) |

### A: Programs and scripts (4)

| Target | Notes | Tested |
|:--|:--|:--|
| `sdust -w` | window size [64] | test_a_sdust.py::test_sdust_program_matches_reference |
| `sdust -t` | score threshold [20] | same |
| `sdust` (no argument) | usage, exit 1 | test_a_sdust.py::test_sdust_usage_stdin_and_gzip |
| `sdust -` (stdin) |  | same |

### A: Behaviours not tied to one name (11)

| Behaviour | Notes | Tested |
|:--|:--|:--|
| FASTA and FASTQ input, multi-line FASTA, lowercase bases, IUPAC codes, U |  | test_a_input.py::test_read_file_forms_map_identically, test_a_input.py::test_target_file_forms_give_the_same_index, test_a_input.py::test_read_case_u_and_iupac, test_a_input.py::test_bseq_read_fasta_and_fastq, test_a_sketch.py::test_case_and_u_give_the_same_minimizers |
| gzip input (target and query) |  | test_a_input.py::test_read_file_forms_map_identically, test_a_input.py::test_target_file_forms_give_the_same_index, test_a_input.py::test_query_and_target_from_stdin (gzip on stdin) |
| several query files in one run (mapped one after another) | main.c loops over files unless --frag | test_a_input.py::test_several_query_files |
| query from stdin (`-`) |  | test_a_input.py::test_query_and_target_from_stdin (also the target) |
| empty query file; empty sequences; reads shorter than k |  | test_a_input.py::test_empty_short_and_unmappable_queries (the SAM line of an empty sequence or name is invalid: a defect of base minimap2, not tested) |
| index built on the fly from FASTA equals the prebuilt .mmi for the same preset |  | test_a_index.py::test_prebuilt_index_equals_on_the_fly, test_a_index.py::test_multi_part_index_mapping |
| `.mmi` round trip: -d output, reload, byte-identical re-dump | .mmi bytes may be pinned | for content (header, names, lengths, sequence, minimizer table) by test_a_index.py::test_mmi_round_trip; a re-dump is not byte-identical in base minimap2 (hash tables are written in table order, which changes on reload), and no document promises byte identity, so bytes are compared only between `-d` and mm_idx_dump of a freshly built index, and pinned for `-d` (test_a_index.py::test_synth_index_bytes) |
| multi-part index (-I small): parts, @SQ warning without --split-prefix |  | test_a_index.py::test_multi_part_index_file, test_a_index.py::test_multi_part_index_mapping (no @SQ, warning) |
| target file that is neither FASTA nor .mmi; missing target file | error, exit 1 | test_a_index.py::test_target_errors (a non-FASTA file is read as an index of no sequence and exit 0: the test requires only that no read maps) |
| duplicate target names | warning from mm_idx_index_name | test_a_index.py::test_duplicate_target_names (duplicate @SQ lines without a warning: a defect of base minimap2) |
| minimizers equal a brute-force (w,k) minimizer scan with minimap2's hash64 | sketch.c; HPC with -H | test_a_sketch.py::test_minimizers_equal_the_reference_scan (exact, 14 (w,k,H) combinations), test_a_sketch.py::test_minimizers_with_tied_kmers, test_a_sketch.py::test_minimizers_around_ambiguous_bases, test_a_index.py::test_mmi_* (index tables against the reference); first-window and before-N edges are a defect of base minimap2, checked only as properties |

## Area B

### B: Command-line options (48)

"Documented": whether the usage text (`minimap2 -h`) and minimap2.1 describe the option.

| Option | Argument | Documented | Notes | Tested |
|:--|:--|:--|:--|:--|
| `-D` / `--no-self` | no | man |  | test_b_ava.py::test_D_drops_self_mappings, test_X_equals_its_expansion; test_b_api.py::test_api_fields_and_flags_match_command_line[-P -D] |
| `-r` | yes | usage, man |  | test_b_chain.py::test_bandwidth_r_and_long_join (exact boundary -r 599,599 vs 600,600; NUM2 long join), test_band_edge; test_b_ksw.py::test_band_limits_global_alignment; test_b_api.py::test_api_fields_and_flags_match_command_line; test_b_ksw.py::test_narrow_band_edge (pinned) |
| `-f` | yes | usage, man |  | test_b_seed.py::test_f_int_sets_occurrence_threshold (-f INT via rl:i), test_occurrence_bounds_U, test_cal_max_occ, test_mapopt_update_occurrence_threshold (-f FLOAT); test_b_api.py::test_api_fields_and_flags_match_command_line; test_b_seed.py::test_many_copy_occurrence_thresholds, test_b_seed.py::test_collect_matches_many_copies (120-copy family) |
| `-g` | yes | usage, man |  | test_b_chain.py::test_max_gap_g_boundary (exact boundary from brute-force anchors); test_chains_never_cross_contigs_strands_or_max_gap; test_b_api.py::test_api_fields_and_flags_match_command_line |
| `-X` | no | usage, man |  | test_b_ava.py::test_X_equals_its_expansion; test_b_presets.py (ava presets) |
| `-s` / `--min-dp-score` | yes | usage, man |  | test_b_align.py::test_peak_score_ms_and_s_boundary (-s ms keeps, ms+1 drops); test_b_api.py::test_api_fields_and_flags_match_command_line |
| `-c` | no | usage, man | base-level alignment; PAF cg:Z (B owns alignment, D owns the cg:Z format) | test_b_align.py (every test: cg:Z, NM, columns 10/11, nn, de recomputed; AS = CIGAR score); test_b_api.py::test_api_fields_and_flags_match_command_line[-c] |
| `-p` | yes | usage, man |  | test_b_hit.py::test_pri_ratio_p (boundary at the brute-force S2/S1); test_b_api.py::test_api_fields_and_flags_match_command_line |
| `-M` / `--mask-level` | yes | man |  | test_b_hit.py::test_mask_level_M (boundary from brute-force chain coordinates); test_b_api.py::test_api_fields_and_flags_match_command_line |
| `-n` / `--min-count` | yes | usage, man |  | test_b_chain.py::test_min_chain_score_and_count_boundaries, test_min_count_and_score_filters; test_b_api.py::test_api_fields_and_flags_match_command_line |
| `-z` | yes | usage, man |  | test_b_align.py::test_zdrop_boundary (exact boundary at the derived Z-drop), test_inversion_rescue (-z INT1,INT2); test_b_ksw.py::test_zdrop_stops_extension; test_b_api.py::test_api_fields_and_flags_match_command_line; test_b_ksw.py::test_zdrop_threshold_is_strict |
| `-A` | yes | usage, man |  | test_b_align.py::test_AS_is_cigar_score_under_custom_scoring; test_b_api.py::test_api_fields_and_flags_match_command_line |
| `-B` | yes | usage, man |  | test_b_align.py::test_AS_is_cigar_score_under_custom_scoring; test_b_api.py::test_api_fields_and_flags_match_command_line |
| `-b` | yes | man |  | test_b_align.py::test_AS_is_cigar_score_under_custom_scoring[transition, transition-equal-B], map-iclr presets; test_b_ksw.py::test_generic_matrix_transitions |
| `-O` | yes | usage, man |  | test_b_align.py::test_AS_is_cigar_score_under_custom_scoring (one value: affine via ksw_extz2; two values: two-piece) |
| `-E` | yes | usage, man |  | test_b_align.py::test_AS_is_cigar_score_under_custom_scoring |
| `-m` / `--min-chain-score` | yes | usage, man |  | test_b_chain.py::test_min_chain_score_and_count_boundaries (s1 from brute force); test_b_ava.py::test_ava_min_chain_score; test_b_api.py::test_api_fields_and_flags_match_command_line |
| `-N` | yes | usage, man |  | test_b_hit.py::test_repeat_reads_secondaries_and_N (-N 2, 1, 0 with warning); test_b_api.py::test_api_fields_and_flags_match_command_line |
| `-P` / `--all-chain` | no | man |  | test_b_hit.py::test_all_chains_P, test_repeat_reads_secondaries_and_N; test_b_ava.py::test_D_drops_self_mappings |
| `-e` | yes | man |  | test_b_seed.py::test_e_samples_repetitive_minimizers, test_collect_matches_against_brute_force; test_b_api.py::test_api_fields_and_flags_match_command_line; test_b_seed.py::test_many_copy_occurrence_thresholds, test_b_seed.py::test_sampled_minimizers_pinned |
| `-U` | yes | man |  | test_b_seed.py::test_occurrence_bounds_U (INT2 > INT1 only: INT2 <= INT1 is a defect of base minimap2, not tested); test_b_seed.py::test_many_copy_occurrence_thresholds |
| `--max-chain-skip` | yes | man |  | test_b_chain.py::test_random_anchor_sets_against_brute_force (off = optimum; 0 loses score), test_chain_score_and_count_substitutions[no-heuristics]; test_b_api.py::test_api_fields_and_flags_match_command_line; test_b_chain.py::test_chaining_heuristics_pinned |
| `--min-dp-len` | yes | - |  | partial; test_b_align.py::test_alignment_options_keep_reads_correct (reads at truth, AS = CIGAR score at 10 and 5000); test_b_api.py::test_api_fields_and_flags_match_command_line; no independently derived exact effect |
| `--no-long-join` | no | man |  | test_b_chain.py::test_bandwidth_r_and_long_join, test_max_gap_g_boundary; test_b_api.py::test_api_fields_and_flags_match_command_line |
| `--end-bonus` | yes | man |  | test_b_align.py::test_end_bonus_boundary (boundary at the brute-force extension drop); test_b_ksw.py::test_extension_matches_reference; test_b_align.py::test_end_bonus_boundary_left |
| `--end-seed-pen` | yes | man |  | test_b_align.py::test_end_seed_pen_boundary (boundary at floor(s - ln g), s brute-force local score) |
| `--for-only` | no | man |  | test_b_ava.py::test_for_only_rev_only; test_b_api.py::test_api_fields_and_flags_match_command_line |
| `--rev-only` | no | man |  | test_b_ava.py::test_for_only_rev_only; test_b_api.py::test_api_fields_and_flags_match_command_line |
| `--heap-sort` | yes | man |  | test_b_ava.py::test_heap_sort_same_output (yes/no identical SAM); test_b_chain.py::test_chain_score_and_count_substitutions[heap]; test_b_api.py::test_api_fields_and_flags_match_command_line |
| `--dual` | yes | man |  | test_b_ava.py::test_ava_ont_reports_true_overlaps (no), test_ava_dual_yes_reports_both_directions (yes) |
| `--max-clip-ratio` | yes | - |  | test_b_align.py::test_max_clip_ratio_boundary (0.333 vs 0.334 on a 1000/1000/1000 read; one-sided clipping kept); test_b_api.py::test_api_fields_and_flags_match_command_line |
| `--min-occ-floor` | yes | man |  | test_b_seed.py::test_occurrence_bounds_U[--min-occ-floor] |
| `--score-N` | yes | man |  | test_b_align.py::test_AS_is_cigar_score_under_custom_scoring[score-N] (3), preset tests (1); --score-N 0 is a defect of base minimap2, not tested |
| `--no-end-flt` | no | man |  | partial; test_b_align.py::test_alignment_options_keep_reads_correct; test_b_api.py::test_api_fields_and_flags_match_command_line; no independently derived exact effect |
| `--hard-mask-level` | no | man |  | test_b_hit.py::test_hard_mask_level; test_b_api.py::test_api_fields_and_flags_match_command_line |
| `--cap-sw-mem` | yes | man | sets max_sw_mat (main.c comment says --cap-sw-mat) | test_b_align.py::test_cap_sw_mem (1: no alignment; 0 and 1g: default output); test_b_api.py::test_api_fields_and_flags_match_command_line |
| `--max-qlen` | yes | man |  | test_b_ava.py::test_max_qlen (L vs L-1, PAF and SAM); test_b_api.py::test_api_fields_and_flags_match_command_line |
| `--max-chain-iter` | yes | man |  | test_b_chain.py::test_max_chain_iter, test_random_anchor_sets_against_brute_force; test_b_api.py::test_api_fields_and_flags_match_command_line; test_b_chain.py::test_chaining_heuristics_pinned |
| `--chain-gap-scale` | yes | man |  | test_b_chain.py::test_chain_score_with_indels, test_gap_cost_formula; test_b_api.py::test_api_fields_and_flags_match_command_line; default (0.8 vs documented 1.0) is a defect of base minimap2, not tested |
| `--alt` | yes | man | listed under indexing options in minimap2.1; changes hit selection (mm_mark_alt) | test_b_hit.py::test_alt_contigs; test_b_seed.py::test_alt_read_api |
| `--alt-drop` | yes | man |  | test_b_hit.py::test_alt_contigs (0 and default 0.15) |
| `--mask-len` | yes | man |  | test_b_hit.py::test_mask_len (boundary at the brute-force uncovered length); test_b_api.py::test_api_fields_and_flags_match_command_line |
| `--rmq` | optional | man | optional argument yes/no; errors with --sr or --splice | test_b_chain.py::test_chain_score_and_count_substitutions[rmq], test_gap_cost_formula (rmq); test_b_presets.py (asm presets, lr:hqae); test_b_api.py::test_api_fields_and_flags_match_command_line |
| `--qstrand` | no | - | errors with -a, -H, --frag or --splice | test_b_ava.py::test_qstrand (coordinates, cg/cs on the query strand, rejected with -a/-H/--frag/--splice); test_b_api.py::test_api_fields_and_flags_match_command_line |
| `--q-occ-frac` | yes | man |  | test_b_seed.py::test_q_occ_frac_drops_query_repeats, test_query_occurrence_filter; test_b_api.py::test_api_fields_and_flags_match_command_line |
| `--chain-skip-scale` | yes | - |  | test_b_chain.py::test_chain_score_with_indels, test_gap_cost_formula; test_b_api.py::test_api_fields_and_flags_match_command_line |
| `--no-hash-name` | no | man |  | test_b_hit.py::test_no_hash_name; test_b_api.py::test_api_fields_and_flags_match_command_line |
| `--rmq-inner` | yes | man |  | test_b_chain.py::test_random_anchor_sets_against_brute_force (1000 helps over 0), test_band_edge, test_chain_score_and_count_substitutions[rmq-inner0]; test_b_api.py::test_api_fields_and_flags_match_command_line; test_b_chain.py::test_chaining_heuristics_pinned |

### B: Presets: mapping parameters (15)

| Preset | Documented | Notes | Tested |
|:--|:--|:--|:--|
| `lr` | - | alias of map-ont (the defaults); undocumented | test_b_presets.py::test_preset_aliases; test_b_api.py::test_preset_mode_flags |
| `map-ont` | usage, man |  | test_b_presets.py::test_preset_equals_documented_options (minimap2.1 expansion, with probe reads sensitive to -g/-r/-z); test_b_presets.py::test_preset_reads_at_truth; test_b_align.py::test_AS_is_cigar_score_under_preset_scoring; pinned output |
| `ava-ont` | usage, man |  | test_b_presets.py::test_preset_equals_documented_options (minimap2.1 expansion, with probe reads sensitive to -g/-r/-z); test_b_ava.py::test_ava_ont_reports_true_overlaps, test_ava_min_chain_score; test_b_align.py; pinned output |
| `map10k` | - | alias of map-pb; undocumented | test_b_presets.py::test_preset_aliases; test_b_api.py::test_preset_mode_flags |
| `map-pb` | usage, man |  | test_b_presets.py::test_preset_equals_documented_options (minimap2.1 expansion, with probe reads sensitive to -g/-r/-z); test_b_presets.py::test_preset_reads_at_truth; test_b_align.py::test_AS_is_cigar_score_under_preset_scoring; pinned output |
| `ava-pb` | usage, man |  | test_b_presets.py::test_preset_equals_documented_options (minimap2.1 expansion, with probe reads sensitive to -g/-r/-z); test_b_align.py; pinned output |
| `lr:hq` | usage, man |  | test_b_presets.py::test_preset_equals_documented_options (minimap2.1 expansion, with probe reads sensitive to -g/-r/-z); test_b_presets.py::test_preset_reads_at_truth; test_b_align.py::test_AS_is_cigar_score_under_preset_scoring; pinned output |
| `map-hifi` | usage, man | also sets -A1 -B4 -O6,26 -E2,1 -s200 | test_b_presets.py::test_preset_equals_documented_options (minimap2.1 expansion, with probe reads sensitive to -g/-r/-z); test_b_presets.py::test_preset_reads_at_truth; test_b_align.py::test_AS_is_cigar_score_under_preset_scoring; pinned output |
| `map-ccs` | - | alias of map-hifi; undocumented | test_b_presets.py::test_preset_aliases; test_b_api.py::test_preset_mode_flags |
| `lr:hqae` | - | high-quality assembly evaluation (options.c comment); undocumented | test_b_presets.py::test_preset_equals_documented_options (expansion from options.c; undocumented), test_preset_reads_at_truth; test_b_align.py::test_AS_is_cigar_score_under_preset_scoring |
| `map-iclr-prerender` | - | undocumented | test_b_presets.py::test_preset_equals_documented_options (expansion from options.c; undocumented), test_preset_reads_at_truth; test_b_align.py::test_AS_is_cigar_score_under_preset_scoring |
| `map-iclr` | usage, man |  | test_b_presets.py::test_preset_equals_documented_options (minimap2.1 expansion, with probe reads sensitive to -g/-r/-z); test_b_presets.py::test_preset_reads_at_truth; test_b_align.py::test_AS_is_cigar_score_under_preset_scoring; pinned output |
| `asm5` | usage, man |  | test_b_presets.py::test_preset_equals_documented_options (minimap2.1 expansion, with probe reads sensitive to -g/-r/-z); test_b_presets.py::test_asm_contigs_at_truth; test_b_align.py::test_AS_is_cigar_score_under_preset_scoring; pinned output |
| `asm10` | usage, man |  | test_b_presets.py::test_preset_equals_documented_options (minimap2.1 expansion, with probe reads sensitive to -g/-r/-z); test_b_presets.py::test_asm_contigs_at_truth; test_b_align.py::test_AS_is_cigar_score_under_preset_scoring; pinned output |
| `asm20` | usage, man |  | test_b_presets.py::test_preset_equals_documented_options (minimap2.1 expansion, with probe reads sensitive to -g/-r/-z); test_b_presets.py::test_asm_contigs_at_truth; test_b_align.py::test_AS_is_cigar_score_under_preset_scoring; pinned output |

### B: minimap.h functions (2)

| Function | Tested |
|:--|:--|
| `mm_mapopt_update` | test_b_seed.py::test_mapopt_update_occurrence_threshold (drv_b_seed mapopt_update) |
| `mm_idx_alt_read` | test_b_seed.py::test_alt_read_api; test_b_hit.py::test_alt_contigs |

### B: minimap.h macros (14)

| Macro | Meaning | Tested |
|:--|:--|:--|
| `MM_F_NO_DIAG` | set by -D, -X, ava presets | -D tests; test_b_api.py::test_preset_mode_flags; test_b_api.py::test_api_fields_and_flags_match_command_line[-P -D] |
| `MM_F_NO_DUAL` | set by --dual=no, -X, ava presets | test_b_ava.py::test_ava_dual_yes_reports_both_directions; test_b_api.py::test_preset_mode_flags |
| `MM_F_CIGAR` | set by -c, -a, --cs, --write-junc; mappy always | test_b_api.py::test_api_fields_and_flags_match_command_line[-c] (p set iff MM_F_CIGAR); all -c tests |
| `MM_F_NO_LJOIN` | set by --no-long-join, -X, ava presets | test_b_api.py::test_api_fields_and_flags_match_command_line[--no-long-join -c]; test_b_chain.py::test_bandwidth_r_and_long_join |
| `MM_F_FOR_ONLY` | set by --for-only | test_b_api.py::test_api_fields_and_flags_match_command_line[--for-only -c]; test_b_ava.py::test_for_only_rev_only |
| `MM_F_REV_ONLY` | set by --rev-only | test_b_api.py::test_api_fields_and_flags_match_command_line[--rev-only]; test_b_ava.py::test_for_only_rev_only |
| `MM_F_HEAP_SORT` | set by --heap-sort, sr and splice:sr presets | test_b_api.py::test_api_fields_and_flags_match_command_line[--heap-sort=yes --no-hash-name]; test_b_api.py::test_preset_mode_flags[sr, splice:sr] |
| `MM_F_ALL_CHAINS` | set by -P, -X, ava presets | test_b_api.py::test_api_fields_and_flags_match_command_line[-P -D]; test_b_api.py::test_preset_mode_flags[ava-*] |
| `MM_F_NO_END_FLT` | set by --no-end-flt | test_b_api.py::test_api_fields_and_flags_match_command_line[--no-end-flt ...] |
| `MM_F_HARD_MLEVEL` | set by --hard-mask-level | test_b_api.py::test_api_fields_and_flags_match_command_line[--hard-mask-level ...]; test_b_hit.py::test_hard_mask_level |
| `MM_F_RMQ` | set by --rmq, asm and lr:hqae presets | test_b_api.py::test_api_fields_and_flags_match_command_line[--rmq ...]; test_b_api.py::test_preset_mode_flags[asm*, lr:hqae] |
| `MM_F_QSTRAND` | set by --qstrand | test_b_api.py::test_api_fields_and_flags_match_command_line[--qstrand -c]; test_b_ava.py::test_qstrand |
| `MM_F_NO_INV` | set by --qstrand | test_b_api.py::test_api_fields_and_flags_match_command_line[--qstrand -c] (set with QSTRAND as --qstrand does) |
| `MM_F_NO_HASH_NAME` | set by --no-hash-name | test_b_api.py::test_api_fields_and_flags_match_command_line[--heap-sort=yes --no-hash-name]; test_b_hit.py::test_no_hash_name |

### B: minimap.h struct fields (79)

Fields of `mm_reg1_t` and `mm_extra_t` are results a C driver can check against the CLI or the simulator's truth; fields of the option structs name the option that sets them.

| Field | Set by / notes | Tested |
|:--|:--|:--|
| `mm_extra_t.capacity` |  | drv_b_api checks capacity >= n_cigar + header words |
| `mm_extra_t.dp_score` |  | test_b_api.py::test_api_fields_and_flags_match_command_line (= AS:i) |
| `mm_extra_t.dp_max` |  | test_b_api.py::test_rank_min_len_and_rank_frac (re-estimated vs dp_max0) |
| `mm_extra_t.dp_max2` |  | test_b_api.py::test_api_result_fields_consistent (<= dp_max0 of primaries) |
| `mm_extra_t.dp_max0` |  | test_b_api.py::test_api_fields_and_flags_match_command_line (= ms:i) |
| `mm_extra_t.n_ambi` |  | test_b_api.py::test_api_fields_and_flags_match_command_line (= nn:i) |
| `mm_extra_t.n_cigar` |  | test_b_api.py::test_api_fields_and_flags_match_command_line (CIGAR printed from n_cigar ops = cg:Z) |
| `mm_extra_t.cigar` |  | test_b_api.py::test_api_fields_and_flags_match_command_line (= cg:Z) |
| `mm_reg1_t.id` |  | test_b_api.py::test_api_fields_and_flags_match_command_line (drv_b_api: id == index; with parent gives tp:A) |
| `mm_reg1_t.cnt` |  | test_b_api.py::test_api_fields_and_flags_match_command_line (= cm:i) |
| `mm_reg1_t.rid` |  | test_b_api.py::test_api_fields_and_flags_match_command_line (= target name) |
| `mm_reg1_t.score` |  | test_b_api.py::test_api_fields_and_flags_match_command_line (= s1:i) |
| `mm_reg1_t.qs` |  | test_b_api.py::test_api_fields_and_flags_match_command_line (= PAF col 3) |
| `mm_reg1_t.qe` |  | test_b_api.py::test_api_fields_and_flags_match_command_line (= PAF col 4) |
| `mm_reg1_t.rs` |  | test_b_api.py::test_api_fields_and_flags_match_command_line (= PAF col 8 (qstrand: converted)) |
| `mm_reg1_t.re` |  | test_b_api.py::test_api_fields_and_flags_match_command_line (= PAF col 9) |
| `mm_reg1_t.parent` |  | test_b_api.py::test_api_fields_and_flags_match_command_line (tp:A (id == parent)) |
| `mm_reg1_t.subsc` |  | test_b_api.py::test_api_fields_and_flags_match_command_line (= s2:i) |
| `mm_reg1_t.as` |  | not covered; internal anchor offset with no documented result |
| `mm_reg1_t.mlen` |  | test_b_api.py::test_api_fields_and_flags_match_command_line (= PAF col 10) |
| `mm_reg1_t.blen` |  | test_b_api.py::test_api_fields_and_flags_match_command_line (= PAF col 11) |
| `mm_reg1_t.n_sub` |  | partial; printed by drv_b_api, enters MAPQ only; no documented value |
| `mm_reg1_t.score0` |  | test_b_api.py::test_api_result_fields_consistent (== score unless split, larger when split) |
| `mm_reg1_t.mapq` |  | test_b_api.py::test_api_fields_and_flags_match_command_line (= PAF col 12) |
| `mm_reg1_t.split` |  | test_b_api.py::test_api_fields_and_flags_match_command_line (= zd:i) |
| `mm_reg1_t.rev` |  | test_b_api.py::test_api_fields_and_flags_match_command_line (= PAF strand) |
| `mm_reg1_t.inv` |  | test_b_api.py::test_api_fields_and_flags_match_command_line (tp:A:I/i) |
| `mm_reg1_t.split_inv` |  | partial; printed; inversions tested through tp:A:I |
| `mm_reg1_t.is_alt` |  | partial; printed; ALT effect tested in test_b_hit::test_alt_contigs |
| `mm_reg1_t.dummy` |  | not covered; unused bits |
| `mm_reg1_t.hash` |  | not covered; tie-breaking hash; its effect is tested through --no-hash-name |
| `mm_reg1_t.div` |  | test_b_api.py::test_api_fields_and_flags_match_command_line (= dv:f (test_b_chain.py::test_chain_score_and_count_substitutions); range checked in test_b_api.py::test_api_result_fields_consistent) |
| `mm_reg1_t.p` |  | test_b_api.py::test_api_fields_and_flags_match_command_line (set iff MM_F_CIGAR (drv_b_api has_p)) |
| `mm_mapopt_t.max_qlen` | --max-qlen | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals --max-qlen) |
| `mm_mapopt_t.bw` | -r | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals -r) |
| `mm_mapopt_t.bw_long` | -r ,NUM2 | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals -r ,NUM2) |
| `mm_mapopt_t.max_gap` | -g | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals -g) |
| `mm_mapopt_t.max_chain_skip` | --max-chain-skip | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals --max-chain-skip) |
| `mm_mapopt_t.max_chain_iter` | --max-chain-iter | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals --max-chain-iter) |
| `mm_mapopt_t.min_cnt` | -n | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals -n) |
| `mm_mapopt_t.min_chain_score` | -m | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals -m) |
| `mm_mapopt_t.chain_gap_scale` | --chain-gap-scale | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals --chain-gap-scale) |
| `mm_mapopt_t.chain_skip_scale` | --chain-skip-scale | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals --chain-skip-scale) |
| `mm_mapopt_t.rmq_size_cap` | no option | not covered; no option and no documented effect; a cap low enough to change output was not constructed |
| `mm_mapopt_t.rmq_inner_dist` | --rmq-inner | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals --rmq-inner) |
| `mm_mapopt_t.rmq_rescue_size` | no option | test_b_api.py::test_rmq_rescue_size_and_ratio |
| `mm_mapopt_t.rmq_rescue_ratio` | no option | test_b_api.py::test_rmq_rescue_size_and_ratio |
| `mm_mapopt_t.mask_level` | -M | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals -M) |
| `mm_mapopt_t.mask_len` | --mask-len | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals --mask-len) |
| `mm_mapopt_t.pri_ratio` | -p | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals -p) |
| `mm_mapopt_t.best_n` | -N | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals -N) |
| `mm_mapopt_t.alt_drop` | --alt-drop | partial; set through --alt-drop only (test_b_hit.py::test_alt_contigs); the API driver has no ALT list |
| `mm_mapopt_t.a` | -A | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals -A) |
| `mm_mapopt_t.b` | -B | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals -B) |
| `mm_mapopt_t.q` | -O | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals -O) |
| `mm_mapopt_t.e` | -E | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals -E) |
| `mm_mapopt_t.q2` | -O ,INT2 | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals -O ,INT2) |
| `mm_mapopt_t.e2` | -E ,INT2 | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals -E ,INT2) |
| `mm_mapopt_t.transition` | -b | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals -b) |
| `mm_mapopt_t.sc_ambi` | --score-N | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals --score-N) |
| `mm_mapopt_t.zdrop` | -z | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals -z) |
| `mm_mapopt_t.zdrop_inv` | -z ,INT2 | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals -z ,INT2) |
| `mm_mapopt_t.end_bonus` | --end-bonus | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals --end-bonus) |
| `mm_mapopt_t.min_dp_max` | -s | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals -s) |
| `mm_mapopt_t.min_ksw_len` | --min-dp-len | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals --min-dp-len) |
| `mm_mapopt_t.anchor_ext_len` | no option | test_b_api.py::test_anchor_ext_len (boundary derived from the brute-force local score) |
| `mm_mapopt_t.anchor_ext_shift` | --end-seed-pen | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals --end-seed-pen (field set directly in test_b_api.py::test_anchor_ext_len)) |
| `mm_mapopt_t.max_clip_ratio` | --max-clip-ratio | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals --max-clip-ratio) |
| `mm_mapopt_t.rank_min_len` | no option | test_b_api.py::test_rank_min_len_and_rank_frac |
| `mm_mapopt_t.rank_frac` | no option | test_b_api.py::test_rank_min_len_and_rank_frac |
| `mm_mapopt_t.mid_occ_frac` | -f FLOAT | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals -f FLOAT) |
| `mm_mapopt_t.q_occ_frac` | --q-occ-frac | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals --q-occ-frac) |
| `mm_mapopt_t.min_mid_occ` | -U, --min-occ-floor | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals -U) |
| `mm_mapopt_t.max_mid_occ` | -U ,INT2 | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals -U ,INT2) |
| `mm_mapopt_t.mid_occ` | -f INT | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals -f INT) |
| `mm_mapopt_t.max_occ` | -f ,INT2 | partial; -f ,INT2 only acts in the short-read re-chaining (area C); API sets it without an observable long-read effect |
| `mm_mapopt_t.max_max_occ` | no option | test_b_api.py::test_max_max_occ; test_b_seed.py::test_collect_matches_against_brute_force[2-8-100] |
| `mm_mapopt_t.occ_dist` | -e | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals -e) |
| `mm_mapopt_t.max_sw_mat` | --cap-sw-mem | test_b_api.py::test_api_fields_and_flags_match_command_line (set directly, output equals --cap-sw-mem) |

### B: Non-static functions (46)

Defined without `static` (prototypes in mmpriv.h, bseq.h, sdust.h, ksw2.h, kalloc.h, kthread.h). A row is ticked when a test exercises the function with a checked result, through the CLI, the API or a C driver.

| File | Function | Notes | Tested |
|:--|:--|:--|:--|
| index.c | `mm_idx_cal_max_occ` |  | test_b_seed.py::test_cal_max_occ; test_b_seed.py::test_many_copy_occurrence_thresholds |
| index.c | `mm_idx_alt_read` |  | test_b_seed.py::test_alt_read_api; test_b_hit.py::test_alt_contigs |
| options.c | `mm_mapopt_update` |  | test_b_seed.py::test_mapopt_update_occurrence_threshold (drv_b_seed mapopt_update) |
| seed.c | `mm_seed_mz_flt` |  | test_b_seed.py::test_query_occurrence_filter, test_q_occ_frac_drops_query_repeats |
| seed.c | `mm_seed_collect_all` |  | test_b_seed.py::test_collect_matches_against_brute_force; tandem flag through test_b_presets.py::test_str_rich_genome_pinned |
| seed.c | `mm_seed_select` |  | test_b_seed.py::test_collect_matches_against_brute_force (dist > 0), test_e_samples_repetitive_minimizers; test_b_seed.py::test_collect_matches_many_copies, test_b_seed.py::test_sampled_minimizers_pinned |
| seed.c | `mm_collect_matches` |  | test_b_seed.py::test_collect_matches_against_brute_force; test_b_seed.py::test_collect_matches_many_copies |
| lchain.c | `mg_chain_backtrack` |  | test_b_chain.py (drv_b_chain: chain extraction, -n/-m filters, random sets) |
| lchain.c | `mg_lchain_dp` |  | test_b_chain.py::test_gap_cost_formula, test_random_anchor_sets_against_brute_force, CLI s1/cm tests; test_b_chain.py::test_chaining_heuristics_pinned |
| lchain.c | `mg_lchain_rmq` |  | test_b_chain.py (drv_b_chain method 1), test_band_edge, test_chain_score_and_count_substitutions[rmq]; test_b_chain.py::test_chaining_heuristics_pinned |
| hit.c | `mm_gen_regs` |  | via CLI: chain coordinates equal the brute-force chain ends (test_b_hit.py::test_mask_level_M), s1/cm (test_b_chain.py) |
| hit.c | `mm_mark_alt` |  | test_b_hit.py::test_alt_contigs |
| hit.c | `mm_split_reg` |  | test_b_align.py::test_zdrop_boundary (zd:i and conserved cm/s1), test_inversion_rescue |
| hit.c | `mm_set_parent` |  | test_b_hit.py::test_mask_level_M, test_hard_mask_level, test_mask_len |
| hit.c | `mm_hit_sort` |  | test_b_hit.py::test_alt_contigs (ranking with the ALT penalty); primaries first in every test |
| hit.c | `mm_set_sam_pri` |  | test_b_align.py::test_chimeric_reads_supplementary_and_SA (one primary, one supplementary); drv_b_api sam_pri checks |
| hit.c | `mm_sync_regs` |  | via -p/-N tests that remove hits (test_b_hit.py), tp:A and SAM flags consistent |
| hit.c | `mm_select_sub` |  | test_b_hit.py::test_pri_ratio_p, test_repeat_reads_secondaries_and_N, test_opposite_strand_secondary_retained |
| hit.c | `mm_filter_strand_retained` |  | partial; test_b_hit.py::test_secondaries_below_pri_ratio_pinned (retained hit kept; the divergence filter itself not triggered) |
| hit.c | `mm_filter_regs` |  | test_b_align.py::test_peak_score_ms_and_s_boundary, test_max_clip_ratio_boundary; test_b_chain.py (-m/-n with -c via API) |
| hit.c | `mm_squeeze_a` |  | via every -c test (alignments from the right anchors) |
| hit.c | `mm_set_mapq2` |  | test_b_hit.py::test_unique_reads_mapq, test_repeat_reads_secondaries_and_N (0 on repeats), test_alt_contigs; test_b_align.py::test_inversion_rescue, test_inversion_mapq_from_flanks; test_b_hit.py::test_mapq_near_identical_copies_pinned |
| align.c | `mm_enlarge_cigar` |  | drv_b_api capacity check on every CIGAR |
| align.c | `mm_event_identity` |  | de:f recomputed by paf_invariant_errors in test_b_align.py tests |
| align.c | `mm_update_dp_max` |  | test_b_api.py::test_rank_min_len_and_rank_frac |
| align.c | `mm_align_skeleton` |  | test_b_align.py (all -c tests) |
| esterr.c | `mm_est_err` |  | test_b_chain.py::test_chain_score_and_count_substitutions (dv:f recomputed) |
| map.c | `mm_map_frag_core` | the mapping core every path goes through | every CLI and API test |
| ksw2_extz2_sse.c | `ksw_extz2_sse41` |  | test_b_ksw.py::test_simd_variants_agree (weak symbol; x86 dispatch build) |
| ksw2_extz2_sse.c | `ksw_extz2_sse2` |  | test_b_ksw.py::test_simd_variants_agree; sse2only configuration |
| ksw2_extz2_sse.c | `ksw_extz2_sse` |  | test_b_ksw.py::test_global_alignment_matches_reference[extz2], test_extension_matches_reference[extz2], test_band_limits_global_alignment, test_gap_left_and_right_alignment |
| ksw2_extd2_sse.c | `ksw_extd2_sse41` |  | test_b_ksw.py::test_simd_variants_agree |
| ksw2_extd2_sse.c | `ksw_extd2_sse2` |  | test_b_ksw.py::test_simd_variants_agree; sse2only configuration |
| ksw2_extd2_sse.c | `ksw_extd2_sse` |  | test_b_ksw.py::test_global_alignment_matches_reference[extd2], test_extension_matches_reference[extd2]; a defect of base minimap2 for narrow bands |
| ksw2_exts2_sse.c | `ksw_exts2_sse41` |  | test_b_ksw.py::test_simd_variants_agree_spliced |
| ksw2_exts2_sse.c | `ksw_exts2_sse2` |  | test_b_ksw.py::test_simd_variants_agree_spliced |
| ksw2_exts2_sse.c | `ksw_exts2_sse` | splice DP; reached through C's splice tests as well | test_b_ksw.py::test_spliced_kernel_places_intron, test_spliced_kernel_long_gap_threshold; test_b_align.py::test_end_seed_pen_boundary |
| ksw2_ll_sse.c | `ksw_ll_qinit` |  | test_b_ksw.py::test_local_alignment_ll_i16 |
| ksw2_ll_sse.c | `ksw_ll_i16` |  | test_b_ksw.py::test_local_alignment_ll_i16; inversion and end-seed tests |
| ksw2_dispatch.c | `__cpuidex` | Intel-derived wrapper with unknown terms: off limits, test only through dispatch | (through dispatch only); off limits; exercised by every kernel call in the dispatch build |
| ksw2_dispatch.c | `ksw_extz2_sse` |  | test_b_ksw.py (entry point in the dispatch build) |
| ksw2_dispatch.c | `ksw_extd2_sse` |  | test_b_ksw.py |
| ksw2_dispatch.c | `ksw_exts2_sse` | splice DP; reached through C's splice tests as well | test_b_ksw.py::test_spliced_kernel_places_intron, test_simd_variants_agree_spliced |
| misc.c | `radix_sort_128x` | macro-generated (KRADIX_SORT_INIT) | test_b_chain.py::test_sort_helpers |
| misc.c | `radix_sort_64` | macro-generated (KRADIX_SORT_INIT) | test_b_chain.py::test_sort_helpers |
| misc.c | `ks_ksmall_uint32_t` | macro-generated (KSORT_INIT_GENERIC); used by mm_idx_cal_max_occ | test_b_chain.py::test_sort_helpers; test_b_seed.py::test_cal_max_occ |

### B: Behaviours not tied to one name (4)

| Behaviour | Notes | Tested |
|:--|:--|:--|
| unique reads map to their true locus; repeats get MAPQ 0 and secondaries |  | test_b_hit.py::test_unique_reads_mapq, test_repeat_reads_secondaries_and_N; test_b_presets.py::test_preset_reads_at_truth |
| AS equals the CIGAR's score under -A/-B/-O/-E for long-read presets | observed on simulated reads; sr differs | test_b_align.py::test_AS_is_cigar_score_under_preset_scoring, test_AS_is_cigar_score_under_custom_scoring |
| inversion rescue: tp:A:I/i alignments |  | test_b_align.py::test_inversion_rescue, test_inversion_mapq_from_flanks |
| chimeric reads: supplementary alignments, SA:Z |  | test_b_align.py::test_chimeric_reads_supplementary_and_SA |

## Area C

### C: Command-line options (20)

"Documented": whether the usage text (`minimap2 -h`) and minimap2.1 describe the option.

| Option | Argument | Documented | Notes | Tested |
|:--|:--|:--|:--|:--|
| `-G` / `--max-intron-len` | yes | usage, man | only takes effect in splice mode (mm_mapopt_max_intron_len) | test_c_splice.py::test_max_intron_len, test_c_api.py::test_max_intron_len_api; order dependence and N > -G is a defect of base minimap2, not tested |
| `-u` | yes | usage, man |  | test_c_splice.py::test_u_option_chooses_strand (b, f, r), test_c_splice.py::test_u_n_ignores_splice_signals (n), test_c_splice.py::test_u_invalid_is_an_error; ts:A under -un is a defect of base minimap2, not tested |
| `-F` | yes | usage, man |  | test_c_pairs.py::test_max_fragment_length (default 800, 500, 1000, pairs at the limit and one base past it) |
| `-C` / `--cost-non-gt-ag` | yes | man |  | test_c_splice.py::test_noncanonical_cost_C_and_J (-J0 cost per site, no effect under -J1) |
| `-J` | yes | usage, man |  | test_c_splice.py::test_noncanonical_cost_C_and_J, test_c_splice.py::test_splice_flank_under_J0, test_c_splice.py::test_junc_bed_decides_a_tie (-J0) |
| `-j` | yes | usage, man |  | test_c_srna.py::test_jump_annotation_extends_read_ends (BED12 and BED6, both read strands, prebuilt index), test_c_srna.py::test_jump_min_match, test_c_srna.py::test_jump_prefers_annotation_over_pass1, test_c_srna.py::test_jump_needs_splice_mode, test_c_srna.py pin; test_c_srna.py::test_jump_with_a_mismatch_near_the_junction |
| `--splice` | no | man |  | test_c_splice.py::test_splice_option_enables_spliced_alignment, test_c_splice.py::test_splice_preset_equivalences, test_c_splice.py::test_max_intron_len |
| `--sr` | optional | man |  | test_c_pairs.py::test_option_values (--sr, =dna, =no, invalid value), test_c_pairs.py::test_splice_sr_single_end_equals_its_definition (=rna) |
| `--frag` | yes | man |  | test_c_pairs.py::test_frag_mode_groups_reads_by_name, test_c_pairs.py::test_three_segment_fragments, test_c_pairs.py::test_option_values |
| `--no-pairing` | no | - | deprecated; same as --pairing=no | test_c_pairs.py::test_pairing_no_aligns_ends_independently (identical to --pairing=no) |
| `--splice-flank` | yes | man |  | partial: test_c_splice.py::test_splice_flank_under_J0 (yes under -J0 and -J1, no under -J0); no effect under -J1 is a defect of base minimap2, not tested |
| `--junc-bed` | yes | man |  | test_c_splice.py::test_junc_bed_decides_a_tie (both gene strands, BED12 and BED6, lines without strand, -J0), test_c_splice.py pins |
| `--junc-bonus` | yes | man |  | test_c_splice.py::test_junc_bed_decides_a_tie (9, 5, 0), test_c_splice.py::test_spsc_scores_decide_a_tie (no effect with --spsc) |
| `--spsc` | yes | man |  | test_c_splice.py::test_spsc_scores_decide_a_tie (both strands), test_c_splice.py::test_spsc_scores_are_capped, test_c_api.py::test_spsc_reader |
| `--junc-pen` | yes | - |  | test_c_splice.py::test_spsc_scores_decide_a_tie (same output as --spsc0), test_c_splice.py::test_junc_bed_decides_a_tie (no effect with --junc-bed) |
| `--pairing` | yes | man |  | test_c_pairs.py::test_pairing_places_repetitive_mates (strong, no), test_c_pairs.py::test_option_values (strong, invalid), test_c_srna.py::test_splice_sr_pairs_against_truth (weak via splice:sr), test_c_pairs.py::test_splice_sr_single_end_equals_its_definition; weak with -x sr is a defect of base minimap2, not tested; test_c_pairs.py::test_pairs_of_repetitive_mates |
| `--jump-min-match` | yes | man |  | test_c_srna.py::test_jump_min_match (5, 8, 9, 13) |
| `--pass1` | yes | man |  | test_c_srna.py::test_pass1_score_threshold (4, 5, 6), test_c_srna.py::test_two_pass_workflow, test_c_srna.py::test_jump_prefers_annotation_over_pass1 |
| `--spsc-scale` | yes | man |  | partial: test_c_splice.py::test_spsc_scores_decide_a_tie (0.5, 0.7, 1), test_c_api.py::test_spsc_reader (0.6, 0.7, none); values above 1 ignored, is a defect of base minimap2, not tested |
| `--spsc0` | yes | man | same field as --junc-pen | test_c_splice.py::test_spsc_scores_decide_a_tie (0, 5, 9; positions not in FILE) |

### C: Presets: mapping parameters (6)

| Preset | Documented | Notes | Tested |
|:--|:--|:--|:--|
| `short` | - | alias of sr; undocumented | test_c_pairs.py::test_option_values (same output as sr), test_c_api.py::test_preset_option_fields (same options) |
| `sr` | usage, man |  | test_c_pairs.py (all tests), test_c_api.py::test_preset_option_fields, test_c_api.py::test_map_frag_api, test_c_pairs.py pin; expansion mismatch (-O12,32, FR orientation) is a defect of base minimap2, not tested |
| `splice` | usage, man |  | test_c_splice.py::test_splice_ont_like_cdna, test_c_splice.py::test_splice_preset_equivalences, test_c_api.py::test_preset_option_fields, test_c_splice.py pins |
| `splice:hq` | usage, man |  | test_c_splice.py::test_u_option_chooses_strand, test_c_splice.py::test_splice_preset_equivalences, test_c_api.py::test_map_splice_api, test_c_splice.py PAF pin |
| `splice:sr` | usage, man |  | test_c_srna.py (all tests), test_c_pairs.py::test_splice_sr_single_end_equals_its_definition, test_c_api.py::test_preset_option_fields, test_c_srna.py pin |
| `cdna` | - | alias of splice; undocumented | test_c_splice.py::test_splice_preset_equivalences, test_c_api.py::test_preset_option_fields |

### C: minimap.h functions (9)

| Function | Tested |
|:--|:--|
| `mm_mapopt_max_intron_len` | test_c_api.py::test_max_intron_len_api, test_c_splice.py::test_max_intron_len |
| `mm_map_frag` | test_c_api.py::test_map_frag_api (pairs, one segment = mm_map, MM_MAX_SEG + 1) |
| `mm_map_file_frag` | test_c_api.py::test_map_file_frag_api (two files vs CLI sr, one file vs CLI splice, missing file) |
| `mm_idx_bed_read` | test_c_api.py::test_bed_junc (BED12, BED6, read_junc 0 and 1), --junc-bed tests (test_c_splice.py) |
| `mm_idx_bed_junc` | test_c_api.py::test_bed_junc |
| `mm_max_spsc_bonus` | test_c_api.py::test_preset_option_fields (pins the current values 20/30/18), test_c_splice.py::test_spsc_scores_are_capped |
| `mm_idx_spsc_read` | test_c_api.py::test_spsc_reader[read] |
| `mm_idx_spsc_read2` | test_c_api.py::test_spsc_reader, --spsc tests (test_c_splice.py) |
| `mm_idx_spsc_get` | test_c_api.py::test_spsc_reader |

### C: minimap.h macros (11)

| Macro | Meaning | Tested |
|:--|:--|:--|
| `MM_F_SPLICE` | set by --splice, splice presets, -u f/b via mm_mapopt_update | test_c_api.py::test_preset_option_fields, test_c_splice.py::test_splice_option_enables_spliced_alignment |
| `MM_F_SPLICE_FOR` | set by -u f/b, splice presets | test_c_api.py::test_preset_option_fields, test_c_splice.py::test_u_option_chooses_strand |
| `MM_F_SPLICE_REV` | set by -u r/b, splice presets | test_c_api.py::test_preset_option_fields, test_c_splice.py::test_u_option_chooses_strand |
| `MM_F_SR` | set by --sr, sr preset | test_c_api.py::test_preset_option_fields, test_c_pairs.py::test_option_values |
| `MM_F_FRAG_MODE` | set by --frag, sr and splice:sr presets | test_c_api.py::test_preset_option_fields, test_c_pairs.py::test_frag_mode_groups_reads_by_name |
| `MM_F_INDEPEND_SEG` | set by --pairing=no, --no-pairing | test_c_pairs.py::test_pairing_no_aligns_ends_independently, test_c_pairs.py::test_pairing_places_repetitive_mates |
| `MM_F_SPLICE_FLANK` | set by --splice-flank, splice presets | partial: test_c_api.py::test_preset_option_fields, test_c_splice.py::test_splice_flank_under_J0 (see --splice-flank) |
| `MM_F_SPLICE_OLD` | set by -J 0 | test_c_splice.py::test_noncanonical_cost_C_and_J, test_c_splice.py::test_splice_flank_under_J0 |
| `MM_F_WEAK_PAIRING` | set by --pairing=weak, splice:sr preset | test_c_api.py::test_preset_option_fields, test_c_srna.py::test_splice_sr_pairs_against_truth |
| `MM_F_SR_RNA` | set by --sr=rna, splice:sr preset | test_c_api.py::test_preset_option_fields, test_c_pairs.py::test_splice_sr_single_end_equals_its_definition |
| `MM_MAX_SEG` | max segments per fragment | test_c_api.py::test_preset_option_fields (255), test_c_api.py::test_map_frag_api (256 segments give no hits) |

### C: minimap.h struct fields (18)

Fields of `mm_reg1_t` and `mm_extra_t` are results a C driver can check against the CLI or the simulator's truth; fields of the option structs name the option that sets them.

| Field | Set by / notes | Tested |
|:--|:--|:--|
| `mm_idx_t.I` | hidden intervals (--junc-bed); C | test_c_api.py::test_bed_junc |
| `mm_idx_t.spsc` | hidden (--spsc); C | test_c_api.py::test_spsc_reader |
| `mm_idx_t.J` | hidden (-j, --pass1); C | test_c_api.py::test_jump_reader |
| `mm_extra_t.trans_strand` |  | test_c_api.py::test_map_splice_api, ts:A in splice and srna tests |
| `mm_reg1_t.proper_frag` |  | test_c_api.py::test_map_frag_api, 0x2 in pairs tests; test_c_pairs.py::test_pairs_of_repetitive_mates |
| `mm_reg1_t.pe_thru` |  | test_c_api.py::test_map_frag_api; test_c_api.py::test_map_frag_api (read-through pairs with adapters) |
| `mm_reg1_t.seg_split` |  | test_c_api.py::test_map_frag_api |
| `mm_reg1_t.seg_id` |  | test_c_api.py::test_map_frag_api |
| `mm_reg1_t.strand_retained` |  | not covered: set by mm_select_sub() between chaining and alignment and filtered by mm_filter_strand_retained(); a retained chain is dropped again after alignment, so the field has no documented result in mm_map() output to assert on |
| `mm_reg1_t.is_spliced` |  | test_c_api.py::test_map_splice_api |
| `mm_mapopt_t.max_gap_ref` | -G (splice), presets | test_c_api.py::test_max_intron_len_api, test_c_api.py::test_preset_option_fields |
| `mm_mapopt_t.max_frag_len` | -F | test_c_api.py::test_preset_option_fields, test_c_pairs.py::test_max_fragment_length |
| `mm_mapopt_t.noncan` | -C | test_c_api.py::test_preset_option_fields, test_c_splice.py::test_noncanonical_cost_C_and_J |
| `mm_mapopt_t.junc_bonus` | --junc-bonus | test_c_api.py::test_preset_option_fields, test_c_splice.py::test_junc_bed_decides_a_tie |
| `mm_mapopt_t.junc_pen` | --junc-pen / --spsc0 | test_c_api.py::test_preset_option_fields, test_c_splice.py::test_spsc_scores_decide_a_tie |
| `mm_mapopt_t.pe_ori` | presets (sr, splice:sr) | test_c_api.py::test_preset_option_fields (FR for sr and splice:sr), every paired test (FR pairs flagged 0x2) |
| `mm_mapopt_t.pe_bonus` | no option | partial: no option and no documented value to check; exercised by every paired run, and flipping its sign in mm_pair() changes proper-pair flags caught by test_c_pairs.py::test_max_fragment_length and test_c_api.py::test_map_frag_api |
| `mm_mapopt_t.jump_min_match` | --jump-min-match | test_c_api.py::test_preset_option_fields (3), test_c_srna.py::test_jump_min_match |

### C: Non-static functions (17)

Defined without `static` (prototypes in mmpriv.h, bseq.h, sdust.h, ksw2.h, kalloc.h, kthread.h). A row is ticked when a test exercises the function with a checked result, through the CLI, the API or a C driver.

| File | Function | Notes | Tested |
|:--|:--|:--|:--|
| index.c | `mm_idx_bed_read` |  | test_c_api.py::test_bed_junc, --junc-bed tests (test_c_splice.py) |
| index.c | `mm_idx_bed_junc` |  | test_c_api.py::test_bed_junc, test_c_splice.py::test_junc_bed_decides_a_tie |
| index.c | `mm_idx_jjump_read` |  | test_c_api.py::test_jump_reader, -j and --pass1 tests |
| index.c | `mm_idx_jump_get` |  | test_c_api.py::test_jump_reader, -j tests (test_c_srna.py) |
| index.c | `mm_idx_spsc_read2` |  | test_c_api.py::test_spsc_reader, --spsc tests (test_c_splice.py) |
| index.c | `mm_idx_spsc_read` |  | test_c_api.py::test_spsc_reader[read] |
| index.c | `mm_idx_spsc_get` |  | test_c_api.py::test_spsc_reader |
| options.c | `mm_mapopt_max_intron_len` |  | test_c_api.py::test_max_intron_len_api |
| options.c | `mm_max_spsc_bonus` |  | test_c_api.py::test_preset_option_fields, test_c_splice.py::test_spsc_scores_are_capped |
| hit.c | `mm_seg_gen` |  | test_c_api.py::test_map_frag_api (seg_split, seg_id), all paired tests (test_c_pairs.py) |
| hit.c | `mm_seg_free` |  | partial: runs after every multi-segment mapping (test_c_api.py::test_map_frag_api, paired tests (test_c_pairs.py)); it only frees memory, so there is no result to check beyond completing |
| jump.c | `mm_jump_split` |  | srna -j, --pass1 and --jump-min-match tests; test_c_srna.py::test_jump_with_a_mismatch_near_the_junction |
| map.c | `mm_map_frag` |  | test_c_api.py::test_map_frag_api |
| map.c | `mm_map_file_frag` |  | test_c_api.py::test_map_file_frag_api |
| pe.c | `mm_select_sub_multi` |  | partial: runs on every paired mapping (test_c_pairs.py); it chooses secondary chains, which `--secondary=no` hides, so only the sr pin (test_c_pairs.py) caught a changed ratio |
| pe.c | `mm_set_pe_thru` |  | test_c_api.py::test_map_frag_api; test_c_api.py::test_map_frag_api (read-through pairs with adapters) |
| pe.c | `mm_pair` |  | test_c_pairs.py::test_max_fragment_length, test_c_pairs.py::test_pairing_places_repetitive_mates, test_c_api.py::test_map_frag_api; test_c_pairs.py::test_pairs_of_repetitive_mates |

### C: Behaviours not tied to one name (3)

| Behaviour | Notes | Tested |
|:--|:--|:--|
| paired-end from two files, interleaved one file, /1 /2 suffixes |  | test_c_pairs.py::test_paired_input_forms_agree, test_c_pairs.py::test_frag_mode_minibatches_keep_mates_together, test_c_pairs.py::test_input_checks, test_c_pairs.py::test_pairs_with_a_split_index |
| frag mode: reads with the same name in one file |  | test_c_pairs.py::test_frag_mode_groups_reads_by_name, test_c_pairs.py::test_three_segment_fragments (middle-segment flags are a defect of base minimap2) |
| spliced reads: introns, ts:A, junction BED bonus |  | test_c_splice.py::test_splice_ont_like_cdna, test_c_splice.py::test_junc_bed_decides_a_tie, test_c_srna.py::test_splice_sr_pairs_against_truth, test_pilot (splice:hq) |

## Area D

### D: Command-line options (33)

"Documented": whether the usage text (`minimap2 -h`) and minimap2.1 describe the option.

| Option | Argument | Documented | Notes | Tested |
|:--|:--|:--|:--|:--|
| `-2` | no | man |  | test_d_invariance.py::test_long_reads_sam, test_d_invariance.py::test_long_reads_paf (output unchanged), test_d_api.py::test_mm_map_file_matches_cli (`-2` = MM_F_2_IO_THREADS), test_d_api.py::test_set_opt_and_mapopt_init (sr, splice:sr); test_d_invariance.py::test_larger_input |
| `-a` / `--sam` | no | usage, man |  | every sam test; test_d_cli.py::test_long_option_names (`--sam`); test_d_api.py::test_mm_map_file_matches_cli (MM_F_OUT_SAM). |
| `-S` | no | - | deprecated: same as --cs=long with a warning | test_d_sam.py::test_cs_variants_and_deprecated_s (warning, long cs). |
| `-K` / `--mb-size` | yes | usage, man |  | inv (output unchanged for -K 1, 1k, 2k, 3k, 5k, 10k), test_d_cli.py::test_long_option_names, test_d_api.py::test_mm_map_file_matches_cli (mini_batch_size 3000). |
| `-t` | yes | usage, man |  | test_d_invariance.py::test_long_reads_sam (1, 2, 3, 4, 8), test_paired_reads, test_spliced_reads; test_d_cli.py::test_usage_shows_current_values_and_preset_first (default 3); test_d_invariance.py::test_larger_input (2000 pairs, 180 long reads; -t 1, 4, 8) |
| `-V` / `--version` | no | usage, man | prints MM_VERSION | test_d_cli.py::test_version, test_d_api.py::test_macros_and_version. |
| `-v` | yes | usage (commented out) | verbosity; commented out of the usage text | test_d_cli.py::test_verbosity (0-3). |
| `-x` | yes | usage, man | preset machinery: applied before other options whatever its position; unknown preset is an error (presets in their own table) | test_d_cli.py::test_usage_shows_current_values_and_preset_first, test_d_cli.py::test_preset_applied_first, test_d_cli.py::test_unknown_preset. |
| `-Q` | no | man |  | test_d_sam.py::test_qualities, test_d_api.py::test_mm_map_file_matches_cli. |
| `-R` | yes | usage, man |  | test_d_sam.py::test_read_group (header, RG:Z, ID length 255/256, errors), test_d_api.py::test_sam_header_writer, test_d_api.py::test_mm_map_file_matches_cli. |
| `-h` / `--help` | no | - | --help: usage to stdout, exit 0; no arguments: usage to stderr, exit 1 | test_d_cli.py::test_help_and_usage, test_d_cli.py::test_usage_shows_current_values_and_preset_first. |
| `-L` | no | usage, man |  | test_d_api.py::test_long_cigar (boundaries 65535/65536 with clips, supplementary, secondary, -Y, --secondary-seq), test_d_sam.py::test_long_cigar_option_keeps_short_cigars. |
| `-y` | no | usage, man |  | test_d_sam.py::test_copy_comment, test_d_paf.py::test_paf_copy_comment, test_d_api.py::test_writers_match_cli (its --paf-no-hit line: a defect of base minimap2). |
| `-Y` | no | usage, man |  | test_d_sam.py::test_soft_clipped_supplementary, test_d_api.py::test_long_cigar. |
| `-o` | yes | usage, man |  | test_d_cli.py::test_output_file. |
| `--seed` | yes | man | seeds the hash that breaks ties between equal hits | test_d_invariance.py::test_seed, test_d_api.py::test_mm_map_file_matches_cli (mm_mapopt_t.seed). |
| `--no-kalloc` | no | man |  | test_d_invariance.py::test_allocator_options, test_d_api.py::test_thread_buffer, test_d_api.py::test_mm_map_matches_cli (no_kalloc=1). |
| `--print-qname` | no | man | debug output to stderr | test_d_cli.py::test_debug_options (output unchanged, every read named on stderr). |
| `--print-seeds` | no | man | debug output to stderr; forces -t 1 | partial: test_d_cli.py::test_debug_options (output unchanged, stderr written); its format and "forces -t 1" are not observable without pinning debug output (tests/README.md rule). |
| `--print-aln-seq` | no | - | debug output to stderr; forces -t 1 | partial: test_d_cli.py::test_debug_options (output unchanged); what it prints is debug output and is not checked |
| `--secondary` | yes | man |  | test_d_sam.py::test_secondary_no_and_n0 (yes, no, invalid value). |
| `--cs` | optional | usage, man |  | test_d_sam.py::test_cs_md, test_d_sam.py::test_cs_variants_and_deprecated_s, test_d_sam.py::test_tags_at_read_ends, test_d_paf.py::test_paf_cs_md, test_d_api.py::test_mm_map_matches_cli. |
| `--MD` | no | usage, man |  | test_d_sam.py::test_cs_md, test_d_sam.py::test_tags_at_read_ends, test_d_paf.py::test_paf_cs_md (with -c; a defect of base minimap2), test_d_api.py::test_mm_map_matches_cli. |
| `--lj-min-ratio` | yes | - | deprecated; prints a warning and is ignored | test_d_cli.py::test_deprecated_lj_min_ratio. |
| `--eqx` | no | usage, man |  | test_d_sam.py::test_eqx, test_d_sam.py::test_tags_at_read_ends. |
| `--paf-no-hit` | no | man |  | test_d_paf.py::test_paf_no_hit. |
| `--sam-hit-only` | no | man |  | test_d_sam.py::test_unmapped_records. |
| `--cap-kalloc` | yes | man |  | test_d_invariance.py::test_allocator_options (1k, 0: output unchanged); its memory effect is a kalloc internal (not pinned). |
| `--print-chains` | no | - | debug output to stderr | partial: test_d_cli.py::test_debug_options (output unchanged); what it prints is debug output and is not checked |
| `--secondary-seq` | no | man |  | test_d_sam.py::test_secondary_seq, test_d_api.py::test_long_cigar. |
| `--ds` | no | usage |  | test_d_sam.py::test_ds_tag, test_d_sam.py::test_tags_at_read_ends, test_d_paf.py::test_paf_ds, test_d_api.py::test_mm_map_matches_cli. |
| `--write-junc` | no | man |  | test_d_paf.py::test_write_junc (all five scored splice signals, MAPQ filter, strand), test_d_api.py::test_mm_map_file_matches_cli. |
| `--dbg-seed-occ` | no | - | debug output to stderr | partial: test_d_cli.py::test_debug_options (output unchanged); what it prints is debug output and is not checked |

### D: Presets: mapping parameters (1)

| Preset | Documented | Notes | Tested |
|:--|:--|:--|:--|
| unknown preset, e.g. `-x foo` or `-x asm7` | - | `[ERROR] unknown preset`, exit 1 | test_d_cli.py::test_unknown_preset, test_d_api.py::test_set_opt_and_mapopt_init |

### D: minimap.h functions (11)

| Function | Tested |
|:--|:--|
| `mm_set_opt` | test_d_api.py::test_set_opt_and_mapopt_init. |
| `mm_check_opt` | test_d_api.py::test_check_opt, test_d_cli.py::test_option_check_errors, test_d_cli.py::test_option_check_boundaries. |
| `mm_tbuf_init` | test_d_api.py::test_thread_buffer, test_d_api.py::test_mm_map_matches_cli. |
| `mm_tbuf_destroy` | test_d_api.py::test_thread_buffer, test_d_api.py::test_mm_map_matches_cli. |
| `mm_tbuf_get_km` | test_d_api.py::test_thread_buffer, test_d_api.py::test_mm_map_matches_cli. |
| `mm_map` | test_d_api.py::test_mm_map_matches_cli, test_d_api.py::test_writers_match_cli, test_d_api.py::test_minimap2_lite. |
| `mm_map_file` | test_d_api.py::test_mm_map_file_matches_cli (11 option sets, return values). |
| `mm_gen_cs` | test_d_api.py::test_mm_map_matches_cli (equal to the CLI tags, returned lengths, buffer reuse, with and without a kalloc pool). |
| `mm_gen_ds` | test_d_api.py::test_mm_map_matches_cli (equal to the CLI tags, returned lengths, buffer reuse, with and without a kalloc pool). |
| `mm_gen_MD` | test_d_api.py::test_mm_map_matches_cli (equal to the CLI tags, returned lengths, buffer reuse, with and without a kalloc pool). |
| `mm_mapopt_init` | test_d_api.py::test_set_opt_and_mapopt_init. |

### D: minimap.h macros (28)

| Macro | Meaning | Tested |
|:--|:--|:--|
| `MM_VERSION` | `--version`, SAM @PG VN | test_d_cli.py::test_version, test_d_sam.py::test_header_lines, test_d_api.py::test_macros_and_version. |
| `MM_F_OUT_SAM` | set by -a | test_d_api.py::test_mm_map_file_matches_cli (each flag gives its option's output); MM_F_LONG_CIGAR, MM_F_SOFTCLIP, MM_F_SECONDARY_SEQ also test_d_api.py::test_long_cigar; MM_F_NO_PRINT_2ND and MM_F_2_IO_THREADS also the presets in test_d_api.py::test_set_opt_and_mapopt_init. |
| `MM_F_NO_QUAL` | set by -Q | test_d_api.py::test_mm_map_file_matches_cli (each flag gives its option's output); MM_F_LONG_CIGAR, MM_F_SOFTCLIP, MM_F_SECONDARY_SEQ also test_d_api.py::test_long_cigar; MM_F_NO_PRINT_2ND and MM_F_2_IO_THREADS also the presets in test_d_api.py::test_set_opt_and_mapopt_init. |
| `MM_F_OUT_CG` | set by -c | test_d_api.py::test_mm_map_file_matches_cli (each flag gives its option's output); MM_F_LONG_CIGAR, MM_F_SOFTCLIP, MM_F_SECONDARY_SEQ also test_d_api.py::test_long_cigar; MM_F_NO_PRINT_2ND and MM_F_2_IO_THREADS also the presets in test_d_api.py::test_set_opt_and_mapopt_init. |
| `MM_F_OUT_CS` | set by --cs, -S | test_d_api.py::test_mm_map_file_matches_cli (each flag gives its option's output); MM_F_LONG_CIGAR, MM_F_SOFTCLIP, MM_F_SECONDARY_SEQ also test_d_api.py::test_long_cigar; MM_F_NO_PRINT_2ND and MM_F_2_IO_THREADS also the presets in test_d_api.py::test_set_opt_and_mapopt_init. |
| `MM_F_OUT_CS_LONG` | set by --cs=long, -S | test_d_api.py::test_mm_map_file_matches_cli (each flag gives its option's output); MM_F_LONG_CIGAR, MM_F_SOFTCLIP, MM_F_SECONDARY_SEQ also test_d_api.py::test_long_cigar; MM_F_NO_PRINT_2ND and MM_F_2_IO_THREADS also the presets in test_d_api.py::test_set_opt_and_mapopt_init. |
| `MM_F_NO_PRINT_2ND` | set by --secondary=no, -N 0, sr and splice:sr presets | test_d_api.py::test_mm_map_file_matches_cli (each flag gives its option's output); MM_F_LONG_CIGAR, MM_F_SOFTCLIP, MM_F_SECONDARY_SEQ also test_d_api.py::test_long_cigar; MM_F_NO_PRINT_2ND and MM_F_2_IO_THREADS also the presets in test_d_api.py::test_set_opt_and_mapopt_init. |
| `MM_F_2_IO_THREADS` | set by -2, sr and splice:sr presets | test_d_api.py::test_mm_map_file_matches_cli (each flag gives its option's output); MM_F_LONG_CIGAR, MM_F_SOFTCLIP, MM_F_SECONDARY_SEQ also test_d_api.py::test_long_cigar; MM_F_NO_PRINT_2ND and MM_F_2_IO_THREADS also the presets in test_d_api.py::test_set_opt_and_mapopt_init. |
| `MM_F_LONG_CIGAR` | set by -L | test_d_api.py::test_mm_map_file_matches_cli (each flag gives its option's output); MM_F_LONG_CIGAR, MM_F_SOFTCLIP, MM_F_SECONDARY_SEQ also test_d_api.py::test_long_cigar; MM_F_NO_PRINT_2ND and MM_F_2_IO_THREADS also the presets in test_d_api.py::test_set_opt_and_mapopt_init. |
| `MM_F_SOFTCLIP` | set by -Y | test_d_api.py::test_mm_map_file_matches_cli (each flag gives its option's output); MM_F_LONG_CIGAR, MM_F_SOFTCLIP, MM_F_SECONDARY_SEQ also test_d_api.py::test_long_cigar; MM_F_NO_PRINT_2ND and MM_F_2_IO_THREADS also the presets in test_d_api.py::test_set_opt_and_mapopt_init. |
| `MM_F_OUT_MD` | set by --MD | test_d_api.py::test_mm_map_file_matches_cli (each flag gives its option's output); MM_F_LONG_CIGAR, MM_F_SOFTCLIP, MM_F_SECONDARY_SEQ also test_d_api.py::test_long_cigar; MM_F_NO_PRINT_2ND and MM_F_2_IO_THREADS also the presets in test_d_api.py::test_set_opt_and_mapopt_init. |
| `MM_F_COPY_COMMENT` | set by -y | test_d_api.py::test_mm_map_file_matches_cli (each flag gives its option's output); MM_F_LONG_CIGAR, MM_F_SOFTCLIP, MM_F_SECONDARY_SEQ also test_d_api.py::test_long_cigar; MM_F_NO_PRINT_2ND and MM_F_2_IO_THREADS also the presets in test_d_api.py::test_set_opt_and_mapopt_init. |
| `MM_F_EQX` | set by --eqx | test_d_api.py::test_mm_map_file_matches_cli (each flag gives its option's output); MM_F_LONG_CIGAR, MM_F_SOFTCLIP, MM_F_SECONDARY_SEQ also test_d_api.py::test_long_cigar; MM_F_NO_PRINT_2ND and MM_F_2_IO_THREADS also the presets in test_d_api.py::test_set_opt_and_mapopt_init. |
| `MM_F_PAF_NO_HIT` | set by --paf-no-hit | test_d_api.py::test_mm_map_file_matches_cli (each flag gives its option's output); MM_F_LONG_CIGAR, MM_F_SOFTCLIP, MM_F_SECONDARY_SEQ also test_d_api.py::test_long_cigar; MM_F_NO_PRINT_2ND and MM_F_2_IO_THREADS also the presets in test_d_api.py::test_set_opt_and_mapopt_init. |
| `MM_F_SAM_HIT_ONLY` | set by --sam-hit-only | test_d_api.py::test_mm_map_file_matches_cli (each flag gives its option's output); MM_F_LONG_CIGAR, MM_F_SOFTCLIP, MM_F_SECONDARY_SEQ also test_d_api.py::test_long_cigar; MM_F_NO_PRINT_2ND and MM_F_2_IO_THREADS also the presets in test_d_api.py::test_set_opt_and_mapopt_init. |
| `MM_F_SECONDARY_SEQ` | set by --secondary-seq | test_d_api.py::test_mm_map_file_matches_cli (each flag gives its option's output); MM_F_LONG_CIGAR, MM_F_SOFTCLIP, MM_F_SECONDARY_SEQ also test_d_api.py::test_long_cigar; MM_F_NO_PRINT_2ND and MM_F_2_IO_THREADS also the presets in test_d_api.py::test_set_opt_and_mapopt_init. |
| `MM_F_OUT_DS` | set by --ds | test_d_api.py::test_mm_map_file_matches_cli (each flag gives its option's output); MM_F_LONG_CIGAR, MM_F_SOFTCLIP, MM_F_SECONDARY_SEQ also test_d_api.py::test_long_cigar; MM_F_NO_PRINT_2ND and MM_F_2_IO_THREADS also the presets in test_d_api.py::test_set_opt_and_mapopt_init. |
| `MM_F_OUT_JUNC` | set by --write-junc | test_d_api.py::test_mm_map_file_matches_cli (each flag gives its option's output); MM_F_LONG_CIGAR, MM_F_SOFTCLIP, MM_F_SECONDARY_SEQ also test_d_api.py::test_long_cigar; MM_F_NO_PRINT_2ND and MM_F_2_IO_THREADS also the presets in test_d_api.py::test_set_opt_and_mapopt_init. |
| `MM_CIGAR_MATCH` | CIGAR operation code in mm_extra_t.cigar | test_d_api.py::test_macros_and_version (BAM codes), test_d_api.py::test_long_cigar (CG:B:I encoding), test_d_api.py::test_mm_map_matches_cli (CIGAR text). |
| `MM_CIGAR_INS` | CIGAR operation code in mm_extra_t.cigar | test_d_api.py::test_macros_and_version (BAM codes), test_d_api.py::test_long_cigar (CG:B:I encoding), test_d_api.py::test_mm_map_matches_cli (CIGAR text) |
| `MM_CIGAR_DEL` | CIGAR operation code in mm_extra_t.cigar | test_d_api.py::test_macros_and_version (BAM codes), test_d_api.py::test_long_cigar (CG:B:I encoding), test_d_api.py::test_mm_map_matches_cli (CIGAR text) |
| `MM_CIGAR_N_SKIP` | CIGAR operation code in mm_extra_t.cigar | test_d_api.py::test_macros_and_version (BAM codes), test_d_api.py::test_long_cigar (CG:B:I encoding), test_d_api.py::test_mm_map_matches_cli (CIGAR text) |
| `MM_CIGAR_SOFTCLIP` | CIGAR operation code in mm_extra_t.cigar | test_d_api.py::test_macros_and_version (BAM codes), test_d_api.py::test_long_cigar (CG:B:I encoding), test_d_api.py::test_mm_map_matches_cli (CIGAR text) |
| `MM_CIGAR_HARDCLIP` | CIGAR operation code in mm_extra_t.cigar | test_d_api.py::test_macros_and_version (BAM codes), test_d_api.py::test_long_cigar (CG:B:I encoding), test_d_api.py::test_mm_map_matches_cli (CIGAR text) |
| `MM_CIGAR_PADDING` | CIGAR operation code in mm_extra_t.cigar | test_d_api.py::test_macros_and_version (BAM codes), test_d_api.py::test_long_cigar (CG:B:I encoding), test_d_api.py::test_mm_map_matches_cli (CIGAR text) |
| `MM_CIGAR_EQ_MATCH` | CIGAR operation code in mm_extra_t.cigar | test_d_api.py::test_macros_and_version (BAM codes), test_d_api.py::test_long_cigar (CG:B:I encoding), test_d_api.py::test_mm_map_matches_cli (CIGAR text) |
| `MM_CIGAR_X_MISMATCH` | CIGAR operation code in mm_extra_t.cigar | test_d_api.py::test_macros_and_version (BAM codes), test_d_api.py::test_long_cigar (CG:B:I encoding), test_d_api.py::test_mm_map_matches_cli (CIGAR text). |
| `MM_CIGAR_STR` | CIGAR operator letters | test_d_api.py::test_macros_and_version (BAM codes), test_d_api.py::test_long_cigar (CG:B:I encoding), test_d_api.py::test_mm_map_matches_cli (CIGAR text). |

### D: minimap.h globals and opaque types (4)

| Name | Notes | Tested |
|:--|:--|:--|
| `mm_verbose` | -v; 3 by default in the CLI | test_d_cli.py::test_verbosity, test_d_api.py::test_check_opt. |
| `mm_dbg_flag` | --print-*, --no-kalloc | test_d_api.py::test_thread_buffer, test_d_api.py::test_mm_map_matches_cli (MM_DBG_NO_KALLOC), test_d_cli.py::test_debug_options. |
| `mm_realtime0` | wall-clock origin for stderr timings | not covered: it only offsets the times in stderr progress lines, which tests/README.md says not to pin. |
| `mm_tbuf_t` / `struct mm_tbuf_s` (`km`, `rep_len`, `frag_gap`) |  | partial: km (test_d_api.py::test_thread_buffer), rep_len = rl:i (test_d_api.py::test_mm_map_matches_cli); frag_gap is an internal chaining bound with no documented value |

### D: minimap.h struct fields (5)

Fields of `mm_reg1_t` and `mm_extra_t` are results a C driver can check against the CLI or the simulator's truth; fields of the option structs name the option that sets them.

| Field | Set by / notes | Tested |
|:--|:--|:--|
| `mm_reg1_t.sam_pri` |  | test_d_api.py::test_mm_map_matches_cli (one per read, the SAM primary), test_d_api.py::test_long_cigar. |
| `mm_mapopt_t.flag` | MM_F_* flags | test_d_api.py::test_mm_map_file_matches_cli, test_d_api.py::test_set_opt_and_mapopt_init. |
| `mm_mapopt_t.seed` | --seed | test_d_api.py::test_mm_map_file_matches_cli (seed 4 = --seed 4), test_d_api.py::test_set_opt_and_mapopt_init (11). |
| `mm_mapopt_t.mini_batch_size` | -K | test_d_api.py::test_mm_map_file_matches_cli, test_d_api.py::test_set_opt_and_mapopt_init. |
| `mm_mapopt_t.cap_kalloc` | --cap-kalloc | test_d_api.py::test_set_opt_and_mapopt_init (500m default, kept by presets), test_d_invariance.py::test_allocator_options (through --cap-kalloc). |

### D: Non-static functions (43)

Defined without `static` (prototypes in mmpriv.h, bseq.h, sdust.h, ksw2.h, kalloc.h, kthread.h). A row is ticked when a test exercises the function with a checked result, through the CLI, the API or a C driver.

| File | Function | Notes | Tested |
|:--|:--|:--|:--|
| options.c | `mm_mapopt_init` | deprecated API (mm_set_opt(0, ...)) | test_d_api.py::test_set_opt_and_mapopt_init. |
| options.c | `mm_set_opt` |  | test_d_api.py::test_set_opt_and_mapopt_init. |
| options.c | `mm_check_opt` |  | test_d_api.py::test_check_opt, test_d_cli.py::test_option_check_errors, test_d_cli.py::test_option_check_boundaries. |
| map.c | `mm_tbuf_init` |  | test_d_api.py::test_thread_buffer, test_d_api.py::test_mm_map_matches_cli. |
| map.c | `mm_tbuf_destroy` |  | test_d_api.py::test_thread_buffer, test_d_api.py::test_mm_map_matches_cli. |
| map.c | `mm_tbuf_get_km` |  | test_d_api.py::test_thread_buffer, test_d_api.py::test_mm_map_matches_cli. |
| map.c | `mm_map` |  | test_d_api.py::test_mm_map_matches_cli, test_d_api.py::test_writers_match_cli, test_d_api.py::test_minimap2_lite. |
| map.c | `mm_map_file` |  | test_d_api.py::test_mm_map_file_matches_cli (11 option sets, return values). |
| format.c | `mm_write_sam_hdr` |  | test_d_api.py::test_sam_header_writer, test_d_sam.py::test_header_lines, test_d_sam.py::test_read_group. |
| format.c | `mm_write_junc` |  | test_d_paf.py::test_write_junc. |
| format.c | `mm_gen_cs_ds_or_MD` |  | test_d_api.py::test_mm_map_matches_cli. |
| format.c | `mm_gen_cs_or_MD` |  | test_d_api.py::test_mm_map_matches_cli. |
| format.c | `mm_gen_cs` |  | test_d_api.py::test_mm_map_matches_cli. |
| format.c | `mm_gen_ds` |  | test_d_api.py::test_mm_map_matches_cli. |
| format.c | `mm_gen_MD` |  | test_d_api.py::test_mm_map_matches_cli. |
| format.c | `mm_write_paf4` |  | test_d_api.py::test_writers_match_cli, all paf tests (paf4 through the CLI, including fragment names). |
| format.c | `mm_write_paf3` |  | test_d_api.py::test_writers_match_cli, all paf tests (paf4 through the CLI, including fragment names). |
| format.c | `mm_write_paf` |  | test_d_api.py::test_writers_match_cli, all paf tests (paf4 through the CLI, including fragment names). |
| format.c | `mm_write_sam3` |  | test_d_api.py::test_writers_match_cli, test_d_api.py::test_long_cigar, all sam tests. |
| format.c | `mm_write_sam2` |  | test_d_api.py::test_writers_match_cli, test_d_api.py::test_long_cigar, all sam tests. |
| format.c | `mm_write_sam` |  | test_d_api.py::test_writers_match_cli, test_d_api.py::test_long_cigar, all sam tests. |
| misc.c | `gettimeofday` | Windows only (_WIN32) | not covered: Windows only (_WIN32). |
| misc.c | `cputime` | two definitions (Windows / POSIX) | (POSIX definition): test_d_api.py::test_misc; the Windows definition is not built here. |
| misc.c | `peakrss` | ru_maxrss units differ by platform: never assert values | partial: test_d_api.py::test_misc checks only that it is non-negative and does not decrease (ru_maxrss is 0 in this sandbox and its units differ by platform). |
| misc.c | `realtime` |  | test_d_api.py::test_misc. |
| misc.c | `mm_err_puts` |  | test_d_api.py::test_misc, test_d_cli.py::test_write_error_exits_nonzero. |
| misc.c | `mm_err_fwrite` |  | test_d_api.py::test_misc (round trip); their error branches are not reachable on a regular file. |
| misc.c | `mm_err_fread` |  | test_d_api.py::test_misc (round trip); their error branches are not reachable on a regular file. |
| kthread.c | `kt_for` |  | test_d_api.py::test_kthread; also every multi-threaded test. |
| kthread.c | `kt_pipeline` |  | test_d_api.py::test_kthread; also every multi-threaded test. |
| kalloc.c | `km_init2` | kalloc internals are not pinned (speed work may change them) | test_d_api.py::test_kalloc (contents, alignment, zeroing, NULL pool; km_stat only as capacity >= live bytes, no internals). |
| kalloc.c | `km_init` |  | test_d_api.py::test_kalloc (contents, alignment, zeroing, NULL pool; km_stat only as capacity >= live bytes, no internals). |
| kalloc.c | `km_destroy` |  | test_d_api.py::test_kalloc (contents, alignment, zeroing, NULL pool; km_stat only as capacity >= live bytes, no internals). |
| kalloc.c | `kfree` |  | test_d_api.py::test_kalloc (contents, alignment, zeroing, NULL pool; km_stat only as capacity >= live bytes, no internals). |
| kalloc.c | `kmalloc` |  | test_d_api.py::test_kalloc (contents, alignment, zeroing, NULL pool; km_stat only as capacity >= live bytes, no internals). |
| kalloc.c | `kcalloc` |  | test_d_api.py::test_kalloc (contents, alignment, zeroing, NULL pool; km_stat only as capacity >= live bytes, no internals). |
| kalloc.c | `krealloc` |  | test_d_api.py::test_kalloc (contents, alignment, zeroing, NULL pool; km_stat only as capacity >= live bytes, no internals). |
| kalloc.c | `krelocate` |  | test_d_api.py::test_kalloc (contents, alignment, zeroing, NULL pool; km_stat only as capacity >= live bytes, no internals). |
| kalloc.c | `km_stat` |  | test_d_api.py::test_kalloc (contents, alignment, zeroing, NULL pool; km_stat only as capacity >= live bytes, no internals). |
| kalloc.c | `km_stat_print` |  | partial: test_d_api.py::test_kalloc runs it; what it prints is kalloc's debug output and is not checked |
| main.c | `liftrlimit` | Linux: raises RLIMIT_AS to the hard limit | (Linux): test_d_cli.py::test_address_space_limit_is_lifted. |
| main.c | `main` | the minimap2 CLI | test_d_cli.py, test_d_sam.py, test_d_paf.py, test_d_invariance.py. |
| example.c | `main` | `minimap2-lite` of `make extra` | test_d_api.py::test_minimap2_lite. |

### D: Programs and scripts (2)

| Target | Notes | Tested |
|:--|:--|:--|
| `minimap2-lite <target> <query>` | example.c; PAF-like lines with cg:Z | test_d_api.py::test_minimap2_lite (FASTA and .mmi target, equal to `minimap2 -c --no-hash-name`). |
| `minimap2-lite` (no argument) | usage, exit 1 | test_d_api.py::test_minimap2_lite. |

### D: Behaviours not tied to one name (7)

| Behaviour | Notes | Tested |
|:--|:--|:--|
| SAM header: @HD, @SQ, @RG (-R), @PG with CL |  | test_d_sam.py::test_header_lines, test_d_sam.py::test_read_group, test_d_api.py::test_sam_header_writer |
| output identical across -t and -K | observed for SE, PE, splice and PAF on simulated data | test_d_invariance.py (SE SAM and PAF, PE sr, splice); test_d_invariance.py::test_larger_input |
| unmapped reads in SAM; --sam-hit-only; --paf-no-hit line format |  | test_d_sam.py::test_unmapped_records, test_d_paf.py::test_paf_no_hit. |
| secondary (SEQ *), supplementary (hard clips; -Y soft) |  | test_d_sam.py::test_secondary_records, test_d_sam.py::test_secondary_seq, test_d_sam.py::test_supplementary_and_sa, test_d_sam.py::test_soft_clipped_supplementary; test_d_sam.py::test_sa_gives_each_alignments_own_mapq |
| -o FILE equals stdout; -o to an unwritable path is an error |  | test_d_cli.py::test_output_file. |
| option errors from mm_check_opt (each returns 1 with a message) |  | test_d_cli.py::test_option_check_errors (33 cases), test_d_cli.py::test_option_check_boundaries (24), test_d_api.py::test_check_opt. |
| libminimap2 used as example.c does (C drivers) |  | test_d_api.py::test_mm_map_matches_cli, test_d_api.py::test_minimap2_lite. |

## Area E

### E: Programs and scripts (34)

| Target | Notes | Tested |
|:--|:--|:--|
| `paftools.js view` |  | needs k8: test_e_paftools.py::test_view_blast_like, test_e_paftools.py::test_view_line_length, test_e_paftools.py::test_view_lastz_cigar_and_maf, usage |
| `paftools.js sam2paf` |  | needs k8: test_e_paftools.py::test_sam2paf (-p -P -U -L, MD->cs, =/X, clips, /1), pilot test_pilot.py::test_paftools_sam2paf_agrees_with_paf_output |
| `paftools.js delta2paf` |  | needs k8: test_e_paftools.py::test_delta2paf |
| `paftools.js splice2bed` |  | needs k8: test_e_paftools.py::test_splice2bed (PAF and SAM, -m, colours); `-n` not tested (uses k8 0.x `Map.put`, absent in k8 1.x) |
| `paftools.js gff2bed` |  | needs k8: test_e_paftools.py::test_gff2bed (GTF, GFF3, -s -j -e -G); `-u`/`-g` not tested |
| `paftools.js gff2junc` |  | needs k8: test_e_paftools.py::test_gff2junc (-f) |
| `paftools.js stat` |  | needs k8: test_e_paftools.py::test_stat_paf (-l -c), test_e_paftools.py::test_stat_sam |
| `paftools.js asmstat` |  | needs k8: test_e_paftools.py::test_asmstat (-l, two assemblies, gzip) |
| `paftools.js asmgene` |  | needs k8: test_e_paftools.py::test_asmgene (-c -i -e); `-a` not tested |
| `paftools.js misjoin` |  | needs k8: test_e_paftools.py::test_misjoin (-l -g -c -e) |
| `paftools.js liftover` |  | needs k8: test_e_paftools.py::test_liftover[liftover] (-l -q -d, stdin, _t5) |
| `paftools.js liftOver` | alias of liftover | needs k8: test_e_paftools.py::test_liftover[liftOver] |
| `paftools.js vcfpair` |  | needs k8: test_e_paftools.py::test_vcfpair (-m -g -s) |
| `paftools.js call` |  | needs k8: test_e_paftools.py::test_call, test_e_paftools.py::test_call_vcf (-l -L -f -s) |
| `paftools.js mapeval` |  | needs k8: test_e_paftools.py::test_mapeval_paf (-r -Q), test_e_paftools.py::test_mapeval_sam_pairs_and_modes (-m 0/1, U line) |
| `paftools.js pafcmp` |  | needs k8: test_e_paftools.py::test_pafcmp (-q); the "additional test alignments" count is not asserted (it increments test_e_paftools.py::test_opt.n_out_high, so it always prints 0) |
| `paftools.js bedcov` |  | needs k8: test_e_paftools.py::test_bedcov (-p -e, BED12); `-d` not tested (calls `Interval.dedup`, which does not exist) |
| `paftools.js mason2fq` |  | needs k8: test_e_paftools.py::test_mason2fq |
| `paftools.js sim2bed` |  | needs k8: test_e_paftools.py::test_sim2bed |
| `paftools.js pbsim2fq` |  | needs k8: test_e_paftools.py::test_pbsim2fq |
| `paftools.js badread2fa` |  | needs k8: test_e_paftools.py::test_badread2fa |
| `paftools.js junceval` |  | needs k8: test_e_paftools.py::test_junceval (-l -p -e) |
| `paftools.js exoneval` |  | needs k8: test_e_paftools.py::test_exoneval (-s) |
| `paftools.js ov-eval` |  | needs k8: test_e_paftools.py::test_ov_eval (-l, stdin) |
| `paftools.js vcfstat` |  | needs k8: test_e_paftools.py::test_vcfstat (stdin) |
| `paftools.js sveval` |  | needs k8: test_e_paftools.py::test_sveval (-e -i) |
| `paftools.js vcfsel` |  | needs k8: test_e_paftools.py::test_vcfsel (-l -L) |
| `paftools.js longcs2seq` |  | needs k8: test_e_paftools.py::test_longcs2seq (-q, short cs error) |
| `paftools.js paf2gff` |  | needs k8: test_e_paftools.py::test_paf2gff (-a, start/stop codons, missing AS) |
| `paftools.js version` |  | needs k8: test_e_paftools.py::test_version_is_the_minimap2_version |
| `paftools.js` (no command) | usage, exit 1 | needs k8: test_e_paftools.py::test_no_command_prints_usage_and_exits_1, test_e_paftools.py::test_subcommand_usage |
| `pafcluster.js -c` |  | needs k8: test_e_pafcluster.py::test_min_coverage |
| `pafcluster.js -d` |  | needs k8: test_e_pafcluster.py::test_max_divergence |
| `pafcluster.js -e` |  | needs k8: test_e_pafcluster.py::test_max_difference (defaults: test_e_pafcluster.py::test_default_thresholds, usage: test_e_pafcluster.py::test_usage) |

### E: mappy and python/minimap2.py (68)

| Name | Notes | Tested |
|:--|:--|:--|
| `mappy.__version__` |  | test_e_mappy.py::test_version_matches_the_sources |
| `Aligner(fn_idx_in=)` |  | test_e_mappy.py::test_aligner_from_fasta, test_e_mappy.py::test_index_file_matches_the_command_line (.mmi) |
| `Aligner(preset=)` |  | test_e_mappy.py::test_preset_index_parameters[*] (13 presets), test_e_mappy.py::test_preset_alignments_match_the_command_line[*] (8 presets) |
| `Aligner(k=)` |  | test_e_mappy.py::test_k_and_w_arguments |
| `Aligner(w=)` |  | test_e_mappy.py::test_k_and_w_arguments |
| `Aligner(min_cnt=)` |  | test_e_mappy.py::test_min_cnt_min_chain_score_min_dp_score |
| `Aligner(min_chain_score=)` |  | same |
| `Aligner(min_dp_score=)` |  | same |
| `Aligner(bw=)` |  | test_e_mappy.py::test_bw_short_reads_with_a_deletion, test_e_mappy.py::test_bw_long_and_long_join |
| `Aligner(bw_long=)` |  | test_e_mappy.py::test_bw_long_and_long_join |
| `Aligner(best_n=)` |  | test_e_mappy.py::test_best_n_limits_secondary_alignments |
| `Aligner(n_threads=)` |  | test_e_mappy.py::test_index_file_matches_the_command_line (index bytes and alignments equal for 1 and 4) |
| `Aligner(fn_idx_out=)` |  | test_e_mappy.py::test_index_file_matches_the_command_line, test_e_mappy.py::test_preset_index_parameters, test_e_mappy.py::test_k_and_w_arguments, test_e_mappy.py::test_seq_argument_builds_a_one_sequence_index (no effect with seq=) |
| `Aligner(max_frag_len=)` |  | test_e_mappy.py::test_max_frag_len_pairs_a_repetitive_mate |
| `Aligner(extra_flags=)` |  | test_e_mappy.py::test_extra_flags_strand_and_eqx, test_e_mappy.py::test_bw_long_and_long_join (MM_F_NO_LJOIN) |
| `Aligner(seq=)` |  | test_e_mappy.py::test_seq_argument_builds_a_one_sequence_index (Aligner.seq() on it not tested: a defect of base minimap2); test_e_mappy.py::test_seq_argument_keeps_the_preset_hpc |
| `Aligner(scoring=)` |  | test_e_mappy.py::test_scoring_and_sc_ambi_decide_end_extension[*] (4, 6 and 7 elements) |
| `Aligner(sc_ambi=)` |  | same |
| `Aligner(max_chain_skip=)` |  | test_e_mappy.py::test_max_chain_skip |
| `Aligner.__bool__` | False when the index failed to load/build | test_e_mappy.py::test_missing_index_is_false, test_e_mappy.py::test_aligner_from_fasta |
| `Aligner.map(seq=)` |  | most tests |
| `Aligner.map(seq2=)` |  | test_e_mappy.py::test_pairs_seq2_and_read_num, test_e_mappy.py::test_max_frag_len_pairs_a_repetitive_mate, test_e_mappy.py::test_thread_buffer_reuse (mate 2 query interval/cs only when full length: a defect of base minimap2); odd-length mates in test_e_mappy.py::test_pairs_seq2_and_read_num |
| `Aligner.map(name=)` |  | test_e_mappy.py::test_name_seeds_the_tie_breaking_hash |
| `Aligner.map(buf=)` |  | test_e_mappy.py::test_thread_buffer_reuse |
| `Aligner.map(cs=)` |  | test_e_mappy.py::test_cs_md_and_ds, test_e_mappy.py::test_cs_md_ds_match_the_command_line, test_e_mappy.py::test_pairs_seq2_and_read_num |
| `Aligner.map(ds=)` |  | test_e_mappy.py::test_cs_md_and_ds, test_e_mappy.py::test_cs_md_ds_match_the_command_line |
| `Aligner.map(MD=)` |  | same |
| `Aligner.map(max_frag_len=)` |  | test_e_mappy.py::test_max_frag_len_pairs_a_repetitive_mate |
| `Aligner.map(extra_flags=)` |  | test_e_mappy.py::test_extra_flags_strand_and_eqx |
| `Aligner.seq(name, start, end)` | None for an unknown name or empty range | test_e_mappy.py::test_aligner_from_fasta (bounds, clipping, unknown name, empty range, N run), test_e_mappy.py::test_index_file_matches_the_command_line |
| `Aligner.k` |  | test_e_mappy.py::test_aligner_from_fasta, test_e_mappy.py::test_preset_index_parameters, test_e_mappy.py::test_k_and_w_arguments |
| `Aligner.w` |  | same |
| `Aligner.n_seq` |  | test_e_mappy.py::test_aligner_from_fasta, test_e_mappy.py::test_seq_argument_builds_a_one_sequence_index |
| `Aligner.seq_names` |  | same, test_e_mappy.py::test_missing_index_is_false |
| `Alignment.ctg` |  | test_e_mappy.py::test_alignments_against_truth_and_recomputed_fields (truth and genome recomputation), test_e_mappy.py::test_preset_alignments_match_the_command_line (equal to the PAF line) |
| `Alignment.ctg_len` |  | test_e_mappy.py::test_alignments_against_truth_and_recomputed_fields (truth and genome recomputation), test_e_mappy.py::test_preset_alignments_match_the_command_line (equal to the PAF line) |
| `Alignment.r_st` |  | test_e_mappy.py::test_alignments_against_truth_and_recomputed_fields (truth and genome recomputation), test_e_mappy.py::test_preset_alignments_match_the_command_line (equal to the PAF line) |
| `Alignment.r_en` |  | test_e_mappy.py::test_alignments_against_truth_and_recomputed_fields (truth and genome recomputation), test_e_mappy.py::test_preset_alignments_match_the_command_line (equal to the PAF line) |
| `Alignment.strand` |  | test_e_mappy.py::test_alignments_against_truth_and_recomputed_fields (truth and genome recomputation), test_e_mappy.py::test_preset_alignments_match_the_command_line (equal to the PAF line), test_e_mappy.py::test_extra_flags_strand_and_eqx |
| `Alignment.trans_strand` |  | test_e_mappy.py::test_trans_strand_of_spliced_reads (truth and ts:A), 0 for genomic reads in test_e_mappy.py::test_alignments_against_truth_and_recomputed_fields |
| `Alignment.blen` |  | test_e_mappy.py::test_alignments_against_truth_and_recomputed_fields (truth and genome recomputation), test_e_mappy.py::test_preset_alignments_match_the_command_line (equal to the PAF line) |
| `Alignment.mlen` |  | test_e_mappy.py::test_alignments_against_truth_and_recomputed_fields (truth and genome recomputation), test_e_mappy.py::test_preset_alignments_match_the_command_line (equal to the PAF line) |
| `Alignment.NM` |  | test_e_mappy.py::test_alignments_against_truth_and_recomputed_fields (truth and genome recomputation), test_e_mappy.py::test_preset_alignments_match_the_command_line (equal to the PAF line) (ambiguous bases included) |
| `Alignment.is_primary` |  | test_e_mappy.py::test_alignments_against_truth_and_recomputed_fields (truth and genome recomputation), test_e_mappy.py::test_preset_alignments_match_the_command_line (equal to the PAF line), test_e_mappy.py::test_best_n_limits_secondary_alignments |
| `Alignment.q_st` |  | test_e_mappy.py::test_alignments_against_truth_and_recomputed_fields (truth and genome recomputation), test_e_mappy.py::test_preset_alignments_match_the_command_line (equal to the PAF line), test_e_mappy.py::test_pairs_seq2_and_read_num (mate 2 only full-length) |
| `Alignment.q_en` |  | same as q_st |
| `Alignment.mapq` |  | test_e_mappy.py::test_alignments_against_truth_and_recomputed_fields (truth and genome recomputation), test_e_mappy.py::test_preset_alignments_match_the_command_line (equal to the PAF line), test_e_mappy.py::test_max_frag_len_pairs_a_repetitive_mate |
| `Alignment.cigar` |  | test_e_mappy.py::test_alignments_against_truth_and_recomputed_fields (truth and genome recomputation), test_e_mappy.py::test_preset_alignments_match_the_command_line (equal to the PAF line), test_e_mappy.py::test_extra_flags_strand_and_eqx (=/X ops) |
| `Alignment.read_num` |  | test_e_mappy.py::test_pairs_seq2_and_read_num, test_e_mappy.py::test_max_frag_len_pairs_a_repetitive_mate |
| `Alignment.cs` |  | test_e_mappy.py::test_cs_md_and_ds (recomputed), test_e_mappy.py::test_cs_md_ds_match_the_command_line |
| `Alignment.ds` |  | test_e_mappy.py::test_cs_md_and_ds (brackets stripped = cs), test_e_mappy.py::test_cs_md_ds_match_the_command_line |
| `Alignment.MD` |  | test_e_mappy.py::test_cs_md_and_ds (recomputed), test_e_mappy.py::test_cs_md_ds_match_the_command_line |
| `Alignment.cigar_str` |  | test_e_mappy.py::test_alignments_against_truth_and_recomputed_fields (truth and genome recomputation), test_e_mappy.py::test_preset_alignments_match_the_command_line (equal to the PAF line), test_e_mappy.py::test_scoring_and_sc_ambi_decide_end_extension |
| `Alignment.__str__` | PAF-like line (q_st onward) with tp, ts, cg, cs/ds/MD | test_e_mappy.py::test_str_is_a_paf_line_from_column_3, test_e_mappy.py::test_trans_strand_of_spliced_reads (ts:A) |
| `ThreadBuffer()` | reusable per-thread buffer passed as map(buf=) | test_e_mappy.py::test_thread_buffer_reuse (single and paired, one buffer for all calls), test_e_mappy.py::test_threads_with_a_buffer_each (four Python threads, one buffer each) |
| `fastx_read(fn, read_comment=False)` | plain or gzip FASTA/FASTQ; '-' is stdin | test_e_mappy.py::test_fastx_read_fasta_fastq_gzip_and_comments, test_e_mappy.py::test_fastx_read_stdin |
| `revcomp(seq)` |  | test_e_mappy.py::test_revcomp |
| `verbose(v=None)` | get/set mm_verbose | test_e_mappy.py::test_verbose_get_and_set |
| `python/minimap2.py -x` |  | test_e_minimap2_py.py::test_preset_k_w_and_tags, test_e_minimap2_py.py::test_band_width |
| `python/minimap2.py -n` |  | test_e_minimap2_py.py::test_min_count, test_e_minimap2_py.py::test_k_and_w_set_the_seed_density |
| `python/minimap2.py -m` | known issue: -m is parsed but not passed on | not covered: parsed but never passed to mappy (a defect of base minimap2) |
| `python/minimap2.py -r` |  | test_e_minimap2_py.py::test_band_width |
| `python/minimap2.py -k` |  | test_e_minimap2_py.py::test_k_and_w_set_the_seed_density, test_e_minimap2_py.py::test_preset_k_w_and_tags |
| `python/minimap2.py -w` |  | test_e_minimap2_py.py::test_k_and_w_set_the_seed_density, test_e_minimap2_py.py::test_preset_k_w_and_tags |
| `python/minimap2.py -c` |  | test_e_minimap2_py.py::test_preset_k_w_and_tags |
| `python/minimap2.py -d` |  | test_e_minimap2_py.py::test_preset_k_w_and_tags |
| `python/minimap2.py -M` |  | test_e_minimap2_py.py::test_preset_k_w_and_tags |
| `python/minimap2.py` usage | fewer than two arguments: usage, exit 1 | test_e_minimap2_py.py::test_usage_with_fewer_than_two_arguments, test_e_minimap2_py.py::test_default_output_is_paf_like |
