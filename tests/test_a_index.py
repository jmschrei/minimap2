"""Area A: index building and the .mmi format (index.c, options.c's index
parameters, main.c's -k -w -H -I -d --idx-no-seq --bucket-bits).

The .mmi files minimap2 writes are decoded by helpers_a.read_mmi (the layout
of mm_idx_dump: MM_IDX_MAGIC, w, k, b, n_seq, flag, names and lengths, one
hash table per bucket, the 4-bit packed sequence) and compared with the
reference minimizer scan of helpers_a, the FASTA they came from and the
documented preset values of minimap2.1. The index API of minimap.h runs
through tests/c/drv_a_index.c.
"""

import os
import re

import pytest

import conftest
import helpers_a as A
import minimap2_helpers as H

# index parameters (w, k, -H) of every preset: minimap2.1 for the documented ones, options.c for the
# undocumented aliases and lr:hqae / map-iclr-prerender
PRESETS = {
	"map-ont": (10, 15, False),  # the default: -k15 -w10
	"lr": (10, 15, False),  # alias of map-ont
	"ava-ont": (5, 15, False),  # -k15 -Xw5
	"map-pb": (10, 19, True),  # -Hk19
	"map10k": (10, 19, True),  # alias of map-pb
	"ava-pb": (5, 19, True),  # -Hk19 -Xw5
	"lr:hq": (19, 19, False),  # -k19 -w19
	"map-hifi": (19, 19, False),  # -xlr:hq ...
	"map-ccs": (19, 19, False),  # alias of map-hifi
	"lr:hqae": (51, 25, False),  # options.c: k 25, w 51
	"map-iclr-prerender": (10, 15, False),  # options.c: k 15
	"map-iclr": (10, 19, False),  # -k19
	"asm5": (19, 19, False),  # -k19 -w19
	"asm10": (19, 19, False),  # -k19 -w19
	"asm20": (10, 19, False),  # -k19 -w10
	"short": (11, 21, False),  # alias of sr
	"sr": (11, 21, False),  # -k21 -w11
	"splice": (5, 15, False),  # -k15 -w5
	"splice:hq": (5, 15, False),  # -xsplice ...
	"splice:sr": (5, 15, False),  # -xsplice:hq ...
	"cdna": (5, 15, False),  # alias of splice
}
ALIASES = {"lr": "map-ont", "map10k": "map-pb", "map-ccs": "map-hifi", "short": "sr", "cdna": "splice"}


@pytest.fixture(scope="module")
def clean(tmp_path_factory):
	"""helpers_a.clean_genome() written as FASTA (60 bp lines)."""
	d = tmp_path_factory.mktemp("clean")
	recs = A.clean_genome()
	return conftest._Data(recs=recs, fasta=A.write_fa(d / "clean.fa", recs, line_width=60), dir=str(d),
		seqs=dict(recs))


def _dump(fa, mmi, args=()):
	r = H.run_mm2(list(args) + ["-t", "2", "-d", str(mmi), str(fa)])
	assert r.returncode == 0, r.stderr[-2000:]
	assert r.stdout == ""
	return str(mmi)


def _check_names_and_seq(part, recs):
	assert part.names == [n for n, _ in recs]
	assert part.lengths == [len(s) for _, s in recs]
	assert part.seq4 == [A.seq_codes(s) for _, s in recs]


def test_mmi_default_index_holds_the_reference_minimizers(clean, tmp_path):
	"""`minimap2 -d` with default options: the magic MMI\\2, w=10 k=15 b=14 flag 0, every name, length
	and base (4-bit codes; lowercase as uppercase, U as T), and a minimizer table equal to the reference
	(w,k)-minimizers of every sequence (none for the 12 bp contig)."""
	mmi = _dump(clean.fasta, tmp_path / "d.mmi")
	with open(mmi, "rb") as f:
		assert f.read(4) == b"MMI\x02"
	parts = A.read_mmi(mmi)
	assert len(parts) == 1
	p = parts[0]
	assert (p.w, p.k, p.b, p.flag) == (10, 15, 14, 0)
	_check_names_and_seq(p, clean.recs)
	assert p.minimizers == A.expected_index([s for _, s in clean.recs], 10, 15)
	assert 5 not in A.stored_minimizers(p)  # tiny12 has no 15-mer
	assert len(A.stored_minimizers(p)[4]) == 1  # short20: 6 k-mers, fewer than w; its smallest one


@pytest.mark.parametrize("args,w,k,hpc", [
	(["-k", "21", "-w", "11"], 11, 21, False),
	(["-k", "11", "-w", "1"], 1, 11, False),  # every k-mer
	(["-w", "255", "-k", "27"], 255, 27, False),
	(["-H"], 10, 15, True),
	(["-H", "-k", "19", "-w", "5"], 5, 19, True),
])
def test_mmi_k_w_H_options(clean, tmp_path, args, w, k, hpc):
	"""-k, -w and -H set the header (flag bit MM_I_HPC for -H) and the minimizers, which equal the
	reference scan with the same (w, k) on the plain or homopolymer-compressed sequence."""
	p = A.read_mmi(_dump(clean.fasta, tmp_path / "o.mmi", args))[0]
	assert (p.w, p.k, p.flag) == (w, k, A.I_HPC if hpc else 0)
	assert p.b == min(14, 2 * k)
	recs = [(n, s) for n, s in clean.recs]
	for n, s in recs:
		assert not A.SketchRef(s, w, k, hpc).has_ties, n
	assert p.minimizers == A.expected_index([s for _, s in recs], w, k, hpc)


def test_mmi_presets(clean, tmp_path):
	"""Every -x preset indexes with its documented (w, k, -H): minimap2.1 for the documented presets and
	options.c for lr:hqae, map-iclr-prerender and the undocumented aliases (lr, map10k, map-ccs, short,
	cdna), whose .mmi files are byte-identical to their targets'. The minimizers equal the reference scan
	for those parameters."""
	by_digest, files = {}, {}
	for preset, (w, k, hpc) in sorted(PRESETS.items()):
		mmi = _dump(clean.fasta, tmp_path / ("%s.mmi" % preset.replace(":", "_")), ["-x", preset])
		files[preset] = mmi
		by_digest.setdefault(A.file_sha256(mmi), []).append(preset)
	for alias, target in ALIASES.items():
		assert A.file_sha256(files[alias]) == A.file_sha256(files[target]), alias
	seqs = [s for _, s in clean.recs]
	for digest, presets in by_digest.items():
		w, k, hpc = PRESETS[presets[0]]
		assert all(PRESETS[p] == (w, k, hpc) for p in presets), presets
		p = A.read_mmi(files[presets[0]])[0]
		assert (p.w, p.k, p.b, p.flag) == (w, k, 14, A.I_HPC if hpc else 0), presets
		_check_names_and_seq(p, clean.recs)
		assert p.minimizers == A.expected_index(seqs, w, k, hpc), presets


def test_preset_is_applied_before_other_options(clean, tmp_path):
	"""The usage text: -x is always applied before other options, so `-k 17 -x sr` indexes with k 17 and
	sr's w 11, and `-H -x map-ont` keeps -H; a later -x replaces the index parameters of an earlier one."""
	p = A.read_mmi(_dump(clean.fasta, tmp_path / "a.mmi", ["-k", "17", "-x", "sr"]))[0]
	assert (p.w, p.k, p.flag) == (11, 17, 0)
	p = A.read_mmi(_dump(clean.fasta, tmp_path / "b.mmi", ["-H", "-x", "map-ont"]))[0]
	assert (p.w, p.k, p.flag) == (10, 15, A.I_HPC)
	p = A.read_mmi(_dump(clean.fasta, tmp_path / "c.mmi", ["-x", "map-pb", "-x", "sr"]))[0]
	assert (p.w, p.k, p.flag) == (11, 21, 0)


def test_bucket_bits(clean, tmp_path):
	"""--bucket-bits sets the number of hash-table buckets (header b), capped at 2k; the minimizer table
	and the sequences are the same for any bucket count, and so is the mapping output."""
	want = A.expected_index([s for _, s in clean.recs], 10, 15)
	for bb in (4, 10, 16):
		p = A.read_mmi(_dump(clean.fasta, tmp_path / ("b%d.mmi" % bb), ["--bucket-bits", str(bb)]))[0]
		assert p.b == bb
		assert p.minimizers == want
		_check_names_and_seq(p, clean.recs)
	p = A.read_mmi(_dump(clean.fasta, tmp_path / "k5.mmi", ["-k", "5", "-w", "3", "--bucket-bits", "14"]))[0]
	assert (p.b, p.k) == (10, 5)
	reads = H.write_fasta(str(tmp_path / "r.fa"), [("q%d" % i, H.mutate(H.Rng(i), clean.seqs["ctgA"][i * 3000:i * 3000 + 2500], 0.02))
		for i in range(6)])
	base = H.map_reads(clean.fasta, reads, ["-a"]).stdout
	for bb in (4, 20):
		assert H.strip_pg(H.map_reads(clean.fasta, reads, ["-a", "--bucket-bits", str(bb)]).stdout) == H.strip_pg(base)


def test_idx_no_seq(clean, tmp_path):
	"""--idx-no-seq: flag MM_I_NO_SEQ, the same names, lengths and minimizers but no sequence section
	(the file is 4 * ceil(total / 8) bytes shorter); mapping with -a needs the sequences and fails with an
	error and exit 1; PAF without -c gives the same output as the full index."""
	full = _dump(clean.fasta, tmp_path / "full.mmi")
	noseq = _dump(clean.fasta, tmp_path / "noseq.mmi", ["--idx-no-seq"])
	pf, pn = A.read_mmi(full)[0], A.read_mmi(noseq)[0]
	assert pn.flag == A.I_NO_SEQ and pn.seq4 is None
	assert (pn.names, pn.lengths, pn.minimizers) == (pf.names, pf.lengths, pf.minimizers)
	total = sum(pf.lengths)
	assert os.path.getsize(full) - os.path.getsize(noseq) == 4 * ((total + 7) // 8)
	reads = H.write_fasta(str(tmp_path / "r.fa"), [("q%d" % i, clean.seqs["ctgA"][i * 4000:i * 4000 + 3000]) for i in range(4)])
	r = H.map_reads(noseq, reads, ["-a"], check=False)
	assert r.returncode == 1 and r.messages
	r = H.map_reads(noseq, reads, ["-c"], check=False)
	assert r.returncode == 1 and r.messages
	assert H.map_reads(noseq, reads).stdout == H.map_reads(full, reads).stdout != ""


@pytest.fixture(scope="module")
def synth_reads(synth, tmp_path_factory):
	"""ONT-like and HiFi-like reads, FR short-read pairs and spliced reads of the standard genomes."""
	d = tmp_path_factory.mktemp("areads")
	ont = H.write_fastq(str(d / "ont.fq"), H.simulate_long_reads(synth.genome, 12, profile="ont", seed=301, length=(1500, 4000)))
	hifi = H.write_fastq(str(d / "hifi.fq"), H.simulate_long_reads(synth.genome, 12, profile="hifi", seed=302, length=(1500, 4000)))
	pairs = H.simulate_pairs(synth.genome, 40, seed=303)
	r1, r2 = H.write_pairs(str(d / "r1.fq"), str(d / "r2.fq"), pairs)
	asm = H.write_fasta(str(d / "asm.fa"), H.simulate_contigs(synth.genome, 2, length=(8000, 15000), seed=304))
	return conftest._Data(ont=ont, hifi=hifi, r1=r1, r2=r2, asm=asm)


@pytest.mark.parametrize("preset,which", [("map-ont", "ont"), ("map-hifi", "hifi"), ("map-pb", "ont"), ("sr", "pair"),
	("asm5", "asm"), ("splice", "hifi")])
def test_prebuilt_index_equals_on_the_fly(synth, synth_reads, index_of, preset, which):
	"""README: an index built with -d and the same preset gives the same alignments as indexing the
	FASTA on the fly: identical SAM apart from @PG."""
	reads = [synth_reads.r1, synth_reads.r2] if which == "pair" else getattr(synth_reads, which)
	a = H.map_reads(synth.fasta, reads, ["-a"], preset=preset)
	b = H.map_reads(index_of(synth.genome, preset), reads, ["-a"], preset=preset)
	assert H.strip_pg(a.stdout) == H.strip_pg(b.stdout)
	assert len(a.sam().records) > 0


def test_prebuilt_index_overrides_k_w_H(synth, synth_reads, index_of):
	"""minimap2.1 (-d): with an index file as the target, -H, -k, -w are overridden by the values stored in
	the index; minimap2 warns and the output equals that without those options."""
	mmi = index_of(synth.genome, "map-ont")
	a = H.map_reads(mmi, synth_reads.ont, ["-a"])
	b = H.map_reads(mmi, synth_reads.ont, ["-a", "-k", "21", "-w", "11", "-H"])
	assert H.strip_pg(a.stdout) == H.strip_pg(b.stdout)
	assert b.messages and not a.messages


def test_dump_while_mapping(clean, tmp_path):
	"""`minimap2 -d idx.mmi ref.fa reads.fa` writes the same index as `-d` alone and the same alignments as
	mapping without -d (the index keeps the sequences even without -a or -c)."""
	reads = H.write_fasta(str(tmp_path / "r.fa"), [("q%d" % i, clean.seqs["ctgA"][i * 4000 + 7:i * 4000 + 3007]) for i in range(4)])
	alone = _dump(clean.fasta, tmp_path / "alone.mmi")
	both = str(tmp_path / "both.mmi")
	r = H.map_reads(clean.fasta, reads, ["-d", both])
	with open(alone, "rb") as f1, open(both, "rb") as f2:
		assert f1.read() == f2.read()
	assert r.stdout == H.map_reads(clean.fasta, reads).stdout != ""


def test_mmi_round_trip(clean, synth, tmp_path, index_of):
	"""An index loaded and dumped again (`minimap2 -d b.mmi a.mmi`, and mm_idx_load/mm_idx_dump through
	the API) holds the same header, names, lengths, sequence and minimizer table. (The hash tables are
	written in table order, so the bytes may differ.)"""
	for fa, src in ((clean.fasta, None), (synth.fasta, index_of(synth.genome, "map-pb"))):
		a = src or _dump(fa, tmp_path / "a.mmi")
		b = _dump(a, tmp_path / "b.mmi")
		pa, pb = A.read_mmi(a)[0], A.read_mmi(b)[0]
		for f in ("w", "k", "b", "flag", "names", "lengths", "seq4", "minimizers"):
			assert getattr(pa, f) == getattr(pb, f), f
	empty = A.write_fa(tmp_path / "empty.txt", [])
	short = H.write_text(str(tmp_path / "short.bin"), "MMI")
	H.assert_c_case("drv_a_index", "dump_load", clean.fasta, tmp_path / "c1.mmi", tmp_path / "c2.mmi", clean.fasta, empty, short)
	pa, pb = A.read_mmi(str(tmp_path / "c1.mmi"))[0], A.read_mmi(str(tmp_path / "c2.mmi"))[0]
	for f in ("w", "k", "b", "flag", "names", "lengths", "seq4", "minimizers"):
		assert getattr(pa, f) == getattr(pb, f), f
	assert pa.minimizers == A.expected_index([s for _, s in clean.recs], 10, 15)
	with open(str(tmp_path / "c1.mmi"), "rb") as f1, open(_dump(clean.fasta, tmp_path / "cli.mmi"), "rb") as f2:
		assert f1.read() == f2.read()  # mm_idx_dump() through the API writes what -d writes


@pytest.mark.pinned
@pytest.mark.thread_invariance
@pytest.mark.parametrize("preset", ["map-ont", "map-pb", "sr", "asm5", "splice"])
def test_synth_index_bytes(synth, tmp_path, preset):
	"""Pins the current output: the SHA-256 of the .mmi that `minimap2 -x PRESET -d` writes for the
	standard synthetic genome, identical at -t 1 and -t 4. The decoded index holds the genome's names,
	lengths and bases, and minimizers that pass helpers_a.index_problems (real k-mers that are window
	minima; every window covered except at the known edges)."""
	digests = []
	for t in (1, 4):
		mmi = str(tmp_path / ("t%d.mmi" % t))
		H.run_mm2(["-x", preset, "-t", str(t), "-d", mmi, synth.fasta], check=True)
		digests.append(A.file_sha256(mmi))
	assert digests[0] == digests[1]
	H.check_pinned("a/synth-%s.mmi.sha256" % preset.replace(":", "_"), digests[0] + "\n")
	w, k, hpc = PRESETS[preset]
	p = A.read_mmi(mmi)[0]
	recs = [(n, synth.genome[n]) for n in synth.genome.names]
	assert (p.w, p.k, p.flag) == (w, k, A.I_HPC if hpc else 0)
	_check_names_and_seq(p, recs)
	assert A.index_problems(p, recs, w, k, hpc) == []


def _part_layout(lengths, batch):
	"""The documented split of -I: parts are consecutive runs of the target sequences, each part ending
	with the sequence that brings it to batch bases or more."""
	parts, cur, tot = [], [], 0
	for i, l in enumerate(lengths):
		cur.append(i)
		tot += l
		if tot >= batch:
			parts.append(cur)
			cur, tot = [], 0
	if cur:
		parts.append(cur)
	return parts


@pytest.fixture(scope="module")
def parted(tmp_path_factory):
	"""Six random contigs of 3, 5, 2, 7, 1 and 4 kbp (no partial sum equals 6000 or 9000)."""
	d = tmp_path_factory.mktemp("parted")
	rng = H.Rng(77)
	recs = [("p%d" % i, rng.dna(l)) for i, l in enumerate((3000, 5000, 2000, 7000, 1000, 4000))]
	return conftest._Data(recs=recs, fasta=A.write_fa(d / "parted.fa", recs), dir=str(d))


@pytest.mark.parametrize("batch", ["6k", "9000", "1"])
def test_multi_part_index_file(parted, tmp_path, batch):
	"""-I NUM with -d writes one index per part back to back in the .mmi: parts are consecutive runs of
	the target sequences, each ending with the sequence that brings it to NUM bases or more (minimap2.1:
	"Load at most NUM target bases"; see known-bugs/a-index-batch-size.md); each part holds the reference
	minimizers of its own sequences with rid counted within the part."""
	num = int(batch[:-1]) * 1000 if batch.endswith("k") else int(batch)
	mmi = _dump(parted.fasta, tmp_path / "mp.mmi", ["-I", batch])
	parts = A.read_mmi(mmi)
	layout = _part_layout([len(s) for _, s in parted.recs], num)
	assert [p.names for p in parts] == [[parted.recs[i][0] for i in ids] for ids in layout]
	for p, ids in zip(parts, layout):
		recs = [parted.recs[i] for i in ids]
		assert (p.w, p.k) == (10, 15)
		_check_names_and_seq(p, recs)
		assert p.minimizers == A.expected_index([s for _, s in recs], 10, 15)
	out = str(tmp_path / "reader.tsv")
	H.assert_c_case("drv_a_index", "reader", mmi, 15, 10, 0, 14, 8000000000, out)
	H.assert_c_case("drv_a_index", "reader", parted.fasta, 10, 15, 0, 14, num, str(tmp_path / "reader2.tsv"), str(tmp_path / "api.mmi"))
	with open(mmi, "rb") as f1, open(str(tmp_path / "api.mmi"), "rb") as f2:
		assert f1.read() == f2.read()  # the reader's fn_out holds every part, as -d writes them
	for f in (out, str(tmp_path / "reader2.tsv")):
		lines = [l.split("\t") for l in H.read_text(f).splitlines()]
		heads = [l for l in lines if l[0] == "part"]
		assert [int(h[1]) for h in heads] == list(range(len(layout)))
		assert [int(h[2]) for h in heads] == [len(ids) for ids in layout]
		assert [h[7] for h in heads] == ["0"] * (len(layout) - 1) + ["1"]  # mm_idx_reader_eof after each part
		assert [(l[1], int(l[2])) for l in lines if l[0] == "seq"] == [(n, len(s)) for n, s in parted.recs]


def test_multi_part_index_mapping(parted, tmp_path):
	"""Mapping with a multi-part index (-I) without --split-prefix: no @SQ lines and a warning that names
	--split-prefix (main.c); minimap2.1 (-I): the queries are read once per part and mapped against it, so
	the SAM holds one block of records per part in query order, each read mapped (at its true place) in the
	block of its contig's part and unmapped in the others. A prebuilt multi-part .mmi gives the same SAM."""
	reads = H.write_fasta(str(tmp_path / "r.fa"), [("from_p%d" % i, s[200:900]) for i, (_, s) in enumerate(parted.recs)])
	r = H.map_reads(parted.fasta, reads, ["-a", "-I", "6k"])
	sam = r.sam(strict=True)
	assert sam.sq == [] and any("split-prefix" in m for m in r.messages)
	layout = _part_layout([len(s) for _, s in parted.recs], 6000)
	n = len(parted.recs)
	assert len(sam.records) == n * len(layout)
	for j, ids in enumerate(layout):
		block = sam.records[j * n:(j + 1) * n]
		assert [x.qname for x in block] == ["from_p%d" % i for i in range(n)]
		for i, x in enumerate(block):
			if i in ids:
				assert (x.flag, x.rname, x.pos0, x.cigar) == (0, "p%d" % i, 200, "700M"), x.line[:80]
			else:
				assert x.is_unmapped, x.line[:80]
	mmi = _dump(parted.fasta, tmp_path / "mp.mmi", ["-I", "6k"])
	r2 = H.map_reads(mmi, reads, ["-a"])
	assert H.strip_pg(r2.stdout) == H.strip_pg(r.stdout)
	paf = H.map_reads(parted.fasta, reads, ["-I", "6k", "-c"]).paf()
	assert [(p.qname, p.tname, p.tstart, p.tend) for p in paf] == \
		[("from_p%d" % i, "p%d" % i, 200, 900) for ids in layout for i in ids]


def test_c_api_reader_single_part(clean, tmp_path):
	"""mm_idx_reader_open/read/eof/close on a FASTA (built with the given options) and on an .mmi (read
	back): one part with index 0, every name and length, the options copied into the reader, is_idx and
	idx_size (the file size) for an index; opening a missing file returns NULL."""
	mmi = _dump(clean.fasta, tmp_path / "x.mmi", ["-k", "19", "-w", "5", "-H"])
	for target, args in ((clean.fasta, (5, 19, 1, 12, 8000000000)), (mmi, (10, 15, 0, 14, 8000000000))):
		out = str(tmp_path / "reader.tsv")
		H.assert_c_case("drv_a_index", "reader", target, *(list(args) + [out]))
		lines = [l.split("\t") for l in H.read_text(out).splitlines()]
		assert lines[0][:8] == ["part", "0", str(len(clean.recs)), "5", "19", "12" if target == clean.fasta else "14", "1", "1"]
		assert [(l[1], int(l[2])) for l in lines[1:]] == [(n, len(s)) for n, s in clean.recs]


def test_c_api_lookup_names_and_sequences(clean, tmp_path):
	"""mm_idx_get finds every minimizer of every sequence (from mm_sketch) at its position in sorted
	position lists that hold nothing else, and nothing for an absent hash; mm_idx_getseq, getseq_rev and
	getseq2 return the bases of any interval of either strand, clip the end and return -1 past the end;
	mm_idx_name2id returns -2 before mm_idx_index_name, then each sequence's index and -1 for an unknown
	name. Run on a FASTA target and on its .mmi."""
	mmi = _dump(clean.fasta, tmp_path / "x.mmi")
	for target in (clean.fasta, mmi):
		H.assert_c_case("drv_a_index", "lookup", target, clean.fasta, 0)


def test_duplicate_target_names(tmp_path):
	"""Two target sequences with one name: mm_idx_index_name() reports the duplicate and mm_idx_name2id()
	returns the first; the command line warns when it builds the name table (here for --alt), and maps to
	both sequences by index (PAF names them as in the FASTA)."""
	rng = H.Rng(88)
	recs = [("chrA", rng.dna(5000)), ("chrB", rng.dna(5000)), ("chrA", rng.dna(6000))]
	fa = A.write_fa(tmp_path / "dup.fa", recs)
	H.assert_c_case("drv_a_index", "lookup", fa, fa, 1)
	reads = H.write_fasta(str(tmp_path / "r.fa"), [("in_first", recs[0][1][100:2100]), ("in_third", recs[2][1][300:2300])])
	alt = H.write_text(str(tmp_path / "alt.txt"), "chrB\n")
	r = H.map_reads(fa, reads, ["-c", "--alt", alt])
	assert r.messages  # the duplicate-name warning
	paf = dict((p.qname, p) for p in r.paf())
	assert (paf["in_first"].tname, paf["in_first"].tlen, paf["in_first"].tstart) == ("chrA", 5000, 100)
	assert (paf["in_third"].tname, paf["in_third"].tlen, paf["in_third"].tstart) == ("chrA", 6000, 300)


def test_c_api_idx_str_init_build_gen(clean, tmp_path):
	"""mm_idx_str: NULL for no sequence; without names it sets MM_I_NO_NAME and keeps no names through a
	dump and load; -H sets MM_I_HPC; bucket bits default to 14 and are capped at 2k; an empty sequence is
	kept with length 0. mm_idx_init sets w < 1 to 1. mm_idx_build (deprecated) and mm_idx_gen (one call per
	part of a batch, with any indexing mini-batch size) build indexes whose sequences and position lists
	hold every minimizer of every sequence and nothing else; mm_idx_build of a missing file returns NULL."""
	H.assert_c_case("drv_a_index", "idx_str")
	H.assert_c_case("drv_a_index", "build", clean.fasta, 10, 15, 0, str(tmp_path / "missing.fa"))
	H.assert_c_case("drv_a_index", "build", clean.fasta, 5, 19, 1, str(tmp_path / "missing.fa"))
	H.assert_c_case("drv_a_index", "gen", clean.fasta, 10, 15, 14, 8000000000, 1, 1 << 18)
	H.assert_c_case("drv_a_index", "gen", clean.fasta, 10, 15, 14, 8000000000, 1, 500)  # many indexing mini-batches
	H.assert_c_case("drv_a_index", "gen", clean.fasta, 7, 17, 10, 21000, 2, 1 << 18)  # parts: ctgA + ctgB (28 kbp), then the rest


def test_c_api_idxopt_defaults():
	"""mm_idxopt_init: -k 15, -w 10 and -I 8G as minimap2.1 documents, no flag, 14 bucket bits and 50M
	bases per indexing mini-batch."""
	H.assert_c_case("drv_a_index", "idxopt")


def test_c_api_idx_stat(clean):
	"""mm_idx_stat reports the k-mer size, w, the HPC flag, the number of sequences, the number of distinct
	minimizers, the singleton percentage, the average occurrence and spacing and the total length, each
	equal to the value from the reference minimizer table (only the numbers are compared)."""
	r = H.run_c_case("drv_a_index", "stat", [clean.fasta])
	assert r.passed, r.report()
	lines = [l for l in r.stderr.splitlines() if "mm_idx_stat" in l]
	nums = set(re.findall(r"(?<![\w.])-?\d+(?:\.\d+)?", " ".join(l.split("]", 1)[1] for l in lines)))
	want = A.expected_index([s for _, s in clean.recs], 10, 15)
	total = sum(len(v) for v in want.values())
	length = sum(len(s) for _, s in clean.recs)
	single = sum(1 for v in want.values() if len(v) == 1)
	for x in ("15", "10", "0", str(len(clean.recs)), str(len(want)), str(length), "%.2f" % (100.0 * single / len(want)),
			"%.3f" % (total / float(len(want))), "%.3f" % (length / float(total))):
		assert x in nums, (x, lines)


def test_target_errors(tmp_path, mt):
	"""A missing target file: an error naming it and exit 1. A target with no FASTA record (text, or an
	empty file) indexes no sequence, so no read maps (exit 0 or an error). -k or -w of 0: exit 1 with an
	error; -k above 28 (the usage text's limit) or -w of 256 or more: a non-zero exit."""
	missing = str(tmp_path / "no-such-target.fa")
	r = H.map_reads(missing, mt.orang, check=False)
	assert r.returncode == 1 and any("no-such-target.fa" in m for m in r.messages)
	junk = H.write_text(str(tmp_path / "junk.txt"), "hello world\nthis is not a sequence file\n")
	empty = H.write_text(str(tmp_path / "empty.fa"), "")
	for t in (junk, empty):
		r = H.map_reads(t, mt.orang, ["-a"], check=False)
		if r.returncode == 0:
			sam = r.sam()
			assert sam.sq == [] and [x.flag for x in sam.records] == [4]
		r = H.map_reads(t, mt.orang, check=False)
		assert r.stdout == ""
	for args in (["-k", "0"], ["-w", "0"]):
		r = H.map_reads(mt.human, mt.orang, args, check=False)
		assert r.returncode == 1 and r.messages and r.stdout == ""
	for args in (["-k", "29"], ["-w", "256"]):
		assert H.map_reads(mt.human, mt.orang, args, check=False).returncode != 0
