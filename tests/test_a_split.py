"""Area A: --split-prefix (splitidx.c, mm_split_merge in map.c).

minimap2.1: with a multi-part index (-I), minimap2 maps the queries against
each part in turn; --split-prefix writes each part's hits to temporary files
and merges them at the end, so the output has one header with every target
sequence and one primary alignment per read. For reads from unique sequence
the merged alignments equal those against a single-part index.

Not tested: --split-prefix with two or more single-end query files aborts
in base minimap2 (known-bugs/a-split-prefix-several-query-files.md).
"""

import os

import pytest

import helpers_a as A
import minimap2_helpers as H


@pytest.fixture(scope="module")
def split_reads(synth, tmp_path_factory):
	d = tmp_path_factory.mktemp("split")
	ont = H.simulate_long_reads(synth.genome, 30, profile="ont", seed=601, length=(1000, 3000))
	uniq = set(r.name for r in ont if synth.genome.is_unique(r.contig, r.pos, r.end))
	pairs = H.simulate_pairs(synth.genome, 50, seed=602, unique=True)
	r1, r2 = H.write_pairs(str(d / "r1.fq"), str(d / "r2.fq"), pairs)
	return dict(ont=H.write_fastq(str(d / "ont.fq"), ont), uniq=uniq, r1=r1, r2=r2, dir=str(d))


def _fields(rec):
	return (rec.flag, rec.rname, rec.pos, rec.cigar, rec.rnext, rec.pnext, rec.tlen)


def test_split_prefix_merges_a_multi_part_index(synth, split_reads, tmp_path):
	"""-I 60k splits the standard genome into two parts (chr1; chr2 and the short contigs). With
	--split-prefix the SAM has @SQ lines for all four contigs in order and one primary per read; every
	primary has the flag, position, CIGAR and mate fields it has against the single-part index, and the
	same MAPQ for reads from unique sequence; the temporary files are removed; PAF is merged the same way."""
	pre = str(tmp_path / "tmp" / "pre")
	os.makedirs(os.path.dirname(pre))
	single = H.map_reads(synth.fasta, split_reads["ont"], ["-a"]).sam()
	r = H.map_reads(synth.fasta, split_reads["ont"], ["-a", "-I", "60k", "--split-prefix", pre])
	split = r.sam()
	assert split.sq == [(n, synth.genome.length(n)) for n in synth.genome.names]
	assert not r.messages
	assert os.listdir(os.path.dirname(pre)) == []
	ps, pm = single.primaries(), split.primaries()
	assert sorted(pm) == sorted(ps)
	for k in ps:
		assert _fields(pm[k]) == _fields(ps[k]), (pm[k].line[:100], ps[k].line[:100])
		if k[0] in split_reads["uniq"]:
			assert pm[k].mapq == ps[k].mapq
	assert H.sam_invariant_errors(split, synth.genome) == []
	paf_single = [p for p in H.map_reads(synth.fasta, split_reads["ont"], ["-c"]).paf() if p.is_primary]
	paf_split = [p for p in H.map_reads(synth.fasta, split_reads["ont"], ["-c", "-I", "60k", "--split-prefix", pre]).paf()
		if p.is_primary]
	cols = lambda p: p.line.split("\t")[:11]
	assert sorted(cols(p) for p in paf_split) == sorted(cols(p) for p in paf_single)
	assert os.listdir(os.path.dirname(pre)) == []


def test_split_prefix_three_parts(tmp_path):
	"""Six contigs (3, 5, 2, 7, 1 and 4 kbp) split by -I 6k into three parts: with --split-prefix each
	read (1% substitutions) is reported once, as a primary on its own contig at its true place, the rid
	of the second and third parts shifted past the earlier parts' sequences, and the SAM equals the one
	from a single-part index."""
	rng = H.Rng(621)
	recs = [("p%d" % i, rng.dna(l)) for i, l in enumerate((3000, 5000, 2000, 7000, 1000, 4000))]
	fa = A.write_fa(tmp_path / "six.fa", recs)
	reads = []
	for i, (n, s) in enumerate(recs):
		for j in range(2):
			st = 100 + 300 * j
			q = H.mutate(rng, s[st:st + 600], 0.01)
			reads.append(("r%d_%d" % (i, j), q if j == 0 else H.revcomp(q), n, st))
	rf = H.write_fasta(str(tmp_path / "r.fa"), [(n, q) for n, q, _, _ in reads])
	single = H.map_reads(fa, rf, ["-a"])
	split = H.map_reads(fa, rf, ["-a", "-I", "6k", "--split-prefix", str(tmp_path / "three")])
	sam = split.sam()
	assert sam.sq == [(n, len(s)) for n, s in recs]
	prim = sam.primaries()
	assert len(prim) == len(reads) and len(sam.records) == len(reads)
	for name, _, ctg, st in reads:
		r = prim[(name, 0)]
		assert (r.rname, r.unclipped_start0, r.is_reverse) == (ctg, st, name.endswith("_1")), r.line[:80]
	assert H.sam_body(split.stdout) == H.sam_body(single.stdout)
	assert not [f for f in os.listdir(str(tmp_path)) if f.startswith("three.")]


def test_split_prefix_paired_reads(synth, split_reads, tmp_path):
	"""Paired short reads (-x sr, two files) with a two-part index and --split-prefix: mates are merged
	as pairs, with the flags, positions, CIGARs, mate fields and TLEN of the single-part mapping and
	consistent pair fields."""
	pre = str(tmp_path / "pe")
	reads = [split_reads["r1"], split_reads["r2"]]
	single = H.map_reads(synth.fasta, reads, ["-a"], preset="sr").sam()
	split = H.map_reads(synth.fasta, reads, ["-a", "-I", "60k", "--split-prefix", pre], preset="sr").sam()
	assert split.sq == single.sq
	ps, pm = single.primaries(), split.primaries()
	assert sorted(pm) == sorted(ps) and len(ps) == 100
	for k in ps:
		assert _fields(pm[k]) == _fields(ps[k]), (pm[k].line[:100], ps[k].line[:100])
	assert H.pair_invariant_errors(split) == []
	assert not [f for f in os.listdir(str(tmp_path)) if f.startswith("pe.")]


def test_split_prefix_with_a_single_part(mt, tmp_path):
	"""With an index of one part --split-prefix changes nothing but the header order: the same records
	and the same header lines."""
	plain = H.map_reads(mt.human, mt.orang, ["-a"]).stdout
	split = H.map_reads(mt.human, mt.orang, ["-a", "--split-prefix", str(tmp_path / "one")]).stdout
	assert H.sam_body(split) == H.sam_body(plain) != ""
	hdr = lambda t: sorted(l for l in t.splitlines() if l.startswith("@") and not l.startswith("@PG"))
	assert hdr(split) == hdr(plain)
	assert H.map_reads(mt.human, mt.orang, ["--split-prefix", str(tmp_path / "two")]).stdout == H.map_reads(mt.human, mt.orang).stdout


def test_split_prefix_rejects_cs_and_md(mt, tmp_path):
	"""minimap2.1 / mm_check_opt: --cs and --MD do not work with --split-prefix: an error and exit 1."""
	for opt in ("--cs", "--MD"):
		r = H.map_reads(mt.human, mt.orang, ["-a", opt, "--split-prefix", str(tmp_path / "x")], check=False)
		assert r.returncode == 1 and r.messages and H.sam_body(r.stdout) == ""


def test_c_api_split_files(tmp_path):
	"""mm_split_init writes a temporary file under PREFIX for each part with its k and sequence names and
	lengths; mm_split_merge_prep reads all parts back into one sequence list (NULL for zero parts or a
	missing file); mm_split_rm_tmp deletes the files, leaving nothing beside the FASTA (the file names
	are not checked)."""
	rng = H.Rng(611)
	recs = [("c%d" % i, rng.dna(l)) for i, l in enumerate((3000, 5000, 2000, 7000, 1000))]
	fa = A.write_fa(tmp_path / "s.fa", recs)
	H.assert_c_case("drv_a_index", "split", fa, 6000, str(tmp_path / "pfx"), 3)
	H.assert_c_case("drv_a_index", "split", fa, 8000000000, str(tmp_path / "one"), 1)
	assert sorted(os.listdir(str(tmp_path))) == ["s.fa"]


@pytest.mark.pinned
def test_split_prefix_output_pin(synth, split_reads, tmp_path):
	"""Pins the current output: the SAM records (SEQ and QUAL, copies of the input, replaced by '*') of
	ONT-like reads mapped with a two-part index (-I 60k) and --split-prefix, the merge step a speed change
	to map.c could alter."""
	r = H.map_reads(synth.fasta, split_reads["ont"], ["-a", "-I", "60k", "--split-prefix", str(tmp_path / "p")])
	body = []
	for line in H.sam_body(r.stdout).splitlines():
		f = line.split("\t")
		f[9] = f[10] = "*"
		body.append("\t".join(f) + "\n")
	H.check_pinned("a/split-prefix-ont.sam", "".join(body))
