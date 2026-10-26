"""Area C: short genomic reads and pairing. Paired input from two files or one
interleaved file with or without /1 /2 suffixes (README), the fragment mode
(`--frag`, reads with the same name in one file, more than two segments),
`-F`, `--pairing`, `--no-pairing`, `--sr`, the `sr` preset and its alias
`short`, proper-pair flags and mate fields (SAM spec), and pairing rescue of
repetitive mates (mm_pair, mm_select_sub_multi and mm_seg_gen through the CLI).

Truth comes from the simulator (minimap2_helpers.simulate_pairs) and from pairs
cut by hand out of the standard genome (conftest.make_synth)."""

import pytest

import minimap2_helpers as H


@pytest.fixture(scope="module")
def sr_pairs(synth, tmp_path_factory):
	"""FR pairs from unique sequence with fragment lengths from 150 to about 1300 bp, plus two pairs at
	each -F limit tested and one base past it; substitutions only, so every mate aligns end to end and
	TLEN equals the fragment length."""
	d = tmp_path_factory.mktemp("sr_pairs")
	pairs = H.simulate_pairs(synth.genome, 150, seed=111, insert_mean=650, insert_sd=250, unique=True, sub_rate=0.005,
		indel_rate=0)
	for i, F in enumerate((500, 501, 800, 801, 1000, 1001)):
		pairs += H.simulate_pairs(synth.genome, 2, seed=120 + i, insert_mean=F, insert_sd=0, unique=True, sub_rate=0.005,
			indel_rate=0, name_prefix="f%d_" % F)
	assert sorted(p.insert for p in pairs if p.name.startswith("f")) == sorted([500, 501, 800, 801, 1000, 1001] * 2)
	f1, f2 = H.write_pairs(str(d / "r1.fq"), str(d / "r2.fq"), pairs)
	return dict(pairs=pairs, f1=f1, f2=f2, dir=d, reads=dict((p.name, [p.r1.seq, p.r2.seq]) for p in pairs))


@pytest.fixture(scope="module")
def sr_default(synth, sr_pairs):
	return H.map_reads(synth.fasta, [sr_pairs["f1"], sr_pairs["f2"]], ["-a"], preset="sr")


def test_paired_input_forms_agree(synth, sr_pairs, sr_default):
	"""README: two files are read in turn and merged into an interleaved stream; adjacent reads with
	the same name (a /[0-9] suffix trimmed) are a pair. Two files and one interleaved file give the
	same SAM records, with or without /1 and /2 suffixes. The suffixes change the read name, which
	seeds the hash that breaks ties, so across the two namings only the placement fields (FLAG to
	TLEN except MAPQ) must agree. SAM drops the suffix from QNAME and sets 0x40 on mates from the first file and
	0x80 on the second; PAF appends /1 and /2 to the read name, even after an existing suffix."""
	d, pairs = sr_pairs["dir"], sr_pairs["pairs"]
	body = H.sam_body(sr_default.stdout)
	s1, s2 = H.write_pairs(str(d / "s1.fq"), str(d / "s2.fq"), pairs, suffix=True)
	il = H.write_interleaved(str(d / "il.fq"), pairs)
	ils = H.write_interleaved(str(d / "ils.fq"), pairs, suffix=True)
	assert H.sam_body(H.map_reads(synth.fasta, [il], ["-a"], preset="sr").stdout) == body
	sfx = H.sam_body(H.map_reads(synth.fasta, [s1, s2], ["-a"], preset="sr").stdout)
	assert H.sam_body(H.map_reads(synth.fasta, [ils], ["-a"], preset="sr").stdout) == sfx

	def place(text):
		return [l.split("\t")[:4] + l.split("\t")[5:9] for l in text.splitlines()]
	assert place(sfx) == place(body)
	sam = sr_default.sam()
	assert sorted(set(r.qname for r in sam.records)) == sorted(p.name for p in pairs)
	first = dict((p.name, p.r1.seq) for p in pairs)
	for r in sam.records:
		seq = r.seq if not r.is_reverse else H.revcomp(r.seq)
		assert r.segment == (1 if seq == first[r.qname] else 2), r.line[:100]
	paf = H.map_reads(synth.fasta, [sr_pairs["f1"], sr_pairs["f2"]], [], preset="sr").paf()
	assert set(p.qname for p in paf) == set("%s/%d" % (p.name, s) for p in pairs for s in (1, 2))
	pafs = H.map_reads(synth.fasta, [s1, s2], [], preset="sr").paf()
	assert set(p.qname for p in pafs) == set("%s/%d/%d" % (p.name, s, s) for p in pairs for s in (1, 2))


@pytest.mark.parametrize("F", [None, "500", "1000"])
def test_max_fragment_length(synth, sr_pairs, sr_default, F):
	"""-F NUM, the maximum fragment length ([800] with -x sr): a pair is flagged 0x2 on both mates
	exactly when its fragment is at most NUM bp. Every mate aligns at its true position and strand
	either way. For proper pairs TLEN is the 5'-to-5' length and its magnitude the fragment length;
	mate fields agree with the mate's primary (SAM spec)."""
	pairs = sr_pairs["pairs"]
	res = sr_default if F is None else H.map_reads(synth.fasta, [sr_pairs["f1"], sr_pairs["f2"]], ["-a", "-F", F], preset="sr")
	sam = res.sam()
	assert H.sam_invariant_errors(sam, synth.genome, sr_pairs["reads"]) == []
	assert H.pair_invariant_errors(sam) == []
	limit = 800 if F is None else int(F)
	prim = sam.primaries()
	n_proper = 0
	for p in pairs:
		a, b = prim[(p.name, 1)], prim[(p.name, 2)]
		for rec, rd in ((a, p.r1), (b, p.r2)):
			assert H.placement_error(rec, rd) == 0 and "S" not in rec.cigar, (rec.line[:120], rd)
		assert a.is_proper_pair == b.is_proper_pair == (p.insert <= limit), (p.insert, a.flag, b.flag)
		if a.is_proper_pair:
			n_proper += 1
			assert a.tlen == H.tlen_5prime(a, b) and abs(a.tlen) == p.insert
	assert 0 < n_proper < len(pairs)


def _cut_pairs(genome, specs, rng):
	"""Pairs cut from the genome: specs = [(name, contig1, start1, contig2, start2)], read 1 forward
	from start1 and read 2 the reverse complement of [start2, start2+150); contig None = random sequence."""
	r1, r2 = [], []
	for name, c1, s1, c2, s2 in specs:
		a = genome.fetch(c1, s1, s1 + 150) if c1 else rng.dna(150)
		b = H.revcomp(genome.fetch(c2, s2, s2 + 150)) if c2 else rng.dna(150)
		r1.append((name, a, "I" * 150))
		r2.append((name, b, "I" * 150))
	return r1, r2


@pytest.fixture(scope="module")
def repeat_pairs(synth, tmp_path_factory):
	"""For each copy of the first repeat family (three exact copies, 600+ bp), two 400 bp fragments:
	one with read 1 inside the copy and read 2 in unique sequence downstream, one with read 2 inside
	and read 1 upstream."""
	d = tmp_path_factory.mktemp("rep")
	g = synth.genome
	reps = [f for f in g.features if f.kind == "repeat" and f.info[0] == 0]
	assert len(reps) == 3 and all(f.end - f.start >= 300 for f in reps)
	specs, truth = [], {}
	for i, f in enumerate(reps):
		fs = f.end - 200  # read 1 = [end-200, end-50) inside the copy, read 2 ends at fs+400
		specs.append(("rep%da" % i, f.contig, fs, f.contig, fs + 250))
		truth["rep%da" % i] = (f.contig, fs, fs + 250, 1)
		fe = f.start + 200  # read 2 = [start+50, start+200) inside the copy
		specs.append(("rep%db" % i, f.contig, fe - 400, f.contig, fe - 150))
		truth["rep%db" % i] = (f.contig, fe - 400, fe - 150, 2)
	r1, r2 = _cut_pairs(g, specs, H.Rng(5))
	f1, f2 = H.write_fastq(str(d / "1.fq"), r1), H.write_fastq(str(d / "2.fq"), r2)
	return dict(f1=f1, f2=f2, truth=truth, reads=dict((a[0], [a[1], b[1]]) for a, b in zip(r1, r2)))


def test_pairing_places_repetitive_mates(synth, repeat_pairs):
	"""A mate inside one of three exact repeat copies cannot be placed alone. With the default strong
	pairing (`-x sr`), the pair's other mate decides: both mates at their true positions, proper pair,
	the repetitive mate with MAPQ above 0. With --pairing=no the two ends are aligned independently: the
	repetitive mate gets MAPQ 0 and no record is flagged 0x2."""
	truth = repeat_pairs["truth"]
	strong = H.map_reads(synth.fasta, [repeat_pairs["f1"], repeat_pairs["f2"]], ["-a"], preset="sr").sam()
	assert H.sam_invariant_errors(strong, synth.genome, repeat_pairs["reads"]) == [] and H.pair_invariant_errors(strong) == []
	prim = strong.primaries()
	for name, (ctg, p1, p2, rep_seg) in truth.items():
		a, b = prim[(name, 1)], prim[(name, 2)]
		assert (a.rname, a.pos0, a.strand, b.rname, b.pos0, b.strand) == (ctg, p1, "+", ctg, p2, "-"), (a.line[:80], b.line[:80])
		assert a.is_proper_pair and b.is_proper_pair and a.mapq > 0 and b.mapq > 0
	indep = H.map_reads(synth.fasta, [repeat_pairs["f1"], repeat_pairs["f2"]], ["-a", "--pairing=no"], preset="sr").sam()
	assert H.pair_invariant_errors(indep) == []
	prim = indep.primaries()
	for name, (ctg, p1, p2, rep_seg) in truth.items():
		assert prim[(name, rep_seg)].mapq == 0 and prim[(name, 3 - rep_seg)].mapq > 0
	assert not any(r.is_proper_pair for r in indep.records)


def test_pairs_of_repetitive_mates(tmp_path):
	"""Both mates inside exact repeats. On one 40 kbp contig, units R and Q (400 bp) have three copies
	each and only R1 and Q1 lie within -F of each other; units S and T have two copies each, and both
	S1-T1 and S2-T2 are 500 bp apart. Alone every mate has MAPQ 0 (--pairing=no). Paired with -x sr:
	a pair from R1 and Q1 has one consistent placement, so both mates are proper, at their true
	positions, with MAPQ above 0, the same for both mates (each mate is equally repetitive and they are
	placed together); a pair from S and T has two equally good placements, so both mates are proper at
	one of them with MAPQ at most 3 (SAM spec: -10 log10 of a 1/2 chance of being wrong)."""
	rng = H.Rng(901)
	ctg = list(rng.dna(40000))
	units = [rng.dna(400) for _ in range(4)]
	for u, starts in zip(units, ((2000, 12000, 33000), (2500, 17000, 37000), (6000, 24000), (6500, 24500))):
		for st in starts:
			ctg[st:st + 400] = list(u)
	g = "".join(ctg)
	fa = H.write_fasta(str(tmp_path / "g.fa"), [("c", g)])
	specs = [("u0", 2150, 2600), ("u1", 2050, 2700), ("u2", 2200, 2550), ("a0", 6150, 6600), ("a1", 6050, 6700)]
	f1 = H.write_fasta(str(tmp_path / "1.fa"), [(n, g[a:a + 150]) for n, a, b in specs])
	f2 = H.write_fasta(str(tmp_path / "2.fa"), [(n, H.revcomp(g[b:b + 150])) for n, a, b in specs])
	indep = H.map_reads(fa, [f1, f2], ["-a", "--pairing=no"], preset="sr").sam()
	assert len(indep.records) == 10 and all(r.mapq == 0 and not r.is_unmapped for r in indep.records)
	sam = H.map_reads(fa, [f1, f2], ["-a"], preset="sr").sam()
	assert H.pair_invariant_errors(sam) == []
	prim = sam.primaries()
	for name, a, b in specs:
		r1, r2 = prim[(name, 1)], prim[(name, 2)]
		assert r1.is_proper_pair and r2.is_proper_pair and r1.cigar == r2.cigar == "150M", (r1.line[:80], r2.line[:80])
		shift = r1.pos0 - a
		if name.startswith("u"):
			assert (r1.pos0, r2.pos0) == (a, b) and r1.mapq == r2.mapq > 0, (name, r1.mapq, r2.mapq)
		else:
			assert shift in (0, 18000) and r2.pos0 - b == shift and r1.mapq <= 3 and r2.mapq <= 3, (name, r1.pos0, r1.mapq, r2.mapq)


def test_pairing_no_aligns_ends_independently(synth, sr_pairs, sr_default):
	"""--pairing=no (minimap2.1): the two ends are aligned independently with no 'properly paired' set;
	--no-pairing is its deprecated spelling and gives identical output. Records keep 0x1, 0x40/0x80 and
	mate fields, and each mate (all from unique sequence) gets the alignment it gets when its file is
	mapped single-end."""
	f1, f2 = sr_pairs["f1"], sr_pairs["f2"]
	a = H.map_reads(synth.fasta, [f1, f2], ["-a", "--pairing=no"], preset="sr")
	b = H.map_reads(synth.fasta, [f1, f2], ["-a", "--no-pairing"], preset="sr")
	assert H.sam_body(a.stdout) == H.sam_body(b.stdout)
	sam = a.sam()
	assert H.pair_invariant_errors(sam) == []
	assert not any(r.is_proper_pair for r in sam.records) and all(r.is_paired for r in sam.records)
	prim = sam.primaries()

	def key(r):
		return (r.rname, r.pos, r.strand, r.cigar, r.mapq, sorted(r.tags.items()))
	for seg, f in ((1, f1), (2, f2)):
		single = H.map_reads(synth.fasta, f, ["-a"], preset="sr").sam()
		for (q, _), r in single.primaries().items():
			assert key(r) == key(prim[(q, seg)]), (r.line[:100], prim[(q, seg)].line[:100])


def test_unmapped_and_split_mates(synth, tmp_path):
	"""SAM mate fields for incomplete pairs: a mapped mate with an unmapped partner carries 0x8 and the
	partner, flag 0x4, takes the mapped mate's RNAME and POS (TLEN 0); mates on different contigs name
	each other's contig in RNEXT, with TLEN 0 and no 0x2; a pair with both ends unmapped has '*' and 0."""
	g = synth.genome
	specs = [("unm", "chr1", 5000, None, 0), ("xctg", "chr1", 6000, "chr2", 6000), ("both", None, 0, None, 0),
		("ok", "chr1", 7000, "chr1", 7200)]
	r1, r2 = _cut_pairs(g, specs, H.Rng(17))
	f1, f2 = H.write_fastq(str(tmp_path / "1.fq"), r1), H.write_fastq(str(tmp_path / "2.fq"), r2)
	sam = H.map_reads(synth.fasta, [f1, f2], ["-a"], preset="sr").sam()
	assert H.pair_invariant_errors(sam) == []
	assert H.sam_invariant_errors(sam, g, dict((a[0], [a[1], b[1]]) for a, b in zip(r1, r2))) == []
	p = sam.primaries()
	m, u = p[("unm", 1)], p[("unm", 2)]
	assert (m.flag, m.rname, m.pos, m.rnext, m.pnext, m.tlen) == (0x1 | 0x8 | 0x40, "chr1", 5001, "=", 5001, 0)
	assert (u.flag, u.rname, u.pos, u.cigar, u.rnext, u.pnext, u.tlen) == (0x1 | 0x4 | 0x80, "chr1", 5001, "*", "=", 5001, 0)
	a, b = p[("xctg", 1)], p[("xctg", 2)]
	assert (a.flag, a.rname, a.pos, a.rnext, a.pnext, a.tlen) == (0x1 | 0x20 | 0x40, "chr1", 6001, "chr2", 6001, 0)
	assert (b.flag, b.rname, b.pos, b.rnext, b.pnext, b.tlen) == (0x1 | 0x10 | 0x80, "chr2", 6001, "chr1", 6001, 0)
	for seg, flag in ((1, 0x1 | 0x4 | 0x8 | 0x40), (2, 0x1 | 0x4 | 0x8 | 0x80)):
		r = p[("both", seg)]
		assert (r.flag, r.rname, r.pos, r.rnext, r.pnext, r.tlen) == (flag, "*", 0, "*", 0, 0)
	a, b = p[("ok", 1)], p[("ok", 2)]
	assert (a.flag, b.flag, a.tlen, b.tlen) == (0x1 | 0x2 | 0x20 | 0x40, 0x1 | 0x2 | 0x10 | 0x80, 350, -350)


def test_three_segment_fragments(synth, tmp_path):
	"""Three query files with --frag=yes make three-segment templates (FI:i, SAMtags); here each
	template is three reads of one strand in fragment order. Per the SAM
	spec: 0x40 alone on the first segment and 0x80 alone on the last (the middle segment's bits are
	known-bugs/d-middle-segment-flags.md); RNEXT/PNEXT and 0x20 give the next segment's primary, and
	the next of the last is the first; TLEN is measured to that next segment. PAF names the segments
	/1, /2, /3."""
	g = synth.genome
	files = [[], [], []]
	truth = {}
	cands = [("chr1" if i % 2 == 0 else "chr2", 1000 + i * 1000) for i in range(80)]
	cands = [x for x in cands if x[1] + 1000 < g.length(x[0]) and g.is_unique(x[0], x[1], x[1] + 800, 100)][:20]
	assert len(cands) == 20
	for i, (c, st) in enumerate(cands):
		# all segments on one strand (the fragment mode pairs segments as given, forward-forward)
		poss, strands = ([st, st + 300, st + 600], ["+"] * 3) if i % 3 else ([st + 600, st + 300, st], ["-"] * 3)
		name = "t%d" % i
		for k in range(3):
			s = g.fetch(c, poss[k], poss[k] + 200, strands[k])
			files[k].append((name, s, "I" * 200))
		truth[name] = (c, poss, strands)
	fs = [H.write_fastq(str(tmp_path / ("s%d.fq" % k)), files[k]) for k in range(3)]
	res = H.map_reads(synth.fasta, fs, ["-a", "--frag=yes"])
	sam = res.sam()
	recs = dict(((r.qname, r.tags["FI"]), r) for r in sam.records if r.is_primary)
	assert len(recs) == 60 and len(sam.records) == 60
	for name, (c, poss, strands) in truth.items():
		for k in range(3):
			r, nxt = recs[(name, k)], recs[(name, (k + 1) % 3)]
			assert (r.rname, r.pos0, r.strand) == (c, poss[k], strands[k])
			assert r.flag & 0x1 and not r.flag & 0x2
			if k != 1:
				assert r.flag & 0xC0 == (0x40 if k == 0 else 0x80)
			assert (r.rnext, r.pnext, bool(r.flag & 0x20)) == ("=", nxt.pos, nxt.is_reverse)
			assert r.tlen == H.tlen_5prime(r, nxt)
	paf = H.map_reads(synth.fasta, fs, ["--frag=yes"]).paf()
	assert set(p.qname for p in paf) == set("%s/%d" % (n, k) for n in truth for k in (1, 2, 3))


def test_frag_mode_groups_reads_by_name(synth, sr_pairs, tmp_path):
	"""In fragment mode, adjacent reads with the same name in one file are one fragment (README; single-
	and paired-end reads can be mixed). In an interleaved file with an unpaired read between pairs, and
	two same-name reads that are not adjacent: `-x sr` pairs only the adjacent mates; `--frag=yes` does
	the same with the default preset; `-x sr --frag=no` treats every read as single-end."""
	pairs = sr_pairs["pairs"][:20]
	recs = []
	for i, p in enumerate(pairs):
		recs += [p.r1, p.r2]
		if i == 9:
			recs.append(("lone", synth.genome.fetch("chr1", 20000, 20150), "I" * 150))
	recs.append(("split", synth.genome.fetch("chr1", 30000, 30150), "I" * 150))
	recs.append(("mid", synth.genome.fetch("chr1", 40000, 40150), "I" * 150))
	recs.append(("split", synth.genome.fetch("chr2", 30000, 30150), "I" * 150))
	fq = H.write_fastq(str(tmp_path / "mix.fq"), recs)
	for args, preset in ((["-a"], "sr"), (["-a", "--frag=yes"], None)):
		sam = H.map_reads(synth.fasta, fq, args, preset=preset).sam()
		assert H.pair_invariant_errors(sam) == []
		for r in sam.records:
			assert r.is_paired == (r.qname not in ("lone", "split", "mid")), r.line[:80]
		assert sorted(r.segment for r in sam.records if r.qname == "split") == [0, 0]
	sam = H.map_reads(synth.fasta, fq, ["-a", "--frag=no"], preset="sr").sam()
	assert len(sam.records) == len(recs) and not any(r.flag & 0xC3 for r in sam.records)


def test_pairs_with_a_split_index(synth, sr_pairs, sr_default, tmp_path):
	"""A multi-part index (-I smaller than the genome) with --split-prefix merges the parts' hits and
	pairs them afterwards (map.c merge_hits). For these pairs from unique sequence every record keeps the
	placement it gets from a single-part index: FLAG (0x2 included), RNAME, POS, CIGAR, RNEXT, PNEXT and
	TLEN."""
	res = H.map_reads(synth.fasta, [sr_pairs["f1"], sr_pairs["f2"]], ["-a", "-I", "60k", "--split-prefix", str(tmp_path / "sp")],
		preset="sr")
	assert H.pair_invariant_errors(res.sam()) == []

	def place(text):
		return [l.split("\t")[:4] + l.split("\t")[5:9] for l in H.sam_body(text).splitlines()]
	assert place(res.stdout) == place(sr_default.stdout)


@pytest.mark.thread_invariance
def test_frag_mode_minibatches_keep_mates_together(synth, sr_pairs, sr_default):
	"""Batches of -K bases end after a pair, not between mates (bseq.c in fragment mode): an interleaved
	file mapped with -K 1k gives the same records as two files with the preset's -K."""
	il = H.write_interleaved(str(sr_pairs["dir"] / "il_k.fq"), sr_pairs["pairs"])
	res = H.map_reads(synth.fasta, il, ["-a", "-K", "1k"], preset="sr")
	assert H.sam_body(res.stdout) == H.sam_body(sr_default.stdout)


def test_input_checks(synth, sr_pairs, tmp_path):
	"""In the sr mode more than two query files is an error (exit 1, message); two files with different
	numbers of records give a warning and only the complete pairs."""
	f1, f2 = sr_pairs["f1"], sr_pairs["f2"]
	r = H.map_reads(synth.fasta, [f1, f2, f1], ["-a"], preset="sr", check=False)
	assert r.returncode == 1 and r.messages and r.stdout == ""
	short = H.write_fastq(str(tmp_path / "short.fq"), [p.r1 for p in sr_pairs["pairs"][:30]])
	r = H.map_reads(synth.fasta, [short, f2], ["-a"], preset="sr")
	assert r.messages
	sam = r.sam()
	assert sorted(set(x.qname for x in sam.records)) == sorted(p.name for p in sr_pairs["pairs"][:30])
	assert len(sam.records) == 60 and H.pair_invariant_errors(sam) == []


def test_option_values(synth, sr_pairs, sr_default):
	"""Warnings for unknown values name the option and fall back to the documented meaning:
	--pairing=bogus to strong (the default), --frag=maybe leaves the mode unchanged, --sr=bogus to dna.
	--sr equals --sr=dna, --sr=no switches the short-read mode off again, and the undocumented preset
	alias `short` equals sr."""
	f = [sr_pairs["f1"], sr_pairs["f2"]]
	body = H.sam_body(sr_default.stdout)
	for opt, name in (("--pairing=bogus", "pairing"), ("--frag=maybe", "frag"), ("--sr=bogus", "sr")):
		r = H.map_reads(synth.fasta, f, ["-a", opt], preset="sr")
		assert any(name in m for m in r.messages) and H.sam_body(r.stdout) == body, opt
	assert H.sam_body(H.map_reads(synth.fasta, f, ["-a", "--pairing=strong"], preset="sr").stdout) == body
	assert H.strip_pg(H.map_reads(synth.fasta, f, ["-a"], preset="short").stdout) == H.strip_pg(sr_default.stdout)
	reads = H.simulate_reads(synth.genome, 100, seed=113, sub_rate=0.02, indel_rate=0.005)
	one = H.write_fastq(str(sr_pairs["dir"] / "indel.fq"), reads)

	def b(*args):
		return H.sam_body(H.map_reads(synth.fasta, one, ["-a"] + list(args)).stdout)
	plain, sr = b(), b("--sr")
	assert sr != plain and b("--sr=dna") == sr and b("--sr", "--sr=no") == plain


def test_splice_sr_single_end_equals_its_definition(synth_tx, tmp_path):
	"""minimap2.1 defines splice:sr as `-xsplice:hq --frag=yes -m25 -s40 -2K100m --heap-sort=yes
	--pairing=weak --sr=rna --min-dp-len=20 --secondary=no`; for single-end reads the two give the
	same output (this exercises --sr=rna and --pairing=weak as options)."""
	reads = H.simulate_spliced_reads(synth_tx.genome, 100, seed=121, profile="short", length=100)
	fq = H.write_fastq(str(tmp_path / "se.fq"), reads)
	a = H.map_reads(synth_tx.fasta, fq, ["-a"], preset="splice:sr")
	b = H.map_reads(synth_tx.fasta, fq, ["-a", "-x", "splice:hq", "--frag=yes", "-m25", "-s40", "-2K100m", "--heap-sort=yes",
		"--pairing=weak", "--sr=rna", "--min-dp-len=20", "--secondary=no"])
	assert H.sam_body(a.stdout) == H.sam_body(b.stdout)


@pytest.mark.pinned
def test_sr_pairs_pinned(synth, tmp_path):
	"""Pins the current output: `-ax sr --cs` (README) on 150 simulated pairs with indels from the
	whole standard genome, repeats and N runs included, SAM without @PG, byte for byte. (No --MD:
	its final 0 is known-bugs/md-trailing-zero.md.)"""
	pairs = H.simulate_pairs(synth.genome, 150, seed=112, insert_mean=350, insert_sd=60, qual=None, n_rate=0.002)
	f1, f2 = H.write_pairs(str(tmp_path / "r1.fa"), str(tmp_path / "r2.fa"), pairs)
	res = H.map_reads(synth.fasta, [f1, f2], ["-a", "--cs"], preset="sr")
	H.check_pinned("c/sr_pairs.sam.gz", H.strip_pg(res.stdout))
