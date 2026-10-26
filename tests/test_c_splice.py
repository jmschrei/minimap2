"""Area C: long-read spliced alignment. The splice presets and their aliases,
`--splice`, `-u`, `-G`, `-C`, `-J`, `--splice-flank`, annotated junctions
(`--junc-bed`, `--junc-bonus`, mm_idx_bed_read/mm_idx_bed_junc) and splice
scores (`--spsc`, `--spsc-scale`, `--spsc0`/`--junc-pen`, mm_idx_spsc_*).

Truth comes from the simulated gene models of the standard transcript genome
and from hand-built loci (helpers_c) whose correct alignment is known by
construction: two equally good junctions that annotation or a splice score
must decide between, a 6 kb intron, and a non-canonical intron. On those loci
the expected AS:i follows from the preset's scoring (`-x splice`: -A1, a
k-long gap costs min{2+k, 32+0*k}, so an intron between two exact 300 bp exons
scores 600 - 32 = 568)."""

import pytest

import helpers_c as C
import minimap2_helpers as H


def _recs(res):
	return [r for r in res.sam().records if not r.is_unmapped]


def _by_name(res):
	return dict((r.qname, r) for r in _recs(res) if r.is_primary)


@pytest.fixture(scope="module")
def hifi_tx(synth_tx, tmp_path_factory):
	"""60 HiFi-like cDNA reads, sense and antisense, from the transcript genome."""
	d = tmp_path_factory.mktemp("hifi_tx")
	reads = H.simulate_spliced_reads(synth_tx.genome, 60, seed=201, profile="hifi")
	return dict(reads=reads, fq=H.write_fastq(str(d / "tx.fq"), reads), truth=dict((r.name, r) for r in reads))


def test_splice_ont_like_cdna(synth_tx, tmp_path):
	"""README `-ax splice ref.fa nanopore-cdna.fa` on ONT-like cDNA (about 7% error), both orientations:
	every read on its true contig and genome strand, at least 75% with exactly the true introns; on
	those, ts:A gives the read's orientation to the transcript (cookbook.md). The PAF of `-cx splice`
	has the same primaries, cg:Z with N for introns and ts:A on spliced lines; tags agree with the genome."""
	g = synth_tx.genome
	reads = H.simulate_spliced_reads(g, 50, seed=202, profile="ont", min_len=500)
	fq = H.write_fastq(str(tmp_path / "ont.fq"), reads)
	res = H.map_reads(synth_tx.fasta, fq, ["-a", "--cs"], preset="splice")
	sam = res.sam()
	assert H.sam_invariant_errors(sam, g, dict((r.name, [r.seq]) for r in reads)) == []
	prim = sam.primaries()
	exact = 0
	for r in reads:
		rec = prim[(r.name, 0)]
		assert not rec.is_unmapped and rec.rname == r.contig and rec.strand == r.strand, (rec.line[:120], r)
		if rec.introns == r.introns and r.introns:
			exact += 1
			assert rec.tags["ts"] == r.ts_strand
	assert exact >= 0.75 * sum(1 for r in reads if r.introns)
	paf = H.map_reads(synth_tx.fasta, fq, ["-c", "--cs"], preset="splice").paf()
	assert H.paf_invariant_errors(paf, g, dict((r.name, r.seq) for r in reads)) == []
	for p in paf:
		if p.is_primary:
			rec = prim[(p.qname, 0)]
			assert (p.tname, p.tstart, p.strand, p.cigar) == (rec.rname, rec.pos0, rec.strand, H.strip_clips(rec.cigar))
			assert p.tags.get("ts") == rec.tags.get("ts")


def test_u_option_chooses_strand(synth_tx, hifi_tx):
	"""-u: f = transcript strand only, r = the reverse, b = both (README: two rounds of alignment,
	keeping the better one and writing it to ts). On HiFi-like spliced reads:
	-uf writes ts:A:+ and -ur ts:A:- on every alignment; a read on the assumed strand gets the same
	CIGAR and AS as with -ub, where ts:A gives its true orientation; a spliced read on the other strand
	scores lower, because its GT-AG signals are not matched."""
	reads, fq = hifi_tx["truth"], hifi_tx["fq"]
	both = _by_name(H.map_reads(synth_tx.fasta, fq, ["-a", "-ub"], preset="splice:hq"))
	assert H.strip_pg(H.map_reads(synth_tx.fasta, fq, ["-a"], preset="splice:hq").stdout) == \
		H.strip_pg(H.map_reads(synth_tx.fasta, fq, ["-a", "-ub"], preset="splice:hq").stdout)
	n_lower = 0
	for u, ts in (("f", "+"), ("r", "-")):
		res = H.map_reads(synth_tx.fasta, fq, ["-a", "-u" + u], preset="splice:hq")
		assert all(r.tags.get("ts") == ts for r in _recs(res))
		one = _by_name(res)
		for name, rd in reads.items():
			b, o = both[name], one[name]
			if rd.ts_strand == ts:
				assert (o.pos, o.cigar, o.tags["AS"]) == (b.pos, b.cigar, b.tags["AS"]), (u, name)
				if b.introns:
					assert b.tags["ts"] == ts
			elif rd.introns and b.introns == rd.introns:
				assert o.tags["AS"] < b.tags["AS"], (u, name)
				n_lower += 1
	assert n_lower >= 20


def test_u_n_ignores_splice_signals(tmp_path):
	"""-un: no attempt to match GT-AG. Across an intron that starts with CA and ends with TC the read
	then aligns with exactly the true intron and AS = 800 matches - 32 (one long gap), while the default
	-ub, which scores splice signals, does not take that path."""
	loc = C.noncanonical_locus(401)
	loc.write(str(tmp_path / "g.fa"))
	fa = H.write_fasta(str(tmp_path / "r.fa"), [("fwd", loc.read), ("rev", H.revcomp(loc.read))])
	cig = "400M%dN400M" % (loc.introns[0][1] - loc.introns[0][0])
	recs = _recs(H.map_reads(loc.fasta, fa, ["-a", "-un"], preset="splice"))
	assert [(r.qname, r.pos0, r.cigar, r.tags["AS"]) for r in recs] == [("fwd", loc.exons[0][0], cig, 768), ("rev", loc.exons[0][0], cig, 768)]
	assert all(r.cigar != cig for r in _recs(H.map_reads(loc.fasta, fa, ["-a"], preset="splice")))


def test_u_invalid_is_an_error(synth_tx, hifi_tx):
	r = H.map_reads(synth_tx.fasta, hifi_tx["fq"], ["-a", "-u", "x"], preset="splice", check=False)
	assert r.returncode == 1 and r.messages and r.stdout == ""


def test_noncanonical_cost_C_and_J(tmp_path):
	"""-C INT: cost of a non-canonical splice site, effective with -J0; -J1 (the default) ignores -C.
	On the CA..TC intron: -J0 -C0 aligns the true intron with AS 768 (no site cost); -J0 -C3 keeps the
	same CIGAR at 3 less per site, 762. Under -J1, -C0, -C20 and no -C give identical output, and -J1
	is the default."""
	loc = C.noncanonical_locus(401)
	loc.write(str(tmp_path / "g.fa"))
	fa = H.write_fasta(str(tmp_path / "r.fa"), [("fwd", loc.read), ("rev", H.revcomp(loc.read))])
	cig = "400M%dN400M" % (loc.introns[0][1] - loc.introns[0][0])

	def run(*args):
		return H.map_reads(loc.fasta, fa, ["-a"] + list(args), preset="splice")
	for c, as_ in (("0", 768), ("3", 762)):
		recs = _recs(run("-J0", "-C" + c))
		assert [(r.pos0, r.cigar, r.tags["AS"]) for r in recs] == [(loc.exons[0][0], cig, as_)] * 2, c
	j1 = H.strip_pg(run().stdout)
	assert H.strip_pg(run("-J1").stdout) == j1
	assert H.strip_pg(run("-C0").stdout) == j1 and H.strip_pg(run("-C20").stdout) == j1
	assert H.strip_pg(run("-J0").stdout) != j1


@pytest.fixture(scope="module")
def ties(tmp_path_factory):
	"""Tie loci on both gene strands (both candidate junctions fully flanked) plus their read files."""
	d = tmp_path_factory.mktemp("ties")
	out = {}
	for strand in "+-":
		loc = C.tie_locus(211, strand=strand)
		loc.write(str(d / ("tie%s.fa" % ("p" if strand == "+" else "m"))))
		loc.reads = H.write_fasta(str(d / ("r%s.fa" % ("p" if strand == "+" else "m"))), [("sense", loc.read if strand == "+" else H.revcomp(loc.read)),
			("anti", H.revcomp(loc.read) if strand == "+" else loc.read)])
		out[strand] = loc
	out["dir"] = d
	return out


def _junctions(res):
	"""{qname: (intron, AS, ts)} of the primary records."""
	return dict((r.qname, (tuple(r.introns), r.tags["AS"], r.tags.get("ts"))) for r in _recs(res) if r.is_primary)


def _tie_base(loc, *args):
	res = H.map_reads(loc.fasta, loc.reads, ["-a"] + list(args), preset="splice")
	j = _junctions(res)
	assert set(j) == {"sense", "anti"}
	assert j["sense"][2] == "+" and j["anti"][2] == "-"
	assert set(v[0] for v in j.values()) <= {(loc.A,), (loc.B,)}
	assert set(v[1] for v in j.values()) == {568}
	return res, j


@pytest.mark.parametrize("strand", ["+", "-"])
def test_junc_bed_decides_a_tie(ties, strand):
	"""--junc-bed: README, an aligned junction matching an annotated one gets a bonus (--junc-bonus,
	9 by default, per splice site). Two GT-AG introns that align the read equally well (AS 568): with
	either one annotated, both reads (sense and antisense) take the annotated intron with AS 568 + 2x9,
	or 568 + 2x5 with --junc-bonus 5. As BED12 (a gene) or as a 6-column junction line the result is
	the same; a junction line without a strand ('.' or no sixth column; minimap2.1: the strand column
	is required) and --junc-bonus 0 give the output of no annotation; --junc-pen, effective with --spsc
	only, changes nothing; -J0 also honours the annotation."""
	loc = ties[strand]
	d = ties["dir"]
	base, _ = _tie_base(loc)
	for which in ("A", "B"):
		iv = getattr(loc, which)
		ex = [(loc.exons[0][0], iv[0]), (iv[1], loc.exons[1][1])]
		g = H.Gene("g" + which, "chr1", strand, ex)
		bed12 = H.write_text(str(d / ("%s%s.bed12" % (which, strand))), H.genes_bed12([g]))
		bed6 = H.write_text(str(d / ("%s%s.bed6" % (which, strand))), H.junctions_bed([g]))
		res = H.map_reads(loc.fasta, loc.reads, ["-a", "--junc-bed", bed12], preset="splice")
		assert _junctions(res) == {"sense": ((iv,), 586, "+"), "anti": ((iv,), 586, "-")}
		assert H.sam_body(H.map_reads(loc.fasta, loc.reads, ["-a", "--junc-bed", bed6], preset="splice").stdout) == H.sam_body(res.stdout)
		assert H.sam_body(H.map_reads(loc.fasta, loc.reads, ["-a", "--junc-bed", bed6, "--junc-pen", "20"], preset="splice").stdout) == H.sam_body(res.stdout)
		j5 = _junctions(H.map_reads(loc.fasta, loc.reads, ["-a", "--junc-bed", bed6, "--junc-bonus", "5"], preset="splice"))
		assert j5 == {"sense": ((iv,), 578, "+"), "anti": ((iv,), 578, "-")}
		j0 = H.map_reads(loc.fasta, loc.reads, ["-a", "--junc-bed", bed6, "--junc-bonus", "0"], preset="splice")
		assert H.sam_body(j0.stdout) == H.sam_body(base.stdout)
		for name, line in (("dot", "chr1\t%d\t%d\tj\t0\t.\n"), ("five", "chr1\t%d\t%d\tj\t0\n")):
			nostrand = H.write_text(str(d / ("%s%s%s.bed" % (which, strand, name))), line % iv)
			res = H.map_reads(loc.fasta, loc.reads, ["-a", "--junc-bed", nostrand], preset="splice")
			assert H.sam_body(res.stdout) == H.sam_body(base.stdout), name
		jj = _junctions(H.map_reads(loc.fasta, loc.reads, ["-a", "-J0", "--junc-bed", bed6], preset="splice"))
		assert jj == {"sense": ((iv,), 586, "+"), "anti": ((iv,), 586, "-")}


def _spsc_file(path, rows):
	return H.write_text(path, "".join("%s\t%d\t%s\t%s\t%d\n" % r for r in rows))


@pytest.mark.parametrize("strand", ["+", "-"])
def test_spsc_scores_decide_a_tie(ties, strand):
	"""--spsc FILE (minimap2.1): contig, offset (bases before the junction), strand, D or A, score; a
	positive score prefers the junction, a negative one disfavours it, scores are scaled by
	--spsc-scale [0.7] and rounded, and --spsc0 [5] (alias --junc-pen) is the penalty for positions not
	in FILE. On the tie: the scored intron wins with AS 568 + 2 x round(score x scale); scored -10
	(-7 per site) it loses to the other while that pays less (2 x --spsc0 of 0 or 5) and wins at
	--spsc0 9; a FILE without this contig costs 2 x --spsc0 per intron;
	--junc-bonus, effective with --junc-bed only, changes nothing. For a '-' gene the donor is the
	intron's right end on the genome."""
	loc = ties[strand]
	d = str(ties["dir"])
	base, _ = _tie_base(loc)
	for which, other in (("A", "B"), ("B", "A")):
		iv, iv2 = getattr(loc, which), getattr(loc, other)
		donor, acceptor = (iv[0], iv[1]) if strand == "+" else (iv[1], iv[0])
		for score, scale, gain in ((10, None, 14), (10, "1", 20), (10, "0.5", 10), (6, "0.7", 8)):
			f = _spsc_file(d + "/sp.tsv", [("chr1", donor, strand, "D", score), ("chr1", acceptor, strand, "A", score)])
			args = ["-a", "--spsc", f] + (["--spsc-scale", scale] if scale else [])
			j = _junctions(H.map_reads(loc.fasta, loc.reads, args, preset="splice"))
			assert j == {"sense": ((iv,), 568 + gain, "+"), "anti": ((iv,), 568 + gain, "-")}, (which, score, scale)
		f = _spsc_file(d + "/neg.tsv", [("chr1", donor, strand, "D", -10), ("chr1", acceptor, strand, "A", -10)])
		for pen, win, as_ in (("0", iv2, 568), ("5", iv2, 558), ("9", iv, 554)):  # -10 scales to -7 per site
			outs = []
			for opt in ("--spsc0", "--junc-pen"):
				res = H.map_reads(loc.fasta, loc.reads, ["-a", "--spsc", f, opt, pen], preset="splice")
				assert _junctions(res) == {"sense": ((win,), as_, "+"), "anti": ((win,), as_, "-")}, (which, opt, pen)
				outs.append(H.sam_body(res.stdout))
			assert outs[0] == outs[1]
	other = _spsc_file(d + "/other.tsv", [("chr2", 100, "+", "D", 3)])
	for pen in ("0", "5", "9"):
		res = H.map_reads(loc.fasta, loc.reads, ["-a", "--spsc", other, "--spsc0", pen], preset="splice")
		j = _junctions(res)
		assert set(v[1] for v in j.values()) == {568 - 2 * int(pen)}
		assert [x[0] for x in j.values()] == [x[0] for x in _junctions(base).values()]
	a = H.map_reads(loc.fasta, loc.reads, ["-a", "--spsc", other], preset="splice")
	b = H.map_reads(loc.fasta, loc.reads, ["-a", "--spsc", other, "--junc-bonus", "20"], preset="splice")
	assert H.sam_body(a.stdout) == H.sam_body(b.stdout)


def test_spsc_scores_are_capped(ties):
	"""Scores beyond what the scoring system allows are capped at mm_max_spsc_bonus(): with
	--spsc-scale 1 a donor and an acceptor scored 100 add 2 x 30 to AS under -x splice (-O2,32) and
	2 x 18 under -x splice:hq (-O6,24), whose unscored AS is 600 - 24 = 576."""
	loc = ties["+"]
	f = _spsc_file(str(ties["dir"]) + "/cap.tsv", [("chr1", loc.A[0], "+", "D", 100), ("chr1", loc.A[1], "+", "A", 100)])
	for preset, as_ in (("splice", 568 + 60), ("splice:hq", 576 + 36)):
		j = _junctions(H.map_reads(loc.fasta, loc.reads, ["-a", "--spsc", f, "--spsc-scale", "1"], preset=preset))
		assert j == {"sense": ((loc.A,), as_, "+"), "anti": ((loc.A,), as_, "-")}, preset


def test_splice_flank_under_J0(tmp_path):
	"""--splice-flank (README): -x splice prefers GT[A/G]..[C/T]AG over GT[C/T]..[A/G]AG;
	--splice-flank=no models GT..AG only. Two GT-AG introns align the read equally well; in locus A
	only intron A has both conserved flanking bases, in locus B only intron B. With the default
	--splice-flank=yes, -J0 and -J1 both take the flanked intron. With -J0 --splice-flank=no the
	flanks are ignored, so loci A and B score like the fully flanked tie locus and take the intron
	it takes. (Under the default -J1, --splice-flank=no has no effect; see known-bugs.)"""
	locs = {}
	for name, fa_, fb in (("A", True, False), ("B", False, True), ("tie", True, True)):
		loc = C.tie_locus(212, flank_a=fa_, flank_b=fb)
		loc.write(str(tmp_path / ("%s.fa" % name)))
		loc.reads = H.write_fasta(str(tmp_path / ("%s.r.fa" % name)), [("sense", loc.read), ("anti", H.revcomp(loc.read))])
		locs[name] = loc

	def intron(loc, *args):
		j = _junctions(H.map_reads(loc.fasta, loc.reads, ["-a"] + list(args), preset="splice"))
		assert len(set(v[0] for v in j.values())) == 1
		return j["sense"][0][0]
	for model in ([], ["-J0"]):
		assert intron(locs["A"], *model) == locs["A"].A
		assert intron(locs["B"], *model) == locs["B"].B
	tie = intron(locs["tie"], "-J0", "--splice-flank=no")
	side = "A" if tie == locs["tie"].A else "B"
	assert intron(locs["A"], "-J0", "--splice-flank=no") == getattr(locs["A"], side)
	assert intron(locs["B"], "-J0", "--splice-flank=no") == getattr(locs["B"], side)


def test_max_intron_len(tmp_path):
	"""-G NUM: maximum intron length, effective in splice mode only. A read across a 6 kb GT-AG intron
	aligns as 400M6000N400M with the preset's -G200k, -G10k or -G7000 (and with `--splice -ub -G10k`),
	and is split into a primary and a supplementary alignment, one per exon, with no N, under -G5999
	or -G2k. Outside splice mode -G changes nothing (map-ont and sr)."""
	loc = C.long_intron_locus(402)
	loc.write(str(tmp_path / "g.fa"))
	fa = H.write_fasta(str(tmp_path / "r.fa"), [("fwd", loc.read), ("rev", H.revcomp(loc.read))])
	(e1s, e1e), (e2s, e2e) = loc.exons
	whole = (e1s, "400M6000N400M")
	for args in (["-x", "splice"], ["-x", "splice", "-G10k"], ["-x", "splice", "-G", "7000"], ["--splice", "-ub", "-G10k"]):
		recs = _recs(H.map_reads(loc.fasta, fa, ["-a"] + args))
		assert [(r.qname, r.pos0, r.cigar) for r in recs] == [("fwd",) + whole, ("rev",) + whole], args
	for G in ("5999", "2k"):
		recs = _recs(H.map_reads(loc.fasta, fa, ["-a", "-G", G], preset="splice"))
		for q in ("fwd", "rev"):
			mine = [r for r in recs if r.qname == q]
			assert len(mine) == 2 and sum(r.is_supplementary for r in mine) == 1, (G, q)
			assert all(not r.introns for r in mine)
			spans = sorted((r.pos0, r.ref_end0) for r in mine)
			assert spans[0][0] == e1s and e1e <= spans[0][1] < e2s and e1e < spans[1][0] <= e2s and spans[1][1] == e2e, spans
	for preset, g in (("map-ont", "1k"), ("sr", "100")):
		a = H.map_reads(loc.fasta, fa, ["-a"], preset=preset)
		b = H.map_reads(loc.fasta, fa, ["-a", "-G", g], preset=preset)
		assert H.sam_body(a.stdout) == H.sam_body(b.stdout), preset


def test_splice_option_enables_spliced_alignment(synth_tx, hifi_tx):
	"""--splice enables the splice mode (minimap2.1) without the preset's other settings: with
	`--splice -ub` and the default scoring, most spliced reads get N operations and every N is a true
	intron, while the default map-ont mode never writes an N operation."""
	g = synth_tx.genome
	true = set((gg.contig, s, e) for gg in g.genes for s, e in gg.introns)
	recs = _recs(H.map_reads(synth_tx.fasta, hifi_tx["fq"], ["-a", "--splice", "-ub"]))
	ivs = [(r.rname,) + iv for r in recs for iv in r.introns]
	assert ivs and all(x in true for x in ivs)
	assert len(set(r.qname for r in recs if r.introns)) >= 0.6 * sum(1 for rd in hifi_tx["reads"] if rd.introns)
	assert all("N" not in r.cigar for r in _recs(H.map_reads(synth_tx.fasta, hifi_tx["fq"], ["-a"])))


def test_splice_preset_equivalences(synth_tx, hifi_tx):
	"""minimap2.1 defines `-x splice` as `-k15 -w5 --splice -g2k -G200k -A1 -B2 -O2,32 -E1,0 -C9 -z200
	-ub --junc-bonus=9 --cap-sw-mem=0 --splice-flank=yes` and `-x splice:hq` as `-xsplice -C5 -O6,24
	-B4`: each gives the same output as its definition. The undocumented alias `cdna` equals splice."""
	fq = hifi_tx["fq"]

	def body(args):
		return H.sam_body(H.map_reads(synth_tx.fasta, fq, ["-a"] + args).stdout)
	splice = body(["-x", "splice"])
	assert splice == body(["-k15", "-w5", "--splice", "-g2k", "-G200k", "-A1", "-B2", "-O2,32", "-E1,0", "-C9", "-z200", "-ub",
		"--junc-bonus=9", "--cap-sw-mem=0", "--splice-flank=yes"])
	assert body(["-x", "splice:hq"]) == body(["-x", "splice", "-C5", "-O6,24", "-B4"])
	assert body(["-x", "cdna"]) == splice


@pytest.mark.pinned
def test_splice_ont_with_annotation_pinned(synth_tx, tmp_path):
	"""Pins the current output: `-ax splice --junc-bed anno.bed12` (README) on 40 ONT-like cDNA reads,
	SAM without @PG, byte for byte."""
	reads = H.simulate_spliced_reads(synth_tx.genome, 40, seed=203, profile="ont", qual=None)
	fa = H.write_fasta(str(tmp_path / "ont.fa"), reads)
	res = H.map_reads(synth_tx.fasta, fa, ["-a", "--junc-bed", synth_tx.bed12], preset="splice")
	H.check_pinned("c/splice_ont_juncbed.sam.gz", H.strip_pg(res.stdout))


@pytest.mark.pinned
def test_splice_hq_paf_pinned(synth_tx, hifi_tx):
	"""Pins the current output: `-cx splice:hq -uf --cs` PAF on the 60 HiFi-like cDNA reads, byte for byte."""
	res = H.map_reads(synth_tx.fasta, hifi_tx["fq"], ["-c", "--cs", "-uf"], preset="splice:hq")
	H.check_pinned("c/splice_hq_uf.paf.gz", res.stdout)
