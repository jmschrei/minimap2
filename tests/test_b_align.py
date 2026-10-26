"""Area B: base-level alignment (align.c, -c) and its options (-A -B -b -O -E --score-N -z -s
--end-bonus --max-clip-ratio --cap-sw-mem --min-dp-len --no-end-flt --end-seed-pen), the AS:i and
ms:i tags, inversion rescue (tp:A:I) and chimeric reads (supplementary alignments, SA:Z).

Expected scores are computed from the CIGAR with the documented scoring (minimap2.1: +A per match,
-B per mismatch, -b per transition, --score-N per column with an ambiguous base, a k-long gap
costing min{O1+k*E1, O2+k*E2})."""

import math
import os

import pytest

import helpers_b as B
import minimap2_helpers as H


@pytest.fixture(scope="module")
def clean(tmp_path_factory):
	"""A 60 + 40 kbp random genome without repeats and its FASTA."""
	d = tmp_path_factory.mktemp("b_align")
	g = H.make_genome(seed=21, contig_lengths=(60000, 40000))

	class D(object):
		pass
	o = D()
	o.genome, o.fasta, o.dir = g, g.write_fasta(str(d / "clean.fa")), str(d)
	return o


@pytest.fixture(scope="module")
def read_sets(clean):
	"""ONT-like reads (about 7% error) and HiFi-like reads (0.4%), both with 0.2% N bases."""
	ont = H.simulate_long_reads(clean.genome, 25, profile="ont", seed=5, length=(1000, 3000), n_rate=0.002)
	hifi = H.simulate_long_reads(clean.genome, 25, profile="hifi", seed=6, length=(1000, 3000), n_rate=0.002)
	return {
		"ont": (H.write_fasta(os.path.join(clean.dir, "ont.fa"), ont), dict((r.name, r.seq) for r in ont), ont),
		"hifi": (H.write_fasta(os.path.join(clean.dir, "hifi.fa"), hifi), dict((r.name, r.seq) for r in hifi), hifi),
	}


# The scoring each preset documents (minimap2.1; map-iclr-prerender and lr:hqae from options.c, which
# minimap2.1 does not list). transition is -b.
PRESET_SCORING = {
	"map-ont": dict(a=2, b=4, q=4, e=2, q2=24, e2=1),
	"lr": dict(a=2, b=4, q=4, e=2, q2=24, e2=1),
	"map-pb": dict(a=2, b=4, q=4, e=2, q2=24, e2=1),
	"map10k": dict(a=2, b=4, q=4, e=2, q2=24, e2=1),
	"lr:hq": dict(a=2, b=4, q=4, e=2, q2=24, e2=1),
	"lr:hqae": dict(a=2, b=4, q=4, e=2, q2=24, e2=1),
	"map-hifi": dict(a=1, b=4, q=6, e=2, q2=26, e2=1),
	"map-ccs": dict(a=1, b=4, q=6, e=2, q2=26, e2=1),
	"map-iclr": dict(a=2, b=6, q=10, e=2, q2=50, e2=1, transition=4),
	"map-iclr-prerender": dict(a=2, b=6, q=10, e=2, q2=50, e2=1, transition=1),
	"asm5": dict(a=1, b=19, q=39, e=3, q2=81, e2=1),
	"asm10": dict(a=1, b=9, q=16, e=2, q2=41, e2=1),
	"asm20": dict(a=1, b=4, q=6, e=2, q2=26, e2=1),
}


def _check_as(clean, paf, seqs, scoring, min_lines):
	assert len(paf) >= min_lines
	assert H.paf_invariant_errors(paf, clean.genome, seqs) == []
	for p in paf:
		s = H.paf_query_part(p, seqs[p.qname])
		want = B.cigar_score_full(clean.genome[p.tname], p.tstart, s, p.cigar, **scoring)
		assert p.tags["AS"] == want, (p.line[:200], want)


@pytest.mark.parametrize("preset", sorted(set(PRESET_SCORING) - {"lr", "map10k", "map-ccs"}))  # aliases: test_b_presets
def test_AS_is_cigar_score_under_preset_scoring(clean, read_sets, preset):
	"""-c with each long-read and assembly preset: AS:i equals the score of the CIGAR under the
	preset's documented -A/-B/-b/-O/-E (and --score-N 1) on reads with N bases; cg:Z, NM, columns 10
	and 11 and nn:i agree with the genome. ONT-like reads for the noisy-read presets, HiFi-like for
	the accurate-read and assembly presets."""
	kind = "hifi" if preset.startswith(("asm", "map-hifi", "map-ccs", "lr:hq")) else "ont"
	fa, seqs, _ = read_sets[kind]
	paf = H.map_reads(clean.fasta, fa, ["-c"], preset=preset).paf()
	_check_as(clean, paf, seqs, PRESET_SCORING[preset], 20)


@pytest.mark.parametrize("args,scoring", [
	(["-A", "3", "-B", "5", "-O", "5,30", "-E", "3,1"], dict(a=3, b=5, q=5, e=3, q2=30, e2=1)),
	(["-A", "1", "-B", "3", "-O", "3", "-E", "1"], dict(a=1, b=3, q=3, e=1, q2=3, e2=1)),
	(["-O", "6", "-E", "2"], dict(a=2, b=4, q=6, e=2, q2=6, e2=2)),
	(["-b", "2"], dict(a=2, b=4, q=4, e=2, q2=24, e2=1, transition=2)),
	(["-B", "6", "-b", "6"], dict(a=2, b=6, q=4, e=2, q2=24, e2=1)),
	(["--score-N", "3"], dict(a=2, b=4, q=4, e=2, q2=24, e2=1, sc_ambi=3)),
], ids=["two-piece", "affine", "affine-default-A", "transition", "transition-equal-B", "score-N"])
def test_AS_is_cigar_score_under_custom_scoring(clean, read_sets, args, scoring):
	"""-A, -B, -b, -O, -E (one value sets both pieces: affine gaps through ksw_extz2; two values the
	two-piece model through ksw_extd2) and --score-N: AS:i equals the CIGAR's score under the given
	scoring on ONT-like reads with N bases."""
	fa, seqs, _ = read_sets["ont"]
	paf = H.map_reads(clean.fasta, fa, ["-c"] + args).paf()
	_check_as(clean, paf, seqs, scoring, 20)


def test_peak_score_ms_and_s_boundary(clean, tmp_path):
	"""ms:i and -s ("minimal peak DP alignment score ... the score of the max scoring segment in the
	alignment"): for reads with substitutions only (gap-free alignments), ms:i is the best sum of
	column scores over a run of columns; -s ms keeps the alignment and -s ms+1 drops it (-c)."""
	reads = H.simulate_reads(clean.genome, 6, length=(300, 900), seed=12, sub_rate=0.06, indel_rate=0, qual=None, name_prefix="pk")
	fa = H.write_fasta(str(tmp_path / "pk.fa"), reads)
	paf = H.map_reads(clean.fasta, fa, ["-c"]).paf()
	seqs = dict((r.name, r.seq) for r in reads)
	ms = {}
	for p in paf:
		assert p.is_primary and "I" not in p.cigar and "D" not in p.cigar, p.line
		s = H.paf_query_part(p, seqs[p.qname])
		ms[p.qname] = B.peak_score(clean.genome[p.tname], p.tstart, s, p.cigar)
		assert p.tags["ms"] == ms[p.qname], p.line
	assert sorted(ms) == sorted(seqs)
	for r in reads[:3]:
		one = H.write_fasta(str(tmp_path / (r.name + ".fa")), [r])
		assert len(H.map_reads(clean.fasta, one, ["-c", "-s", ms[r.name]]).paf()) == 1
		assert H.map_reads(clean.fasta, one, ["-c", "-s", ms[r.name] + 1]).stdout == ""


def test_end_bonus_boundary(clean, tmp_path):
	"""--end-bonus INT ("score bonus when alignment extends to the end of the query"): a read of 900
	reference bases whose last 7 bases match nothing nearby. The right extension stops at the best
	score unless reaching the read end scores more once the bonus is added: with the drop D from the
	best extension score to the best score ending at the read end (brute-force extension DP),
	--end-bonus D leaves the 7 bases clipped and D+1 aligns them."""
	ref = clean.genome["chr2"]
	start = 15000
	seq = list(ref[start:start + 900])
	for i in range(893, 900):  # avoid the reference base and its neighbours, so no shifted match either
		bad = set(ref[start + i - 1:start + i + 2].upper())
		seq[i] = [c for c in "ACGT" if c not in bad][0]
	seq = "".join(seq)
	tail_q = seq[800:]
	tail_t = ref[start + 800:start + 1100]
	mx, mqe = B.affine_ext(tail_q, tail_t, 2, 4, 4, 2, 24, 1)
	D = mx - mqe
	assert 7 < D < 30
	fa = H.write_fasta(str(tmp_path / "eb.fa"), [("eb", seq)])
	for bonus, qend in ((D, 893), (D + 1, 900), (0, 893), (200, 900)):
		p = H.map_reads(clean.fasta, fa, ["-c", "--end-bonus", bonus]).paf()
		assert len(p) == 1 and p[0].qstart == 0 and p[0].qend == qend, (bonus, D, p)


def test_zdrop_boundary(clean, tmp_path):
	"""-z INT: a read of 2000 reference bases, 300 bases with 50% substitutions and 2000 more. Let D
	be the largest Z-drop along its unbroken alignment (from -z 10000): the running score's fall below
	its maximum minus E times the change of diagonal. -z D keeps one alignment; -z D-1 truncates it
	into two alignments of the flanks (zd:i:1 on the left one, zd:i:2 on the right one) whose cm:i and
	s1:i add up to the unbroken chain's."""
	ref = clean.genome["chr1"]
	rng = H.Rng(3)
	seq = ref[40000:42000] + H.mutate(rng, ref[42000:42300], 0.5) + ref[42300:44300]
	fa = H.write_fasta(str(tmp_path / "div.fa"), [("div", seq)])
	whole = H.map_reads(clean.fasta, fa, ["-c", "-z", "10000"]).paf()
	assert len(whole) == 1 and (whole[0].qstart, whole[0].qend) == (0, 4300)
	p = whole[0]
	D = B.max_zdrop(ref, p.tstart, seq[p.qstart:p.qend], p.cigar)
	assert 100 < D < 400
	assert len(H.map_reads(clean.fasta, fa, ["-c", "-z", D]).paf()) == 1
	split = H.map_reads(clean.fasta, fa, ["-c", "-z", D - 1]).paf()
	assert sorted((x.qstart < 100, x.qend > 4200, x.tags.get("zd")) for x in split) == [(False, True, 2), (True, False, 1)], split
	for x in split:
		assert x.tstart - x.qstart == 40000 and x.cigar == "%dM" % (x.qend - x.qstart)
	# splitting the chain at the drop splits its minimizers and its score between the two parts
	assert sum(x.tags["cm"] for x in split) == p.tags["cm"] and sum(x.tags["s1"] for x in split) == p.tags["s1"]


@pytest.mark.parametrize("length,args", [(300, []), (600, []), (200, ["-z", "200,100"])])
def test_inversion_rescue(clean, tmp_path, length, args):
	"""Inversion rescue (-z INT1,INT2; tp:A:I): a read of 2000 reference bases, the reverse complement
	of the next `length` bases, and 2000 more. The flanks align on + with zd:i:1/zd:i:2, and the
	inverted piece aligns on - exactly at its origin as an inversion (tp:A:I) with -c, with a mapping
	quality no higher than its flanks'."""
	ref = clean.genome["chr1"]
	seq = ref[30000:32000] + H.revcomp(ref[32000:32000 + length]) + ref[32000 + length:34000 + length]
	fa = H.write_fasta(str(tmp_path / "inv.fa"), [("inv", seq)])
	paf = H.map_reads(clean.fasta, fa, ["-c"] + args).paf()
	inv = [p for p in paf if p.tags["tp"] == "I"]
	assert len(inv) == 1, paf
	i = inv[0]
	assert (i.strand, i.qstart, i.qend, i.tstart, i.tend, i.cigar) == ("-", 2000, 2000 + length, 32000, 32000 + length, "%dM" % length), i.line
	flanks = sorted((p.qstart, p.qend, p.tstart, p.strand, p.tags.get("zd")) for p in paf if p.tags["tp"] == "P")
	assert flanks == [(0, 2000, 30000, "+", 1), (2000 + length, 4000 + length, 32000 + length, "+", 2)], flanks
	assert i.mapq <= min(p.mapq for p in paf if p.tags["tp"] == "P")
	assert H.paf_invariant_errors(paf, clean.genome, {"inv": seq}) == []


def test_chimeric_reads_supplementary_and_SA(clean, tmp_path):
	"""Chimeric reads (two pieces from loci at least 10 kbp apart, HiFi-like errors): in SAM each read
	has one primary and one supplementary record (0x800), each placing its piece at the simulated
	locus and strand (unclipped coordinates within 10 bp); the supplementary is hard-clipped; each
	record's SA:Z names the other record's contig, position and strand."""
	reads = H.simulate_chimeric_reads(clean.genome, 10, profile="hifi", length=(2000, 5000), seed=9)
	fa = H.write_fasta(str(tmp_path / "chim.fa"), reads)
	sam = H.map_reads(clean.fasta, fa, ["-a"]).sam()
	assert H.sam_invariant_errors(sam, clean.genome, dict((r.name, [r.seq]) for r in reads)) == []
	by = sam.by_qname()
	for r in reads:
		recs = by[r.name]
		assert sorted(x.is_supplementary for x in recs) == [False, True], [x.line for x in recs]
		for x in recs:
			assert not x.is_secondary and x.mapq > 0
			seg = [s for s in r.segments if s.contig == x.rname and (s.strand == "-") == x.is_reverse]
			lc = x.left_clip
			ok = False
			for s in seg:
				a, b = (s.qstart, s.qend) if s.strand == "+" else (len(r.seq) - s.qend, len(r.seq) - s.qstart)
				ok = ok or (abs(x.pos0 - s.pos) <= 10 and abs(lc - a) <= 10)
			assert ok, (x.line[:200], r.segments)
			if x.is_supplementary:
				assert "H" in x.cigar and "S" not in x.cigar
			other = [y for y in recs if y is not x][0]
			sa = x.tags["SA"].rstrip(";").split(",")
			assert (sa[0], int(sa[1]), sa[2]) == (other.rname, other.pos, "-" if other.is_reverse else "+"), (x.tags["SA"], other.line[:100])


def test_max_clip_ratio_boundary(clean, tmp_path):
	"""--max-clip-ratio: an alignment is dropped when both its unaligned ends exceed RATIO times the
	read length. A read of 1000 random bases, 1000 reference bases and 1000 random bases (both ends
	1000 of 3000) is dropped at 0.333 and kept at 0.334 and by default (-c); a read clipped at one end
	only (1000 random + 2000 reference bases) is kept at 0.3."""
	rng = H.Rng(7)
	seq = rng.dna(1000) + clean.genome["chr1"][5000:6000] + rng.dna(1000)
	fa = H.write_fasta(str(tmp_path / "clip.fa"), [("clip", seq)])
	assert H.map_reads(clean.fasta, fa, ["-c", "--max-clip-ratio", "0.333"]).stdout == ""
	for args in (["--max-clip-ratio", "0.334"], []):
		p = H.map_reads(clean.fasta, fa, ["-c"] + args).paf()
		assert [(x.qstart, x.qend, x.tstart, x.tend) for x in p] == [(1000, 2000, 5000, 6000)], (args, p)
	one = H.write_fasta(str(tmp_path / "clip1.fa"), [("clip1", rng.dna(1000) + clean.genome["chr1"][5000:7000])])
	p = H.map_reads(clean.fasta, one, ["-c", "--max-clip-ratio", "0.3"]).paf()
	assert [(x.qstart, x.qend, x.tstart, x.tend) for x in p] == [(1000, 3000, 5000, 7000)], p


def test_cap_sw_mem(clean, read_sets):
	"""--cap-sw-mem NUM ("skip alignment if the DP matrix size is above NUM; 0 to disable"): at 1 every
	DP is skipped and no read keeps an alignment; 0 and a cap above every matrix give the same output
	as the default."""
	fa, seqs, _ = read_sets["ont"]
	base = H.strip_pg(H.map_reads(clean.fasta, fa, ["-c"]).stdout)
	assert len(base.splitlines()) >= 20
	assert H.map_reads(clean.fasta, fa, ["-c", "--cap-sw-mem", "1"]).stdout == ""
	for cap in ("0", "1g"):
		assert H.map_reads(clean.fasta, fa, ["-c", "--cap-sw-mem", cap]).stdout == base


@pytest.mark.parametrize("args", [["--min-dp-len", "10"], ["--min-dp-len", "5000"], ["--no-end-flt"], ["--max-chain-skip", "0"]])
def test_alignment_options_keep_reads_correct(clean, read_sets, args):
	"""--min-dp-len (gap-filling DP only between anchors at least that far apart; a huge value
	aligns from the first to the last anchor in one DP) and --no-end-flt (keep the chain-end seeds):
	every ONT-like read still aligns at its simulated locus within 30 bp at both ends, and AS:i equals
	the CIGAR's score."""
	fa, seqs, reads = read_sets["ont"]
	paf = H.map_reads(clean.fasta, fa, ["-c"] + args).paf()
	_check_as(clean, paf, seqs, PRESET_SCORING["map-ont"], 20)
	prim = dict((p.qname, p) for p in paf if p.is_primary)
	for r in reads:
		err = H.paf_placement_error(prim[r.name], r)
		assert err is not None and err <= 30, (args, prim[r.name].line[:150], r)


def test_end_seed_pen_boundary(tmp_path):
	"""--end-seed-pen INT (splice mode: "drop a terminal anchor if s < log(g)+INT, where s is the local
	alignment score around the anchor and g the length of the terminal gap in the chain"): a cDNA read
	whose 19 bp first exon holds one anchor, 1500 bp of intron away from the next. With s the
	brute-force local score of the anchor's k-mer plus 20 bases each side and g the distance to the
	next anchor, INT = floor(s - ln g) keeps the exon (19M1500N600M) and INT + 1 drops it (soft-clipped,
	600M)."""
	rng = H.Rng(1)
	ex1, intron, ex2 = rng.dna(19), "GT" + rng.dna(1496) + "AG", rng.dna(600)
	ex1, ex2 = ex1[:-1] + "C", "T" + ex2[1:]
	ref = rng.dna(3000) + ex1 + intron + ex2 + rng.dna(3000)
	read = ex1 + ex2
	fa = H.write_fasta(str(tmp_path / "sp.fa"), [("ref", ref)])
	rq = H.write_fasta(str(tmp_path / "spr.fa"), [("tx", read)])
	idx = B.MinimizerIndex({"ref": ref}, 5, 15)
	anc = sorted([a for a in B.anchors(read, idx) if a.strand == "+"], key=lambda a: a.x)
	assert anc[0].y < 19 <= anc[1].y
	g = anc[1].x - anc[0].x
	qs, qe = max(0, anc[0].y + 1 - 15 - 20), anc[0].y + 1 + 20
	rs, re = anc[0].x + 1 - 15 - 20, anc[0].x + 1 + 20
	s = B.local_score(read[qs:qe], ref[rs:re], 1, 2, 2, 1)
	keep = int(math.floor(s - math.log(g)))
	assert 15 < math.log(g) + keep + 1  # the anchor (span 15) is examined at INT + 1
	for pen, want in ((keep, (0, "19M1500N600M")), (keep + 1, (19, "600M")), (6, (0, "19M1500N600M"))):
		p = H.map_reads(fa, rq, ["-c", "--end-seed-pen", pen], preset="splice").paf()
		assert [(x.qstart, x.cigar) for x in p] == [want], (pen, s, g, p)


@pytest.fixture(scope="module")
def special(synth, tmp_path_factory):
	"""Reads exercising the rarer alignment paths on the standard synthetic genome: an inversion, a
	600 bp deletion, a 300 bp stretch of 50% substitutions, 400 N bases, three chimeric reads, and
	ONT-like reads overlapping repeats, inverted repeats, N runs and low-complexity tracts."""
	d = tmp_path_factory.mktemp("b_special")
	g = synth.genome
	ref = g["chr1"]
	rng = H.Rng(501)
	recs = [("inv", ref[30000:32000] + H.revcomp(ref[32000:32400]) + ref[32400:34400]),
		("del600", ref[20000:22000] + ref[22600:24600]),
		("div", ref[40000:42000] + H.mutate(rng, ref[42000:42300], 0.5) + ref[42300:44300]),
		("nrun", ref[60000:61500] + "N" * 400 + ref[61900:63400])]
	recs += [(c.name, c.seq) for c in H.simulate_chimeric_reads(g, 3, profile="hifi", length=(2000, 4000), seed=502)]
	feats = [f for f in g.features if f.kind in ("repeat", "inversion", "N", "lowcomp")]
	for i, f in enumerate(feats):
		s = max(0, f.start - 800)
		piece = g[f.contig][s:f.end + 800]
		recs.append(("feat%d" % i, H.mutate(rng, piece, 0.02) if i % 2 else H.revcomp(piece)))
	return H.write_fasta(str(d / "special.fa"), recs)


@pytest.mark.pinned
@pytest.mark.parametrize("name,preset,extra", [("special-map-ont.sam.gz", "map-ont", ["-a", "--cs", "-Y"]),
	("special-map-hifi.paf.gz", "map-hifi", ["-c", "--cs"]), ("special-asm20.paf.gz", "asm20", ["-c"])])
def test_special_reads_pinned(synth, special, name, preset, extra):
	"""Pins the current output of base minimap2 on reads with inversions, a long deletion, a divergent
	stretch, an N run, chimeric pieces and genome features (@PG removed)."""
	H.check_pinned("b/" + name, H.strip_pg(H.map_reads(synth.fasta, special, extra, preset=preset).stdout))
