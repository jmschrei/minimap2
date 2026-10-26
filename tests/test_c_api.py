"""Area C: the splice and short-read parts of the libminimap2 API through the
C driver tests/c/drv_c_api.c: preset option values against minimap2.1,
mm_mapopt_max_intron_len, the junction and splice-score readers
(mm_idx_bed_read, mm_idx_bed_junc, mm_idx_spsc_read, mm_idx_spsc_read2,
mm_idx_spsc_get, mm_idx_jjump_read, mm_idx_jump_get), mm_map in splice mode,
mm_map_frag on read pairs (proper_frag, seg_split, seg_id, pe_thru) and
mm_map_file_frag against the command line."""

import pytest

import helpers_c as C
import minimap2_helpers as H

DRV = "drv_c_api"


def test_preset_option_fields():
	"""mm_set_opt() for sr/short, splice/cdna, splice:hq and splice:sr sets what minimap2.1 lists for
	each preset (flags MM_F_SR, MM_F_FRAG_MODE, MM_F_SPLICE*, MM_F_WEAK_PAIRING, MM_F_SR_RNA; -F 800,
	-G 200k, -C, --junc-bonus 9, --spsc0 5, --jump-min-match 3, forward-reverse pe_ori for the paired
	presets); the aliases equal their presets; MM_MAX_SEG is 255; a misspelt splice preset is an error."""
	H.assert_c_case(DRV, "preset_fields")


def test_max_intron_len_api():
	"""mm_mapopt_max_intron_len() sets the maximum reference gap and both band widths in splice mode
	only, and ignores values that are not positive."""
	H.assert_c_case(DRV, "max_intron_len")


@pytest.fixture(scope="module")
def tx_files(synth_tx, tmp_path_factory):
	d = tmp_path_factory.mktemp("api")
	g = synth_tx.genome
	extra = "chrUn\t100\t900\tjunk\t0\t+\n"  # unknown contig: skipped
	bed12 = H.write_text(str(d / "anno.bed12"), H.genes_bed12(g.genes) + extra)
	bed6 = H.write_text(str(d / "junc.bed6"), H.junctions_bed(g.genes) + H.junctions_bed(g.genes[:2]) + extra)
	return dict(dir=d, g=g, fa=synth_tx.fasta, bed12=bed12, bed6=bed6)


@pytest.mark.parametrize("bed,read_junc", [("bed12", 1), ("bed6", 1), ("bed12", 0)])
def test_bed_junc(tx_files, bed, read_junc):
	"""mm_idx_bed_read() then mm_idx_bed_junc(): in a window, bits 1 and 2 mark the first and last base
	of each '+' intron lying wholly inside it, 8 and 4 those of '-' introns (the annotation used by
	--junc-bed); every other base is 0. A BED12 gene contributes the gaps between its blocks, a BED6
	line its own interval; with read_junc=0 a BED12 line is one interval from chromStart to chromEnd.
	Duplicate lines and lines on unknown contigs change nothing. Before any BED is read, and for a contig
	out of range, the call returns -1 with the window zeroed."""
	g = tx_files["g"]
	if read_junc:
		ivs = [(s, e, gg.strand, gg.contig) for gg in g.genes for s, e in gg.introns]
	else:
		ivs = [(gg.start, gg.end, gg.strand, gg.contig) for gg in g.genes]
	for ctg in ("chr1", "chr2"):
		mine = [(s, e, st) for s, e, st, c in ivs if c == ctg]
		windows = [(0, g.length(ctg))]
		a = mine[1]  # a window that cuts through the second interval
		windows.append((a[0] + 1, g.length(ctg) - 1))
		windows.append((max(0, a[0] - 50), a[1] + 50))
		for st, en in windows:
			bits = C.intron_bits(mine, st, en)
			H.assert_c_case(DRV, "bed_junc", tx_files["fa"], tx_files[bed], read_junc, ctg, st, en, C.encode_bits(bits))


def _scaled(score, scale, max_sc):
	if 0.0 < scale < 1.0:  # minimap2.1: scaled by --spsc-scale and rounded to the nearest integer
		score = int(score * scale + 0.499) if score > 0 else int(score * scale - 0.499)
	return max(-max_sc, min(max_sc, score))


@pytest.mark.parametrize("scale", ["read", "0.7", "0.6"])
def test_spsc_reader(tx_files, scale):
	"""mm_idx_spsc_read2() reads --spsc files (minimap2.1: contig, offset, + or -, D or A, score): the
	score is scaled and rounded, capped at +-max_sc, and stored per contig and strand; reading back with
	mm_idx_spsc_get() gives it at its offset and nothing (0xff) elsewhere, the larger score where an
	offset appears twice. Offsets 0 and the contig length, lines with an unknown contig, a type other
	than D or A, or fewer than five fields are ignored. mm_idx_spsc_read() is the same without scaling.
	mm_idx_spsc_get() returns the window length (clamped to the contig) or -1 without scores or for a
	contig out of range."""
	g = tx_files["g"]
	n1, n2 = g.length("chr1"), g.length("chr2")
	rows = [("chr1", 1200, "+", "D", 10), ("chr1", 1900, "+", "A", -10), ("chr1", 2500, "-", "A", 3),
		("chr1", 2600, "-", "D", -4), ("chr2", 700, "+", "D", 100), ("chr2", 800, "-", "A", -100),
		("chr1", 3000, "+", "D", 2), ("chr1", 3000, "+", "D", 9),  # duplicate offset: the larger wins
		("chr1", 0, "+", "D", 7), ("chr1", n1, "+", "A", 7), ("chr2", n2, "-", "D", 7),  # contig ends: ignored
		("chrUn", 500, "+", "D", 7), ("chr1", 4000, "+", "X", 7)]  # unknown contig and type: ignored
	text = "".join("%s\t%d\t%s\t%s\t%d\n" % r for r in rows) + "chr1\t4100\t+\tD\n"  # four fields: ignored
	f = H.write_text(str(tx_files["dir"] / ("sc_%s.tsv" % scale)), text)
	max_sc = 30
	sc = 1.0 if scale == "read" else float(scale)
	want = {}
	for ctg, pos, strand, typ, score in rows:
		if ctg == "chrUn" or typ not in "DA" or pos <= 0 or pos >= g.length(ctg):
			continue
		v = _scaled(score, sc, max_sc)
		k = (ctg, strand, pos)
		if k not in want or v > want[k][1]:
			want[k] = (typ, v)
	exp = ["%s:%s:%d:%s:%d" % (k[0], k[1], k[2], t, v) for k, (t, v) in sorted(want.items())]
	H.assert_c_case(DRV, "spsc", tx_files["fa"], f, max_sc, scale, *exp)


def test_jump_reader(tx_files):
	"""mm_idx_jjump_read() and mm_idx_jump_get(), behind -j and --pass1: every intron of the annotation
	gives a jump from its start to its end and one from its end to its start; mm_idx_jump_get() returns
	those whose origin lies in (st, en], sorted by origin. Junctions from -j carry MM_JUNC_ANNO (1),
	from --pass1 MM_JUNC_MISC (2), both when a junction is in both; --pass1 rows scoring below 5 are
	left out (minimap2.1)."""
	g = tx_files["g"]
	d = tx_files["dir"]
	anno = [gg for gg in g.genes[:4]]
	a_ivs = set((gg.contig, s, e) for gg in anno for s, e in gg.introns)
	p_rows = []
	for i, gg in enumerate(g.genes[2:]):
		for j, (s, e) in enumerate(gg.introns):
			p_rows.append((gg.contig, s, e, 6 if (i + j) % 3 else 4, gg.strand))
	p1 = H.write_text(str(d / "pass1.bed"), "".join("%s\t%d\t%d\tr\t%d\t%s\n" % r for r in p_rows))
	a1 = H.write_text(str(d / "anno4.bed12"), H.genes_bed12(anno))
	p_ivs = set((c, s, e) for c, s, e, sc, _ in p_rows if sc >= 5)
	for ctg in ("chr1", "chr2"):
		flags = {}
		for c, s, e in a_ivs | p_ivs:
			if c != ctg:
				continue
			f = (1 if (c, s, e) in a_ivs else 0) | (2 if (c, s, e) in p_ivs else 0)
			flags[(s, e)] = f
			flags[(e, s)] = f
		for st, en in ((0, g.length(ctg)), (20001, 30001), (30001, 47001)):
			exp = ",".join("%d:%d:%d" % (o, o2, f) for (o, o2), f in sorted(flags.items()) if st < o <= en)
			H.assert_c_case(DRV, "jump_get", tx_files["fa"], a1, p1, ctg, st, en, exp)
			H.assert_c_case(DRV, "jump_get", tx_files["fa"], a1, "-", ctg, st, en,
				",".join("%d:%d:1" % (o, o2) for (o, o2), f in sorted(flags.items()) if st < o <= en and f & 1))


def test_map_splice_api(synth_tx, tmp_path):
	"""mm_map() with the splice:hq preset on HiFi-like cDNA reads: each primary hit lies on the true
	contig and strand, is_spliced is set exactly when its CIGAR has N, at least 57 of 60 have exactly
	the true introns, and on those mm_extra_t.trans_strand is 1 for a read on the transcript strand
	and 2 for the opposite (ts:A in cookbook.md)."""
	reads = H.simulate_spliced_reads(synth_tx.genome, 60, seed=501, profile="hifi", qual=None)
	fa = H.write_fasta(str(tmp_path / "r.fa"), reads)
	tsv = C.spliced_truth_tsv(str(tmp_path / "truth.tsv"), reads)
	H.assert_c_case(DRV, "map_splice", synth_tx.fasta, "splice:hq", fa, tsv, 57)


def test_map_frag_api(synth, tmp_path):
	"""mm_map_frag() on forward-reverse pairs (mate 2 reverse-complemented first, as map.c does for
	the sr preset): one primary per mate at its true contig, strand and position, proper_frag on both;
	every hit has seg_split set and seg_id = its mate's index; pe_thru marks exactly the pairs whose
	fragment is no longer than the reads (both mates cover the same bases). mm_map_frag() with one
	segment gives the hits of mm_map(); more than MM_MAX_SEG segments give no hits."""
	g = synth.genome
	pairs = H.simulate_pairs(g, 30, seed=502, insert_mean=320, insert_sd=40, min_insert=200, unique=True, indel_rate=0, qual=None)
	thru = H.simulate_pairs(g, 10, seed=503, insert_mean=150, insert_sd=0, unique=True, indel_rate=0, qual=None, name_prefix="t")
	allp = pairs + thru
	f1, f2 = H.write_pairs(str(tmp_path / "1.fa"), str(tmp_path / "2.fa"), allp)
	tsv = str(tmp_path / "truth.tsv")
	with open(tsv, "w") as f:
		for p in allp:
			f.write("%s\t%s\t%d\t%s\t%d\t%s\t%d\n" % (p.name, p.contig, p.r1.pos, p.r1.strand, p.r2.pos, p.r2.strand, p.insert))
	assert all(p.insert == 150 for p in thru) and all(p.insert >= 200 for p in pairs)
	H.assert_c_case(DRV, "map_frag", synth.fasta, f1, f2, tsv, 0)


def _driver_sam(res):
	return "".join(l + "\n" for l in res.stdout.splitlines() if l and not l.startswith(("#", "ok ", "not ok ")))


def test_map_file_frag_api(synth, synth_tx, tmp_path):
	"""mm_map_file_frag() with two query files writes the same SAM records as `minimap2 -ax sr r1 r2`,
	and with one file the same as `minimap2 -ax splice reads`; a query file that cannot be opened
	makes it return -1."""
	pairs = H.simulate_pairs(synth.genome, 60, seed=504)
	f1, f2 = H.write_pairs(str(tmp_path / "1.fq"), str(tmp_path / "2.fq"), pairs)
	res = H.assert_c_case(DRV, "map_file_frag", synth.fasta, "sr", f1, f2)
	cli = H.map_reads(synth.fasta, [f1, f2], ["-a"], preset="sr")
	assert _driver_sam(res) == H.sam_body(cli.stdout) and len(cli.sam().records) >= 120
	reads = H.simulate_spliced_reads(synth_tx.genome, 20, seed=505, profile="hifi")
	fq = H.write_fastq(str(tmp_path / "tx.fq"), reads)
	res = H.assert_c_case(DRV, "map_file_frag", synth_tx.fasta, "splice", fq)
	assert _driver_sam(res) == H.sam_body(H.map_reads(synth_tx.fasta, fq, ["-a"], preset="splice").stdout)
