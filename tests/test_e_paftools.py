"""misc/paftools.js, every subcommand, on small hand-made inputs (area E).

Each input is a few lines built so that the expected output can be worked out
by hand from what the subcommand is documented to do (misc/README.md, the
usage texts, the PAF/SAM/VCF/BED/GTF/MAF formats). The scripts need the k8
JavaScript shell; these tests skip when k8 is not on PATH.
"""

import gzip

import pytest

import minimap2_helpers as H

pytestmark = pytest.mark.k8


def _w(d, name, lines):
	"""Write lines (a list of str, or of lists joined with tabs) to d/name, newline-terminated; gzip for .gz."""
	text = "".join(("\t".join(str(x) for x in l) if isinstance(l, (list, tuple)) else l) + "\n" for l in lines)
	p = str(d / name)
	if p.endswith(".gz"):
		with open(p, "wb") as f:
			f.write(gzip.compress(text.encode(), mtime=0))
	else:
		with open(p, "w") as f:
			f.write(text)
	return p


def _pt(*args, **kw):
	"""Run `k8 paftools.js <args>`; check=False to allow a non-zero exit."""
	return H.run_paftools([str(a) for a in args], **kw)


def _tsv(text):
	return [l.split("\t") for l in text.splitlines()]


# ---------------------------------------------------------------------------
# main, version
# ---------------------------------------------------------------------------


def test_no_command_prints_usage_and_exits_1(k8):
	"""`paftools.js` alone prints the command list to stdout and exits 1; an unknown command fails."""
	r = _pt(check=False)
	assert r.returncode == 1
	assert r.stdout.startswith("Usage: paftools.js <command> [arguments]")
	for cmd in ("view", "splice2bed", "sam2paf", "delta2paf", "gff2bed", "gff2junc", "longcs2seq", "stat", "asmstat",
			"asmgene", "misjoin", "liftover", "call", "bedcov", "vcfstat", "sveval", "version", "mapeval", "pafcmp",
			"mason2fq", "sim2bed", "pbsim2fq", "badread2fa", "junceval", "exoneval", "ov-eval"):
		assert "\n  %s " % cmd in r.stdout, cmd
	bad = _pt("no-such-command", check=False)
	assert bad.returncode != 0 and "no-such-command" in bad.stderr


def test_version_is_the_minimap2_version(k8, mm2_bin):
	"""`paftools.js version` prints the release it ships with, the same string as `minimap2 --version`."""
	assert _pt("version").stdout == H.run_mm2(["--version"], check=True).stdout


# ---------------------------------------------------------------------------
# Conversion: view, sam2paf, delta2paf, splice2bed, gff2bed, gff2junc, longcs2seq, paf2gff
# ---------------------------------------------------------------------------

# q1[2,18) aligned to t1[5,21): 5 matches, a/g mismatch, insertion c, 4 matches, deletion t, 5 matches
VIEW_PAF = ["q1", 20, 2, 18, "+", "t1", 30, 5, 21, 14, 17, 60, "tp:A:P", "AS:i:9", "cg:Z:6M1I4M1D5M", "cs:Z::5*ag+c:4-t:5"]


def _pad(x, n, right):
	s = str(x)
	return s.ljust(n) if right else s.rjust(n)


def _blast_block(rs, qs, strand, slen, elen, ref, mid, qry):
	if strand == "+":
		st, en = qs + slen[1] + 1, qs + elen[1]
	else:
		st, en = qs - slen[1], qs - elen[1] + 1
	return [" ".join(["Ref+:", _pad(rs + slen[0] + 1, 10, False), ref, _pad(rs + elen[0], 10, True)]),
		" " * 17 + mid,
		" ".join(["Qry" + strand + ":", _pad(st, 10, False), qry, _pad(en, 10, True)])]


def test_view_blast_like(k8, tmp_path):
	"""view (default -f aln): a '>' line with the PAF line minus cg/cs plus mismatch (mm), insertion and deletion
	open (oi, od) and extension (ei, ed) counts, then Ref/match/Qry rows with 1-based start and end positions
	(Qry counts down on the reverse strand) and '//'. A line without cs is skipped."""
	fwd = list(VIEW_PAF)
	rev = list(VIEW_PAF)
	rev[0], rev[4] = "q2", "-"
	nocs = ["q3", 20, 0, 10, "+", "t1", 30, 0, 10, 10, 10, 60, "tp:A:P", "cg:Z:10M"]
	# multi-base gaps: ei and ed count bases, oi and od events
	gaps = ["q4", 20, 0, 10, "+", "t1", 30, 0, 11, 8, 13, 60, "tp:A:P", "cs:Z::3+gg:2-tta:3"]
	p = _w(tmp_path, "v.paf", [fwd, nocs, rev, gaps])
	out = _pt("view", p).stdout.splitlines()
	ref, mid, qry = "=====a-====t=====", "|||||  |||| |||||", "=====gc====-====="
	head = "\tq1\t20\t2\t18\t+\tt1\t30\t5\t21\t14\t17\t60\ttp:A:P\tAS:i:9"
	tail = "\tmm:i:1\toi:i:1\tei:i:1\tod:i:1\ted:i:1"
	want = [">" + head[1:] + tail] + _blast_block(5, 2, "+", [0, 0], [16, 16], ref, mid, qry) + ["//"]
	want += [">" + head[1:].replace("q1", "q2").replace("\t+\t", "\t-\t") + tail]
	want += _blast_block(5, 18, "-", [0, 0], [16, 16], ref, mid, qry) + ["//"]
	want += [">q4\t20\t0\t10\t+\tt1\t30\t0\t11\t8\t13\t60\ttp:A:P\tmm:i:0\toi:i:1\tei:i:2\tod:i:1\ted:i:3"]
	want += _blast_block(0, 0, "+", [0, 0], [11, 10], "===--==tta===", "|||  ||   |||", "===gg==---===") + ["//"]
	assert out == want
	assert want[1] == "Ref+:          6 =====a-====t===== 21        "
	assert want[8] == "Qry-:         18 =====gc====-===== 3         "


def test_view_line_length(k8, tmp_path):
	"""view -l 10 wraps the alignment every 10 columns; blocks are separated by an empty line and each block
	carries the positions of its own first and last base. With -l 5 the first run of matches fills a line
	exactly, so the mismatch after it starts the next line."""
	p = _w(tmp_path, "v.paf", [VIEW_PAF])
	out = _pt("view", "-l", "10", p).stdout.splitlines()
	want = [out[0]]
	want += _blast_block(5, 2, "+", [0, 0], [9, 10], "=====a-===", "|||||  |||", "=====gc===")
	want += [""] + _blast_block(5, 2, "+", [9, 10], [16, 16], "=t=====", "| |||||", "=-=====") + ["//"]
	assert out == want
	out = _pt("view", "-l", "5", p).stdout.splitlines()
	want = [out[0]] + _blast_block(5, 2, "+", [0, 0], [5, 5], "=====", "|||||", "=====")
	want += [""] + _blast_block(5, 2, "+", [5, 5], [9, 10], "a-===", "  |||", "gc===")
	want += [""] + _blast_block(5, 2, "+", [9, 10], [14, 14], "=t===", "| |||", "=-===")
	want += [""] + _blast_block(5, 2, "+", [14, 14], [16, 16], "==", "||", "==") + ["//"]
	assert out == want


def test_view_lastz_cigar_and_maf(k8, tmp_path):
	"""view -f lastz-cigar writes `cigar: qname qs qe strand tname ts te + score` and the cg operations as
	`op len` pairs; -f maf (which needs the long cs form) writes an `a score` line and one `s` line per
	sequence with start, aligned length, strand and length (query start counted on its own strand)."""
	p = _w(tmp_path, "v.paf", [VIEW_PAF])
	assert _pt("view", "-f", "lastz-cigar", p).stdout == "cigar: q1 2 18 + t1 5 21 + 9 M 6 I 1 M 4 D 1 M 5\n"
	long_cs = "cs:Z:=ACGTA*ag+c=TTGC-t=GGCAT"
	fwd = VIEW_PAF[:-1] + [long_cs]
	rev = ["q2", 20, 1, 17, "-", "t1", 30, 5, 21, 14, 17, 60, "tp:A:P", "AS:i:9", "cg:Z:6M1I4M1D5M", long_cs]
	p2 = _w(tmp_path, "long.paf", [fwd, rev])
	out = _pt("view", "-f", "maf", p2).stdout.splitlines()
	ref, qry = "ACGTAa-TTGCtGGCAT", "ACGTAgcTTGC-GGCAT"
	assert out == ["##maf version=1", "",
		"a 9", "s t1          5         16 +         30 " + ref, "s q1          2         16 +         20 " + qry, "",
		"a 9", "s t1          5         16 +         30 " + ref, "s q2          3         16 -         20 " + qry, ""]
	bad = _pt("view", "-f", "maf", p, check=False)
	assert bad.returncode != 0 and "--cs=long" in bad.stderr


# A SAM record with soft clips, a mismatch, an insertion and a deletion (SEQ in reference orientation):
# t1:5 ACGTAC|GATC C AGTCA ; read TT ACGTTC G GATC - AGTCA GGG
SAM2PAF_SAM = [
	"@HD\tVN:1.6",
	"@SQ\tSN:t1\tLN:100",
	["r1", 0, "t1", 6, 60, "2S6M1I4M1D5M3S", "*", 0, 0, "TTACGTTCGGATCAGTCAGGG", "*", "NM:i:3", "MD:Z:4A5^C5"],
	["r1", 2048, "t1", 81, 60, "16H5M", "*", 0, 0, "ACGTA", "*", "NM:i:0"],
	["r2", 16, "t1", 31, 7, "3H10M2H", "*", 0, 0, "ACGTACGTAC", "*", "NM:i:0"],
	["r2", 256, "t1", 61, 0, "10M", "*", 0, 0, "*", "*", "NM:i:1"],
	["r3", 4, "*", 0, 0, "*", "*", 0, 0, "ACGT", "*"],
	["r4", 65, "t1", 91, 60, "4M", "*", 0, 0, "ACGT", "*"],
	["r5", 0, "t1", 11, 60, "3=1X2=", "*", 0, 0, "ACGTAC", "*", "NM:i:1"],
]
SAM2PAF_OUT = {
	"r1": ["r1", "21", "2", "18", "+", "t1", "100", "5", "21", "14", "17", "60", "tp:A:P", "NM:i:3", "mm:i:1", "gn:i:2",
		"go:i:2", "cg:Z:6M1I4M1D5M"],
	"r1s": ["r1", "21", "16", "21", "+", "t1", "100", "80", "85", "5", "5", "60", "tp:A:P", "NM:i:0", "mm:i:0", "gn:i:0",
		"go:i:0", "cg:Z:5M"],
	"r2": ["r2", "15", "2", "12", "-", "t1", "100", "30", "40", "10", "10", "7", "tp:A:P", "NM:i:0", "mm:i:0", "gn:i:0",
		"go:i:0", "cg:Z:10M"],
	"r2s": ["r2", "10", "0", "10", "+", "t1", "100", "60", "70", "9", "10", "0", "tp:A:S", "NM:i:1", "mm:i:1", "gn:i:0",
		"go:i:0", "cg:Z:10M"],
	"r4": ["r4/1", "4", "0", "4", "+", "t1", "100", "90", "94", "4", "4", "60", "tp:A:P", "gn:i:0", "go:i:0", "cg:Z:4M"],
	"r5": ["r5", "6", "0", "6", "+", "t1", "100", "10", "16", "5", "6", "60", "tp:A:P", "NM:i:1", "mm:i:1", "gn:i:0",
		"go:i:0", "cg:Z:3=1X2="],
}


def test_sam2paf(k8, tmp_path):
	"""sam2paf: query length and interval from the clips (counted from the other end on the reverse strand),
	target length from @SQ, columns 10/11 from NM and the CIGAR (=/X counted directly), tp:A:S for
	secondaries, /1 on a first mate, cg:Z without clips, and cs:Z rebuilt from MD and SEQ (here
	:4*at:1+g:4-c:5). Unmapped records are dropped; -U writes them as `name len 0 0 * * 0 0 0 0 0 0`;
	-p drops secondaries, -P also supplementaries; -L writes the long cs form."""
	p = _w(tmp_path, "in.sam", SAM2PAF_SAM)
	out = _tsv(_pt("sam2paf", p).stdout)
	o = SAM2PAF_OUT
	assert [l[:-1] if l[0] == "r1" and l[2] == "2" else l for l in out] == [o["r1"], o["r1s"], o["r2"], o["r2s"], o["r4"], o["r5"]]
	assert out[0][-1].lower() == "cs:z::4*at:1+g:4-c:5"
	unmapped = _tsv(_pt("sam2paf", "-U", p).stdout)
	assert ["r3", "4", "0", "0", "*", "*", "0", "0", "0", "0", "0", "0"] in unmapped and len(unmapped) == 7
	assert [l[0] + ":" + l[2] for l in _tsv(_pt("sam2paf", "-p", p).stdout)] == ["r1:2", "r1:16", "r2:2", "r4/1:0", "r5:0"]
	assert [l[0] + ":" + l[2] for l in _tsv(_pt("sam2paf", "-P", p).stdout)] == ["r1:2", "r2:2", "r4/1:0", "r5:0"]
	lng = _tsv(_pt("sam2paf", "-L", "-", stdin=open(p).read()).stdout)
	assert lng[0][:-1] == o["r1"] and lng[0][-1].lower() == "cs:z:=acgt*at=c+g=gatc-c=agtca"


def test_delta2paf(k8, tmp_path):
	"""delta2paf turns MUMmer delta records into PAF: 1-based inclusive coordinates become 0-based half-open,
	decreasing query coordinates mean the reverse strand, each positive delta d is d-1 matches then a
	deletion, each negative one |d|-1 matches then an insertion, and column 10 is the block length minus
	the error count (also NM:i). A run of 1s after a delta extends the same gap (here 3M2D7M)."""
	p = _w(tmp_path, "in.delta", ["/r.fa /q.fa", "NUCMER", ">t1 q1 1000 800", "101 200 51 150 2 2 0", "5", "-10", "0",
		"201 300 400 301 0 0 0", "0", ">t2 q2 500 400", "11 20 1 10 1 1 0", "0", "31 42 21 30 2 2 0", "4", "1", "0"])
	assert _tsv(_pt("delta2paf", p).stdout) == [
		["q1", "800", "50", "150", "+", "t1", "1000", "100", "200", "99", "101", "0", "NM:i:2", "cg:Z:4M1D9M1I86M"],
		["q1", "800", "300", "400", "-", "t1", "1000", "200", "300", "100", "100", "0", "NM:i:0", "cg:Z:100M"],
		["q2", "400", "0", "10", "+", "t2", "500", "10", "20", "9", "10", "0", "NM:i:1", "cg:Z:10M"],
		["q2", "400", "20", "30", "+", "t2", "500", "30", "42", "10", "12", "0", "NM:i:2", "cg:Z:3M2D7M"]]
	r = _pt("delta2paf", check=False)
	assert r.returncode == 1 and r.stdout.startswith("Usage: paftools.js delta2paf")


def test_splice2bed(k8, tmp_path):
	"""splice2bed writes BED12 per alignment: blocks split at N operations (M and D extend a block), score
	1000*matches/block length from PAF (1000 from SAM), thick start/end = start/end, colour blue for a
	read's single primary, red when it has several, and secondaries (no s2 in PAF, 0x100 in SAM) only with
	-m, in green. A paired SAM read gets /1 or /2."""
	paf = [
		["q1", 300, 0, 300, "+", "chrA", 10000, 1000, 2350, 290, 300, 60, "tp:A:P", "s1:i:200", "s2:i:0", "cg:Z:100M1000N150M50D50M"],
		["q1", 300, 0, 300, "+", "chrA", 10000, 5000, 5300, 280, 300, 0, "tp:A:S", "s1:i:100", "cg:Z:300M"],
		["q3", 200, 0, 100, "-", "chrB", 9000, 100, 200, 100, 100, 60, "tp:A:P", "s2:i:0", "cg:Z:100M"],
		["q3", 200, 100, 200, "-", "chrB", 9000, 500, 700, 100, 100, 60, "tp:A:P", "s2:i:0", "cg:Z:50M100N50M"],
	]
	p = _w(tmp_path, "s.paf", paf)
	blue, red, green = "0,128,255", "255,0,0", "0,192,0"
	want = [["chrA", "1000", "2350", "q1", "966", "+", "1000", "2350", blue, "2", "100,250,", "0,1100,"],
		["chrB", "100", "200", "q3", "1000", "-", "100", "200", red, "1", "100,", "0,"],
		["chrB", "500", "700", "q3", "1000", "-", "500", "700", red, "2", "50,50,", "0,150,"]]
	assert _tsv(_pt("splice2bed", p).stdout) == want
	sec = ["chrA", "5000", "5300", "q1", "933", "+", "5000", "5300", green, "1", "300,", "0,"]
	assert _tsv(_pt("splice2bed", "-m", p).stdout) == [want[0], sec] + want[1:]
	sam = _w(tmp_path, "s.sam", ["@SQ\tSN:chrA\tLN:10000",
		["s1", 0, "chrA", 501, 60, "10S50M200N40M", "*", 0, 0, "A" * 100, "*"],
		["s2", 83, "chrA", 1001, 60, "30M", "=", 1201, 0, "C" * 30, "*"],
		["s3", 4, "*", 0, 0, "*", "*", 0, 0, "ACGT", "*"]])
	assert _tsv(_pt("splice2bed", sam).stdout) == [
		["chrA", "500", "790", "s1", "1000", "+", "500", "790", blue, "2", "50,40,", "0,250,"],
		["chrA", "1000", "1030", "s2/1", "1000", "-", "1000", "1030", blue, "1", "30,", "0,"]]


GTF = [
	'chr1\tsrc\tgene\t101\t1000\t.\t+\t.\tgene_id "G1"; gene_type "protein_coding"; gene_name "Abc";',
	'chr1\tsrc\ttranscript\t101\t1000\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; transcript_type "protein_coding"; gene_name "Abc"; transcript_name "Abc-201"; tag "Ensembl_canonical";',
	'chr1\tsrc\texon\t101\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; transcript_type "protein_coding"; gene_name "Abc"; transcript_name "Abc-201"; tag "Ensembl_canonical";',
	'chr1\tsrc\tCDS\t151\t200\t.\t+\t0\tgene_id "G1"; transcript_id "T1"; transcript_type "protein_coding"; gene_name "Abc"; transcript_name "Abc-201"; tag "Ensembl_canonical";',
	'chr1\tsrc\texon\t501\t1000\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; transcript_type "protein_coding"; gene_name "Abc"; transcript_name "Abc-201"; tag "Ensembl_canonical";',
	'chr1\tsrc\tCDS\t501\t800\t.\t+\t1\tgene_id "G1"; transcript_id "T1"; transcript_type "protein_coding"; gene_name "Abc"; transcript_name "Abc-201"; tag "Ensembl_canonical";',
	'chr2\tsrc\texon\t2001\t2100\t.\t-\t.\tgene_id "G2"; transcript_id "T2"; transcript_type "lincRNA"; gene_name "Xyz"; transcript_name "Xyz-201";',
	'chr2\tsrc\texon\t1001\t1050\t.\t-\t.\tgene_id "G2"; transcript_id "T2"; transcript_type "lincRNA"; gene_name "Xyz"; transcript_name "Xyz-201";',
]


def test_gff2bed(k8, tmp_path):
	"""gff2bed: one BED12 line per transcript from its exons (sorted), name transcript_id|type|gene, thick
	start/end from the CDS (the transcript when there is none), colour by type; -s names
	transcript_name|type; -j writes one line per intron instead; -e keeps transcripts tagged
	Ensembl_canonical; -G writes gene lines; GFF3 attributes (key=value) work the same way."""
	p = _w(tmp_path, "a.gtf", ["#comment"] + GTF)
	t1 = ["chr1", "100", "1000", "T1|protein_coding|Abc", "1000", "+", "150", "800", "0,128,255", "2", "100,500,", "0,400,"]
	t2 = ["chr2", "1000", "2100", "T2|lincRNA|Xyz", "1000", "-", "1000", "2100", "0,192,0", "2", "50,100,", "0,1000,"]
	assert _tsv(_pt("gff2bed", p).stdout) == [t1, t2]
	assert [l[3] for l in _tsv(_pt("gff2bed", "-s", p).stdout)] == ["Abc-201|protein_coding", "Xyz-201|lincRNA"]
	assert _tsv(_pt("gff2bed", "-j", p).stdout) == [["chr1", "200", "500", "T1|protein_coding|Abc", "1000", "+"],
		["chr2", "1050", "2000", "T2|lincRNA|Xyz", "1000", "-"]]
	assert _tsv(_pt("gff2bed", "-e", p).stdout) == [t1]
	assert _tsv(_pt("gff2bed", "-G", p).stdout) == [["chr1", "100", "1000", "G1|protein_coding|Abc", "1000", "+"]]
	gff3 = _w(tmp_path, "a.gff3", ["##gff-version 3", "chr3\tsrc\tmRNA\t1\t500\t.\t+\t.\tID=tx3;Parent=g3",
		"chr3\tsrc\texon\t11\t100\t.\t+\t.\tParent=tx3;transcript_id=tx3;gene_name=Foo;transcript_biotype=mRNA",
		"chr3\tsrc\texon\t201\t300\t.\t+\t.\tParent=tx3;transcript_id=tx3;gene_name=Foo;transcript_biotype=mRNA"])
	assert _tsv(_pt("gff2bed", gff3).stdout) == [["chr3", "10", "300", "tx3|mRNA|Foo", "1000", "+", "10", "300",
		"0,128,255", "2", "90,100,", "0,190,"]]


def test_gff2junc(k8, tmp_path):
	"""gff2junc writes the gaps between consecutive CDS (or -f feature) records of each Parent as BED6
	`chr end_of_previous start_of_next parent 0 strand`, whatever order the records come in."""
	p = _w(tmp_path, "a.gff3", ["##gff-version 3", "chr1\tsrc\tmRNA\t1\t1000\t.\t+\t.\tID=t1",
		"chr1\tsrc\tCDS\t101\t200\t.\t+\t0\tParent=t1", "chr1\tsrc\tCDS\t401\t500\t.\t+\t2\tParent=t1",
		"chr1\tsrc\tCDS\t701\t800\t.\t+\t1\tParent=t1", "chr2\tsrc\tCDS\t1001\t1100\t.\t-\t0\tParent=t2",
		"chr2\tsrc\tCDS\t501\t600\t.\t-\t2\tParent=t2", "chr2\tsrc\texon\t1\t50\t.\t-\t.\tParent=t3",
		"chr2\tsrc\texon\t101\t150\t.\t-\t.\tParent=t3"])
	assert _tsv(_pt("gff2junc", p).stdout) == [["chr1", "200", "400", "t1", "0", "+"], ["chr1", "500", "700", "t1", "0", "+"],
		["chr2", "600", "1000", "t2", "0", "-"]]
	assert _tsv(_pt("gff2junc", "-f", "exon", p).stdout) == [["chr2", "50", "100", "t3", "0", "-"]]


def test_longcs2seq(k8, tmp_path):
	"""longcs2seq rebuilds the aligned target (default) or query (-q) sequence from a long-form cs, named
	name_start_end; lines without cs are skipped and a short-form cs is an error."""
	p = _w(tmp_path, "l.paf", [VIEW_PAF[:-1] + ["cs:Z:=ACGTA*ag+c=TTGC-t=GGCAT"], VIEW_PAF[:12]])
	assert _pt("longcs2seq", p).stdout == ">t1_5_21\nACGTAATTGCTGGCAT\n"
	assert _pt("longcs2seq", "-q", p).stdout == ">q1_2_18\nACGTAGCTTGCGGCAT\n"
	assert _pt("longcs2seq", _w(tmp_path, "s.paf", [VIEW_PAF]), check=False).returncode != 0


def test_paf2gff(k8, tmp_path):
	"""paf2gff writes a GFF3 transcript (score AS, identity = matches/block length, aa_start/aa_end) and one CDS
	per block between N operations, 1-based; with -a CIGAR lengths are in amino acids (x3) and da:i:0/do:i:0
	add start and stop codons, the stop extending the transcript; a second line of a query is #2."""
	nt = ["q1", 300, 0, 300, "+", "chr1", 10000, 1000, 1600, 290, 300, 60, "AS:i:500", "cg:Z:100M200N100M100N100M"]
	p = _w(tmp_path, "n.paf", [nt, nt[:7] + [5000, 5300, 300, 300, 60, "AS:i:600", "cg:Z:300M"]])
	attr = "transcript_id=q1#%d;transcript_type=protein_coding"
	assert _tsv(_pt("paf2gff", p).stdout) == [
		["chr1", "paf2gff", "transcript", "1001", "1600", "500", "+", ".", attr % 1 + ";identity=0.9667;aa_start=0;aa_end=0"],
		["chr1", "paf2gff", "CDS", "1001", "1100", "0", "+", "0", attr % 1 + ";frameshift=0"],
		["chr1", "paf2gff", "CDS", "1301", "1400", "0", "+", "0", attr % 1 + ";frameshift=0"],
		["chr1", "paf2gff", "CDS", "1501", "1600", "0", "+", "0", attr % 1 + ";frameshift=0"],
		["chr1", "paf2gff", "transcript", "5001", "5300", "600", "+", ".", attr % 2 + ";identity=1.0000;aa_start=0;aa_end=0"],
		["chr1", "paf2gff", "CDS", "5001", "5300", "0", "+", "0", attr % 2 + ";frameshift=0"]]
	aa = _w(tmp_path, "a.paf", [["q2", 100, 0, 100, "+", "chr2", 5000, 100, 400, 95, 100, 0, "AS:i:300", "np:i:98",
		"da:i:0", "do:i:0", "cg:Z:100M"]])
	a2 = "transcript_id=q2#1;transcript_type=protein_coding"
	assert _tsv(_pt("paf2gff", "-a", aa).stdout) == [
		["chr2", "paf2gff", "transcript", "101", "403", "300", "+", ".", a2 + ";identity=0.9500;positive=2.9400;aa_start=0;"
			"aa_end=0;dist_start_codon=0;dist_stop_codon=0"],
		["chr2", "paf2gff", "start_codon", "101", "103", "0", "+", ".", a2 + ";frameshift=0"],
		["chr2", "paf2gff", "CDS", "101", "400", "0", "+", "0", a2 + ";frameshift=0"],
		["chr2", "paf2gff", "stop_codon", "401", "403", "0", "+", ".", a2 + ";frameshift=0"]]
	# on the reverse strand the stop codon lies before the alignment on the reference and the start codon at its end
	rev = _w(tmp_path, "r.paf", [["q3", 100, 0, 100, "-", "chr2", 5000, 100, 400, 95, 100, 0, "AS:i:300", "np:i:98",
		"da:i:0", "do:i:0", "cg:Z:100M"]])
	a3 = "transcript_id=q3#1;transcript_type=protein_coding"
	assert _tsv(_pt("paf2gff", "-a", rev).stdout) == [
		["chr2", "paf2gff", "transcript", "98", "400", "300", "-", ".", a3 + ";identity=0.9500;positive=2.9400;aa_start=0;"
			"aa_end=0;dist_start_codon=0;dist_stop_codon=0"],
		["chr2", "paf2gff", "stop_codon", "98", "100", "0", "-", ".", a3 + ";frameshift=0"],
		["chr2", "paf2gff", "CDS", "101", "400", "0", "-", "0", a3 + ";frameshift=0"],
		["chr2", "paf2gff", "start_codon", "398", "400", "0", "-", ".", a3 + ";frameshift=0"]]
	assert _pt("paf2gff", _w(tmp_path, "noas.paf", [nt[:12] + [nt[13]]]), check=False).returncode != 0


# ---------------------------------------------------------------------------
# Statistics: stat, asmstat, asmgene, misjoin, bedcov, vcfstat
# ---------------------------------------------------------------------------

STAT_PAF = [
	["q1", 1000, 0, 1000, "+", "t1", 10000, 100, 1095, 980, 1015, 60, "tp:A:P", "NM:i:20", "s1:i:900", "s2:i:0", "cg:Z:500M10I300M5D190M"],
	["q1", 1000, 0, 1000, "+", "t1", 10000, 3000, 4000, 900, 1000, 0, "tp:A:S", "NM:i:50", "s1:i:500", "cg:Z:1000M"],
	["q2", 2000, 100, 1900, "-", "t1", 10000, 2000, 3740, 1735, 1800, 60, "tp:A:P", "NM:i:65", "s1:i:1700", "s2:i:10", "cg:Z:800M60I940M"],
	["q2", 2000, 1800, 2000, "+", "t1", 10000, 4000, 5700, 200, 1700, 60, "tp:A:P", "NM:i:1500", "s1:i:150", "s2:i:0", "cg:Z:100M1500D100M"],
	["q3", 500, 0, 0, "*", "*", 0, 0, 0, 0, 0, 0],
]


def _stat_lines(n_seq, n_pri, n_2nd, l_tot, l_cov, n_sub, ins, dels):
	bins = ["[0,50)", "[50,100)", "[100,300)", "[300,400)", "[400,1000)", "[1000,inf)"]
	return (["Number of mapped sequences: %d" % n_seq, "Number of primary alignments: %d" % n_pri,
		"Number of secondary alignments: %d" % n_2nd, "Number of primary alignments with >65535 CIGAR operations: 0",
		"Number of bases in mapped sequences: %d" % l_tot, "Number of mapped bases: %d" % l_cov,
		"Number of substitutions: %d" % n_sub] + ["Number of insertions in %s: %d" % (b, n) for b, n in zip(bins, ins)] +
		["Number of deletions in %s: %d" % (b, n) for b, n in zip(bins, dels)])


def test_stat_paf(k8, tmp_path):
	"""stat on PAF: lines without s2 are secondary; mapped bases are the union of a query's aligned intervals;
	substitutions are NM minus gap lengths; gaps are binned by length. -l N lists gaps of at least N bp
	(query position and target position after the gap); -c prints per alignment: name, query length, mapQ,
	unaligned length, NM, aligned columns plus gaps, mismatches plus gap opens, aligned columns plus gap
	opens."""
	p = _w(tmp_path, "s.paf", STAT_PAF)
	assert _pt("stat", p).stdout.splitlines() == _stat_lines(2, 3, 1, 3000, 2900, 10, [1, 1, 0, 0, 0, 0], [1, 0, 0, 0, 0, 1])
	assert _tsv(_pt("stat", "-l", "50", p).stdout) == [["q2", "860", "-", "t1", "2800", "I", "60"],
		["q2", "100", "+", "t1", "5600", "D", "1500"]]
	assert _tsv(_pt("stat", "-c", p).stdout) == [["q1", "1000", "60", "0", "20", "1005", "7", "992"],
		["q2", "2000", "60", "200", "65", "1800", "6", "1741"], ["q2", "2000", "60", "1800", "1500", "1700", "1", "201"]]
	# gaps on both sides of every bin boundary; -l 400 lists the gaps of 400 bp and longer
	gaps = [49, "I", 99, "D", 50, "I", 100, "D", 300, "I", 399, "D", 999, "I", 400, "D", 1000, "D"]
	cg = "100M" + "".join("%d%s100M" % (gaps[i], gaps[i + 1]) for i in range(0, len(gaps), 2))
	b = _w(tmp_path, "b.paf", [["q9", 2398, 0, 2398, "+", "t9", 5000, 0, 2998, 1000, 4396, 60, "tp:A:P", "NM:i:3396",
		"s1:i:900", "s2:i:0", "cg:Z:" + cg]])
	assert _pt("stat", b).stdout.splitlines() == _stat_lines(1, 1, 0, 2398, 2398, 0, [1, 1, 0, 1, 1, 0], [0, 1, 1, 1, 1, 1])
	assert _tsv(_pt("stat", "-l", "400", b).stdout) == [["q9", "2098", "+", "t9", "1298", "I", "999"],
		["q9", "2198", "+", "t9", "1798", "D", "400"], ["q9", "2298", "+", "t9", "2898", "D", "1000"]]


def test_stat_sam(k8, tmp_path):
	"""stat on SAM: 0x100 records are secondary, the query interval comes from the clips (from the other end on
	the reverse strand) and the read length from CIGAR plus clips; read from stdin with '-'. r2's primary
	covers read bases [10,100) (its trailing soft clip is the start of the read) and its supplementary
	[0,10), so all 100 bases of r2 are mapped."""
	sam = ["@SQ\tSN:t1\tLN:1000", ["r1", 0, "t1", 1, 60, "5S95M", "*", 0, 0, "A" * 100, "*", "NM:i:2"],
		["r1", 256, "t1", 501, 0, "95M5S", "*", 0, 0, "*", "*", "NM:i:3"],
		["r2", 16, "t1", 201, 60, "50M10D40M10S", "*", 0, 0, "C" * 100, "*", "NM:i:13"],
		["r2", 2048, "t1", 701, 60, "10M90H", "*", 0, 0, "ACGTACGTAC", "*", "NM:i:0"],
		["r3", 4, "*", 0, 0, "*", "*", 0, 0, "ACGT", "*"]]
	text = "".join("\t".join(str(x) for x in l) + "\n" if isinstance(l, list) else l + "\n" for l in sam)
	assert _pt("stat", "-", stdin=text).stdout.splitlines() == _stat_lines(2, 3, 1, 200, 195, 5, [0] * 6, [1, 0, 0, 0, 0, 0])


def test_asmstat(k8, tmp_path):
	"""asmstat on a two-contig assembly (12 kbp contig split across t1 and t2, 4 kbp contig with 3 kbp on t1)
	against a 20 kbp reference (.fai), with -l 1000: assembly length, reference bases covered and their
	percentage, duplicated reference fraction, query coverage, NG75/NG50 of contigs, NGA50 and auNGA of
	aligned blocks, breaks and break points; one column per assembly, named after the file without
	.paf(.gz)."""
	fai = _w(tmp_path, "ref.fa.fai", [["t1", 10000, 4, 60, 61], ["t2", 10000, 10173, 60, 61]])
	lines = [["c1", 12000, 0, 6000, "+", "t1", 10000, 0, 6000, 6000, 6000, 60, "tp:A:P", "NM:i:0", "cg:Z:6000M"],
		["c1", 12000, 6000, 12000, "+", "t2", 10000, 2000, 8000, 6000, 6000, 60, "tp:A:P", "NM:i:0", "cg:Z:6000M"],
		["c2", 4000, 0, 3000, "-", "t1", 10000, 6000, 9000, 3000, 3000, 60, "tp:A:P", "NM:i:0", "cg:Z:3000M"]]
	a1 = _w(tmp_path, "asm1.paf", lines)
	a2 = _w(tmp_path, "asm2.paf.gz", lines)
	out = _tsv(_pt("asmstat", "-l", "1000", fai, a1, a2).stdout)
	vals = [("Length", "16000"), ("l_cov", "15000"), ("Rcov", "75.00%"), ("Rdup", "0.00%"), ("Qcov", "93.75%"),
		("NG75", "4000"), ("NG50", "12000"), ("NGA50", "6000"), ("AUNGA", "4050"), ("#breaks", "1"),
		("bp(1000,0)", "1"), ("bp(1000,10k)", "1")]
	assert out == [["Metric", a1[:-4], a2[:-7]]] + [[k, v, v] for k, v in vals]


ASMGENE_REF = [
	["g1", 1000, 0, 1000, "+", "chr1", 100000, 10000, 12000, 1000, 1000, 60],
	["g2", 1000, 0, 1000, "+", "chr2", 100000, 20000, 22000, 1000, 1000, 60],
	["g3", 1000, 0, 1000, "+", "chr3", 100000, 30000, 32000, 1000, 1000, 60],
	["g4", 1000, 0, 1000, "+", "chr4", 100000, 40000, 42000, 1000, 1000, 60],
	["g4", 1000, 0, 1000, "+", "chr4", 100000, 50000, 52000, 1000, 1000, 60],
	["g5", 1000, 0, 1000, "+", "chr5", 100000, 50000, 52000, 1000, 1000, 60],
	["g6", 1000, 0, 1000, "+", "chr6", 100000, 60000, 62000, 1000, 1000, 60],
	["g7", 1000, 0, 1000, "+", "chr7", 100000, 70000, 72000, 1000, 1000, 60],
	["g9", 1000, 0, 1000, "+", "chr9", 100000, 90000, 92000, 1000, 1000, 60],
]
ASMGENE_ASM = [
	["g1", 1000, 0, 1000, "+", "a1", 100000, 10000, 12000, 1000, 1000, 60],
	["g2", 1000, 0, 1000, "+", "a2", 100000, 20000, 22000, 1000, 1000, 60],
	["g2", 1000, 0, 1000, "+", "a2", 100000, 40000, 42000, 1000, 1000, 60],
	["g3", 1000, 0, 600, "+", "a3", 100000, 30000, 31000, 600, 600, 60],
	["g4", 1000, 0, 1000, "+", "a4", 100000, 40000, 42000, 1000, 1000, 60],
	["g4", 1000, 0, 1000, "+", "a4", 100000, 50000, 52000, 1000, 1000, 60],
	["g6", 1000, 0, 500, "+", "a6", 100000, 60000, 61000, 500, 500, 60],
	["g6", 1000, 500, 1000, "+", "a6", 100000, 63000, 64000, 500, 500, 60],
	["g7", 1000, 0, 200, "+", "a7", 100000, 70000, 70500, 200, 200, 60],
	["g9", 1000, 0, 1000, "+", "a9", 100000, 90000, 92000, 900, 1000, 60],
]


def test_asmgene(k8, tmp_path):
	"""asmgene counts, per file, genes that are single-copy in the first (reference) file: full single copy,
	full duplicated, fragmented (covered only by several partial alignments), partial >=50%, >=10%, <10%
	(missing or below -i identity included), and for multi-copy genes how many stay multi-copy and their
	copies. -c sets the full-length cutoff, -i the identity cutoff, -e lists the problem genes."""
	ref = _w(tmp_path, "ref.paf", ASMGENE_REF)
	asm = _w(tmp_path, "asm.paf", ASMGENE_ASM)
	rows = ["full_sgl", "full_dup", "frag", "part50+", "part10+", "part10-", "dup_cnt", "dup_sum"]

	def table(*opts):
		out = _tsv(_pt("asmgene", *(list(opts) + [ref, asm])).stdout)
		return out[0], dict((l[1], (int(l[2]), int(l[3]))) for l in out[1:] if l[0] == "X"), [l for l in out if l[0] not in "HX"]
	head, x, err = table()
	assert head == ["H", "Metric", ref[:-4], asm[:-4]] and sorted(x) == sorted(rows) and err == []
	assert [x[r] for r in rows] == [(7, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 2), (1, 1), (2, 2)]
	_, x, _ = table("-c", "0.5")
	assert [x[r] for r in rows] == [(7, 2), (0, 2), (0, 0), (0, 0), (0, 1), (0, 2), (1, 1), (2, 2)]
	_, x, _ = table("-i", "0.85")
	assert [x[r] for r in rows] == [(7, 2), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (1, 1), (2, 2)]
	_, _, err = table("-e")
	assert sorted((l[0], l[2]) for l in err) == [("0", "g9"), ("1", "g7"), ("5", "g3"), ("D", "g2"), ("F", "g6"), ("M", "g5")]
	assert all(l[1] == asm[:-4] for l in err)


MISJOIN_PAF = [
	["c1", 100000, 0, 30000, "+", "chrA", 1000000, 100000, 130000, 30000, 30000, 60],
	["c1", 100000, 30000, 60000, "+", "chrB", 1000000, 200000, 230000, 30000, 30000, 60],
	["c1", 100000, 60000, 60500, "+", "chrC", 1000000, 1000, 1500, 500, 500, 60],
	["c2", 100000, 0, 30000, "+", "chrA", 1000000, 300000, 330000, 30000, 30000, 60],
	["c2", 100000, 30000, 60000, "+", "chrA", 1000000, 400000, 430000, 30000, 30000, 60],
	["c3", 100000, 0, 30000, "+", "chrA", 1000000, 500000, 530000, 30000, 30000, 60],
	["c3", 100000, 30000, 60000, "-", "chrA", 1000000, 530000, 560000, 30000, 30000, 60],
	["c3", 100000, 60000, 90000, "+", "chrA", 1000000, 560000, 590000, 30000, 30000, 60],
	["c4", 100000, 0, 30000, "+", "chrA", 1000000, 700000, 730000, 30000, 30000, 60],
	["c4", 100000, 30000, 60000, "-", "chrA", 1000000, 730000, 760000, 30000, 30000, 60],
	["c5", 100000, 0, 90000, "+", "chrA", 1000000, 800000, 890000, 90000, 90000, 60],
	# reverse strand without a gap: the second query block lies just before the first on the reference
	["c6", 100000, 0, 30000, "-", "chrA", 1000000, 930000, 960000, 30000, 30000, 60],
	["c6", 100000, 30000, 60000, "-", "chrA", 1000000, 900000, 930000, 30000, 30000, 60],
	# forward then two contiguous reverse blocks: an inversion that runs to the contig end
	["c7", 100000, 0, 30000, "+", "chrA", 1000000, 600000, 630000, 30000, 30000, 60],
	["c7", 100000, 30000, 60000, "-", "chrA", 1000000, 660000, 690000, 30000, 30000, 60],
	["c7", 100000, 60000, 90000, "-", "chrA", 1000000, 630000, 660000, 30000, 30000, 60],
]


def test_misjoin(k8, tmp_path):
	"""misjoin, with blocks of at least -l (here 1k) and gaps above -g (1k): consecutive blocks of a contig on
	different chromosomes are an inter-chromosomal misjoin, on the same strand far apart a gap (on the
	reverse strand the reference runs backwards, so c6 has none), a reversed block between two forward ones
	a closed inversion, a reversed run that reaches the last block (c4, c7) an inversion at the contig end.
	The second number counts those touching -c centromeres; -e prints the misjoins."""
	p = _w(tmp_path, "m.paf", MISJOIN_PAF)
	assert _pt("misjoin", "-l", "1k", "-g", "1k", p).stdout.splitlines() == ["# inter-chromosomal misjoins: 1,0",
		"# intra-chromosomal gaps: 1,0", "# candidate inversions in the middle: 1,0", "# candidate inversions at contig ends: 2,0"]
	cen = _w(tmp_path, "cen.bed", [["chrB", 195000, 235000]])
	assert _pt("misjoin", "-l", "1k", "-g", "1k", "-c", cen, p).stdout.splitlines()[0] == "# inter-chromosomal misjoins: 1,1"
	assert _pt("misjoin", p).stdout.splitlines() == ["# inter-chromosomal misjoins: 0,0", "# intra-chromosomal gaps: 0,0",
		"# candidate inversions in the middle: 0,0", "# candidate inversions at contig ends: 0,0"]
	err = [l for l in _tsv(_pt("misjoin", "-l", "1k", "-g", "1k", "-e", p).stdout) if l[0] != "C" and not l[0].startswith("#")]
	assert [(l[0], l[1], l[6]) for l in err] == [("J", "c1", "chrA"), ("J", "c1", "chrB"), ("G", "c2", "chrA"),
		("G", "c2", "chrA"), ("M", "c3", "chrA"), ("M", "c3", "chrA"), ("M", "c3", "chrA")]
	assert err[0][1:] == [str(x) for x in MISJOIN_PAF[0]]


def test_bedcov(k8, tmp_path):
	"""bedcov counts the bases of the target BED (second file; BED12 by blocks) covered by the merged regions
	of the first: -p prints `F` per target (name columns, length, covered), the totals go to stderr with a
	percentage; -e drops targets overlapping its BED."""
	reg = _w(tmp_path, "reg.bed", [["chr1", 100, 200], ["chr1", 150, 300], ["chr1", 500, 600], ["chr2", 0, 50]])
	tgt = _w(tmp_path, "tgt.bed", [["chr1", 0, 1000], ["chr1", 250, 550], ["chr2", 10, 20], ["chr3", 0, 100],
		["chr1", 50, 650, "tx", 0, "+", 50, 650, 0, 2, "100,100,", "0,450,"]])
	# tx blocks [50,150) and [500,600): 50 + 100 bases covered by the merged regions [100,300) and [500,600)
	r = _pt("bedcov", "-p", reg, tgt)
	assert _tsv(r.stdout) == [["F", "chr1", "0", "1000", "1000", "300"], ["F", "chr1", "250", "550", "300", "100"],
		["F", "chr2", "10", "20", "10", "10"], ["F", "chr1", "50", "650", "tx", "200", "150"]]
	assert "# target bases: 1610" in r.stderr and "# target bases overlapping regions: 560 (34.78%)" in r.stderr
	r = _pt("bedcov", reg, tgt)
	assert r.stdout == "" and "560 (34.78%)" in r.stderr
	excl = _w(tmp_path, "excl.bed", [["chr1", 900, 950]])
	r = _pt("bedcov", "-p", "-e", excl, reg, tgt)
	assert _tsv(r.stdout)[0] == ["F", "chr1", "0", "1000", "0", "0"]
	assert "# target bases: 610" in r.stderr and "260 (42.62%)" in r.stderr


VCF_HEAD = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1"]


def test_vcfstat(k8, tmp_path):
	"""vcfstat counts substitutions (with ts/tv), insertions and deletions binned by length, over every ALT
	allele; symbolic alleles are skipped; the VCF may come on stdin."""
	recs = [["chr1", 100, ".", "A", "G"], ["chr1", 200, ".", "C", "A"], ["chr1", 300, ".", "T", "TA"],
		["chr1", 400, ".", "GAC", "G"], ["chr1", 500, ".", "A", "<DEL>"], ["chr1", 600, ".", "A", "C,AT"],
		["chr1", 700, ".", "A", "A" + "C" * 60], ["chr1", 800, ".", "A" + "G" * 7000, "A"]]
	text = "".join(l + "\n" for l in VCF_HEAD) + "".join("\t".join(map(str, r + [30, "PASS", ".", "GT", "1/1"])) + "\n" for r in recs)
	assert _pt("vcfstat", stdin=text).stdout.splitlines() == ["# substitutions: 3", "ts/tv: 0.500", "# insertions: 3",
		"# 1bp insertions: 2", "# 2bp insertions: 0", "# [3,50) insertions: 0", "# [50,1000) insertions: 1",
		"# [1000,7000) insertions: 0", "# >=7000 insertions: 0", "# deletions: 2", "# 1bp deletions: 0",
		"# 2bp deletions: 1", "# [3,50) deletions: 0", "# [50,1000) deletions: 0", "# [1000,7000) deletions: 0",
		"# >=7000 deletions: 1"]


# ---------------------------------------------------------------------------
# liftover, call, vcfpair, sveval, vcfsel
# ---------------------------------------------------------------------------

LIFT_PAF = [
	["q1", 1000, 100, 900, "+", "t1", 5000, 1000, 1802, 790, 802, 60, "tp:A:P", "cg:Z:300M2D500M"],
	["q1", 1000, 100, 900, "+", "t1", 5000, 3000, 3800, 790, 800, 60, "tp:A:S", "cg:Z:800M"],
	["q2", 1000, 100, 900, "-", "t1", 5000, 3000, 3800, 800, 800, 60, "tp:A:P", "cg:Z:800M"],
	["q3", 1000, 0, 1000, "+", "t2", 5000, 0, 1000, 1000, 1000, 3, "tp:A:P", "cg:Z:1000M"],
	["q4", 1000, 0, 1000, "+", "t3", 5000, 0, 990, 990, 1000, 60, "tp:A:P", "cg:Z:500M10I490M"],
]


@pytest.mark.parametrize("cmd", ["liftover", "liftOver"])
def test_liftover(k8, tmp_path, cmd):
	"""liftover maps query-coordinate BED intervals through primary alignments of at least -l bp (default 50
	kbp) and mapQ -q (5): bases inside the alignment move with the CIGAR (deletions shift later bases), a
	reverse-strand alignment flips the interval, an end outside the alignment is clipped to it and marked
	_t5/_t3, a base inside an insertion maps to the insertion point; -d drops alignments above that
	gap-compressed divergence (mismatches plus gap opens over that plus matches: 11/801 for q1, 1/991 for
	q4). liftOver is the same command."""
	paf = _w(tmp_path, "l.paf", LIFT_PAF)
	bed = _w(tmp_path, "q.bed", [["q1", 50, 150], ["q1", 200, 250], ["q1", 450, 460], ["q2", 200, 250], ["q3", 10, 20],
		["q4", 495, 505]])
	assert _pt(cmd, paf, bed).stdout == ""
	want = [["t1", "1000", "1050", "q1_50_150_t5", "0", "+"], ["t1", "1100", "1150", "q1_200_250", "0", "+"],
		["t1", "1352", "1362", "q1_450_460", "0", "+"], ["t1", "3650", "3700", "q2_200_250", "0", "-"]]
	q3, q4 = ["t2", "10", "20", "q3_10_20", "0", "+"], ["t3", "495", "501", "q4_495_505", "0", "+"]
	assert _tsv(_pt(cmd, "-l", "500", paf, bed).stdout) == want + [q4]
	assert _tsv(_pt(cmd, "-l", "500", "-q", "2", paf, bed).stdout) == want + [q3, q4]
	assert _tsv(_pt(cmd, "-l", "500", "-d", "0.001", paf, bed).stdout) == want[3:]
	assert _tsv(_pt(cmd, "-l", "500", "-d", "0.014", paf, bed).stdout) == want + [q4]
	assert _tsv(_pt(cmd, "-l", "500", "-d", "0.0137", paf, bed).stdout) == want[3:] + [q4]
	assert _tsv(_pt(cmd, "-l", "500", "-d", "0.05", "-", bed, stdin=open(paf).read()).stdout) == want + [q4]


CALL_PAF = [
	["qA", 20001, 0, 20001, "+", "chr1", 100000, 1000, 21000, 19996, 20002, 60, "tp:A:P", "cs:Z::100*ag:100+tt:100-c:19698"],
	["qB", 25000, 0, 25000, "+", "chr1", 100000, 15000, 40000, 24998, 25000, 60, "tp:A:P", "cs:Z::1000*ct:13999*ga:9999"],
	["qC", 1000, 0, 1000, "-", "chr2", 50000, 5000, 6000, 996, 1002, 60, "tp:A:P", "cs:Z::10*ac:100*tc:100+ga:100-tt:686"],
]


def test_call(k8, tmp_path):
	"""call on coordinate-sorted asm-to-ref PAF with cs (misc/README.md): R lines give regions covered by
	exactly one contig, V lines give each variant as chr, start, end, depth, mapQ, REF, ALT, query, query
	start, end, strand (query positions on the forward query, so they count down along a reverse-strand
	contig); depth counts the alignments over the variant. The summary on stderr counts depth-1 variants.
	Alignments of at least -l count toward depth and the single-coverage regions; only those of at least
	-L give variants and R lines."""
	p = _w(tmp_path, "c.paf", CALL_PAF)
	r = _pt("call", "-l", "1000", "-L", "1000", p)
	assert _tsv(r.stdout) == [
		["R", "chr1", "1000", "15000"],
		["V", "chr1", "1100", "1101", "1", "60", "a", "g", "qA", "100", "101", "+"],
		["V", "chr1", "1201", "1201", "1", "60", "-", "tt", "qA", "201", "203", "+"],
		["V", "chr1", "1301", "1302", "1", "60", "c", "-", "qA", "303", "303", "+"],
		["R", "chr1", "21000", "40000"],
		["V", "chr1", "16000", "16001", "2", "60", "c", "t", "qB", "1000", "1001", "+"],
		["V", "chr1", "30000", "30001", "1", "60", "g", "a", "qB", "15000", "15001", "+"],
		["R", "chr2", "5000", "6000"],
		["V", "chr2", "5010", "5011", "1", "60", "a", "c", "qC", "989", "990", "-"],
		["V", "chr2", "5111", "5112", "1", "60", "t", "c", "qC", "888", "889", "-"],
		["V", "chr2", "5212", "5212", "1", "60", "-", "ga", "qC", "786", "788", "-"],
		["V", "chr2", "5312", "5314", "1", "60", "tt", "-", "qC", "686", "686", "-"]]
	err = r.stderr.splitlines()
	assert "34000 reference bases covered by exactly one contig" in err
	assert "4 substitutions; ts/tv = 3.000" in err and "1 1bp deletions" in err and "2 2bp insertions" in err
	assert "1 2bp deletions" in err and "0 1bp insertions" in err
	assert _pt("call", p).stdout == ""  # default -L 50000 / -l 10000: every alignment is too short
	# -L 20500: qA (20002 columns) and qC still count toward depth and coverage but give no variants or R lines
	r = _pt("call", "-l", "1000", "-L", "20500", p)
	assert _tsv(r.stdout) == [["R", "chr1", "21000", "40000"],
		["V", "chr1", "16000", "16001", "2", "60", "c", "t", "qB", "1000", "1001", "+"],
		["V", "chr1", "30000", "30001", "1", "60", "g", "a", "qB", "15000", "15001", "+"]]
	assert "19000 reference bases covered by exactly one contig" in r.stderr.splitlines()


def test_call_vcf(k8, tmp_path):
	"""call -f ref.fa writes VCF: a header with one ##contig per reference sequence and the -s sample name,
	then depth-1 variants with 1-based POS, indels anchored on the preceding reference base, and QNAME,
	QSTART (1-based) and QSTRAND in INFO."""
	rng = H.Rng(77)
	chr1, chr2 = rng.dna(100000), rng.dna(50000)
	ref = H.write_fasta(str(tmp_path / "ref.fa"), [("chr1", chr1), ("chr2", chr2)], line_width=60)
	p = _w(tmp_path, "c.paf", CALL_PAF[:2])
	out = _pt("call", "-l", "1000", "-L", "1000", "-f", ref, "-s", "asm1", p).stdout.splitlines()
	assert out[:3] == ["##fileformat=VCFv4.1", "##contig=<ID=chr1,length=100000>", "##contig=<ID=chr2,length=50000>"]
	assert out[3:7] == ['##INFO=<ID=QNAME,Number=1,Type=String,Description="Query name">',
		'##INFO=<ID=QSTART,Number=1,Type=Integer,Description="Query start">',
		'##INFO=<ID=QSTRAND,Number=1,Type=String,Description="Query strand">',
		'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">']
	assert out[7] == "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tasm1"
	a, b = chr1[1200].upper(), chr1[1300].upper()
	assert _tsv("\n".join(out[8:])) == [
		["chr1", "1101", ".", "A", "G", "60", ".", "QNAME=qA;QSTART=101;QSTRAND=+", "GT", "1/1"],
		["chr1", "1201", ".", a, a + "TT", "60", ".", "QNAME=qA;QSTART=202;QSTRAND=+", "GT", "1/1"],
		["chr1", "1301", ".", b + "C", b, "60", ".", "QNAME=qA;QSTART=304;QSTRAND=+", "GT", "1/1"],
		["chr1", "30001", ".", "G", "A", "60", ".", "QNAME=qB;QSTART=15001;QSTRAND=+", "GT", "1/1"]]


def test_vcfpair(k8, tmp_path):
	"""vcfpair merges a two-haplotype VCF into one phased sample: GT h1|h2 (an uncalled or heterozygous
	haplotype becomes '.' and adds GAP1/GAP2 or HET1/HET2 to FILTER), AD summed, QUAL 30; only autosomes and
	X (and Y with -m) are kept; ##source lines are dropped and FILTER definitions added; -s names the
	sample. With -m -g 37, an uncalled haplotype on haploid X (outside the PARs) or Y is not filtered;
	-m without -g is an error."""
	head = ["##fileformat=VCFv4.2", "##source=test", "##contig=<ID=chr1,length=1000>", "##contig=<ID=chrUn_x,length=10>",
		"##contig=<ID=chrX,length=156040895>", "##contig=<ID=chrY,length=57227415>",
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tH1\tH2"]
	recs = [["chr1", 100, ".", "A", "G", 50, "PASS", ".", "GT:AD", "1/1:0,10", "1/1:0,12"],
		["chr1", 200, ".", "C", "T", 50, "PASS", ".", "GT:AD", "0/0:10,0", "1/1:0,8"],
		["chr1", 300, ".", "G", "GA", 50, "PASS", ".", "GT:AD", "./.:0,0", "1/1:0,9"],
		["chr1", 400, ".", "T", "C", 50, "PASS", ".", "GT:AD", "0/1:5,5", "1/1:0,7"],
		["chrUn_x", 5, ".", "A", "C", 50, "PASS", ".", "GT:AD", "1/1:0,3", "1/1:0,3"],
		["chrX", 1000, ".", "A", "G", 50, "PASS", ".", "GT:AD", "./.:0,0", "1/1:0,5"],
		["chrX", 3000000, ".", "A", "G", 50, "PASS", ".", "GT:AD", "./.:0,0", "1/1:0,5"],
		["chrY", 500, ".", "C", "T", 50, "PASS", ".", "GT:AD", "1/1:0,4", "./.:0,0"]]
	p = _w(tmp_path, "pair.vcf", head + recs)
	filt = ['##FILTER=<ID=HET1,Description="Heterozygous in the first haplotype">',
		'##FILTER=<ID=HET2,Description="Heterozygous in the second haplotype">',
		'##FILTER=<ID=GAP1,Description="Uncalled in the first haplotype">',
		'##FILTER=<ID=GAP2,Description="Uncalled in the second haplotype">']
	out = _pt("vcfpair", p).stdout.splitlines()
	assert out[:7] == ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=1000>", "##contig=<ID=chrX,length=156040895>",
		head[6]] + filt[:3]
	assert out[7:9] == [filt[3], "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsyndip"]
	body = [["chr1", "100", ".", "A", "G", "30", ".", ".", "GT:AD", "1|1:0,22"],
		["chr1", "200", ".", "C", "T", "30", ".", ".", "GT:AD", "0|1:10,8"],
		["chr1", "300", ".", "G", "GA", "30", "GAP1", ".", "GT:AD", ".|1:0,9"],
		["chr1", "400", ".", "T", "C", "30", "HET1", ".", "GT:AD", ".|1:5,12"],
		["chrX", "1000", ".", "A", "G", "30", "GAP1", ".", "GT:AD", ".|1:0,5"],
		["chrX", "3000000", ".", "A", "G", "30", "GAP1", ".", "GT:AD", ".|1:0,5"]]
	assert _tsv("\n".join(out[9:])) == body
	male = _pt("vcfpair", "-m", "-g", "37", "-s", "NA24385", p).stdout.splitlines()
	assert "##contig=<ID=chrY,length=57227415>" in male and male[9].endswith("\tNA24385")
	assert _tsv("\n".join(male[10:])) == body[:5] + [body[5][:6] + ["."] + body[5][7:],
		["chrY", "500", ".", "C", "T", "30", ".", ".", "GT:AD", "1|.:0,4"]]
	assert _pt("vcfpair", "-m", p, check=False).returncode != 0


SV_BASE = [["chr1", 1000, ".", "A", "<DEL>", 30, "PASS", "SVTYPE=DEL;SVLEN=-100;END=1100"],
	["chr1", 5000, ".", "A", "A" + "C" * 200, 30, "PASS", "."],
	["chr1", 20000, ".", "A", "<DEL>", 30, "PASS", "SVTYPE=DEL;SVLEN=-300;END=20300"],
	["chr1", 30000, ".", "A", "<INV>", 30, "PASS", "SVTYPE=INV;END=31000"],
	["chr1", 40000, ".", "A", "<DEL>", 30, "LowQual", "SVTYPE=DEL;SVLEN=-500"],
	["chr1", 50000, ".", "A", "<DEL>", 30, "PASS", "SVTYPE=DEL;SVLEN=-40"],
	["chr1", 70000, ".", "A", "<DEL>", 30, "PASS", "SVTYPE=DEL;SVLEN=-1000;END=71000"]]
SV_CALL = [["chr1", 1050, ".", "A", "<DEL>", 30, "PASS", "SVTYPE=DEL;SVLEN=-110;END=1160"],
	["chr1", 5100, ".", "A", "A" + "G" * 150, 30, "PASS", "."],
	["chr1", 60000, ".", "A", "<DEL>", 30, "PASS", "SVTYPE=DEL;SVLEN=-80;END=60080"],
	["chr1", 71700, ".", "A", "<DEL>", 30, "PASS", "SVTYPE=DEL;SVLEN=-1000;END=72700"]]


def test_sveval(k8, tmp_path):
	"""sveval compares PASS SVs of at least -i (50) bp: a base SV is found when a call of the same type lies
	within the window (500 bp plus half its length) and its length differs by less than 50 bp or by less
	than half; SN = found/base, PC = matched calls/calls, F1 their mean. Inversions, breakends and filtered
	records are ignored; -e lists the false negatives and positives. The 1 kbp deletions at 70000 and
	71700 are 700 bp apart, inside the 500 + 500 bp window."""
	base = _w(tmp_path, "base.vcf", VCF_HEAD[:1] + ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"] + SV_BASE)
	call = _w(tmp_path, "call.vcf", VCF_HEAD[:1] + ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"] + SV_CALL)
	assert _tsv(_pt("sveval", base, call).stdout) == [["SN", "4", "3", "0.750000"], ["PC", "4", "3", "0.750000"],
		["F1", "0.750000"]]
	out = _tsv(_pt("sveval", "-e", base, call).stdout)
	assert out[:2] == [["FN", "chr1", "19999", "20300", "-300"], ["FP", "chr1", "59999", "60080", "-80"]]
	assert _tsv(_pt("sveval", "-i", "120", base, call).stdout) == [["SN", "3", "2", "0.666667"], ["PC", "2", "2", "1.000000"],
		["F1", "0.833333"]]


def test_vcfsel(k8, tmp_path):
	"""vcfsel keeps header lines and the records whose largest allele length change is at least -l and whose
	smallest is at most -L (SVLEN counts for symbolic alleles)."""
	recs = [["chr1", 100, ".", "A", "G", 30, "PASS", "."], ["chr1", 200, ".", "A", "A" + "T" * 60, 30, "PASS", "."],
		["chr1", 300, ".", "A", "<DEL>", 30, "PASS", "SVLEN=-500;END=800"], ["chr1", 400, ".", "A", "C,ACGT", 30, "PASS", "."]]
	p = _w(tmp_path, "s.vcf", VCF_HEAD[:1] + ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"] + recs)

	def pos(*opts):
		lines = _pt("vcfsel", *(list(opts) + [p])).stdout.splitlines()
		assert lines[:2] == [VCF_HEAD[0], "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
		return [int(l.split("\t")[1]) for l in lines[2:]]
	assert pos() == [100, 200, 300, 400]
	assert pos("-l", "50") == [200, 300]
	assert pos("-L", "10") == [100, 400]
	assert pos("-l", "50", "-L", "100") == [200]


# ---------------------------------------------------------------------------
# Evaluation: mapeval, pafcmp, junceval, exoneval, ov-eval
# ---------------------------------------------------------------------------

MAPEVAL_PAF = [
	["S1_1!chr1!1000!2000!+", 1000, 0, 1000, "+", "chr1", 100000, 1000, 2000, 990, 1000, 60, "tp:A:P", "s1:i:900", "s2:i:0"],
	["S1_2!chr1!5000!6000!-", 1000, 0, 1000, "-", "chr1", 100000, 5100, 6100, 990, 1000, 70, "tp:A:P", "s1:i:900", "s2:i:0"],
	["S1_3!chr2!100!900!+", 800, 0, 800, "+", "chr1", 100000, 100, 900, 800, 800, 30, "tp:A:P", "s1:i:700", "s2:i:600"],
	["S1_4!chr1!3000!4000!+", 1000, 0, 1000, "-", "chr1", 100000, 3000, 4000, 990, 1000, 10, "tp:A:P", "s1:i:900", "s2:i:800"],
	["S1_5!chr1!7000!8000!+", 1000, 0, 1000, "+", "chr1", 100000, 7000, 8000, 990, 1000, 10, "tp:A:P", "s1:i:900", "s2:i:800"],
	["S1_6!chr1!9000!9500!+", 500, 0, 500, "+", "chr1", 100000, 20000, 20500, 490, 500, 0, "tp:A:S", "s1:i:400"],
]


def test_mapeval_paf(k8, tmp_path):
	"""mapeval on PAF of reads named by pbsim2fq (name!chr!start!end!strand): a mapping is correct on the true
	chromosome and strand with overlap/union above -r (0.1). Q lines, from the highest mapQ down, give the
	threshold, mapped reads, wrong ones, cumulative error rate and cumulative mapped reads, with a new line
	wherever errors appear; mapQ is capped at 60, secondaries (s1 without s2) are ignored. -Q prints wrong
	mappings at or above that mapQ."""
	p = _w(tmp_path, "m.paf", MAPEVAL_PAF)
	assert _tsv(_pt("mapeval", p).stdout) == [["Q", "60", "2", "0", "0.000000000", "2"],
		["Q", "30", "1", "1", "0.333333333", "3"], ["Q", "10", "2", "1", "0.400000000", "5"]]
	assert _tsv(_pt("mapeval", "-r", "0.9", p).stdout) == [["Q", "60", "2", "1", "0.500000000", "2"],
		["Q", "30", "1", "1", "0.666666667", "3"], ["Q", "10", "2", "1", "0.600000000", "5"]]
	out = _tsv(_pt("mapeval", "-Q", "20", p).stdout)
	assert out[0] == ["E", "S1_3!chr2!100!900!+", "chr1", "100", "900", "+", "30", "800"] and out[1][0] == "Q"


def test_mapeval_sam_pairs_and_modes(k8, tmp_path):
	"""mapeval on SAM with mason2fq pair names (name!chr!s1_s2!e1_e2!strands/N, mate from the flag): each mate
	is checked against its own interval, unmapped records give the U line. -m 0 judges a read by its longest
	alignment (matched bases), -m 1 by its first, -m 2 by all its primary alignments: wrong at the highest
	mapQ of its wrong ones if any is wrong (S1_8: mapQ 0), otherwise correct at its highest mapQ."""
	sam = ["@SQ\tSN:chr1\tLN:100000",
		["p1!chr1!100_400!200_500!+-", 99, "chr1", 101, 60, "100M", "=", 401, 400, "A" * 100, "*"],
		["p1!chr1!100_400!200_500!+-", 147, "chr1", 401, 60, "100M", "=", 101, -400, "A" * 100, "*"],
		["p2!chr1!1000_1300!1100_1400!+-", 77, "*", 0, 0, "*", "*", 0, 0, "A" * 100, "*"]]
	p = _w(tmp_path, "m.sam", sam)
	assert _tsv(_pt("mapeval", p).stdout) == [["Q", "60", "2", "0", "0.000000000", "2"], ["U", "1"]]
	two = _w(tmp_path, "two.paf", [["S1_7!chr1!1000!2000!+", 1000, 0, 100, "+", "chr1", 100000, 50000, 50100, 100, 100, 60,
		"tp:A:P", "s1:i:90", "s2:i:0"], ["S1_7!chr1!1000!2000!+", 1000, 100, 1000, "+", "chr1", 100000, 1100, 2000, 900, 900,
		5, "tp:A:P", "s1:i:800", "s2:i:0"]])
	assert _tsv(_pt("mapeval", "-m", "0", two).stdout) == [["Q", "5", "1", "0", "0.000000000", "1"]]
	assert _tsv(_pt("mapeval", "-m", "1", two).stdout) == [["Q", "60", "1", "1", "1.000000000", "1"]]
	three = _w(tmp_path, "three.paf", [open(two).read().rstrip("\n"),
		["S1_8!chr1!3000!4000!+", 1000, 0, 1000, "+", "chr1", 100000, 3000, 4000, 1000, 1000, 20, "tp:A:P", "s1:i:900", "s2:i:0"],
		["S1_8!chr1!3000!4000!+", 1000, 0, 1000, "+", "chr1", 100000, 80000, 81000, 990, 1000, 0, "tp:A:P", "s1:i:890",
			"s2:i:0"]])
	assert _tsv(_pt("mapeval", "-m", "2", three).stdout) == [["Q", "60", "1", "1", "1.000000000", "1"],
		["Q", "0", "1", "1", "1.000000000", "2"]]


def test_pafcmp(k8, tmp_path):
	"""pafcmp compares a test PAF with a base PAF for queries of at least 5 kbp with one non-secondary base
	line: a test mapping overlapping the base one by half the union is a hit, otherwise W (if the base
	has mapQ >= -q); base queries absent from the test are M; X lines summarize."""
	def l(q, ql, chrom, st, en, mapq, tp="P"):
		return [q, ql, 0, ql, "+", chrom, 1000000, st, en, ql, ql, mapq, "tp:A:" + tp]
	base = _w(tmp_path, "base.paf", [l("r1", 10000, "chr1", 1000, 11000, 60), l("r2", 10000, "chr1", 50000, 60000, 60),
		l("r3", 10000, "chr2", 100, 10100, 60), l("r4", 3000, "chr1", 0, 3000, 60), l("r5", 10000, "chr1", 0, 5000, 60),
		l("r5", 10000, "chr1", 8000, 13000, 60), l("r6", 10000, "chr3", 0, 10000, 5), l("r6", 10000, "chr4", 0, 10000, 0, "S"),
		l("r8", 10000, "chr1", 100000, 110000, 60)])
	# r1 overlaps its base mapping by 9000 of 11000 bases (a hit), r8 by 2000 of 18000 (wrong)
	test = _w(tmp_path, "test.paf", [l("r1", 10000, "chr1", 2000, 12000, 60), l("r2", 10000, "chr1", 500000, 510000, 60),
		l("r6", 10000, "chr3", 400000, 410000, 60), l("r7", 10000, "chr1", 0, 10000, 60),
		l("r8", 10000, "chr1", 108000, 118000, 60)])
	out = _tsv(_pt("pafcmp", base, test).stdout)
	assert out[:3] == [["W", "r2", "chr1", "50000", "60000", "60", "chr1", "500000", "510000", "60"],
		["W", "r8", "chr1", "100000", "110000", "60", "chr1", "108000", "118000", "60"], ["M", "r3", "chr2", "100", "10100", "60"]]
	assert [x[1] for x in out[3:7]] == ["4 base alignments with mapQ>=10", "1 base alignments correctly mapped by test",
		"2 wrong test alignment", "1 base alignments missing"]
	out = _tsv(_pt("pafcmp", "-q", "1", base, test).stdout)
	assert [x[0] + x[1] for x in out[:4]] == ["Wr2", "Wr6", "Wr8", "Mr3"] and out[4][1] == "5 base alignments with mapQ>=1"


JUNC_GTF = [
	'chr1\tsrc\texon\t101\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";',
	'chr1\tsrc\texon\t301\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1";',
	'chr1\tsrc\texon\t601\t700\t.\t+\t.\tgene_id "G1"; transcript_id "T1";',
	'chr1\tsrc\texon\t101\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T2";',
	'chr1\tsrc\texon\t501\t700\t.\t+\t.\tgene_id "G1"; transcript_id "T2";',
]
JUNC_SAM = [
	"@SQ\tSN:chr1\tLN:10000", "@SQ\tSN:chr2\tLN:10000",
	["q1", 0, "chr1", 151, 60, "50M100N100M200N50M", "*", 0, 0, "*", "*"],
	["q1", 256, "chr1", 151, 0, "50M100N100M200N50M", "*", 0, 0, "*", "*"],
	["q2", 0, "chr1", 151, 60, "50M298N50M", "*", 0, 0, "*", "*"],
	["q3", 0, "chr1", 301, 60, "100M", "*", 0, 0, "*", "*"],
	["q4", 0, "chr2", 101, 60, "10M50N10M", "*", 0, 0, "*", "*"],
	["q5", 4, "*", 0, 0, "*", "*", 0, 0, "*", "*"],
	["q6", 0, "chr1", 1001, 60, "10M100N10M", "*", 0, 0, "*", "*"],
]


def test_junceval(k8, tmp_path):
	"""junceval compares the introns of primary alignments with the annotated ones (GTF exons; -l tolerance):
	counts of unmapped and mapped reads, alignments, unspliced alignments, predicted introns, introns
	overlapping no annotated intron, and exact matches. -p lists each overlapping or novel intron with the
	annotated introns it overlaps (C correct, P partial, N novel); -e only the wrong ones."""
	gtf = _w(tmp_path, "a.gtf", JUNC_GTF)
	sam = _w(tmp_path, "a.sam", JUNC_SAM)
	summary = ["# unmapped reads: 1", "# mapped reads: 5", "# primary alignments: 5", "# singletons: 1",
		"# predicted introns: 5", "# non-overlapping introns: 2"]
	assert _pt("junceval", gtf, sam).stdout.splitlines() == summary + ["# correct introns: 2 (40.00%)"]
	assert _pt("junceval", "-l", "2", gtf, sam).stdout.splitlines() == summary + ["# correct introns: 3 (60.00%)"]
	assert _tsv(_pt("junceval", "-p", gtf, sam).stdout) == [
		["C", "q1", "1", "chr1", "200", "300", "[(200,300), (200,500)]"],
		["C", "q1", "2", "chr1", "400", "600", "[(200,500), (400,600)]"],
		["P", "q2", "1", "chr1", "200", "498", "[(200,300), (200,500), (400,600)]"],
		["N", "q6", "1", "chr1", "1010", "1110"]]
	assert [l[0] for l in _tsv(_pt("junceval", "-e", gtf, sam).stdout)] == ["P", "N"]


def test_exoneval(k8, tmp_path):
	"""exoneval compares the exons (blocks between N operations) of primary alignments with annotated exons:
	predicted exons, exons overlapping none, exact matches; -s adds base-level sensitivity (annotated bases
	covered) and specificity (predicted bases annotated) over merged exons."""
	gtf = _w(tmp_path, "a.gtf", JUNC_GTF)
	sam = _w(tmp_path, "a.sam", JUNC_SAM)
	summary = ["# unmapped reads: 1", "# mapped reads: 5", "# primary alignments: 5", "# predicted exons: 10",
		"# non-overlapping exons: 4", "# correct exons: 2 (20.00%)"]
	assert _pt("exoneval", gtf, sam).stdout.splitlines() == summary
	assert _pt("exoneval", "-s", gtf, sam).stdout.splitlines() == summary + ["Base Sn: 248 / 400 = 62.00%",
		"Base Sp: 248 / 270 = 91.85%"]


def test_ov_eval(k8, tmp_path):
	"""ov-eval infers read overlaps of at least -l (2000) bp from sorted primary read-to-reference mappings
	(mapQ >= 10, mapped over >= 95% of the read) and reports how many the overlapper's PAF misses."""
	def m(q, st, mapq=60, tp="P"):
		return [q, 10000, 0, 10000, "+", "chr1", 1000000, st, st + 10000, 9900, 10000, mapq, "tp:A:" + tp]
	ref = _w(tmp_path, "ref.paf", [m("r1", 0), m("r2", 5000), m("r3", 14000), m("r4", 16000), m("r5", 17000, 5),
		m("r6", 18000, 60, "S")])
	ov = _w(tmp_path, "ov.paf", [["r2", 10000, 0, 5000, "+", "r1", 10000, 5000, 10000, 4900, 5000, 0]])
	assert _pt("ov-eval", ref, ov).stdout.splitlines() == ["2 overlaps inferred from the reference mapping",
		"1 missed by the read overlapper", "50.00% sensitivity"]
	assert _pt("ov-eval", "-l", "6000", "-", ov, stdin=open(ref).read()).stdout.splitlines() == [
		"1 overlaps inferred from the reference mapping", "1 missed by the read overlapper", "0.00% sensitivity"]


# ---------------------------------------------------------------------------
# Simulated-read converters: mason2fq, sim2bed, pbsim2fq, badread2fa
# ---------------------------------------------------------------------------


def test_mason2fq(k8, tmp_path):
	"""mason2fq turns mason2's SAM into FASTQ named name!chr!start!end!strand (pairs: starts, ends and strands of
	both mates joined, /1 and /2), with the XE:XS:XI counts as comment; reverse-strand reads are turned back
	into their sequenced orientation (sequence reverse-complemented, qualities reversed)."""
	sam = _w(tmp_path, "m.sam", ["@SQ\tSN:chr1\tLN:1000",
		["simulated.1", 99, "chr1", 101, 60, "10M", "=", 201, 110, "ACGTACGTAC", "IIIIIIIIII", "XE:i:1", "XS:i:0", "XI:i:0"],
		["simulated.1", 147, "chr1", 201, 60, "10M", "=", 101, -110, "GGGGCCCCAA", "ABCDEFGHIJ", "XE:i:2", "XS:i:1", "XI:i:0"],
		["simulated.2", 0, "chr1", 301, 60, "3M2D2M", "*", 0, 0, "ACGTT", "IIIII", "XE:i:0"]])
	assert _pt("mason2fq", sam).stdout.splitlines() == [
		"@1!chr1!100_200!110_210!+-/1 1:0:0", "ACGTACGTAC", "+", "IIIIIIIIII",
		"@1!chr1!100_200!110_210!+-/2 2:1:0", "TTGGGGCCCC", "+", "JIHGFEDCBA",
		"@2!chr1!300!307!+ 0:0:0", "ACGTT", "+", "IIIII"]


def test_sim2bed(k8, tmp_path):
	"""sim2bed turns simulated read names into BED6: a mason2 pair name gives the interval and strand of its own
	mate, a pbsim/badread name its interval (start and end ordered); other lines are skipped."""
	names = ["1!chr1!100_200!110_210!+-/1", "1!chr1!100_200!110_210!+-/2", "S1_5!chr2!5000!4000!-", "junk"]
	p = _w(tmp_path, "names.txt", names)
	assert _tsv(_pt("sim2bed", p).stdout) == [["chr1", "100", "110", names[0], "0", "+"],
		["chr1", "200", "210", names[1], "0", "-"], ["chr2", "4000", "5000", names[2], "0", "-"]]


def test_pbsim2fq(k8, tmp_path):
	"""pbsim2fq turns PBSIM MAF into FASTA named read!chr!start!end!strand, the chromosome taken from the .fai
	by the read's S<n>_ prefix and the interval from the reference row; gaps are removed, reverse-strand
	reads reverse-complemented, and reads containing NN dropped."""
	fai = _w(tmp_path, "ref.fa.fai", [["chr1", 1000, 6, 60, 61], ["chr2", 500, 1030, 60, 61]])
	maf = _w(tmp_path, "r.maf", ["a score=0", "s ref 100 20 + 1000 ACGTACGTAC-GTACGTACGT",
		"s S1_1 0 21 + 21 ACGTACGTACAGTACGTACGT", "", "a score=0", "s ref 0 10 + 500 ACGTTTTTCC",
		"s S2_1 0 10 - 10 ACGTTTTTCC", "", "a score=0", "s ref 20 10 + 1000 ACGTACGTAC", "s S1_2 0 10 + 10 ACGNNCGTAC", ""])
	assert _pt("pbsim2fq", fai, maf).stdout.splitlines() == [">S1_1!chr1!100!120!+", "ACGTACGTACAGTACGTACGT",
		">S2_1!chr2!0!10!-", "GGAAAAACGT"]


def test_badread2fa(k8, tmp_path):
	"""badread2fa turns Badread FASTQ into FASTA named S<n>!chr!start!end!strand with ri:f read identity,
	numbering every read; reverse-strand intervals are converted to forward coordinates with the .fai
	length; chimeras and junk reads are dropped and counted on stderr."""
	fai = _w(tmp_path, "ref.fa.fai", [["chr1", 10000], ["chr2", 5000]])
	fq = _w(tmp_path, "b.fq", ["@read1 chr1,+strand,100-600 length=500 error-free_length=500 read_identity=95.5%", "ACGT", "+", "IIII",
		"@read2 chr2,-strand,1000-1500 length=500 error-free_length=500 read_identity=90%", "GGCC", "+", "IIII",
		"@read3 chr1,+strand,1-100 chimera length=99 read_identity=80%", "TTTT", "+", "IIII",
		"@read4 random_seq length=4 read_identity=50%", "AAAA", "+", "IIII",
		"@read5 chr1,+strand,0-50 length=50 error-free_length=50 read_identity=99.1%", "CCAA", "+", "IIII"])
	r = _pt("badread2fa", fai, fq)
	assert r.stdout.splitlines() == [">S1!chr1!100!600!+\tri:f:95.5", "ACGT", ">S2!chr2!3500!4000!-\tri:f:90", "GGCC",
		">S5!chr1!0!50!+\tri:f:99.1", "CCAA"]
	assert "discarded 2 reads" in r.stderr


# ---------------------------------------------------------------------------
# Usage texts
# ---------------------------------------------------------------------------

USAGE_EXIT1 = ["view", "sam2paf", "delta2paf", "splice2bed", "gff2bed", "stat", "asmstat", "asmgene", "liftover", "vcfpair",
	"call", "bedcov", "mason2fq", "sim2bed", "pbsim2fq", "badread2fa", "junceval", "exoneval", "ov-eval", "mapeval"]
USAGE_EXIT0 = ["misjoin", "sveval", "vcfsel", "pafcmp", "longcs2seq", "paf2gff", "gff2junc"]


@pytest.mark.parametrize("cmd", USAGE_EXIT1 + USAGE_EXIT0)
def test_subcommand_usage(k8, cmd):
	"""Each subcommand without its required arguments prints a usage line naming itself (mapeval on stderr)
	and exits 1, or 0 for the subcommands that return after printing it."""
	r = _pt(cmd, check=False)
	text = r.stderr if cmd == "mapeval" else r.stdout
	assert "Usage:" in text and cmd in text
	assert r.returncode == (1 if cmd in USAGE_EXIT1 else 0)
