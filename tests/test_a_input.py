"""Area A: reading sequences (bseq.c) and the command line's input handling:
FASTA and FASTQ, multi-line records, gzip, lowercase, U and IUPAC bases,
several query files, stdin, empty files, empty and short sequences.

The reader functions run through tests/c/drv_a_bseq.c and are compared with
the records the test wrote; the command line is checked by mapping the same
reads written in different forms.
"""

import gzip
import os

import pytest

import conftest
import helpers_a as A
import minimap2_helpers as H


@pytest.fixture(scope="module")
def ref(tmp_path_factory):
	"""A 30 kbp random target and 12 reads (1-3 kbp, about 1% substitutions) from it, half reverse
	complemented."""
	d = tmp_path_factory.mktemp("ainput")
	rng = H.Rng(201)
	seq = rng.dna(30000)
	fa = A.write_fa(d / "ref.fa", [("tgt", seq)])
	reads = []
	for i in range(12):
		L = rng.randint(1000, 3000)
		st = rng.below(len(seq) - L)
		s = H.mutate(rng, seq[st:st + L], 0.01)
		reads.append(("read%d" % i, s if i % 2 == 0 else H.revcomp(s), "".join(chr(35 + rng.below(30)) for _ in range(L))))
	return conftest._Data(seq=seq, fasta=fa, reads=reads, dir=str(d))


def _sam_without_qual(text):
	out = []
	for l in H.strip_pg(text).splitlines():
		f = l.split("\t")
		if not l.startswith("@"):
			f[10] = "*"
		out.append("\t".join(f))
	return out


def test_read_file_forms_map_identically(ref, tmp_path):
	"""The same reads as one-line FASTA, multi-line FASTA (60 and 7 bases per line), FASTQ, multi-line
	FASTQ and gzip FASTA and FASTQ: identical PAF, and identical SAM apart from QUAL ('*' for FASTA)."""
	recs = [(n, s) for n, s, _ in ref.reads]
	fq_text = "".join("@%s\n%s\n+\n%s\n" % (n, "\n".join(s[i:i + 50] for i in range(0, len(s), 50)), "\n".join(q[i:i + 50]
		for i in range(0, len(q), 50))) for n, s, q in ref.reads)
	files = [
		H.write_fasta(str(tmp_path / "a.fa"), recs),
		H.write_fasta(str(tmp_path / "b.fa"), recs, line_width=60),
		H.write_fasta(str(tmp_path / "c.fa"), recs, line_width=7),
		H.write_fastq(str(tmp_path / "d.fq"), ref.reads),
		H.write_text(str(tmp_path / "e.fq"), fq_text),
		H.write_fasta(str(tmp_path / "f.fa.gz"), recs, line_width=60),
		H.write_fastq(str(tmp_path / "g.fq.gz"), ref.reads),
	]
	with open(files[5], "rb") as f:
		assert f.read(2) == b"\x1f\x8b"
	paf = [H.map_reads(ref.fasta, f, ["-c"]).stdout for f in files]
	assert all(p == paf[0] for p in paf) and len(H.parse_paf(paf[0])) >= len(recs)
	sam = [_sam_without_qual(H.map_reads(ref.fasta, f, ["-a"]).stdout) for f in files]
	assert all(s == sam[0] for s in sam)
	fq_sam = H.map_reads(ref.fasta, files[6], ["-a"]).sam()
	quals = dict((n, q) for n, _, q in ref.reads)
	for r in fq_sam.records:  # FASTQ qualities reach the output, reversed with the read
		if r.is_primary:
			assert r.qual == (quals[r.qname] if not r.is_reverse else quals[r.qname][::-1])


def test_target_file_forms_give_the_same_index(ref, tmp_path):
	"""A target as FASTA with any line width, gzip FASTA, FASTQ, in lowercase, or with U for T (RNA)
	gives a byte-identical .mmi (minimap2 stores bases as 2-bit codes plus N) and the same alignments."""
	seq = ref.seq
	forms = [
		A.write_fa(tmp_path / "t1.fa", [("tgt", seq)], line_width=60),
		A.write_fa(tmp_path / "t2.fa", [("tgt", seq)], line_width=13),
		A.write_fa(tmp_path / "t3.fa.gz", [("tgt", seq)], line_width=80),
		H.write_fastq(str(tmp_path / "t4.fq"), [("tgt", seq, "I" * len(seq))]),
		A.write_fa(tmp_path / "t5.fa", [("tgt", seq.lower())]),
		A.write_fa(tmp_path / "t6.fa", [("tgt", seq[:15000].replace("T", "U") + seq[15000:].lower().replace("t", "u"))]),
	]
	digests = []
	for i, f in enumerate(forms):
		mmi = str(tmp_path / ("t%d.mmi" % i))
		H.run_mm2(["-t", "1", "-d", mmi, f], check=True)
		digests.append(A.file_sha256(mmi))
	assert len(set(digests)) == 1
	p = A.read_mmi(mmi)[0]
	assert p.names == ["tgt"] and p.seq4 == [A.seq_codes(seq)]
	reads = H.write_fasta(str(tmp_path / "r.fa"), [(n, s) for n, s, _ in ref.reads])
	out = [H.strip_pg(H.map_reads(f, reads, ["-a"]).stdout) for f in forms]
	assert all(o == out[0] for o in out)


def test_read_case_u_and_iupac(ref, tmp_path):
	"""A read in lowercase or with U for T maps exactly as in uppercase DNA (same PAF -c --cs columns);
	in SAM its SEQ is the read with U written as T. A read with IUPAC codes (R, Y, K, M, S, W, B, D, H, V,
	N) in place of bases maps to the same place, and nn:i counts those bases."""
	n, s, _ = ref.reads[0]
	iupac = list(s)
	pos = list(range(300, 300 + 11 * 37, 37))
	for p, c in zip(pos, "RYKMSWBDHVN"):
		iupac[p] = c
	recs = [("upper", s), ("lower", s.lower()), ("rna", s.replace("T", "U")), ("rna_lower", s.lower().replace("t", "u")),
		("iupac", "".join(iupac))]
	fa = H.write_fasta(str(tmp_path / "v.fa"), recs)
	paf = dict((p.qname, p) for p in H.map_reads(ref.fasta, fa, ["-c", "--cs"]).paf() if p.is_primary)
	cols = lambda p: p.line.split("\t", 1)[1]
	for name in ("lower", "rna", "rna_lower"):
		assert cols(paf[name]) == cols(paf["upper"]), name
	u, x = paf["upper"], paf["iupac"]
	assert (x.tname, x.strand, x.qstart, x.qend, x.tstart, x.tend) == (u.tname, u.strand, u.qstart, u.qend, u.tstart, u.tend)
	assert x.tags["nn"] == 11 and u.tags["nn"] == 0
	sam = H.map_reads(ref.fasta, fa, ["-a"]).sam()
	for r in sam.primaries().values():
		want = dict(recs)[r.qname].replace("U", "T").replace("u", "t")
		assert (H.revcomp(r.seq) if r.is_reverse else r.seq).upper() == want.upper()


def test_several_query_files(ref, tmp_path):
	"""Several query files in one run (main.c maps them one after another): one SAM header, then the
	records of each file in turn, the same records as separate runs; also for PAF, and for files of
	different formats."""
	a = H.write_fasta(str(tmp_path / "a.fa"), [(n, s) for n, s, _ in ref.reads[:5]])
	b = H.write_fastq(str(tmp_path / "b.fq.gz"), ref.reads[5:9])
	c = H.write_fasta(str(tmp_path / "c.fa"), [(n, s) for n, s, _ in ref.reads[9:]])
	both = H.map_reads(ref.fasta, [a, b, c], ["-a"])
	parts = [H.map_reads(ref.fasta, f, ["-a"]) for f in (a, b, c)]
	assert H.strip_pg(both.stdout) == H.strip_pg(parts[0].stdout) + "".join(H.sam_body(p.stdout) for p in parts[1:])
	assert [r.qname for r in both.sam().records if r.is_primary] == [n for n, _, _ in ref.reads]
	assert H.map_reads(ref.fasta, [a, b, c]).stdout == "".join(H.map_reads(ref.fasta, f).stdout for f in (a, b, c))


def test_query_and_target_from_stdin(ref, tmp_path):
	"""'-' reads the query, or the target, from standard input, plain or gzip, with the same result as
	the file."""
	reads = H.write_fastq(str(tmp_path / "r.fq"), ref.reads)
	want = H.strip_pg(H.map_reads(ref.fasta, reads, ["-a"]).stdout)
	with open(reads, "rb") as f:
		data = f.read()
	assert H.strip_pg(H.map_reads(ref.fasta, "-", ["-a"], stdin=data).stdout) == want
	assert H.strip_pg(H.map_reads(ref.fasta, "-", ["-a"], stdin=gzip.compress(data)).stdout) == want
	with open(ref.fasta, "rb") as f:
		tdata = f.read()
	assert H.strip_pg(H.map_reads("-", reads, ["-a"], stdin=tdata).stdout) == want
	mmi = str(tmp_path / "stdin.mmi")
	H.run_mm2(["-t", "1", "-d", mmi, "-"], stdin=tdata, check=True)
	assert A.read_mmi(mmi)[0].names == ["tgt"]


def test_empty_short_and_unmappable_queries(ref, tmp_path):
	"""An empty query file gives a SAM header only (or no PAF) and exit 0. Each query record gets one SAM
	line in input order: an empty sequence, a sequence shorter than k, an all-N sequence and a record
	with no name are unmapped (flag 4) next to a mapped read; PAF holds only the mapped read, and
	--paf-no-hit lists the others with their lengths. A record with no name draws a warning."""
	empty = H.write_text(str(tmp_path / "empty.fa"), "")
	r = H.map_reads(ref.fasta, empty, ["-a"])
	assert H.sam_body(r.stdout) == "" and r.sam().sq == [("tgt", 30000)]
	assert H.map_reads(ref.fasta, empty).stdout == ""
	n, s, _ = ref.reads[2]
	text = ">empty\n\n>short\nACGTACGTAC\n>allN\n%s\n>%s\n%s\n>\nGATTACAGATTACAGATTACAGATTACA\n" % ("N" * 500, n, s)
	q = H.write_text(str(tmp_path / "q.fa"), text)
	r = H.map_reads(ref.fasta, q, ["-a"])
	sam = H.parse_sam(r.stdout, strict=False)  # strict parsing rejects the empty SEQ and QNAME (known-bugs/a-empty-records-in-sam.md)
	prim = [x for x in sam.records if not x.is_secondary and not x.is_supplementary]
	assert [x.qname for x in prim] == ["empty", "short", "allN", n, ""]
	assert [x.flag for x in prim if x.qname != n] == [4, 4, 4, 4]
	assert not prim[3].is_unmapped
	assert r.messages  # the empty name
	paf = H.map_reads(ref.fasta, q).paf()
	assert [p.qname for p in paf if p.is_primary] == [n]
	nohit = H.map_reads(ref.fasta, q, ["--paf-no-hit"]).paf()
	assert [(p.qname, p.qlen) for p in nohit if p.is_unmapped] == [("empty", 0), ("short", 10), ("allN", 500), ("", 28)]


# --- the reader through the C driver -------------------------------------------------------------


def _read_driver(path, chunk, with_qual=1, with_comment=1, frag=0, fn=3, stdin=None, tmp=None):
	out = os.path.join(tmp, "recs.tsv")
	r = H.run_c_case("drv_a_bseq", "read", [path, chunk, with_qual, with_comment, frag, fn, out], stdin=stdin)
	assert r.passed, r.report()
	batches = []
	for line in H.read_text(out).splitlines():
		b, name, l, seq, qual, comment = line.split("\t")
		if int(b) == len(batches):
			batches.append([])
		batches[-1].append((name, seq, None if qual == "*" else qual, None if comment == "*" else comment))
		assert int(l) == len(seq)
	return batches, r


def _check_batches(batches, records, chunk, frag=False):
	"""The batches hold the records in order; each batch except the last reaches chunk bases. Without
	frag mode a batch ends with the record that reaches chunk bases; in frag mode a batch never splits
	records whose names differ only by a trailing /N."""
	assert [x for b in batches for x in b] == records
	sizes = [sum(len(x[1]) for x in b) for b in batches]
	assert all(s >= chunk for s in sizes[:-1]), sizes
	if not frag:
		for b in batches:
			assert sum(len(x[1]) for x in b[:-1]) < chunk
	else:
		stem = lambda n: n[:-2] if len(n) >= 3 and n[-2] == "/" and n[-1].isdigit() else n
		for b1, b2 in zip(batches, batches[1:]):
			assert stem(b1[-1][0]) != stem(b2[0][0])


def test_bseq_read_fasta_and_fastq(tmp_path):
	"""mm_bseq_read3 returns each record's name (up to the first space), sequence (U and u turned into T
	and t, any other byte kept), quality (FASTQ, with_qual) and comment (with_comment), in file order and
	in batches of about chunk_size bases; mm_bseq_read2 and mm_bseq_read return no comments; the reader
	takes multi-line FASTA and FASTQ, gzip and standard input, and mm_bseq_eof is true at the end."""
	rng = H.Rng(211)
	recs = []
	for i in range(14):
		s = rng.dna(rng.randint(1, 300))
		if i == 3:
			s = s.lower()
		if i == 4:
			s = s[:50] + "UuNnRYacgt-*" + s[50:]
		if i == 5:
			s = ""
		recs.append(("r%d" % i, s, "".join(chr(33 + rng.below(40)) for _ in range(len(s))), "c%d x=%d y" % (i, i) if i % 3 else None))
	fa = H.write_text(str(tmp_path / "a.fa"), "".join(">%s%s\n%s\n" % (n, " " + c if c else "", "\n".join(s[j:j + 60]
		for j in range(0, len(s), 60))) for n, s, _, c in recs))
	fq = H.write_text(str(tmp_path / "a.fq"), "".join("@%s%s\n%s\n+\n%s\n" % (n, " " + c if c else "", s, q) for n, s, q, c in recs))
	fqgz = H.write_fastq(str(tmp_path / "a.fq.gz"), [(n, s, q, c) for n, s, q, c in recs])
	conv = lambda s: s.replace("U", "T").replace("u", "t")
	want_fa = [(n, conv(s), None, c) for n, s, _, c in recs]
	want_fq = [(n, conv(s), q if s else None, c) for n, s, q, c in recs]
	for chunk in (1, 250, 10 ** 9):
		b, _ = _read_driver(fa, chunk, tmp=str(tmp_path))
		_check_batches(b, want_fa, chunk)
		for f in (fq, fqgz):
			b, _ = _read_driver(f, chunk, tmp=str(tmp_path))
			_check_batches(b, want_fq, chunk)
	no_c = [(n, s, q, None) for n, s, q, _ in want_fq]
	no_q = [(n, s, None, c) for n, s, _, c in want_fq]
	for fn in (2, 1):
		b, _ = _read_driver(fq, 250, fn=fn, tmp=str(tmp_path))
		_check_batches(b, no_c, 250)
	b, _ = _read_driver(fq, 250, with_qual=0, tmp=str(tmp_path))
	_check_batches(b, no_q, 250)
	with open(fq, "rb") as f:
		b, _ = _read_driver("-", 100, stdin=f.read(), tmp=str(tmp_path))
	_check_batches(b, want_fq, 100)
	nameless = H.write_text(str(tmp_path / "nameless.fa"), ">\nACGT\n>ok\nACGT\n")
	b, r = _read_driver(nameless, 1, tmp=str(tmp_path))
	assert [x[0] for bb in b for x in bb] == ["", "ok"] and r.stderr.strip()
	H.assert_c_case("drv_a_bseq", "open_missing", str(tmp_path / "missing.fa"))


def test_bseq_read_frag_mode_keeps_mates_together(tmp_path):
	"""In fragment mode (frag_mode, as -x sr uses for an interleaved file) a batch is not ended between
	reads whose names differ only by /1 and /2, even at chunk size 1; without it, chunk size 1 gives one
	record per batch."""
	rng = H.Rng(221)
	recs = []
	for i in range(9):
		for m in (1, 2):
			recs.append(("p%d/%d" % (i, m), rng.dna(rng.randint(20, 80)), None, None))
		if i % 3 == 0:
			recs.append(("single%d" % i, rng.dna(30), None, None))
	f = H.write_fasta(str(tmp_path / "il.fa"), [(n, s) for n, s, _, _ in recs])
	for chunk in (1, 100, 10 ** 6):
		b, _ = _read_driver(f, chunk, frag=1, tmp=str(tmp_path))
		_check_batches(b, recs, chunk, frag=True)
	b, _ = _read_driver(f, 1, frag=0, tmp=str(tmp_path))
	assert [len(x) for x in b] == [1] * len(recs)


def test_bseq_read_frag_two_files(tmp_path):
	"""mm_bseq_read_frag2 and mm_bseq_read_frag read one record from each file in turn (read 1, read 2,
	read 1, ...) in batches of whole pairs that end with the pair reaching chunk_size; with fewer records
	in one file they stop at its end and warn."""
	rng = H.Rng(231)
	r1 = [("q%d" % i, rng.dna(rng.randint(30, 90)), None, "a%d b" % i) for i in range(10)]
	r2 = [("q%d" % i, rng.dna(rng.randint(30, 90)), None, None) for i in range(10)]
	f1 = H.write_text(str(tmp_path / "r1.fa"), "".join(">%s %s\n%s\n" % (n, c, s) for n, s, _, c in r1))
	f2 = H.write_fasta(str(tmp_path / "r2.fa"), [(n, s) for n, s, _, _ in r2])
	f2short = H.write_fasta(str(tmp_path / "r2s.fa"), [(n, s) for n, s, _, _ in r2[:7]])
	inter = [x for a, b in zip(r1, r2) for x in (a, b)]
	for fn, keep_c in ((2, True), (1, False)):
		for chunk in (1, 200, 10 ** 6):
			out = str(tmp_path / "f.tsv")
			r = H.run_c_case("drv_a_bseq", "frag", [f1, f2, chunk, 1, 1, fn, out])
			assert r.passed, r.report()
			batches = []
			for line in H.read_text(out).splitlines():
				b, name, l, seq, qual, comment = line.split("\t")
				if int(b) == len(batches):
					batches.append([])
				batches[-1].append((name, seq, None, None if comment == "*" else comment))
			want = [(n, s, None, c if keep_c else None) for n, s, _, c in inter]
			assert [x for b in batches for x in b] == want
			sizes = [sum(len(x[1]) for x in b) for b in batches]
			assert all(s >= chunk for s in sizes[:-1])
			assert all(sum(len(x[1]) for x in b[:-2]) < chunk for b in batches)
	out = str(tmp_path / "s.tsv")
	r = H.run_c_case("drv_a_bseq", "frag", [f1, f2short, 10 ** 6, 1, 1, 2, out])
	assert r.passed, r.report()
	assert len(H.read_text(out).splitlines()) == 14 and r.stderr.strip()
