"""SAM output (format.c): header lines, read groups, unmapped, secondary and supplementary
records, SEQ and QUAL, SA:Z, -Q, -Y, -y, --eqx, --secondary, --secondary-seq, -N 0,
--sam-hit-only, -L, cs:Z, MD:Z and ds:Z, and templates of three segments.

Expected values come from the SAM specification, minimap2.1 and the usage text, and from the
genome and the simulated reads; one run of an option is compared with the default run where the
documentation says the option changes only one thing.
"""

import re

import pytest

import helpers_d as D
import minimap2_helpers as H


@pytest.fixture(scope="module")
def mixed(synth, tmp_path_factory):
	"""(Mixed reads, FASTQ path with comments, FASTA path, {qname: [seq]})."""
	m = D.mixed_reads(synth.genome)
	d = tmp_path_factory.mktemp("d_sam")
	fq = D.write_mixed(d / "mixed.fq", m)
	fa = D.write_mixed(d / "mixed.fa", m, comments=False, fmt="fasta")
	return m, fq, fa, dict((n, [s]) for n, s in m.seqs.items())


@pytest.fixture(scope="module")
def base_sam(synth, mixed):
	"""`minimap2 -a` on the mixed reads (default preset)."""
	return H.map_reads(synth.fasta, mixed[1], ["-a"])


def _records(text):
	return H.parse_sam(text).records


def test_header_lines(synth, mixed, base_sam):
	"""@HD VN:1.6 SO:unsorted GO:query; one @SQ per contig in FASTA order with its length; no @RG
	without -R; one @PG (ID and PN minimap2, VN the version, CL `minimap2` and the arguments as
	given); every header line comes before the first record."""
	text = base_sam.stdout
	head = [l for l in D.lines(text) if l.startswith("@")]
	assert D.lines(text)[:len(head)] == head
	sam = base_sam.sam()
	assert head[0] == "@HD\tVN:1.6\tSO:unsorted\tGO:query"
	assert sam.sq == [(n, synth.genome.length(n)) for n in synth.genome.names]
	assert [l for l in head if l.startswith("@SQ")] == ["@SQ\tSN:%s\tLN:%d" % x for x in sam.sq]
	assert sam.rg == [] and sam.co == []
	assert len(sam.pg) == 1
	pg = sam.pg[0]
	assert pg["ID"] == "minimap2" and pg["PN"] == "minimap2" and pg["VN"] == D.mm_version()
	assert pg["CL"] == "minimap2 " + " ".join(base_sam.argv[1:])
	assert head[-1].startswith("@PG\t") and len(head) == 2 + len(sam.sq)


def test_read_group(synth, mixed, base_sam):
	"""-R (minimap2.1, usage: a line like '@RG\\tID:foo\\tSM:bar'): the header gets the @RG line
	with its escaped TABs turned into TABs, and every record, unmapped and secondary ones included,
	carries RG:Z with the ID; nothing else changes. PAF has no read groups and is unchanged. An
	ID of 255 characters is accepted, 256 is an error ("@RG:ID is longer than 255 characters"),
	as are a line not starting with @RG, a literal TAB, and a line without ID: exit 1 with a
	message and no records."""
	fa, fq = base_sam.argv[-2], mixed[1]
	r = H.map_reads(fa, fq, ["-a", "-R", "@RG\\tID:grp1\\tSM:sample 1\\tLB:lib"])
	sam = r.sam()
	assert sam.rg == [{"ID": "grp1", "SM": "sample 1", "LB": "lib"}]
	assert [l for l in D.lines(r.stdout) if l.startswith("@RG")] == ["@RG\tID:grp1\tSM:sample 1\tLB:lib"]
	assert len(sam.records) == len(_records(base_sam.stdout))
	for g, b in zip(sam.records, _records(base_sam.stdout)):
		assert g.tags["RG"] == "grp1" and D.drop_tag(g.line, "RG") == b.line
	head = [l for l in D.lines(r.stdout) if l.startswith("@") and not l.startswith(("@RG", "@PG"))]
	assert head == [l for l in D.lines(base_sam.stdout) if l.startswith("@") and not l.startswith("@PG")]
	paf = H.map_reads(fa, fq, ["-R", "@RG\\tID:grp1"])
	assert paf.stdout == H.map_reads(fa, fq).stdout
	long_id = "x" * 255
	r = H.map_reads(fa, fq, ["-a", "-R", "@RG\\tID:" + long_id + "\\tSM:s"])
	assert all(x.tags["RG"] == long_id for x in r.sam().records)
	for bad in ("@RG\\tID:" + "x" * 256, "RG\\tID:a", "@RG\tID:a", "@RG\\tSM:a"):
		e = H.map_reads(fa, fq, ["-a", "-R", bad], check=False)
		assert e.returncode == 1 and e.messages and D.body_lines(e.stdout) == [], bad


def test_every_read_has_one_primary_and_consistent_records(synth, mixed, base_sam):
	"""Every input read has exactly one primary record, in input order; SEQ is the read (reverse
	complemented on the reverse strand), the CIGAR spans the read, and NM equals the edit distance
	recomputed from the genome (sam_invariant_errors)."""
	m, _, _, reads = mixed
	sam = base_sam.sam()
	assert H.sam_invariant_errors(sam, synth.genome, reads) == []
	prim = [r.qname for r in sam.records if r.is_primary]
	assert prim == m.names()
	for r in sam.records:
		if r.is_unmapped:
			continue
		assert r.query_length == len(m.seqs[r.qname])
		assert r.tags["tp"] in ("P", "S", "I", "i")


def test_unmapped_records(mixed, base_sam):
	"""A read without a hit is one record with FLAG 4, RNAME '*', POS 0, MAPQ 0, CIGAR '*', RNEXT
	'*', PNEXT 0, TLEN 0 and the read's SEQ and QUAL as given (SAM spec, unmapped reads);
	--sam-hit-only drops exactly these records and leaves every other line as it was."""
	m = mixed[0]
	recs = _records(base_sam.stdout)
	un = [r for r in recs if r.is_unmapped]
	assert sorted(r.qname for r in un) == sorted(m.names("random"))
	for r in un:
		assert (r.flag, r.rname, r.pos, r.mapq, r.cigar, r.rnext, r.pnext, r.tlen) == (4, "*", 0, 0, "*", "*", 0, 0)
		assert r.seq == m.seqs[r.qname] and r.qual == m.quals[r.qname]
	only = H.map_reads(base_sam.argv[-2], mixed[1], ["-a", "--sam-hit-only"])
	want = [l for l in D.body_lines(base_sam.stdout) if not (int(l.split("\t")[1]) & 4)]
	assert D.body_lines(only.stdout) == want


def test_secondary_records(synth, mixed, base_sam):
	"""Reads copied from one unit of an exact repeat get secondary records (FLAG 0x100, tp:A:S) at
	the other copies, with SEQ and QUAL '*', soft clips only and a CIGAR spanning the whole read;
	the primary and each secondary lie on one of the repeat's copies."""
	m = mixed[0]
	by = H.parse_sam(base_sam.stdout).by_qname()
	for name in m.names("repeat"):
		recs = by[name]
		sec = [r for r in recs if r.is_secondary]
		assert len(sec) >= 1 and sum(1 for r in recs if r.is_primary) == 1
		copies = set((c, s) for c, s, e, _ in m.spans[name])
		for r in recs:
			assert (r.rname, r.pos0) in copies, r.line[:120]
		for r in sec:
			assert r.seq == "*" and r.qual == "*" and r.tags["tp"] == "S" and not (r.flag & 0x800)
			assert "H" not in r.cigar and r.query_length == len(m.seqs[name])


def test_secondary_no_and_n0(synth, mixed, base_sam):
	"""--secondary=no prints the default output without its secondary records (minimap2.1);
	--secondary=yes is the default; a value other than yes/no prints a warning naming
	--secondary; -N 0 warns and gives the --secondary=no output (its own warning says it is
	changed to '-N 5 --secondary=no')."""
	fa, fq = base_sam.argv[-2], mixed[1]
	want = [l for l in D.body_lines(base_sam.stdout) if not (int(l.split("\t")[1]) & 0x100)]
	assert len(want) < len(D.body_lines(base_sam.stdout))
	no = H.map_reads(fa, fq, ["-a", "--secondary=no"])
	assert D.body_lines(no.stdout) == want and no.messages == []
	yes = H.map_reads(fa, fq, ["-a", "--secondary=yes"])
	assert D.body_lines(yes.stdout) == D.body_lines(base_sam.stdout)
	bad = H.map_reads(fa, fq, ["-a", "--secondary=maybe"])
	assert any("--secondary" in l for l in bad.messages)
	assert D.body_lines(bad.stdout) == D.body_lines(base_sam.stdout)
	n0 = H.map_reads(fa, fq, ["-a", "-N", "0"])
	assert n0.messages and D.body_lines(n0.stdout) == want


def test_secondary_seq(synth, mixed, base_sam):
	"""--secondary-seq (minimap2.1): secondary records carry the aligned part of the read (in
	alignment orientation) and its qualities, with hard clips; every other record is unchanged."""
	m = mixed[0]
	r = H.map_reads(base_sam.argv[-2], mixed[1], ["-a", "--secondary-seq"])
	got, base = _records(r.stdout), _records(base_sam.stdout)
	assert len(got) == len(base)
	n_sec = 0
	for g, b in zip(got, base):
		if not b.is_secondary:
			assert g.line == b.line
			continue
		n_sec += 1
		f = g.line.split("\t")
		fb = b.line.split("\t")
		assert f[:5] + f[6:9] + f[11:] == fb[:5] + fb[6:9] + fb[11:]
		assert g.cigar.replace("H", "S") == b.cigar
		h5, h3 = D.hard_clips(g.cigar)
		seq, qual = m.seqs[g.qname], m.quals[g.qname]
		if g.is_reverse:
			seq, qual = H.revcomp(seq), qual[::-1]
		assert g.seq == seq[h5:len(seq) - h3] and g.qual == qual[h5:len(qual) - h3]
		assert g.query_length == len(seq)
	assert n_sec > 0
	assert H.sam_invariant_errors(H.parse_sam(r.stdout), synth.genome, mixed[3]) == []


def _approx_cigar_ok(entry, rec, read_len):
	"""The SA:Z CIGAR of an alignment is approximate (minimap2.1): it must have the record's
	clipping at each end as soft clips, the record's reference span and the whole read."""
	cig = entry[3]
	lc, rc = H.cigar_clips(cig)
	return (lc, rc) == (rec.left_clip, rec.right_clip) and "H" not in cig and \
		H.cigar_ref_span(cig) == rec.ref_span and H.cigar_query_span(cig) == read_len


def test_supplementary_and_sa(synth, mixed, base_sam):
	"""Chimeric reads give one primary and one or more supplementary records (FLAG 0x800); a
	supplementary record is hard-clipped and carries the aligned part of the read; each piece's
	records lie at its true contig and strand. Every primary and supplementary record of a split
	read lists all the others in SA:Z (rname, pos, strand, mapQ and NM as in that record; the
	approximate CIGAR has its clips and spans), and records of reads that are not split have no SA."""
	m = mixed[0]
	by = H.parse_sam(base_sam.stdout).by_qname()
	for name in m.names("chimeric"):
		recs = [r for r in by[name] if not r.is_secondary]
		sup = [r for r in recs if r.is_supplementary]
		assert len(recs) >= 2 and len(sup) == len(recs) - 1, [r.line[:100] for r in recs]
		seq = m.seqs[name]
		for r in sup:
			h5, h3 = D.hard_clips(r.cigar)
			assert "S" not in r.cigar and (h5 > 0 or h3 > 0)
			assert r.seq == D.record_seq(r, seq)
		for r in recs:
			others = [o for o in recs if o is not r]
			entries = _sa_sorted(D.sa_entries(r))
			assert len(entries) == len(others)
			for e, o in zip(entries, _sa_sorted([(o.rname, o.pos, o.strand, None, o.mapq, o.tags["NM"]) for o in others])):
				assert (e[0], e[1], e[2], e[4], e[5]) == (o[0], o[1], o[2], o[4], o[5])
			for o in others:
				e = [x for x in D.sa_entries(r) if (x[0], x[1], x[2]) == (o.rname, o.pos, o.strand)][0]
				assert _approx_cigar_ok(e, o, len(seq)), (e, o.cigar)
		truth = set((c, st) for c, s, e, st in m.spans[name])
		assert set((r.rname, r.strand) for r in recs) <= truth
		if "N" in seq:
			assert all(r.tags["nn"] > 0 for r in recs)
	for name in m.names("unique") + m.names("repeat"):
		assert all("SA" not in r.tags for r in by[name])


def test_sa_gives_each_alignments_own_mapq(synth, tmp_path):
	"""SAM spec, SA:Z: each entry is rname,pos,strand,CIGAR,mapQ,NM of another alignment of the read.
	Chimeric reads of 2.5 kbp of unique sequence followed by a whole copy of an exact repeat (either
	strand): the unique piece is the primary with a high MAPQ, the repeat piece a supplementary
	alignment with a low one, and each record's SA:Z entry carries the MAPQ (and position and
	strand) of the other record, not its own."""
	g = synth.genome
	reps = [f for f in g.features if f.kind == "repeat"][:4]
	starts = [st for st in range(1000, 70000, 997) if g.is_unique("chr1", st, st + 2500, 200)]
	reads = []
	for i, f in enumerate(reps):
		st = [x for x in starts if abs(x - f.start) > 5000][i]
		reads.append(("cr%d" % i, g.fetch("chr1", st, st + 2500) + g.fetch(f.contig, f.start, f.end, "+-"[i % 2])))
	sam = H.map_reads(synth.fasta, H.write_fasta(str(tmp_path / "cr.fa"), reads), ["-a"]).sam()
	n_diff = 0
	for name, _ in reads:
		recs = [r for r in sam.by_qname()[name] if not r.is_secondary]
		assert len(recs) == 2 and sum(r.is_supplementary for r in recs) == 1, [r.line[:80] for r in recs]
		for r, o in (recs, recs[::-1]):
			e = D.sa_entries(r)
			assert len(e) == 1 and (e[0][0], e[0][1], e[0][2], e[0][4], e[0][5]) == (o.rname, o.pos, o.strand, o.mapq, o.tags["NM"])
		n_diff += recs[0].mapq != recs[1].mapq
	assert n_diff == len(reads)


def _sa_sorted(entries):
	return sorted(entries, key=lambda e: (e[0], e[1], e[2]))


def test_soft_clipped_supplementary(synth, mixed, base_sam):
	"""-Y (minimap2.1): supplementary records use soft clipping and carry the whole read and its
	qualities; primary and unmapped records are as without -Y. Secondary records keep their
	fields except SEQ and QUAL, which may be '*' or the whole read (SAM spec; base minimap2 writes
	the whole read under -Y)."""
	m = mixed[0]
	r = H.map_reads(base_sam.argv[-2], mixed[1], ["-a", "-Y"])
	got, base = _records(r.stdout), _records(base_sam.stdout)
	assert len(got) == len(base)
	n = 0
	for g, b in zip(got, base):
		if b.is_secondary:
			fg, fb = g.line.split("\t"), b.line.split("\t")
			assert fg[:9] + fg[11:] == fb[:9] + fb[11:]
			assert g.seq in ("*", D.record_seq(g, m.seqs[g.qname]))
			continue
		if not b.is_supplementary:
			assert g.line == b.line
			continue
		n += 1
		assert g.cigar == b.cigar.replace("H", "S") and "H" not in g.cigar
		seq, qual = m.seqs[g.qname], m.quals[g.qname]
		assert g.seq == (H.revcomp(seq) if g.is_reverse else seq)
		assert g.qual == (qual[::-1] if g.is_reverse else qual)
		fg, fb = g.line.split("\t"), b.line.split("\t")
		assert fg[:5] + fg[6:9] + fg[11:] == fb[:5] + fb[6:9] + fb[11:]
	assert n == len(m.names("chimeric"))
	assert H.sam_invariant_errors(H.parse_sam(r.stdout), synth.genome, mixed[3]) == []


def test_qualities(synth, mixed, base_sam):
	"""QUAL is the input quality string, reversed on the reverse strand (and cut like SEQ on
	hard-clipped records); -Q writes '*' and changes nothing else; FASTA input gives '*'."""
	m = mixed[0]
	recs = _records(base_sam.stdout)
	for r in recs:
		if r.seq == "*":
			continue
		q = m.quals[r.qname]
		if r.is_reverse:
			q = q[::-1]
		h5, h3 = D.hard_clips(r.cigar)
		assert r.qual == q[h5:len(q) - h3]
	noq = _records(H.map_reads(base_sam.argv[-2], mixed[1], ["-a", "-Q"]).stdout)
	assert len(noq) == len(recs)
	for g, b in zip(noq, recs):
		f, fb = g.line.split("\t"), b.line.split("\t")
		assert f[10] == "*" and f[:10] + f[11:] == fb[:10] + fb[11:]
	fa = _records(H.map_reads(base_sam.argv[-2], mixed[2], ["-a"]).stdout)
	assert [r.qual for r in fa] == ["*"] * len(fa)
	assert [r.line.split("\t")[:10] for r in fa] == [r.line.split("\t")[:10] for r in recs]


def test_copy_comment(synth, mixed, base_sam):
	"""-y (minimap2.1): the FASTQ comment is appended to every record of its read, after the
	other optional fields, as written (a comment holding two TAB-separated tags gives two
	fields); reads without a comment and every other field are unchanged."""
	m = mixed[0]
	r = H.map_reads(base_sam.argv[-2], mixed[1], ["-a", "-y"])
	got, base = r.stdout.split("\n"), base_sam.stdout.split("\n")
	got = [l for l in got if l and not l.startswith("@")]
	base = [l for l in base if l and not l.startswith("@")]
	assert len(got) == len(base)
	n = 0
	for g, b in zip(got, base):
		c = m.comments[b.split("\t")[0]]
		if c is None:
			assert g == b
		else:
			assert g == b + "\t" + c
			n += 1
	assert n > 5
	H.parse_sam(r.stdout)


def test_eqx(synth, mixed, base_sam):
	"""--eqx (minimap2.1): CIGARs use = for runs of equal bases and X for mismatches instead of
	M, exactly where the read and the genome agree or differ; merging them back gives the default
	CIGAR, and nothing else changes. With -c the PAF cg:Z uses the same operators."""
	m = mixed[0]
	r = H.map_reads(base_sam.argv[-2], mixed[1], ["-a", "--eqx"])
	got, base = _records(r.stdout), _records(base_sam.stdout)
	assert len(got) == len(base)
	n_x = 0
	for g, b in zip(got, base):
		fg, fb = g.line.split("\t"), b.line.split("\t")
		assert fg[:5] + fg[6:] == fb[:5] + fb[6:]
		if b.is_unmapped:
			assert g.cigar == "*"
			continue
		assert "M" not in g.cigar and D.merge_eqx(g.cigar) == b.cigar
		s = D.record_seq(g, m.seqs[g.qname])
		assert g.cigar == _with_hard(g.cigar, D.eqx_cigar(synth.genome[g.rname], g.pos0, s, D.merge_eqx(g.cigar)))
		n_x += g.cigar.count("X")
	assert n_x > 50
	paf = H.map_reads(base_sam.argv[-2], mixed[1], ["-c", "--eqx"]).paf()
	sam_cig = [H.strip_clips(x.cigar) for x in got if not x.is_unmapped]
	assert [p.cigar for p in paf] == sam_cig


def _with_hard(orig, inner):
	"""inner (a CIGAR without hard clips) with the hard clips of orig put back."""
	ops = H.parse_cigar(orig)
	pre = "%d%s" % ops[0] if ops[0][1] == "H" else ""
	post = "%d%s" % ops[-1] if ops[-1][1] == "H" else ""
	return pre + inner + post


def test_long_cigar_option_keeps_short_cigars(mixed, base_sam):
	"""-L moves a CIGAR to the CG tag only above 65535 operations (minimap2.1); on ordinary reads
	the output is unchanged. (The long case is in test_d_api.py through a C driver.)"""
	r = H.map_reads(base_sam.argv[-2], mixed[1], ["-a", "-L"])
	assert D.body_lines(r.stdout) == D.body_lines(base_sam.stdout)


@pytest.mark.parametrize("opt,tag,long_form", [
	("--cs", "cs", False), ("--cs=short", "cs", False), ("--cs=long", "cs", True), ("--MD", "MD", None)])
def test_cs_md(synth, mixed, base_sam, opt, tag, long_form):
	"""--cs[=short|long] and --MD (minimap2.1): every mapped record carries the tag, equal to its
	recomputation from the genome (cs:Z as the man page defines it; MD as the SAM spec, without
	the final 0 base minimap2 omits); the rest of each record is the default output."""
	r = H.map_reads(base_sam.argv[-2], mixed[1], ["-a", opt])
	sam = H.parse_sam(r.stdout)
	assert H.sam_invariant_errors(sam, synth.genome, mixed[3]) == []
	base = _records(base_sam.stdout)
	assert len(sam.records) == len(base)
	for g, b in zip(sam.records, base):
		if b.is_unmapped:
			assert g.line == b.line
			continue
		v = g.tags[tag]
		if long_form is not None:
			assert v.startswith("=") == long_form or v[0] in "*+-~"
			assert (":" in v) == (not long_form)
		assert D.drop_tag(g.line, tag) == b.line


def test_cs_variants_and_deprecated_s(mixed, base_sam):
	"""--cs=none writes no cs; an unknown --cs value warns naming --cs and writes the short form;
	the deprecated -S warns and writes the long form (TARGETS.md: same as --cs=long)."""
	fa, fq = base_sam.argv[-2], mixed[1]
	none = H.map_reads(fa, fq, ["-a", "--cs=none"])
	assert D.body_lines(none.stdout) == D.body_lines(base_sam.stdout)
	short = H.map_reads(fa, fq, ["-a", "--cs"])
	bad = H.map_reads(fa, fq, ["-a", "--cs=medium"])
	assert any("--cs" in l for l in bad.messages)
	assert D.body_lines(bad.stdout) == D.body_lines(short.stdout)
	long_ = H.map_reads(fa, fq, ["-a", "--cs=long"])
	s = H.map_reads(fa, fq, ["-a", "-S"])
	assert s.messages and any("-S" in l for l in s.messages)
	assert D.body_lines(s.stdout) == D.body_lines(long_.stdout)


def test_ds_tag(synth, mixed, base_sam):
	"""--ds: ds:Z is cs:Z (short form) with the ambiguous bases of each indel in square brackets
	(NEWS: "indicates uncertainty in INDEL positions. It is an extension to cs"); without the
	brackets it equals the --cs run's cs:Z, and it equals the value recomputed from the genome
	(helpers_d.expected_ds). Indels in random sequence often border repeated bases, so brackets
	appear."""
	m = mixed[0]
	fa, fq = base_sam.argv[-2], mixed[1]
	ds = _records(H.map_reads(fa, fq, ["-a", "--ds"]).stdout)
	cs = _records(H.map_reads(fa, fq, ["-a", "--cs"]).stdout)
	assert len(ds) == len(cs)
	n_br = 0
	for d, c in zip(ds, cs):
		if c.is_unmapped:
			continue
		v = d.tags["ds"]
		assert D.strip_ds(v) == c.tags["cs"]
		s = D.record_seq(d, m.seqs[d.qname])
		assert v == D.expected_ds(synth.genome[d.rname], d.pos0, s, d.cigar), d.qname
		n_br += v.count("[")
	assert n_br > 20


@pytest.fixture(scope="module")
def end_reads(synth, tmp_path_factory):
	"""Error-free 120-bp reads from unique sequence with substitutions placed at the read ends (first,
	second, last or second-to-last base, or both ends); -x sr's end bonus extends alignments over
	them, so tags start or end with a mismatch or a one-base match."""
	reads = H.simulate_reads(synth.genome, 24, length=120, seed=611, sub_rate=0.0, indel_rate=0.0, unique=True, unique_pad=50)
	recs = []
	for i, r in enumerate(reads):
		s = list(r.fwd_seq)
		n = len(s)
		for p in [(1,), (n - 2,), (0,), (n - 1,), (1, n - 2), (2, n - 3)][i % 6]:
			s[p] = "A" if s[p] != "A" else "C"
		fwd = "".join(s)
		recs.append((r.name, fwd if r.strand == "+" else H.revcomp(fwd), r.qual))
	fq = H.write_fastq(str(tmp_path_factory.mktemp("d_ends") / "ends.fq"), recs)
	return fq, dict((n, [q]) for n, q, _ in recs)


@pytest.mark.parametrize("opt,tag", [("--MD", "MD"), ("--cs", "cs"), ("--cs=long", "cs"), ("--ds", "ds"), ("--eqx", None)])
def test_tags_at_read_ends(synth, end_reads, opt, tag):
	"""MD:Z, cs:Z (short and long), ds:Z and --eqx CIGARs of alignments that start or end with a
	mismatch or a one-base match equal their recomputation from the genome (the first or last
	number of MD may be 0 or 1; MD omits a final 0, known-bugs/md-trailing-zero.md)."""
	fq, reads = end_reads
	sam = H.map_reads(synth.fasta, fq, ["-a", opt], preset="sr").sam()
	long_cs = opt == "--cs=long"
	assert H.sam_invariant_errors(sam, synth.genome, reads) == []
	recs = [r for r in sam.records if r.is_primary]
	assert len(recs) == len(reads) and all(not r.is_unmapped and "S" not in r.cigar for r in recs)
	vals = [r.tags[tag] for r in recs] if tag else [r.cigar for r in recs]
	for r, v in zip(recs, vals):
		s = D.record_seq(r, reads[r.qname][0])
		ref = synth.genome[r.rname]
		if tag == "ds":
			assert v == D.expected_ds(ref, r.pos0, s, r.cigar)
		elif long_cs:
			assert v == H.compute_cs(ref, r.pos0, s, r.cigar, long_form=True)
		elif tag is None:
			assert v == D.eqx_cigar(ref, r.pos0, s, D.merge_eqx(v))
	if tag == "MD":
		assert any(re.search(r"[A-Z]1$", v) for v in vals) and any(v.startswith("0") for v in vals)
		assert any(re.match(r"^1[A-Z]", v) for v in vals) and any(re.search(r"[A-Z]$", v) for v in vals)
	elif tag is None:
		assert any(v.startswith("1X") for v in vals) and any(v.endswith("X1=") for v in vals)
	else:
		assert any(v.startswith("*") for v in vals) and any(v[-3] == "*" or v.endswith("*" + v[-2:]) for v in vals)


@pytest.fixture(scope="module")
def three_seg(synth, tmp_path_factory):
	"""Templates of three segments: --frag=yes groups adjacent reads with the same name."""
	g = synth.genome
	reads = H.simulate_reads(g, 9, length=150, seed=541, unique=True)
	rng = H.Rng(542)
	recs = []
	for i in range(3):
		for j in range(3):
			r = reads[3 * i + j]
			recs.append(("tpl%d" % i, r.seq, r.qual))
	recs += [("tplu", reads[0].seq, reads[0].qual), ("tplu", rng.dna(150), "I" * 150), ("tplu", reads[1].seq, reads[1].qual)]
	fq = H.write_fastq(str(tmp_path_factory.mktemp("d_tpl") / "tpl.fq"), recs)
	truth = {}
	for i in range(3):
		truth["tpl%d" % i] = reads[3 * i:3 * i + 3]
	truth["tplu"] = [reads[0], None, reads[1]]
	return fq, truth


def test_three_segment_template(synth, three_seg):
	"""With --frag=yes, adjacent reads of one name form a template of three segments (SAM spec
	1.4): every record has FLAG 0x1 and FI:i (the segment index); the first segment has 0x40 and
	not 0x80, the last 0x80 and not 0x40 (the middle one: known-bugs/d-middle-segment-flags.md);
	RNEXT/PNEXT give the next segment's primary position, the first segment being next after the
	last; 0x8 and 0x20 follow that next segment; an unmapped segment is placed at the position of
	a mapped segment of its template; mapped segments lie at the truth. PAF names the segments
	name/1, name/2 and name/3."""
	fq, truth = three_seg
	sam = H.map_reads(synth.fasta, fq, ["-a", "--frag=yes"]).sam()
	by = sam.by_qname()
	assert sorted(by) == sorted(truth)
	for name, recs in by.items():
		prim = [r for r in recs if r.is_primary]
		assert len(prim) == 3
		assert [r.tags["FI"] for r in prim] == [0, 1, 2]
		for i, r in enumerate(prim):
			assert r.flag & 0x1
			if i != 1:  # the middle segment: known-bugs/d-middle-segment-flags.md
				assert bool(r.flag & 0x40) == (i == 0) and bool(r.flag & 0x80) == (i == 2)
			nxt = prim[(i + 1) % 3]
			if nxt.is_unmapped:
				assert r.flag & 0x8 and (r.rnext, r.pnext) == ("=", r.pos)
			else:
				assert not (r.flag & 0x8) and bool(r.flag & 0x20) == nxt.is_reverse
				assert r.pnext == nxt.pos and r.rnext == ("=" if nxt.rname == r.rname else nxt.rname)
			t = truth[name][i]
			if t is None:
				assert r.is_unmapped and (r.rname, r.pos) in [(o.rname, o.pos) for o in prim if not o.is_unmapped]
			else:
				err = H.placement_error(r, t)
				assert err is not None and err <= 5, (r.line[:120], t)
	paf = H.map_reads(synth.fasta, fq, ["--frag=yes"]).paf()
	assert sorted(set(p.qname for p in paf)) == sorted("%s/%d" % (n, k) for n in truth for k in (1, 2, 3) if not (n == "tplu" and k == 2))


@pytest.mark.pinned
def test_pinned_sam_output(mixed, base_sam):
	"""Pins the current output: `-a -y --cs --secondary-seq -Y -Q` on the mixed reads (unique,
	repeat, chimeric, low-complexity, ambiguous and unmapped reads), without @PG."""
	r = H.map_reads(base_sam.argv[-2], mixed[1], ["-a", "-y", "--cs", "--secondary-seq", "-Y", "-Q"])
	H.check_pinned("d/mixed-a-y-cs-secseq-Y-Q.sam.gz", H.strip_pg(r.stdout))
