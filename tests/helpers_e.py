"""Helpers for area E: the mappy binding, python/minimap2.py and the misc/ scripts.

hit_key / paf_key       one tuple per alignment, comparable between a mappy Alignment and a PAF line
cli_hits                minimap2 -c --secondary=yes on a read file -> {qname: [paf_key, ...]}
mappy_hits              Aligner.map() over records -> {name: [hit_key, ...]}
hit_problems            differences between two such dicts, as readable strings
ref_oriented            the aligned part of a read in reference orientation
mm_flags                the MM_F_* constants of minimap.h in the source tree under test
mmi_header              (magic, w, k, b, n_seq, flag) of an .mmi file
mappy_env               environment that lets a child Python import the cached mappy build
"""

import os
import re
import struct

import minimap2_build
import minimap2_helpers as H


def hit_key(h):
	"""The fields of a mappy Alignment that a `minimap2 -c` PAF line also carries."""
	return (h.ctg, h.ctg_len, h.r_st, h.r_en, h.q_st, h.q_en, h.strand, h.mapq, h.is_primary, h.cigar_str,
		h.NM, h.mlen, h.blen)


def paf_key(p):
	"""hit_key() of the alignment a PAF line describes (strand as +1/-1, NM from NM:i)."""
	return (p.tname, p.tlen, p.tstart, p.tend, p.qstart, p.qend, 1 if p.strand == "+" else -1, p.mapq, p.is_primary,
		p.cigar, p.tags.get("NM"), p.nmatch, p.alen)


def cli_hits(target, reads, args=(), preset=None, threads=None):
	"""Run `minimap2 [-x preset] -c --secondary=yes [args]` and group PAF lines by query name, in output order.

	mappy always computes CIGARs (-c) and returns secondary alignments, so the
	command line must print them too (presets such as sr turn them off)."""
	res = H.map_reads(target, reads, ["-c", "--secondary=yes"] + list(args), preset=preset, threads=threads)
	out = {}
	for p in res.paf():
		out.setdefault(p.qname, []).append(paf_key(p))
	return out


def mappy_hits(aligner, records, use_name=False, **map_kw):
	"""{name: [hit_key, ...]} from aligner.map(seq) for (name, seq, ...) records or SimReads.

	With use_name=True the read name is passed as map(name=...), which seeds the
	tie-breaking hash as the command line does; without it mappy behaves like
	`minimap2 --no-hash-name`."""
	out = {}
	for r in records:
		name, seq = (r.name, r.seq) if hasattr(r, "seq") else (r[0], r[1])
		kw = dict(map_kw)
		if use_name:
			kw["name"] = name
		hits = [hit_key(h) for h in aligner.map(seq, **kw)]
		if hits:
			out[name] = hits
	return out


def hit_problems(got, want, names=None, max_report=8):
	"""Readable differences between two {name: [key, ...]} dicts (empty when they agree)."""
	names = sorted(set(got) | set(want)) if names is None else names
	out = []
	for n in names:
		if got.get(n, []) != want.get(n, []):
			out.append("%s:\n  got  %r\n  want %r" % (n, got.get(n, []), want.get(n, [])))
	if len(out) > max_report:
		out = out[:max_report] + ["... %d more" % (len(out) - max_report)]
	return out


def ref_oriented(seq, h):
	"""The bases of `seq` that alignment h covers, in reference orientation."""
	s = seq[h.q_st:h.q_en]
	return s if h.strand > 0 else H.revcomp(s)


_FLAGS = {}


def mm_flags():
	"""{name: value} of the MM_F_* and MM_I_* constants in minimap.h of the sources under test."""
	src = minimap2_build.src_dir()
	if src not in _FLAGS:
		with open(os.path.join(src, "minimap.h")) as f:
			text = f.read()
		d = {}
		for m in re.finditer(r"^#define\s+(MM_[FI]_\w+)\s+\(?(0x[0-9a-fA-F]+|\d+)(LL)?\)?", text, re.M):
			d[m.group(1)] = int(m.group(2), 0)
		_FLAGS[src] = d
	return _FLAGS[src]


def mmi_header(path):
	"""(magic, w, k, b, n_seq, flag) from the first 24 bytes of an .mmi file (MM_IDX_MAGIC then five uint32)."""
	with open(path, "rb") as f:
		head = f.read(24)
	return (head[:4],) + struct.unpack("<5I", head[4:24])


def mappy_env():
	"""Environment for a child Python process that imports the cached mappy build."""
	lib = minimap2_build.mappy_dir()
	return {"PYTHONPATH": lib + os.pathsep + os.environ.get("PYTHONPATH", "")}


def source_text(rel):
	"""Text of a file in the source tree under test (e.g. 'setup.py')."""
	with open(os.path.join(minimap2_build.src_dir(), rel)) as f:
		return f.read()
