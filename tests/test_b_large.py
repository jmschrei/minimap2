"""A moderately large case, mapped with one thread and with eight.

helpers_b.large_genome() is 1.5 Mbp with segmental duplications, inverted repeats, low-complexity
tracts, N runs and a 150-copy family whose minimizers exceed the default occurrence threshold, so
high-occurrence seeds are filtered and sampled (-f, -e, --q-occ-frac) on reads that overlap it.
ONT-like and HiFi-like reads are mapped at -t 1, at -t 8 and at -t 8 with a small -K (several
mini-batches in flight); the outputs must be byte-identical, reads from unique sequence must land at
their simulated loci, every record must agree with the genome, and the -t 1 output is pinned."""

import pytest

import helpers_b as B
import minimap2_helpers as H


@pytest.fixture(scope="module")
def large(tmp_path_factory):
	d = tmp_path_factory.mktemp("b_large")
	g = B.large_genome()
	ont = H.simulate_long_reads(g, 60, profile="ont", seed=93)
	hifi = H.simulate_long_reads(g, 40, profile="hifi", seed=94)
	return dict(genome=g, fasta=g.write_fasta(str(d / "large.fa")), ont=ont, hifi=hifi,
		ont_fa=H.write_fasta(str(d / "ont.fa"), ont), hifi_fa=H.write_fasta(str(d / "hifi.fa"), hifi))


def _three_ways(target, reads, args, preset):
	"""Output at -t 1, -t 8 and -t 8 -K 50k (@PG removed); fails unless all three are identical."""
	outs = [H.strip_pg(H.map_reads(target, reads, args + extra, preset=preset).stdout)
		for extra in (["-t", "1"], ["-t", "8"], ["-t", "8", "-K", "50k"])]
	assert outs[0] == outs[1], "-t 1 and -t 8 differ"
	assert outs[0] == outs[2], "-t 1 and -t 8 -K 50k differ"
	return outs[0]


def _family_reads(g, reads):
	return [r for r in reads if g.features_overlapping(r.contig, r.pos, r.end, kinds=("repeat",)) and r.contig == "chr3"]


@pytest.mark.pinned
@pytest.mark.thread_invariance
def test_large_map_ont_sam(large):
	"""Pins the current output (SAM, -a --cs, SEQ and QUAL replaced by '*') of 60 ONT-like reads under
	map-ont, identical at -t 1, -t 8 and -t 8 -K 50k. Independently: one primary per read, every record
	consistent with the genome (CIGAR span, NM, cs), reads from unique sequence placed within 30 bp of
	their simulated locus, and some reads overlapping the 150-copy family report a repetitive length."""
	g, reads = large["genome"], large["ont"]
	sam_text = _three_ways(large["fasta"], large["ont_fa"], ["-a", "--cs"], "map-ont")
	sam = H.parse_sam(sam_text)
	assert H.sam_invariant_errors(sam, g, dict((r.name, [r.seq]) for r in reads)) == []
	prim = sam.primaries()
	assert sorted(k[0] for k in prim) == sorted(r.name for r in reads)
	n = 0
	for r in reads:
		if g.is_unique(r.contig, r.pos, r.end):
			n += 1
			err = H.placement_error(prim[(r.name, 0)], r)
			assert err is not None and err <= 30, (prim[(r.name, 0)].line[:150], r)
	assert n >= 30
	assert len(_family_reads(g, reads)) >= 4
	body = []
	for line in H.sam_body(sam_text).splitlines():
		f = line.split("\t")
		f[9] = f[10] = "*"
		body.append("\t".join(f) + "\n")
	H.check_pinned("b/large-map-ont.sam.gz", "".join(body))


@pytest.mark.pinned
@pytest.mark.thread_invariance
def test_large_map_hifi_paf(large):
	"""Pins the current output (PAF, -c --cs) of 40 HiFi-like reads under map-hifi, identical at -t 1,
	-t 8 and -t 8 -K 50k. Independently: every line consistent with the genome (cg spans, columns 10 and
	11, NM, cs, de), reads from unique sequence placed within 10 bp of their simulated locus, and the
	reads overlapping the 150-copy family report a repetitive length (rl:i) on their primary."""
	g, reads = large["genome"], large["hifi"]
	paf = H.parse_paf(_three_ways(large["fasta"], large["hifi_fa"], ["-c", "--cs"], "map-hifi"))
	assert H.paf_invariant_errors(paf, g, dict((r.name, r.seq) for r in reads)) == []
	prim = {}  # the longest tp:A:P line of each read
	for p in paf:
		if p.is_primary and (p.qname not in prim or p.qend - p.qstart > prim[p.qname].qend - prim[p.qname].qstart):
			prim[p.qname] = p
	assert sorted(prim) == sorted(r.name for r in reads)
	n = 0
	for r in reads:
		if g.is_unique(r.contig, r.pos, r.end):
			n += 1
			err = H.paf_placement_error(prim[r.name], r)
			assert err is not None and err <= 10, (prim[r.name].line[:150], r)
	assert n >= 20
	assert any(prim[r.name].tags["rl"] > 0 for r in _family_reads(g, reads) if r.name in prim)
	H.check_pinned("b/large-map-hifi.paf.gz", "".join(p.line + "\n" for p in paf))
