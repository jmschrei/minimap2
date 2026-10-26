"""Area B: the mapping parameters of the long-read and assembly presets (-x lr map-ont ava-ont map10k
map-pb ava-pb lr:hq map-hifi map-ccs lr:hqae map-iclr-prerender map-iclr asm5 asm10 asm20).

Each preset must equal the options minimap2.1 lists for it (lr:hqae and map-iclr-prerender, which
minimap2.1 does not list, against their definition in options.c), its aliases must be identical, and
its reads must land at their simulated loci. Pinned outputs on fixed simulated inputs guard the
mapping code paths of each preset byte for byte."""

import os

import pytest

import minimap2_helpers as H


# minimap2.1's expansion of each preset; options apply on top of the defaults (map-ont)
EXPANSION = {
	"map-ont": "",
	"lr:hq": "-k19 -w19 -U50,500 -g10k",
	"map-hifi": "-k19 -w19 -U50,500 -g10k -A1 -B4 -O6,26 -E2,1 -s200",
	"map-pb": "-Hk19",
	"map-iclr": "-k19 -B6 -b4 -O10,50",
	"asm5": "-k19 -w19 -U50,500 --rmq -r1k,100k -g10k -A1 -B19 -O39,81 -E3,1 -s200 -z200 -N50",
	"asm10": "-k19 -w19 -U50,500 --rmq -r1k,100k -g10k -A1 -B9 -O16,41 -E2,1 -s200 -z200 -N50",
	"asm20": "-k19 -w10 -U50,500 --rmq -r1k,100k -g10k -A1 -B4 -O6,26 -E2,1 -s200 -z200 -N50",
	"ava-pb": "-Hk19 -Xw5 -e0 -m100",
	"ava-ont": "-k15 -Xw5 -e0 -m100 -r2k",
	# not in minimap2.1; from options.c
	"lr:hqae": "-k25 -w51 --rmq -U50,500 --rmq-inner 5000 -e200 -N100 --chain-gap-scale 5",
	"map-iclr-prerender": "-k15 -B6 -b1 -O10,50",
}

READS = {"map-ont": "ont", "lr:hq": "hifi", "map-hifi": "hifi", "map-pb": "clr", "map-iclr": "hifi", "asm5": "ctg",
	"asm10": "ctg", "asm20": "ctg", "ava-pb": "clr", "ava-ont": "ont", "lr:hqae": "hifi_u", "map-iclr-prerender": "hifi"}


@pytest.fixture(scope="module")
def sim(synth, tmp_path_factory):
	"""Reads on the standard synthetic genome: ONT-like, HiFi-like and CLR-like (25 each, 1-4 kbp,
	anywhere, so some fall in repeats, inverted repeats, N runs and low-complexity tracts), 15
	HiFi-like reads from unique sequence (for lr:hqae, which is slow on low-complexity reads) and
	three assembly-like contigs of 15-25 kbp, each with one inverted piece."""
	d = tmp_path_factory.mktemp("b_presets")
	g = synth.genome
	sets = {
		"ont": H.simulate_long_reads(g, 25, profile="ont", seed=201, length=(1000, 4000)),
		"hifi": H.simulate_long_reads(g, 25, profile="hifi", seed=202, length=(1000, 4000)),
		"clr": H.simulate_long_reads(g, 25, profile="clr", seed=203, length=(1000, 4000)),
		"hifi_u": H.simulate_long_reads(g, 15, profile="hifi", seed=205, length=(1000, 4000), unique=True),
		"ctg": H.simulate_contigs(g, 3, length=(15000, 25000), seed=204, inversions=1, inversion_len=(800, 1500)),
	}
	out = {}
	for name, reads in sets.items():
		out[name] = (H.write_fasta(str(d / (name + ".fa")), reads), reads)
	return out


@pytest.mark.parametrize("preset", sorted(EXPANSION))
def test_preset_equals_documented_options(synth, sim, preset):
	"""-x PRESET gives byte-identical PAF with -c (chaining tags, CIGAR, alignment scores) to the
	default preset plus the options minimap2.1 lists for it; for the all-vs-all presets, PAF without -c
	of the reads against themselves."""
	fa, _ = sim[READS[preset]]
	target = fa if preset.startswith("ava") else synth.fasta
	modes = [[]] if preset.startswith("ava") else [["-c"]]
	for extra in modes:
		a = H.map_reads(target, fa, extra, preset=preset)
		b = H.map_reads(target, fa, EXPANSION[preset].split() + extra)
		assert a.stdout.strip(), (preset, extra)
		assert H.strip_pg(a.stdout) == H.strip_pg(b.stdout), (preset, extra)


@pytest.mark.parametrize("alias,preset", [("lr", "map-ont"), ("map10k", "map-pb"), ("map-ccs", "map-hifi")])
def test_preset_aliases(synth, sim, alias, preset):
	"""lr, map10k and map-ccs are aliases (options.c): their SAM output equals map-ont's, map-pb's and
	map-hifi's (@PG aside); map-ont is also the default (no -x)."""
	fa, _ = sim[READS[preset]]
	a = H.map_reads(synth.fasta, fa, ["-a"], preset=alias)
	b = H.map_reads(synth.fasta, fa, ["-a"], preset=preset)
	assert H.strip_pg(a.stdout) == H.strip_pg(b.stdout) and a.stdout.count("\n") > 25
	if alias == "lr":
		assert H.strip_pg(H.map_reads(synth.fasta, fa, ["-a"]).stdout) == H.strip_pg(a.stdout)


@pytest.mark.parametrize("preset", ["map-ont", "map-pb", "lr:hq", "map-hifi", "lr:hqae", "map-iclr", "map-iclr-prerender"])
def test_preset_reads_at_truth(synth, sim, preset):
	"""Reads from unique sequence (no repeat, inverted repeat, N run or low-complexity tract within
	the read) map primarily to their simulated contig and strand with both ends within 30 bp (-c), and
	every SAM record is consistent with the genome (CIGAR, SEQ, NM, MD)."""
	fa, reads = sim[READS[preset]]
	paf = H.map_reads(synth.fasta, fa, ["-c"], preset=preset).paf()
	prim = dict((p.qname, p) for p in paf if p.is_primary)
	n = 0
	for r in reads:
		if not synth.genome.is_unique(r.contig, r.pos, r.end):
			continue
		n += 1
		err = H.paf_placement_error(prim[r.name], r)
		assert err is not None and err <= 30, (preset, prim[r.name].line[:150], r)
	assert n >= 10
	sam = H.map_reads(synth.fasta, fa, ["-a", "--MD"], preset=preset).sam()
	assert H.sam_invariant_errors(sam, synth.genome, dict((r.name, [r.seq]) for r in reads)) == []


@pytest.mark.parametrize("preset", ["asm5", "asm10", "asm20"])
def test_asm_contigs_at_truth(synth, sim, preset):
	"""Assembly presets on contigs with one inverted piece: every simulated piece longer than 1 kbp is
	covered by an alignment on its contig and strand whose ends lie within 100 bp of the piece's (-c;
	asm5's gap penalties can clip a piece's end at a nearby indel), and the alignments agree with the
	genome."""
	fa, contigs = sim["ctg"]
	paf = H.map_reads(synth.fasta, fa, ["-c"], preset=preset).paf()
	assert H.paf_invariant_errors(paf, synth.genome, dict((c.name, c.seq) for c in contigs)) == []
	for c in contigs:
		for s in c.segments:
			if s.end - s.pos < 1000:
				continue
			hit = [p for p in paf if p.qname == c.name and p.tname == s.contig and p.strand == s.strand
				and abs(p.tstart - s.pos) <= 100 and abs(p.tend - s.end) <= 100]
			assert hit, (preset, c.name, s, [p.line[:120] for p in paf if p.qname == c.name])


# Pinned outputs: (pin name, preset, read set, extra arguments)
PINS = [
	("map-ont.sam.gz", "map-ont", "ont", ["-a", "--cs"]),
	("map-hifi.paf.gz", "map-hifi", "hifi", ["-c", "--cs"]),
	("map-pb.sam.gz", "map-pb", "clr", ["-a"]),
	("lr-hq.paf.gz", "lr:hq", "hifi", ["-c", "--MD"]),
	("lr-hqae.paf.gz", "lr:hqae", "hifi", ["-c"]),  # includes the slow low-complexity reads
	("map-iclr.sam.gz", "map-iclr", "hifi", ["-a"]),
	("asm5.paf.gz", "asm5", "ctg", ["-c", "--cs"]),
	("asm20.sam.gz", "asm20", "ctg", ["-a"]),
	("map-ont-nocigar.paf.gz", "map-ont", "ont", []),
]


@pytest.mark.pinned
@pytest.mark.parametrize("name,preset,reads,extra", PINS, ids=[p[0] for p in PINS])
def test_preset_output_pinned(synth, sim, name, preset, reads, extra):
	"""Pins the current output of base minimap2 for each preset on the simulated reads (@PG removed)."""
	fa, _ = sim[reads]
	H.check_pinned("b/" + name, H.strip_pg(H.map_reads(synth.fasta, fa, extra, preset=preset).stdout))


@pytest.mark.pinned
@pytest.mark.parametrize("preset", ["ava-ont", "ava-pb"])
def test_ava_output_pinned(sim, preset):
	"""Pins the current all-vs-all output of base minimap2 (reads against themselves)."""
	fa, _ = sim[READS[preset]]
	H.check_pinned("b/%s.paf.gz" % preset, H.map_reads(fa, fa, [], preset=preset).stdout)
