'use strict';

// check/genotype.js, worker side (genotyping spec §5.6-§5.8, §7.3): prepare()'s cores and haplotypes, the semi-global alignment,
// the exact-core rule and its gap handling, copy merging and the ortholog test, the megablast fallback over recorded rows
// (check/blast.js megablast arguments and parser), the three-primer prediction, and the §2.13 results block rebuilt from recorded
// assembly windows. Offline and pure; sequences come from test/primers/fixtures/check_core/genotype/.

const test = require('node:test');
const fs = require('fs');
const path = require('path');
const should = require('should');

const stubs = require('../fixtures/check_core/genotype/stubs');
const P = require('../fixtures/check_core/genotype/products');
const W = require('../fixtures/verdicts/fake_world');
const genotype = require('../../../api/helpers/primers/check/genotype');
const blast = require('../../../api/helpers/primers/check/blast');
const classify = require('../../../api/helpers/primers/check/classify');
const variation = require('../../../api/helpers/primers/variation/normalize');

const FIXTURES = path.join(__dirname, '..', 'fixtures', 'check_core', 'genotype');
const EXAMPLE = require(path.join(FIXTURES, 'results_2_13.json'));
const WINDOWS = require(path.join(FIXTURES, 'assembly_windows.json')).windows;
const g = stubs.bases;
const CFG = stubs.makeCfg();
const PARAMS = classify.DEFAULT_PARAMS;
const I = genotype._internal;

// prepare()'s variant block for a variant alone (no sets), over the recorded window.
function variantOnly(vcf) {
  const len = stubs.regionLengthOf('sorghum_bicolor', vcf.region);
  const win = genotype.windowFor(vcf, len);
  const start = Math.max(win.start, stubs.WINDOW.start);
  const end = Math.min(win.end, stubs.WINDOW.end);
  return genotype.prepare({ variant: vcf, sets: [] }, [], variation.sequenceWindow(g(start, end), start, len), { cfg: CFG, links: [] }).variant;
}

const V871 = variantOnly(stubs.VARIANTS.rs871475760);

// The reference segment [segStart, segEnd] aligned to a genome string whose first base is genome coordinate 1.
function readOn(variant, segStart, segEnd, genomeSeq) {
  const al = genotype.semiGlobal(g(segStart, segEnd), genomeSeq);
  return I.readAlignment(al, segStart, variant, function (k) { return 1 + k; });
}

// The reference [lo, hi] with bases from..to replaced by repl.
function edited(lo, hi, from, to, repl) {
  return g(lo, from - 1) + repl + g(to + 1, hi);
}

function megablastRows(systemName) {
  return fs.readFileSync(path.join(FIXTURES, 'megablast_rs871475760.tsv'), 'utf8').split('\n')
    .filter(function (l) { return l && l[0] !== '#' && l.split('\t')[0] === systemName; })
    .map(function (l) { return blast.parseMegablastLine(l.slice(l.indexOf('\t') + 1)); });
}

// One synthetic amplicons.finalize product at 1:100-200.
function amp(o) {
  return Object.assign({ region: '1', start: 100, end: 200, size: 101, strand: 1, orientation: 'LR', likelihood: 'likely', left_mm: 0, right_mm: 0,
    left_3p_mm: 0, right_3p_mm: 0, left_mm_pos: [], right_mm_pos: [], terminal_mismatch: false, approx: false }, o);
}
const FORWARD = { orientation: 'forward', deliberate_mismatch: { as_ref: null, as_alt: null }, in_shift_tract: { as_ref: false, as_alt: false } };
const COPY = { region: '1', start: 50, end: 250 };
const side = function (amplicons, unlikely) { return { amplicons: amplicons || [], unlikely: unlikely || [] }; };

// ---- prepare: cores, flanks and haplotypes ------------------------------------------------------------------------------------

test('prepare: the five cores, flanks and haplotype spans of the §5.6 table; the public variant block of §2.13', function () {
  const rows = [
    [stubs.VARIANTS.rs871475760, 0, [11108, 11110, 'TCT', 'TAT'], 15, [11093, 11125, 33, 33]],
    [{ region: '1', position: 11193, ref: 'C', alt: 'T' }, 0, [11192, 11194, 'TCT', 'TTT'], 15, [11177, 11209, 33, 33]],
    [stubs.VARIANTS.rs5413864115, 2, [11281, 11286, 'TCAAAG', 'TCAAG'], 15, [11266, 11301, 36, 35]],
    [stubs.VARIANTS.tmp_1_11502_C_CGT, 0, [11501, 11503, 'CCA', 'CCGTA'], 15, [11486, 11518, 33, 35]],
    [stubs.VARIANTS.rs5413863413, 13, [13734, 13752, 'ATTGGTGGTGGTGGTGGTA', 'ATTGGTGGTGGTGGTA'], 22, [13712, 13774, 63, 60]]
  ];
  rows.forEach(function (row) {
    const v = variantOnly(row[0]);
    const label = v.key;
    should(v.shift).equal(row[1], label);
    should(v.core).eql({ start: row[2][0], end: row[2][1], ref: row[2][2], alt: row[2][3] }, label);
    should(v.flank).equal(row[3], label);
    should([v.haplotypes.start, v.haplotypes.end, v.haplotypes.ref.length, v.haplotypes.alt.length]).eql(row[4], label);
  });
  should(V871.haplotypes).match({ ref: 'ATAGTCATACTCTATTCTGAATTTCTCGCTAGT', alt: 'ATAGTCATACTCTATTATGAATTTCTCGCTAGT' });
  should(genotype.emptyResults(P.prepared(stubs.BODY_2_11, CFG)).variant).eql(EXAMPLE.variant);
  should(genotype.megablastQuery(P.prepared(stubs.BODY_2_11, CFG))).eql({ region: '1', start: 10894, end: 11324 });
});

// ---- alignment and the core --------------------------------------------------------------------------------------------------

test('semiGlobal: the query is fully consumed, the target ends are free, and a 1 bp indel in the flank leaves the core mapped', function () {
  const query = g(11080, 11140);
  const deletion = 'GGGGGGG' + edited(11080, 11140, 11095, 11095, '') + 'CCCCC';
  let al = genotype.semiGlobal(query, deletion);
  should([al.cost, al.ops.length, al.ops.replace(/M/g, ''), al.qFirst, al.qLast, al.qToT[0], al.qToT[60]]).eql([1, 61, 'I', 0, 60, 7, 66]);
  let read = I.readAlignment(al, 11080, V871, function (k) { return 1 + k; });
  should(read).match({ observed: 'TCT', call: 'ref', flank_edits: 1, variant_position: 7 + 29 - 1 + 1, aligned_length: 61, identity: 98.36, gap_compressed_identity: 98.36 });

  const insertion = 'TTT' + edited(11080, 11140, 11095, 11094, 'A') + 'AAAA';
  al = genotype.semiGlobal(query, insertion);
  should([al.cost, al.ops.length, al.ops.replace(/M/g, '')]).eql([1, 62, 'D']);
  read = I.readAlignment(al, 11080, V871, function (k) { return 1 + k; });
  should(read).match({ observed: 'TCT', call: 'ref', flank_edits: 1, aligned_length: 62 });
  // 6 matches; the I run and the D run count as one column each: 6 / 8
  should(genotype.gapCompressedIdentity('MMMIIIMMDM')).equal(75);
  should(genotype.gapCompressedIdentity({ ops: 'MMMM' })).equal(100);
});

test('the core read through gaps (§5.6 step 3): an end on a gap moves outward; with no aligned base outward the anchor is missing', function () {
  // rs871475760 with S1's R_s 1:11018-11182 on the reference 1:10978-11222 with one base deleted
  should(readOn(V871, 11018, 11182, edited(10978, 11222, 11110, 11110, ''))).match({ call: 'other', observed: 'TCG' });
  should(readOn(V871, 11018, 11182, edited(10978, 11222, 11108, 11108, ''))).match({ call: 'ref', observed: 'TCT' });
  const v115 = variantOnly(stubs.VARIANTS.rs5413864115);
  should(readOn(v115, 11207, 11394, edited(11167, 11434, 11281, 11281, ''))).match({ call: 'ref', observed: 'TCAAAG' });
  const v413 = variantOnly(stubs.VARIANTS.rs5413863413);
  should(readOn(v413, 13682, 13884, edited(13642, 13924, 13752, 13752, ''))).match({ call: 'other', observed: 'ATTGGTGGTGGTGGTGGTG' });
  // ten query bases before the genome's four: the core, on query bases 2-4, lies in the leading gap run
  const al = genotype.semiGlobal('GGGGGGGGGGACGT', 'ACGT');
  should(al.ops).equal('IIIIIIIIIIMMMM');
  const v = { position: 2, core: { start: 2, end: 4, ref: 'GGG', alt: 'GAG' } };
  should(I.readAlignment(al, 1, v, function (k) { return 1 + k; })).match({ call: 'missing', observed: null, flank_edits: null, variant_position: null, locus: null });
  should(I.readCore(al, 1, 3)).equal(null);
});

test('exact-core rule at the repeats: rs5413864115 TCAAAG ref, TCAAG alt, TCAG and TCAAAAG other; rs5413863413 one TGG unit deleted alt, two deleted or one inserted other', function () {
  const v115 = variantOnly(stubs.VARIANTS.rs5413864115);
  const tract115 = function (repl) { return readOn(v115, 11207, 11394, edited(11167, 11434, 11283, 11285, repl)); };
  should([['AAA'], ['AA'], ['A'], ['AAAA']].map(function (r) { const x = tract115(r[0]); return [x.call, x.observed]; })).eql([
    ['ref', 'TCAAAG'], ['alt', 'TCAAG'], ['other', 'TCAG'], ['other', 'TCAAAAG']
  ]);
  const v413 = variantOnly(stubs.VARIANTS.rs5413863413);
  const tract = g(13736, 13751);
  const tract413 = function (repl) { return readOn(v413, 13682, 13884, edited(13642, 13924, 13736, 13751, repl)); };
  should([tract, tract.slice(3), tract.slice(6), tract.slice(0, 3) + tract].map(function (r) { const x = tract413(r); return [x.call, x.observed]; })).eql([
    ['ref', 'ATTGGTGGTGGTGGTGGTA'], ['alt', 'ATTGGTGGTGGTGGTA'], ['other', 'ATTGGTGGTGGTA'], ['other', 'ATTGGTGGTGGTGGTGGTGGTA']
  ]);
});

test('third lengths of the insertion tmp_1_11502_C_CGT (GTGT, G) and a third SNP allele at rs871475760 are other', function () {
  const v = variantOnly(stubs.VARIANTS.tmp_1_11502_C_CGT);
  const ins = function (bases) { const x = readOn(v, 11408, 11579, g(11368, 11502) + bases + g(11503, 11619)); return [x.call, x.observed]; };
  should([ins(''), ins('GT'), ins('GTGT'), ins('G')]).eql([['ref', 'CCA'], ['alt', 'CCGTA'], ['other', 'CCGTGTA'], ['other', 'CCGA']]);
  should(readOn(V871, 11018, 11182, edited(10978, 11222, 11109, 11109, 'T'))).match({ call: 'other', observed: 'TTT', flank_edits: 0 });
});

// ---- copies and orthologs ---------------------------------------------------------------------------------------------------

test('ortholog test: a 35 bp deletion gives identity 94.9 % but gap-compressed 99.85 %, an ortholog; +21 % size is a paralog; ortholog: true wins', function () {
  const seg = g(10500, 11185);
  const al = genotype.semiGlobal(seg, 'GGGGG' + seg.slice(0, 300) + seg.slice(335) + 'CCCCC');
  should([al.cost, al.ops.length, genotype.gapCompressedIdentity(al)]).eql([35, 686, 99.85]);
  const read = I.readAlignment(al, 10500, V871, function (k) { return 1 + k; });
  should(read).match({ identity: 94.9, gap_compressed_identity: 99.85, call: 'ref', flank_edits: 35 });
  const copyOf = function (size, ortholog, extra) {
    return genotype.mergeCopies([Object.assign({}, read, { region: '1', start: 1, end: size, strand: 1, anchor: true, ortholog: ortholog, size: size, ref_size: 100 }, extra || {})])[0];
  };
  should(I.isOrthologous(copyOf(100, null), CFG.check)).be.true();
  should(I.isOrthologous(copyOf(120, null), CFG.check)).be.true();
  should(I.isOrthologous(copyOf(121, null), CFG.check)).be.false();
  should(I.isOrthologous(copyOf(79, false), CFG.check)).be.false();
  should(I.isOrthologous(copyOf(121, true, { matches: 10, gc_columns: 100 }), CFG.check)).be.true();
  should(copyOf(100, true).copy).match({ ortholog: true, anchors: 1, source: 'amplicon' });
});

test('a truncated anchor that aligns over 43 % of R_s is removed by the identity test (step 6), never read as an allele', function () {
  const al = genotype.semiGlobal(g(11018, 11182), g(10980, 11060));
  const read = I.readAlignment(al, 11018, V871, function (k) { return 10980 + k; });
  should([read.identity, read.aligned_length]).eql([43.03, 165]);
  const group = genotype.mergeCopies([Object.assign(read, { region: '1', start: 11068, end: 11132, strand: 1, anchor: true, ortholog: null, size: 65, ref_size: 65 })])[0];
  should(I.isOrthologous(group, CFG.check)).be.false();
});

test('copy merging and every primer call of §2.13: the recorded windows of bicolorv5, pi180348 (two copies) and pi329250 rebuild results.genotyping exactly', async function () {
  const prepared = P.prepared(stubs.BODY_2_11, CFG);
  const call = function (name, display, isReference, products, deps) {
    return genotype.callGenome({ prepared: prepared, system_name: name, display_name: display, is_reference: isReference, products: products, cfg: CFG.check, params: PARAMS }, deps);
  };
  const win = function (name, start) {
    const w = WINDOWS.find(function (x) { return x.system_name === name && x.start === start; });
    return P.genome(w.region, w.seq, w.start);
  };
  const oneCopy = async function (name, display, genome, isReference) {
    return call(name, display, isReference, prepared.sets.map(function (set) { return P.productsForSet(genome, set, 1); }), P.fetchFrom([genome]));
  };
  const reference = await oneCopy('sorghum_bicolor', 'Sb bicolor BTx623 v3', P.genome('1', g(9000, 15500), 9000), true);
  const v5 = await oneCopy('sorghum_bicolorv5', 'Sb bicolor BTx623 v5', win('sorghum_bicolorv5', 31300), false);
  const pi329250 = await oneCopy('sorghum_pi329250', 'Sb verticilliflorum PI329250 ADAR', win('sorghum_pi329250', 203700), false);
  // pi180348: the inverted copy (1:14800-15400, products RL) and the direct one (1:38100-38700); both windows are region 1
  const minus = win('sorghum_pi180348', 14800);
  const plus = win('sorghum_pi180348', 38100);
  const products = prepared.sets.map(function (set) {
    return P.mergeProducts(P.productsForSet(minus, set, -1), P.productsForSet(plus, set, 1));
  });
  const both = {
    fetch: async function (region, start, end) {
      const w = [minus, plus].find(function (x) { return start >= x.start && end <= x.end; });
      if (region !== '1' || !w) throw new Error('no bases for ' + region + ':' + start + '-' + end);
      return w.seq.slice(start - w.start, end - w.start + 1);
    },
    regionLength: async function () { return 84410038; }
  };
  const pi180348 = await call('sorghum_pi180348', 'Sb bicolor PI180348 Juar (IS 12876)', false, products, both);

  should(v5.genome.copies).eql(EXAMPLE.genomes[1].copies);
  should(pi180348.genome.copies.map(function (c) { return [c.start, c.end, c.strand, c.aligned_length, c.identity, c.flank_edits]; })).eql([
    [15028, 15132, -1, 192, 98.44, 2], [38342, 38446, 1, 193, 98.45, 2]
  ]);
  should([v5.megablast, pi180348.megablast, v5.anchors, pi180348.anchors, v5.failed_reads]).eql([null, null, 2, 4, 0]);
  const block = genotype.emptyResults(prepared);
  genotype.writeResults(block, {
    reference: reference,
    genomes: [v5, pi180348, pi329250],
    specificity: prepared.sets.map(function () { return { ref_verdict: 'specific', alt_verdict: 'specific', off_target_count: 0 }; })
  });
  should(block).eql(EXAMPLE);
  // pi180348's KASP common primer reports likely from the ALT pair's product, although the REF pair's products come first
  should(block.sets[0].genomes[1].common_primer).eql({ status: 'match', likelihood: 'likely', mm_pos: [], residual_mm_pos: [] });
});

test('shift tracts (§4.12, §5.7): REF ref; ALT unknown with shift_tract_uncertain; third repeat lengths are other, whatever the primers predict', async function () {
  const rows = [
    ['rs5413864115', [['REF', g(11283, 11285), 'ref', 'ref', true], ['ALT', g(11284, 11285), 'alt', 'unknown', null],
      ['two-A deletion', 'A', 'other', 'alt', false], ['AAAA', 'AAAA', 'other', 'ref', false]], [11283, 11285], 11282],
    ['rs5413863413', [['REF', g(13736, 13751), 'ref', 'ref', true], ['ALT', g(13739, 13751), 'alt', 'unknown', null],
      ['two units deleted', g(13742, 13751), 'other', 'none', true], ['one unit inserted', g(13736, 13738) + g(13736, 13751), 'other', 'ref', false]], [13736, 13751], 13735]
  ];
  for (const row of rows) {
    const body = stubs.body(stubs.VARIANTS[row[0]], stubs.SETS[row[0]]);
    const prepared = P.prepared(body, CFG);
    prepared.sets.forEach(function (s) { should(s.in_shift_tract).eql({ as_ref: true, as_alt: false }, row[0] + ' ' + s.id); });
    for (const gen of row[1]) {
      const seq = W.randomSeq(1000, 91) + g(row[3] - 350, row[2][0] - 1) + gen[1] + g(row[2][1] + 1, row[3] + 350) + W.randomSeq(1000, 92);
      const genome = P.genome('x', seq, 1);
      for (const set of prepared.sets) {
        const products = prepared.sets.map(function (s) { return P.productsForSet(genome, s, 1); });
        const out = await genotype.callGenome({ prepared: prepared, system_name: 'x', products: products, cfg: CFG.check, params: PARAMS }, P.fetchFrom([genome]));
        const c = out.sets[prepared.sets.indexOf(set)];
        const label = row[0] + ' ' + set.id + ' ' + gen[0];
        should([out.genome.allele, c.predicted, c.agrees]).eql([gen[2], gen[3], gen[4]], label);
        if (gen[0] === 'ALT') should([c.ref_primer.status, c.reasons]).eql(['uncertain', ['shift_tract_uncertain']], label);
        if (gen[0] === 'REF') should([c.ref_primer.status, c.alt_primer.status, c.common_primer.status]).eql(['match', 'terminal_mismatch', 'match'], label);
        if (gen[0] === 'two-A deletion') should(c.strength).equal('weak', label);
      }
    }
  }
});

// ---- the megablast fallback ---------------------------------------------------------------------------------------------------

test('megablast rows: the recorded pi180348 HSPs are two ALT copies; pi154987 has four; pi536008 (a scaffold end) and is36143 keep none', function () {
  const prepared = P.prepared(stubs.BODY_2_11, CFG);
  const query = genotype.megablastQuery(prepared);
  const res = function (name) { return genotype.megablastCopies(megablastRows(name), query, prepared, CFG.check); };
  const pi180348 = res('sorghum_pi180348');
  should([pi180348.hsps, pi180348.kept, pi180348.covering]).eql([5, 2, 2]);
  should(pi180348.copies.map(function (c) { return c.copy; })).eql([
    { region: '1', start: 14876, end: 15306, strand: -1, variant_position: 15091, identity: 99.3, gap_compressed_identity: 99.3, aligned_length: 431,
      observed: 'TAT', flank_edits: 2, call: 'alt', anchors: 0, ortholog: null, source: 'megablast' },
    { region: '1', start: 38167, end: 38598, strand: 1, variant_position: 38383, identity: 99.31, gap_compressed_identity: 99.31, aligned_length: 432,
      observed: 'TAT', flank_edits: 2, call: 'alt', anchors: 0, ortholog: null, source: 'megablast' }
  ]);
  should(res('sorghum_pi154987').copies.map(function (c) { return [c.copy.region, c.copy.strand, c.copy.call]; })).eql([
    ['scaffold_1532', -1, 'alt'], ['scaffold_1532', 1, 'alt'], ['1', 1, 'alt'], ['1', -1, 'alt']
  ]);
  should(res('sorghum_bicolorv5').copies.map(function (c) { return [c.copy.start, c.copy.variant_position, c.copy.call, c.copy.identity]; })).eql([[31349, 31564, 'ref', 100]]);
  // pi536008's only locus HSP covers query bases 160-431 (272 of 431, below the 0.8 cover); is36143 has only 94-98 bp repeat hits
  should([res('sorghum_pi536008').kept, res('sorghum_is36143').kept]).eql([0, 0]);
  should(megablastRows('sorghum_pi536008')[0]).match({ sseqid: 'scaffold2100', sstart: 3908, send: 3637, strand: -1, qstart: 160, qend: 431, length: 272 });
});

test('callGenome without an orthologous copy: megablast calls it (source megablast); budget and failure leave it missing with their reasons', async function () {
  const prepared = P.prepared(stubs.BODY_2_11, CFG);
  const input = { prepared: prepared, system_name: 'sorghum_pi180348', display_name: 'PI180348', products: [], cfg: CFG.check, params: PARAMS };
  const called = await genotype.callGenome(input, { megablast: async function () { return { status: 'ok', rows: megablastRows('sorghum_pi180348'), query: genotype.megablastQuery(prepared) }; } });
  should(called.genome).match({ allele: 'alt', observed: 'TAT', source: 'megablast', orthologous_copies: 2, paralog_copies: 0, reason: null });
  should(called.sets.map(function (c) { return [c.predicted, c.agrees]; })).eql([['none', false], ['none', false]]);
  should(called.megablast).equal('ok');
  const outcome = async function (mb) { return (await genotype.callGenome(input, mb ? { megablast: async function () { return mb; } } : {})).genome; };
  should(await outcome({ status: 'budget' })).match({ allele: 'missing', source: null, copies: [], reason: 'fallback_budget' });
  should(await outcome({ status: 'failed' })).match({ allele: 'missing', reason: 'fallback_failed' });
  should(await outcome({ status: 'ok', rows: megablastRows('sorghum_is36143'), query: genotype.megablastQuery(prepared) })).match({ allele: 'missing', reason: 'no_orthologous_copy' });
  should(await outcome({ status: 'ok', rows: megablastRows('sorghum_pi536008'), query: genotype.megablastQuery(prepared) })).match({ allele: 'missing', reason: 'no_orthologous_copy' });
  should(await outcome(null)).match({ allele: 'missing', reason: 'no_orthologous_copy' });
});

test('check/blast.js megablast: the §5.6 arguments, the query FASTA and the 11-column parser', function () {
  should(blast.buildMegablastArgs({ db: '/scratch/x/Prefix.dna.toplevel.nal' })).eql(['-task', 'megablast', '-db', '/scratch/x/Prefix.dna.toplevel', '-query', '-',
    '-dust', 'no', '-soft_masking', 'false', '-evalue', '1e-20', '-max_target_seqs', '50', '-max_hsps', '10', '-num_threads', '1',
    '-outfmt', '6 sseqid sstart send sstrand pident length bitscore qstart qend qseq sseq']);
  should(blast.GENOTYPE_OUTFMT).equal('6 sseqid sstart send sstrand pident length bitscore qstart qend qseq sseq');
  should(function () { blast.buildMegablastArgs({ db: '-remote' }); }).throw(TypeError);
  should(function () { blast.buildMegablastArgs({ db: 'a b' }); }).throw(TypeError);
  should(function () { blast.buildMegablastArgs({ db: 'x', threads: 0 }); }).throw(TypeError);
  should(blast.megablastQueryFasta('acgtN')).equal('>genotype\nACGTN\n');
  should(function () { blast.megablastQueryFasta('ACGT\n>q1\nAC'); }).throw(TypeError);
  should(function () { blast.megablastQueryFasta(''); }).throw(TypeError);
  should(blast.parseMegablastLine('1\t38167\t38170\tplus\t75.000\t4\t1.19e+03\t10\t12\tAC-G\tACTG')).eql({
    sseqid: '1', sstart: 38167, send: 38170, strand: 1, pident: 75, length: 4, bitscore: 1190, qstart: 10, qend: 12, qseq: 'AC-G', sseq: 'ACTG'
  });
  should(blast.parseMegablastLine('')).equal(null);
  should(blast.parseMegablastLine('# comment')).equal(null);
  [
    '1\t1\t4\tplus\t100\t4\t8\t1\t4\tACGT',
    '1\t1\t4\tboth\t100\t4\t8\t1\t4\tACGT\tACGT',
    '1\t4\t1\tplus\t100\t4\t8\t1\t4\tACGT\tACGT',
    '1\t1\t4\tplus\t100\t4\t8\t1\t4\tACGT\tACG',
    '1\t1\t4\tplus\tn/a\t4\t8\t1\t4\tACGT\tACGT',
    '\t1\t4\tplus\t100\t4\t8\t1\t4\tACGT\tACGT'
  ].forEach(function (line) { should(function () { blast.parseMegablastLine(line); }).throw({ code: 'BLAST_PARSE' }); });
});

// ---- prediction rows and results ----------------------------------------------------------------------------------------------

test('predict: the rows the fake world does not reach (approx alignment, third allele, off-locus dyes alone, a declared -3, a blocked common site)', function () {
  const call = function (products, allele, copies, set) {
    return genotype.predict(set || FORWARD, products, { allele: allele, copies: copies === undefined ? [COPY] : copies }, PARAMS);
  };
  // row 4: an approximate allele-specific alignment
  let c = call({ ref: side([amp({ left_mm_pos: null, left_mm: 1 })]), alt: side([amp({ left_mm_pos: [1], left_mm: 1, likelihood: 'likely_weak' })]) }, 'ref');
  should([c.ref_primer.status, c.predicted, c.reasons, c.agrees]).eql(['unknown', 'unknown', ['approx_alignment'], null]);
  // row 7 on a third allele
  c = call({ ref: side([amp({})]), alt: side([amp({})]) }, 'other');
  should([c.predicted, c.strength, c.reasons, c.agrees]).eql(['both', 'normal', ['third_allele'], false]);
  // no orthologous copy: an amplifying ALT-pair paralog alone lights its dye
  c = call({ ref: side([], [amp({ left_mm_pos: [2, 1], left_mm: 2, likelihood: 'unlikely' })]), alt: side([amp({})]) }, 'missing', []);
  should([c.predicted, c.reasons, c.off_locus_products, c.agrees, c.ref_primer.status]).eql(['alt', ['alt_signal_off_locus', 'no_orthologous_copy'], 1, false, 'no_product']);
  // a deliberate -3 on both allele-specific primers (§5.7): [3] -> residual [] match; [3,1] on an unlikely product -> blocked
  const minus3 = Object.assign({}, FORWARD, { deliberate_mismatch: { as_ref: 3, as_alt: 3 } });
  c = call({ ref: side([amp({ left_mm_pos: [3], left_mm: 1 })]), alt: side([], [amp({ left_mm_pos: [3, 1], left_mm: 2, likelihood: 'unlikely' })]) }, 'ref', undefined, minus3);
  should([c.ref_primer, c.alt_primer, c.predicted]).eql([
    { status: 'match', likelihood: 'likely', mm_pos: [3], residual_mm_pos: [] },
    { status: 'blocked', likelihood: 'unlikely', mm_pos: [3, 1], residual_mm_pos: [1] },
    'ref'
  ]);
  // a common site with 2 mismatches in its 3' window is blocked on its own: no_call even though the product likelihood is unlikely
  c = call({ ref: side([], [amp({ right_mm_pos: [2, 1], right_mm: 2, likelihood: 'unlikely' })]), alt: side() }, 'ref');
  should([c.common_primer.status, c.predicted, c.reasons, c.agrees]).eql(['blocked', 'no_call', ['common_primer_3p_mismatch'], null]);
  // unavailable
  should(genotype.predict(FORWARD, {}, { allele: 'unavailable', copies: [] }, PARAMS)).match({ predicted: 'unknown', agrees: null, off_locus_products: 0 });
  // LL/RR products are never read
  c = call({ ref: side([amp({ orientation: 'LL' })]), alt: side() }, 'ref');
  should([c.ref_primer.status, c.predicted]).eql(['no_product', 'none']);
});

test('writeResults: summaries count pan-genome genomes only and add up; the control warns on a weak or off-locus prediction and fails otherwise', function () {
  const prepared = P.prepared(stubs.BODY_2_11, CFG);
  const entry = function (name, allele, predicted, extra) {
    return {
      genome: { system_name: name, display_name: name, is_reference: name === 'ref', allele: allele, observed: null, source: null, copies: [], orthologous_copies: 0, paralog_copies: 0, reason: null },
      sets: prepared.sets.map(function () {
        return Object.assign({ system_name: name, ref_primer: {}, alt_primer: {}, common_primer: {}, predicted: predicted, strength: null, agrees: null, reasons: [], off_locus_products: 0 }, extra || {});
      })
    };
  };
  const block = genotype.emptyResults(prepared);
  should(block.sets.map(function (s) { return [s.specificity, s.control, s.reference, s.genomes]; })).eql([[null, null, null, []], [null, null, null, []]]);
  genotype.writeResults(block, {
    reference: entry('ref', 'ref', 'ref', { strength: 'weak', reasons: ['common_primer_weak'], agrees: true }),
    genomes: [entry('a', 'alt', 'alt', { agrees: true }), entry('b', 'other', 'both', { agrees: false }), genotype.unavailableEntry(prepared, { system_name: 'c' }, 'blast_error')],
    specificity: null
  });
  should(block.summary).eql({ genomes_total: 3, ref: 0, alt: 1, other: 1, ambiguous: 0, missing: 0, unavailable: 1 });
  should(block.genomes.map(function (x) { return x.system_name; })).eql(['ref', 'a', 'b', 'c']);
  should(block.sets[0].summary).eql({ genomes_total: 3, predicted_ref: 0, predicted_alt: 1, both: 1, none: 0, no_call: 0, unknown: 1, weak: 0, agree: 1, disagree: 1, not_comparable: 1 });
  should(block.sets[0].control).eql({ status: 'warn', allele: 'ref', reasons: ['weak', 'common_primer_weak'] });
  should(block.sets[0].specificity).equal(null);
  should(I.referenceControl({ allele: 'ref' }, { predicted: 'ref', strength: 'normal', off_locus_products: 1, reasons: ['ref_signal_off_locus'] })).eql({ status: 'warn', allele: 'ref', reasons: ['off_locus_products', 'ref_signal_off_locus'] });
  should(I.referenceControl({ allele: 'other' }, { predicted: 'none', strength: null, off_locus_products: 0, reasons: [] })).eql({ status: 'fail', allele: 'other', reasons: ['allele_not_ref', 'prediction_not_ref'] });
  should(genotype.unavailableEntry(prepared, { system_name: 'c', display_name: 'C', is_reference: false }, 'no_such_reason').genome).match({ allele: 'unavailable', reason: 'call_failed', display_name: 'C' });
});

test('fillSequences reads a segment and haplotypes that ran past the submit-time window', async function () {
  const prepared = P.prepared(stubs.BODY_2_11, CFG);
  const want = JSON.parse(JSON.stringify([prepared.sets[0].segment, prepared.variant.core, prepared.variant.haplotypes]));
  prepared.sets[0].segment.sequence = null;
  prepared.variant.haplotypes.ref = null;
  prepared.variant.haplotypes.alt = null;
  const calls = [];
  await genotype.fillSequences(prepared, async function (region, start, end) { calls.push([region, start, end]); return g(start, end).toLowerCase(); });
  should([prepared.sets[0].segment, prepared.variant.core, prepared.variant.haplotypes]).eql(want);
  should(calls).eql([['1', 11018, 11182], ['1', 11093, 11125]]);
});
