'use strict';

const test = require('node:test');
const should = require('should');
const fs = require('fs');
const path = require('path');

const normalize = require('../../../api/helpers/primers/variation/normalize');

// Spec §3.4, §3.5, §3.9, §4.17 and the variation_normalize rows of §7.3. Offline, from recorded fixtures
// (captured 2026-09-15, read-only):
//   design/sorghum_bicolor_1_10500-12100.plus.txt, _13400-14000.plus.txt: template.seq of
//     POST 127.0.0.1:50111/sorghum_v11/primers/design {"mode":"region","system_name":"sorghum_bicolor",
//     "region":{"region":"1","start":S,"end":E,"strand":1},"template_only":true}
//     (byte-identical to the spec's fix-pass copies, and to design/sorghum_bicolor_1_11080-15099.plus.txt)
//   variation/overlap_1_<start>-<end>.json: GET https://data.gramene.org/pansite-ensembl-115/overlap/region/
//     sorghum_bicolor/1:<start>-<end>?feature=variation, raw bodies (the §4.3 example template windows, and 13600-13820)
//   variation/variation_<id>.json: GET .../variation/sorghum_bicolor/<id>; variation_not_found.json is rs0000000001 (HTTP 400)

const FIX = path.join(__dirname, '..', 'fixtures');
const CHR1_LENGTH = 80884392;

function windowFixture(start, end) {
  const seq = fs.readFileSync(path.join(FIX, 'design', 'sorghum_bicolor_1_' + start + '-' + end + '.plus.txt'), 'utf8').trim();
  return normalize.sequenceWindow(seq, start, CHR1_LENGTH);
}

function variationFixture(name) {
  return JSON.parse(fs.readFileSync(path.join(FIX, 'variation', name), 'utf8'));
}

const G1 = windowFixture(10500, 12100);
const G2 = windowFixture(13400, 14000);

// The four windows overlap, so most records arrive two to four times, as from overlapping client chunks.
const RECORDS = [].concat.apply([], ['10709-11509', '10793-11593', '10883-11685', '11102-11903'].map(function (w) {
  return variationFixture('overlap_1_' + w + '.json');
}));
const RECORDS_13 = variationFixture('overlap_1_13600-13820.json');

function record(id, records) {
  const r = (records || RECORDS).find(function (x) { return x.id === id; });
  should.exist(r, id);
  return r;
}

function entriesOf(records, genome, opts) {
  return normalize.recordsToEntries(records, genome, Object.assign({ region: '1' }, opts)).entries;
}

function only(records, genome, opts) {
  const list = entriesOf(records, genome, opts);
  list.should.have.length(1);
  return list[0];
}

let allCache = null;
function all() {
  if (!allCache) allCache = entriesOf(RECORDS, G1);
  return allCache;
}

function byKey(key) {
  const e = all().find(function (x) { return x.key === key; });
  should.exist(e, key);
  return e;
}

function site(position, refBase, altBase, altMapsTo) {
  return { position: position, ref_base: refBase, alt_base: altBase, alt_maps_to: altMapsTo };
}

function throwsWith(fn, status, code, details) {
  let err = null;
  try {
    fn();
  } catch (e) {
    err = e;
  }
  should.exist(err, 'expected ' + code);
  err.should.have.property('status', status);
  err.should.have.property('code', code);
  if (details) err.details.should.eql(details);
  return err;
}

// ALT-haplotype span [qStart, qEnd] -> reference blocks through altToRef (as §4.9 genomic.blocks).
function altBlocks(v, qStart, qEnd) {
  const blocks = [];
  let inserted = 0;
  for (let q = qStart; q <= qEnd; q++) {
    const r = normalize.altToRef(v, q);
    if (r === null) {
      inserted += 1;
      continue;
    }
    const last = blocks[blocks.length - 1];
    if (last && last.end === r - 1) last.end = r;
    else blocks.push({ start: r, end: r });
  }
  return { blocks: blocks, inserted_bases: inserted };
}

// ---- the §3.4 table ------------------------------------------------------------------------------------

test('rs871475760 1:11109 C/A: snv, no shift or tract, zone [11109,11109], both discriminating {11109,C,A,11109}', function () {
  const e = only([record('rs871475760')], G1);
  e.key.should.equal('1:11109:C:A');
  e.kind.should.equal('snv');
  e.shift.should.equal(0);
  should(normalize.tractOf(e.vcf, e.shift)).be.null();
  e.zone.should.eql({ start: 11109, end: 11109 });
  e.discriminating.should.eql({ forward: site(11109, 'C', 'A', 11109), reverse: site(11109, 'C', 'A', 11109) });
  e.label.should.equal('1:11109 C/A');
  e.ref_verified.should.be.true();
});

test('rs871475760 by id: the §2.9 variant block (without the design-only requested_id and submission_sequence)', function () {
  const lookup = variationFixture('variation_rs871475760.json');
  lookup.synonyms.should.eql(['.', 'tmp_1_11109_C_A']);
  const e = only([record('rs871475760')], G1, { requested_id: 'rs871475760', synonyms: lookup.synonyms });
  e.should.eql({
    key: '1:11109:C:A',
    ids: ['rs871475760'],
    synonyms: ['tmp_1_11109_C_A'],
    label: '1:11109 C/A',
    kind: 'snv',
    region: '1',
    vcf: { position: 11109, ref: 'C', alt: 'A' },
    minimal: { start: 11109, end: 11109, ref: 'C', alt: 'A' },
    alleles: ['C', 'A'],
    multiallelic: null,
    shift: 0,
    zone: { start: 11109, end: 11109 },
    discriminating: {
      forward: { position: 11109, ref_base: 'C', alt_base: 'A', alt_maps_to: 11109 },
      reverse: { position: 11109, ref_base: 'C', alt_base: 'A', alt_maps_to: 11109 }
    },
    records: [{ id: 'rs871475760', source: 'EVA', ems: false }],
    ems: false,
    consequence: 'downstream_gene_variant',
    ref_verified: true,
    designable: true,
    issues: []
  });
});

test('synonyms: Ensembl "." is dropped, only the lookup fills them, and a merged id is never a synonym', function () {
  byKey('1:11109:C:A').synonyms.should.eql([]);
  const del = only([record('rs5413864115')], G1, {
    requested_id: 'rs5413864115', synonyms: variationFixture('variation_rs5413864115.json').synonyms
  });
  del.synonyms.should.eql(['tmp_1_11282_CA_C']);
  const ins = only([record('tmp_1_11502_C_CGT'), record('rs5413863549')], G1,
    { requested_id: 'tmp_1_11502_C_CGT', synonyms: ['.', 'rs5413863549', 'x.y', 'x.y'] });
  ins.synonyms.should.eql(['x.y']);
});

test('tmp_1_11193_C_T 1:11193 C/T: key 1:11193:C:T, EMS from source EMS_PMID38100514_Jiao', function () {
  const e = only([record('tmp_1_11193_C_T')], G1);
  e.key.should.equal('1:11193:C:T');
  e.shift.should.equal(0);
  e.zone.should.eql({ start: 11193, end: 11193 });
  e.discriminating.should.eql({ forward: site(11193, 'C', 'T', 11193), reverse: site(11193, 'C', 'T', 11193) });
  e.records.should.eql([{ id: 'tmp_1_11193_C_T', source: 'EMS_PMID38100514_Jiao', ems: true }]);
  e.ems.should.be.true();
  // the pattern is injected, never read from config
  only([record('tmp_1_11318_C_T')], G1).ems.should.be.false();
  only([record('tmp_1_11318_C_T')], G1, { ems_source_pattern: '^BAP_' }).ems.should.be.true();
});

test('rs5413864115 1:11283 A/- in AAA: key 1:11282:CA:C, shift 2, tract [11283,11285], zone [11282,11286]', function () {
  const e = only([record('rs5413864115')], G1);
  e.key.should.equal('1:11282:CA:C');
  e.kind.should.equal('deletion');
  e.shift.should.equal(2);
  normalize.tractOf(e.vcf, e.shift).should.eql({ start: 11283, end: 11285 });
  e.zone.should.eql({ start: 11282, end: 11286 });
  e.discriminating.forward.should.eql(site(11285, 'A', 'G', 11286));
  e.discriminating.reverse.should.eql(site(11283, 'A', 'C', 11282));
  e.minimal.should.eql({ start: 11283, end: 11283, ref: 'A', alt: '-' });
  e.label.should.equal('1:11283-11283 A/-');
  e.alleles.should.eql(['A', '-']);
});

test('insertion conventions: SAP (start 11503, end 11502) and EVA (start = end = 11503) both give 1:11502:C:CGT', function () {
  const sap = record('tmp_1_11502_C_CGT');
  const eva = record('rs5413863549');
  [sap.start, sap.end, eva.start, eva.end].should.eql([11503, 11502, 11503, 11503]);
  only([sap], G1).key.should.equal('1:11502:C:CGT');
  only([eva], G1).key.should.equal('1:11502:C:CGT');

  const e = only([eva, sap], G1);
  e.ids.should.eql(['rs5413863549', 'tmp_1_11502_C_CGT']);
  only([eva, sap], G1, { requested_id: 'tmp_1_11502_C_CGT' }).ids.should.eql(['tmp_1_11502_C_CGT', 'rs5413863549']);
  e.kind.should.equal('insertion');
  e.shift.should.equal(0);
  should(normalize.tractOf(e.vcf, e.shift)).be.null();
  e.discriminating.forward.should.eql(site(11503, 'A', 'G', null));
  e.discriminating.reverse.should.eql(site(11502, 'C', 'T', null));
  e.zone.should.eql({ start: 11502, end: 11503 });
  e.minimal.should.eql({ start: 11503, end: 11502, ref: '-', alt: 'GT' });
  e.label.should.equal('1:11502^11503 -/GT');
});

test('tmp_1_11502_C_CGT by id: the §2.4 entry, duplicate ids merged with the requested id first', function () {
  const list = entriesOf(RECORDS, G1, {
    requested_id: 'tmp_1_11502_C_CGT', synonyms: variationFixture('variation_tmp_1_11502_C_CGT.json').synonyms
  }).filter(function (e) { return e.ids.indexOf('tmp_1_11502_C_CGT') >= 0; });
  list.should.eql([{
    key: '1:11502:C:CGT', ids: ['tmp_1_11502_C_CGT', 'rs5413863549'], synonyms: [], label: '1:11502^11503 -/GT',
    kind: 'insertion', region: '1',
    vcf: { position: 11502, ref: 'C', alt: 'CGT' }, minimal: { start: 11503, end: 11502, ref: '-', alt: 'GT' },
    alleles: ['-', 'GT'], multiallelic: null, shift: 0, zone: { start: 11502, end: 11503 },
    discriminating: { forward: { position: 11503, ref_base: 'A', alt_base: 'G', alt_maps_to: null },
      reverse: { position: 11502, ref_base: 'C', alt_base: 'T', alt_maps_to: null } },
    records: [{ id: 'tmp_1_11502_C_CGT', source: 'SAP_PMID35653240_Boatwri', ems: false },
      { id: 'rs5413863549', source: 'EVA', ems: false }],
    ems: false, consequence: '3_prime_UTR_variant', ref_verified: true, designable: true, issues: []
  }]);
});

test('rs5413863413 1:13736-13738 TGG/- in five TGG units: key 1:13735:TTGG:T, shift 13, zone [13735,13752]', function () {
  const e = only([record('rs5413863413', RECORDS_13)], G2);
  e.key.should.equal('1:13735:TTGG:T');
  e.kind.should.equal('deletion');
  e.shift.should.equal(13);
  normalize.tractOf(e.vcf, e.shift).should.eql({ start: 13736, end: 13751 });
  e.discriminating.forward.should.eql(site(13749, 'G', 'A', 13752));
  e.discriminating.reverse.should.eql(site(13738, 'G', 'T', 13735));
  e.zone.should.eql({ start: 13735, end: 13752 });
  e.minimal.should.eql({ start: 13736, end: 13738, ref: 'TGG', alt: '-' });
  e.label.should.equal('1:13736-13738 TGG/-');

  // The BAP insertion of one unit at the same repeat is another event with its own key (and the same zone).
  const ins = only([record('tmp_1_13735_T_TTGG', RECORDS_13)], G2);
  ins.key.should.equal('1:13735:T:TTGG');
  ins.shift.should.equal(16);
  normalize.tractOf(ins.vcf, ins.shift).should.eql({ start: 13736, end: 13751 });
  ins.discriminating.should.eql({ forward: site(13752, 'A', 'G', 13749), reverse: site(13735, 'T', 'G', null) });
  ins.zone.should.eql({ start: 13735, end: 13752 });
});

// ---- multi-allelic sites, '*' and merging --------------------------------------------------------------

test('rs5413863494 C/T/G gives two entries with other_alts; rs5413863901 C/T/* gives one entry plus STAR_ALLELE', function () {
  const multi = entriesOf([record('rs5413863494')], G1);
  multi.map(function (e) { return [e.key, e.multiallelic, e.designable]; }).should.eql([
    ['1:10718:C:G', { alleles: ['C', 'T', 'G'], other_alts: ['T'] }, true],
    ['1:10718:C:T', { alleles: ['C', 'T', 'G'], other_alts: ['G'] }, true]
  ]);
  multi[0].alleles.should.eql(['C', 'T', 'G']);

  const res = normalize.recordsToEntries([record('rs5413863901')], G1, { region: '1' });
  res.skipped.should.eql({ count: 0, reasons: {} });
  res.entries.should.have.length(1);
  const star = res.entries[0];
  star.key.should.equal('1:11318:C:T');
  star.multiallelic.should.eql({ alleles: ['C', 'T', '*'], other_alts: ['*'] });
  star.issues.map(function (i) { return [i.code, i.details]; }).should.eql([['STAR_ALLELE', { ids: ['rs5413863901'] }]]);
  star.designable.should.be.true();

  // the BAP record with the same key merges into it, and the site keeps the three-allele list
  const merged = only([record('tmp_1_11318_C_T'), record('rs5413863901')], G1);
  merged.ids.should.eql(['rs5413863901', 'tmp_1_11318_C_T']);
  merged.alleles.should.eql(['C', 'T', '*']);
  merged.issues.map(function (i) { return i.code; }).should.eql(['STAR_ALLELE']);
  only([record('tmp_1_11318_C_T'), record('rs5413863901')], G1, { requested_id: 'tmp_1_11318_C_T' })
    .multiallelic.should.eql({ alleles: ['C', 'T', '*'], other_alts: ['*'] });
});

test('the four overlapping windows: each record counted once, ids never repeated, no record skipped', function () {
  RECORDS.should.have.length(237);
  const res = normalize.recordsToEntries(RECORDS, G1, { region: '1' });
  res.skipped.should.eql({ count: 0, reasons: {} });
  res.entries.should.have.length(89);
  const kinds = {};
  res.entries.forEach(function (e) {
    kinds[e.kind] = (kinds[e.kind] || 0) + 1;
    e.ids.should.eql(Array.from(new Set(e.ids)));
  });
  kinds.should.eql({ snv: 79, insertion: 4, deletion: 6 });
  byKey('1:11109:C:A').ids.should.eql(['rs871475760']);
  byKey('1:11553:A:AC').ids.should.eql(['rs5413863374', 'tmp_1_11553_A_AC']);
});

test('ids: the requested id first, then rs ids by number, then the rest by code unit', function () {
  const r = record('rs871475760');
  const recs = ['tmp_b', 'rs10', 'rs9', 'Abc'].map(function (id) { return Object.assign({}, r, { id: id }); });
  only(recs, G1).ids.should.eql(['rs9', 'rs10', 'Abc', 'tmp_b']);
  only(recs, G1, { requested_id: 'tmp_b' }).ids.should.eql(['tmp_b', 'rs9', 'rs10', 'Abc']);
});

// ---- the §2.3 listing --------------------------------------------------------------------------------

test('window 1:11180-11290: the four entries of §2.3, deep-equal and in order', function () {
  const listed = all().filter(function (e) { return normalize.inWindow(e.minimal, 11180, 11290); });
  listed.should.eql([
    { key: '1:11182:A:G', ids: ['rs873026643'], synonyms: [], label: '1:11182 A/G', kind: 'snv', region: '1',
      vcf: { position: 11182, ref: 'A', alt: 'G' }, minimal: { start: 11182, end: 11182, ref: 'A', alt: 'G' },
      alleles: ['A', 'G'], multiallelic: null, shift: 0, zone: { start: 11182, end: 11182 },
      discriminating: { forward: { position: 11182, ref_base: 'A', alt_base: 'G', alt_maps_to: 11182 },
        reverse: { position: 11182, ref_base: 'A', alt_base: 'G', alt_maps_to: 11182 } },
      records: [{ id: 'rs873026643', source: 'EVA', ems: false }], ems: false,
      consequence: '3_prime_UTR_variant', ref_verified: true, designable: true, issues: [] },
    { key: '1:11193:C:T', ids: ['tmp_1_11193_C_T'], synonyms: [], label: '1:11193 C/T', kind: 'snv', region: '1',
      vcf: { position: 11193, ref: 'C', alt: 'T' }, minimal: { start: 11193, end: 11193, ref: 'C', alt: 'T' },
      alleles: ['C', 'T'], multiallelic: null, shift: 0, zone: { start: 11193, end: 11193 },
      discriminating: { forward: { position: 11193, ref_base: 'C', alt_base: 'T', alt_maps_to: 11193 },
        reverse: { position: 11193, ref_base: 'C', alt_base: 'T', alt_maps_to: 11193 } },
      records: [{ id: 'tmp_1_11193_C_T', source: 'EMS_PMID38100514_Jiao', ems: true }], ems: true,
      consequence: '3_prime_UTR_variant', ref_verified: true, designable: true, issues: [] },
    { key: '1:11203:C:T', ids: ['tmp_1_11203_C_T'], synonyms: [], label: '1:11203 C/T', kind: 'snv', region: '1',
      vcf: { position: 11203, ref: 'C', alt: 'T' }, minimal: { start: 11203, end: 11203, ref: 'C', alt: 'T' },
      alleles: ['C', 'T'], multiallelic: null, shift: 0, zone: { start: 11203, end: 11203 },
      discriminating: { forward: { position: 11203, ref_base: 'C', alt_base: 'T', alt_maps_to: 11203 },
        reverse: { position: 11203, ref_base: 'C', alt_base: 'T', alt_maps_to: 11203 } },
      records: [{ id: 'tmp_1_11203_C_T', source: 'EMS_PMID29378822_Addo-Qu', ems: true }], ems: true,
      consequence: '3_prime_UTR_variant', ref_verified: true, designable: true, issues: [] },
    { key: '1:11282:CA:C', ids: ['rs5413864115'], synonyms: [], label: '1:11283-11283 A/-', kind: 'deletion', region: '1',
      vcf: { position: 11282, ref: 'CA', alt: 'C' }, minimal: { start: 11283, end: 11283, ref: 'A', alt: '-' },
      alleles: ['A', '-'], multiallelic: null, shift: 2, zone: { start: 11282, end: 11286 },
      discriminating: { forward: { position: 11285, ref_base: 'A', alt_base: 'G', alt_maps_to: 11286 },
        reverse: { position: 11283, ref_base: 'A', alt_base: 'C', alt_maps_to: 11282 } },
      records: [{ id: 'rs5413864115', source: 'EVA', ems: false }], ems: false,
      consequence: '3_prime_UTR_variant', ref_verified: true, designable: true, issues: [] }
  ]);
});

test('inWindow: minimal spans overlap the window; an insertion is in when either flanking coordinate is', function () {
  const ins = { start: 11503, end: 11502, ref: '-', alt: 'GT' };
  normalize.inWindow(ins, 11400, 11502).should.be.true();
  normalize.inWindow(ins, 11503, 11600).should.be.true();
  normalize.inWindow(ins, 11504, 11600).should.be.false();
  normalize.inWindow(ins, 11400, 11501).should.be.false();
  const del = { start: 11283, end: 11283, ref: 'A', alt: '-' };
  normalize.inWindow(del, 11283, 11283).should.be.true();
  normalize.inWindow(del, 11284, 11290).should.be.false();
  normalize.inWindow({ start: 13736, end: 13738, ref: 'TGG', alt: '-' }, 13738, 13800).should.be.true();
});

// ---- manual variants (§3.9) ---------------------------------------------------------------------------

test('parseManual: VCF and Ensembl styles become pseudo-records (insertion: end = start - 1)', function () {
  normalize.parseManual({ region: '1', position: 11282, ref: 'ca', alt: 'C' }).should.eql(
    { id: null, source: null, seq_region_name: '1', start: 11282, end: 11283, alleles: ['CA', 'C'] });
  normalize.parseManual({ region: '1', position: 11283, ref: 'A', alt: '-' }).should.eql(
    { id: null, source: null, seq_region_name: '1', start: 11283, end: 11283, alleles: ['A', '-'] });
  normalize.parseManual({ region: '1', position: 11503, ref: '-', alt: 'gt' }).should.eql(
    { id: null, source: null, seq_region_name: '1', start: 11503, end: 11502, alleles: ['-', 'GT'] });
});

test('parseManual: §2.7 rule 3 and the allele limits', function () {
  const base = { region: '1', position: 11109 };
  function manual(extra) {
    return function () { normalize.parseManual(Object.assign({}, base, extra)); };
  }
  throwsWith(manual({ ref: '-', alt: '-' }), 400, 'INVALID_VARIANT', { reason: 'alleles' });
  throwsWith(manual({ ref: 'c', alt: 'C' }), 400, 'INVALID_VARIANT', { reason: 'alleles' });
  throwsWith(manual({ ref: '-', alt: 'XY' }), 400, 'INVALID_VARIANT', { reason: 'alleles' });
  throwsWith(manual({ ref: 'C', alt: 'N' }), 400, 'UNSUPPORTED_ALLELE', { allele: 'N' });
  throwsWith(manual({ ref: 'C', alt: '*' }), 400, 'UNSUPPORTED_ALLELE', { allele: '*' });
  throwsWith(manual({ ref: 'C'.repeat(51), alt: 'C' }), 400, 'INVALID_VARIANT', { reason: 'allele_too_long' });
  normalize.parseManual(Object.assign({}, base, { ref: 'C'.repeat(10), alt: 'C' }), { max_allele_length: 10 }).end.should.equal(11118);
  throwsWith(function () {
    normalize.parseManual(Object.assign({}, base, { ref: 'C'.repeat(11), alt: 'C' }), { max_allele_length: 10 });
  }, 400, 'INVALID_VARIANT', { reason: 'allele_too_long' });
  throwsWith(manual({ ref: 'C', alt: 'A', position: 0 }), 400, 'INVALID_VARIANT', { reason: 'id_or_manual' });
  throwsWith(function () { normalize.parseManual({ position: 11109, ref: 'C', alt: 'A' }); }, 400, 'INVALID_VARIANT', { reason: 'id_or_manual' });
});

test('manual {11283,"A","-"} and {11282,"CA","C"} give the same entry; with the overlap it is rs5413864115', function () {
  const ensemblStyle = normalize.parseManual({ region: '1', position: 11283, ref: 'A', alt: '-' });
  const vcfStyle = normalize.parseManual({ region: '1', position: 11282, ref: 'CA', alt: 'C' });
  const a = only([ensemblStyle], G1);
  a.should.eql(only([vcfStyle], G1));
  a.key.should.equal('1:11282:CA:C');
  a.ids.should.eql([]);
  a.records.should.eql([]);
  a.ems.should.be.false();
  a.alleles.should.eql(['A', '-']);
  a.shift.should.equal(2);
  a.zone.should.eql({ start: 11282, end: 11286 });

  // §3.9: on a genome with variation data the template-window overlap fills ids and records for the same key
  const listing = byKey('1:11282:CA:C');
  entriesOf(RECORDS.concat([ensemblStyle]), G1).find(function (e) { return e.key === listing.key; }).should.eql(listing);
  entriesOf([vcfStyle].concat(RECORDS), G1).find(function (e) { return e.key === listing.key; }).should.eql(listing);

  const insEnsembl = normalize.parseManual({ region: '1', position: 11503, ref: '-', alt: 'GT' });
  const insVcf = normalize.parseManual({ region: '1', position: 11502, ref: 'C', alt: 'CGT' });
  only([insEnsembl], G1).should.eql(only([insVcf], G1));
  only([insEnsembl], G1).key.should.equal('1:11502:C:CGT');

  const snv = normalize.parseManual({ region: '1', position: 11109, ref: 'C', alt: 'A' });
  entriesOf(RECORDS.concat([snv]), G1).find(function (e) { return e.key === '1:11109:C:A'; }).ids.should.eql(['rs871475760']);
});

test('manual ref A at 1:11109: REF_MISMATCH {given: "A", genome: "C"}, not verified, not designable', function () {
  const e = only([normalize.parseManual({ region: '1', position: 11109, ref: 'A', alt: 'T' })], G1);
  e.ref_verified.should.be.false();
  e.designable.should.be.false();
  e.issues.should.have.length(1);
  e.issues[0].code.should.equal('REF_MISMATCH');
  e.issues[0].details.should.eql({ region: '1', position: 11109, given: 'A', genome: 'C' });
  // the given allele still yields geometry, so a listing never breaks on a bad record
  e.zone.should.eql({ start: 11109, end: 11109 });
});

// ---- left normalization, shift and the helpers --------------------------------------------------------

test('leftNormalize: right-shifted and VCF-style descriptions of one event give one key', function () {
  normalize.leftNormalize({ position: 11284, ref: 'AA', alt: 'A' }, G1).should.eql({ position: 11282, ref: 'CA', alt: 'C' });
  normalize.leftNormalize({ position: 11285, ref: 'A', alt: '' }, G1).should.eql({ position: 11282, ref: 'CA', alt: 'C' });
  normalize.leftNormalize({ position: 11281, ref: 'tca', alt: 'tc' }, G1).should.eql({ position: 11282, ref: 'CA', alt: 'C' });
  normalize.leftNormalize({ position: 13748, ref: 'TGGT', alt: 'T' }, G2).should.eql({ position: 13735, ref: 'TTGG', alt: 'T' });
  normalize.leftNormalize({ position: 11502, ref: 'C', alt: 'CGT' }, G1).should.eql({ position: 11502, ref: 'C', alt: 'CGT' });
  normalize.leftNormalize({ position: 11109, ref: 'C', alt: 'A' }, G1).should.eql({ position: 11109, ref: 'C', alt: 'A' });
  (function () { normalize.leftNormalize({ position: 11109, ref: 'C', alt: 'c' }, G1); }).should.throw(/differ/);
  only([{ id: 'x', seq_region_name: '1', start: 11285, end: 11285, alleles: ['A', '-'] }], G1).key.should.equal('1:11282:CA:C');
});

test('kindOf, minimalOf and labelOf for snv, mnv, insertion, deletion and complex', function () {
  normalize.kindOf('C', 'A').should.equal('snv');
  normalize.kindOf('AC', 'GT').should.equal('mnv');
  normalize.kindOf('C', 'CGT').should.equal('insertion');
  normalize.kindOf('TTGG', 'T').should.equal('deletion');
  normalize.kindOf('CC', 'T').should.equal('complex');
  normalize.minimalOf({ position: 13735, ref: 'TTGG', alt: 'T' }).should.eql({ start: 13736, end: 13738, ref: 'TGG', alt: '-' });
  normalize.minimalOf({ position: 8, ref: 'AC', alt: 'GT' }).should.eql({ start: 8, end: 9, ref: 'AC', alt: 'GT' });
  normalize.labelOf('c', { start: 8, end: 9, ref: 'AC', alt: 'GT' }).should.equal('c:8-9 AC/GT');
  normalize.labelOf('c', { start: 9, end: 10, ref: 'CC', alt: 'T' }).should.equal('c:9-10 CC/T');

  // synthetic region c: CATTTTGACCGTA
  const c = normalize.sequenceWindow('CATTTTGACCGTA', 1, 13);
  const mnv = only([{ id: 'mnv', start: 8, end: 9, alleles: ['AC', 'GT'] }], c, { region: 'c' });
  [mnv.key, mnv.kind, mnv.label].should.eql(['c:8:AC:GT', 'mnv', 'c:8-9 AC/GT']);
  mnv.discriminating.should.eql({ forward: site(8, 'A', 'G', 8), reverse: site(9, 'C', 'T', 9) });
  mnv.zone.should.eql({ start: 8, end: 9 });
  // A complex event has no anchor base, so its ALT bases map to no reference coordinate (§3.5 altToRef).
  const cpx = only([{ id: 'cpx', start: 9, end: 10, alleles: ['CC', 'T'] }], c, { region: 'c' });
  [cpx.key, cpx.kind, cpx.shift].should.eql(['c:9:CC:T', 'complex', 0]);
  cpx.discriminating.should.eql({ forward: site(9, 'C', 'T', null), reverse: site(10, 'C', 'T', null) });
  cpx.zone.should.eql({ start: 9, end: 10 });
});

test('altToRef reproduces the §4.9 ALT primer blocks', function () {
  // rs5413864115 forward ALT: same template span 375-403 on alt_seq (template 1:10883-11685)
  const del = { position: 11282, ref: 'CA', alt: 'C' };
  altBlocks(del, 11257, 11285).should.eql({ blocks: [{ start: 11257, end: 11282 }, { start: 11284, end: 11286 }], inserted_bases: 0 });
  // reverse ALT: the REF span 1:11283-11307 shifted by delta -1
  altBlocks(del, 11282, 11306).should.eql({ blocks: [{ start: 11282, end: 11282 }, { start: 11284, end: 11307 }], inserted_bases: 0 });
  // tmp_1_11502_C_CGT reverse ALT: the REF span 1:11502-11529 shifted by delta +2
  const ins = { position: 11502, ref: 'C', alt: 'CGT' };
  altBlocks(ins, 11504, 11531).should.eql({ blocks: [{ start: 11503, end: 11529 }], inserted_bases: 1 });
  [11501, 11502, 11503, 11504, 11505].map(function (q) { return normalize.altToRef(ins, q); }).should.eql([11501, 11502, null, null, 11503]);
  altBlocks({ position: 11109, ref: 'C', alt: 'A' }, 11109, 11132).should.eql({ blocks: [{ start: 11109, end: 11132 }], inserted_bases: 0 });
});

test('haplotypes: the §5.6 cores keep the anchor base, unlike the shift tract', function () {
  function core(genome, v, shift, K) {
    const start = v.position - 1;
    const end = v.position + v.ref.length - 1 + shift + 1;
    const c = normalize.haplotypes(v, genome, start, end);
    const h = normalize.haplotypes(v, genome, start - K, end + K);
    return [start, end, c.ref, c.alt, h.start, h.end, h.ref.length, h.alt.length];
  }
  const snv = { position: 11109, ref: 'C', alt: 'A' };
  core(G1, snv, 0, 15).should.eql([11108, 11110, 'TCT', 'TAT', 11093, 11125, 33, 33]);
  normalize.haplotypes(snv, G1, 11093, 11125).should.eql({
    start: 11093, end: 11125, ref: 'ATAGTCATACTCTATTCTGAATTTCTCGCTAGT', alt: 'ATAGTCATACTCTATTATGAATTTCTCGCTAGT'
  });
  core(G1, { position: 11193, ref: 'C', alt: 'T' }, 0, 15).should.eql([11192, 11194, 'TCT', 'TTT', 11177, 11209, 33, 33]);
  const del = { position: 11282, ref: 'CA', alt: 'C' };
  core(G1, del, 2, 15).should.eql([11281, 11286, 'TCAAAG', 'TCAAG', 11266, 11301, 36, 35]);
  const tract = normalize.tractOf(del, 2);
  [tract.start - 2, tract.end + 1].should.eql([11281, 11286]);
  core(G1, { position: 11502, ref: 'C', alt: 'CGT' }, 0, 15).should.eql([11501, 11503, 'CCA', 'CCGTA', 11486, 11518, 33, 35]);
  core(G2, { position: 13735, ref: 'TTGG', alt: 'T' }, 13, 22)
    .should.eql([13734, 13752, 'ATTGGTGGTGGTGGTGGTA', 'ATTGGTGGTGGTGGTA', 13712, 13774, 63, 60]);
  (function () { normalize.haplotypes(del, G1, 11283, 11300); }).should.throw(/cover/);
});

test('shift above max_shift: REPEAT_TOO_LONG, no zone; a repeat running to the region end has no discriminating base', function () {
  // region c: G + 30 A + CGT
  const run = normalize.sequenceWindow('G' + 'A'.repeat(30) + 'CGT', 1, 34);
  const v = { position: 1, ref: 'GA', alt: 'G' };
  normalize.shiftOf(v, run).should.equal(29);
  normalize.shiftOf(v, run, { max_shift: 5 }).should.equal(6);
  normalize.shiftOf({ position: 1, ref: 'G', alt: 'A' }, run).should.equal(0);
  const ok = only([{ id: 'r', start: 31, end: 31, alleles: ['A', '-'] }], run, { region: 'c' });
  [ok.key, ok.shift, ok.zone].should.eql(['c:1:GA:G', 29, { start: 1, end: 32 }]);
  ok.discriminating.should.eql({ forward: site(31, 'A', 'C', 32), reverse: site(2, 'A', 'G', 1) });

  const capped = only([{ id: 'r', start: 31, end: 31, alleles: ['A', '-'] }], run, { region: 'c', max_shift: 5 });
  [capped.shift, capped.zone, capped.discriminating, capped.designable].should.eql([6, null, null, false]);
  capped.issues.map(function (i) { return [i.code, i.details]; }).should.eql([['REPEAT_TOO_LONG', { region: 'c', position: 1, max: 5 }]]);

  const edge = normalize.sequenceWindow('GC' + 'A'.repeat(10), 1, 12);
  const err = (function () {
    try {
      normalize.discriminating({ position: 2, ref: 'CA', alt: 'C' }, edge);
    } catch (e) {
      return e;
    }
    return null;
  })();
  should.exist(err);
  err.code.should.equal('NO_DISCRIMINATING_BASE');
  const atEnd = only([{ id: 'e', start: 5, end: 5, alleles: ['A', '-'] }], edge, { region: 'c' });
  [atEnd.key, atEnd.shift, atEnd.zone, atEnd.designable].should.eql(['c:2:CA:C', 9, null, false]);
  atEnd.issues[0].code.should.equal('REPEAT_TOO_LONG');
});

test('allele issues, skipped records and the window accessor', function () {
  const c = normalize.sequenceWindow('CATTTTGACCGTA', 1, 13);
  const res = normalize.recordsToEntries([
    { id: 'edge', start: 1, end: 0, alleles: ['-', 'A'] },
    { id: 'same', start: 2, end: 2, alleles: ['A', 'a'] },
    { id: 'star', start: 2, end: 2, alleles: ['*', 'A'] },
    { id: 'one', start: 2, end: 2, alleles: ['A'] },
    { id: 'coord', start: '2', end: 2, alleles: ['A', 'C'] },
    { id: 'dup', start: 2, end: 2, alleles: ['A', 'G'] },
    { id: 'dup', start: 2, end: 2, alleles: ['A', 'G'] },
    { id: 'nn', start: 11, end: 11, alleles: ['G', 'N'] },
    { id: 'long', start: 3, end: 3, alleles: ['T', 'T'.repeat(51)] }
  ], c, { region: 'c' });
  res.skipped.should.eql({ count: 5, reasons: { region_edge: 1, identical_alleles: 1, star_reference: 1, invalid_alleles: 1, invalid_coordinates: 1 } });
  res.entries.map(function (e) {
    return [e.key, e.ids, e.designable, e.issues.map(function (i) { return [i.code, i.details]; })];
  }).should.eql([
    ['c:2:A:ATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT', ['long'], false, [['ALLELE_TOO_LONG', { length: 51, max: 50 }]]],
    ['c:2:A:G', ['dup'], true, []],
    ['c:11:G:N', ['nn'], false, [['UNSUPPORTED_ALLELE', { allele: 'N' }]]]
  ]);
  normalize.recordsToEntries([{ id: 'x', start: 2, end: 2, alleles: ['A', 'C'] }], c).skipped
    .should.eql({ count: 1, reasons: { invalid_region: 1 } });

  const win = (function () {
    try {
      G1.bases(10499, 10500);
    } catch (e) {
      return e;
    }
    return null;
  })();
  should.exist(win);
  win.should.be.instanceOf(normalize.SequenceWindowError);
  [win.code, win.needed].should.eql(['SEQUENCE_WINDOW', { start: 10499, end: 10500 }]);
  (function () { c.bases(0, 1); }).should.throw(/outside the region/);
  G1.bases(11282, 11286).should.equal('CAAAG');
  G1.bases(11500, 11499).should.equal('');
});

test('§3.5 geometry holds for every recorded entry: zone = tract ± 1 for each shiftable indel', function () {
  const entries = all().concat(entriesOf(RECORDS_13, G2));
  let shiftable = 0;
  entries.forEach(function (e) {
    const v = e.vcf;
    const d = e.discriminating;
    const z = e.zone;
    const vEnd = v.position + v.ref.length - 1;
    (z.start <= v.position && z.end >= vEnd).should.be.true(e.key);
    [d.forward.position, d.reverse.position, d.forward.alt_maps_to, d.reverse.alt_maps_to].forEach(function (p) {
      if (p !== null) (p >= z.start && p <= z.end).should.be.true(e.key);
    });
    const tract = normalize.tractOf(v, e.shift);
    if (!tract) return;
    shiftable += 1;
    z.should.eql({ start: tract.start - 1, end: tract.end + 1 });
    if (e.kind === 'deletion') {
      [d.forward.position, d.reverse.position].should.eql([tract.end - v.ref.length + 2, vEnd]);
      [d.forward.alt_maps_to, d.reverse.alt_maps_to].should.eql([tract.end + 1, v.position]);
    } else {
      [d.forward.position, d.reverse.position].should.eql([tract.end + 1, v.position]);
    }
  });
  shiftable.should.equal(8);
});

// ---- submission strings (§4.17) -----------------------------------------------------------------------

test('submissionSequence: the four §4.17 strings, with IUPAC-coded neighbours and nothing omitted', function () {
  const expected = {
    '1:11109:C:A': 'YCCTCAAAAAGCTTCTCTAAGTGGTTATCCGAATATAGTCATACTCTATT[C/A]TGAATTTCTCGCTAGTCAAAGATAACAAAAATAGCATATTCTGGATTTCT',
    '1:11193:C:T': 'GCATATTCTGGATTTCTCWGTAGACAAAGATAGATAACARAAATAGCTCT[C/T]TAGAGTATACACAATATATTAAAAAGTTGTTAGAGAGTGAAAATATATAG',
    '1:11282:CA:C': 'AAAATATATAGAAAACAATTTTATACAGATGATTTTCCAAATGATGATTC[A/]AAGTGTGAAATTTGRAAAGWCTCTTRGASATGMTYTAAGTGGAAGGAACA',
    '1:11502:C:CGT': 'ATGTTAGGATCTTTGCAACCCWGTGTTGCGTGCAATCTCGGTATCTCRCC[/GT]ATATGAYGTTAGGWTTTCTTWTCCTGCAACTGCCAAGAGAAYAAATATAT'
  };
  Object.keys(expected).forEach(function (key) {
    normalize.submissionSequence(byKey(key), G1, all()).should.eql({ sequence: expected[key], omitted_ids: [] });
  });

  // rs873774986 C/T at 1:11059 is the leading Y; EMS tmp_1_11069_G_A stays the reference G
  const s871 = normalize.submissionSequence(byKey('1:11109:C:A'), G1, all()).sequence;
  [s871[0], s871[10]].should.eql(['Y', 'G']);
  byKey('1:11069:G:A').ems.should.be.true();
  // the right flank starts after ']' at minimal.end + 1
  function rightFlankBase(key, pos) {
    const s = normalize.submissionSequence(byKey(key), G1, all()).sequence;
    return s[s.indexOf(']') + 1 + pos - (byKey(key).minimal.end + 1)];
  }
  // rs5413863901 C/T/* at 1:11318 counts as biallelic: Y, and no omission
  rightFlankBase('1:11282:CA:C', 11318).should.equal('Y');
  // EMS tmp_1_11203_C_T stays the reference C in the tmp_1_11193_C_T string
  rightFlankBase('1:11193:C:T', 11203).should.equal('C');
});

test('submissionSequence: indel and multi-allelic neighbours stay reference and are listed for SUBMISSION_NEIGHBOURS_OMITTED', function () {
  // rs5413863234 G/A at 1:11500: the tmp_1_11502_C_CGT insertion sits in its right flank
  normalize.submissionSequence(byKey('1:11500:G:A'), G1, all()).should.eql({
    sequence: 'ATGATGTTAGGATCTTTGCAACCCWGTGTTGCGTGCAATCTCGGTATCTC[G/A]CCATATGAYGTTAGGWTTTCTTWTCCTGCAACTGCCAAGAGAAYAAATAT',
    omitted_ids: ['rs5413863549', 'tmp_1_11502_C_CGT']
  });
  // rs5413863688 C/T at 1:10710: rs5413863494 C/T/G at 10718 stays C
  normalize.submissionSequence(byKey('1:10710:C:T'), G1, all()).should.eql({
    sequence: 'TATTTTCTTTTTACTTACGAGTTCTTAGGAATTCATCGAGTCAACCTCGG[C/T]GAGTTCTCAAGTTCTRCRTGAATACCTCTTGAYCGTGGCATCCRGGCGTA',
    omitted_ids: ['rs5413863494']
  });
  // two biallelic records with different alts at one position make a multi-allelic site
  const second = Object.assign({}, record('rs873774986'), { id: 'other_alt', alleles: ['C', 'G'] });
  const neighbours = entriesOf(RECORDS.concat([second]), G1);
  const r = normalize.submissionSequence(byKey('1:11109:C:A'), G1, neighbours);
  r.sequence[0].should.equal('C');
  // listed by position, then key: 1:11059:C:G before 1:11059:C:T
  r.omitted_ids.should.eql(['other_alt', 'rs873774986']);
  // EMS neighbours never count, and the target itself is never a neighbour
  normalize.submissionSequence(byKey('1:11109:C:A'), G1, [byKey('1:11109:C:A'), byKey('1:11069:G:A')]).sequence
    .should.equal('CCCTCAAAAAGCTTCTCTAAGTGGTTATCCGAATATAGTCATACTCTATT[C/A]TGAATTTCTCGCTAGTCAAAGATAACAAAAATAGCATATTCTGGATTTCT');
});

test('IUPAC: the six two-base codes', function () {
  normalize.IUPAC.should.eql({ AC: 'M', AG: 'R', AT: 'W', CG: 'S', CT: 'Y', GT: 'K' });
  Object.isFrozen(normalize.IUPAC).should.be.true();
});
