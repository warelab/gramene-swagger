'use strict';

// variation/index.js (genotyping spec §2.3, §2.4, §2.7 rules 7-14, §3.2, §3.6-§3.9, §4.15), offline: the real client
// over a fake fetch that serves the recorded Ensembl bodies (M3's fixtures; a lookup without a recorded body is built
// from the recorded overlap record of that id), a fake sequence module serving the recorded FASTA windows (N elsewhere)
// that counts its reads, and a stub assembly resolver.

require('../../../api/helpers/primers/node_compat');

const test = require('node:test');
const should = require('should');
const fs = require('fs');
const path = require('path');

const config = require('../../../api/helpers/primers/config');
const { PrimerHttpError } = require('../../../api/helpers/primers/errors');
const normalize = require('../../../api/helpers/primers/variation/normalize');
const { createVariationClient } = require('../../../api/helpers/primers/variation/client');
const variation = require('../../../api/helpers/primers/variation');

const FIX = path.join(__dirname, '..', 'fixtures');
const CHR1_LENGTH = 80884392;
const REGION2_LENGTH = 5000;
const SILENT = { log: function () {}, info: function () {}, warn: function () {}, error: function () {} };

function readText(rel) {
  return fs.readFileSync(path.join(FIX, rel), 'utf8');
}

const PIECES = [
  { start: 10500, seq: readText('design/sorghum_bicolor_1_10500-12100.plus.txt').trim() },
  { start: 13400, seq: readText('design/sorghum_bicolor_1_13400-14000.plus.txt').trim() }
];

function baseAt(p) {
  for (const piece of PIECES) {
    if (p >= piece.start && p < piece.start + piece.seq.length) return piece.seq[p - piece.start];
  }
  return 'N';
}

// Every recorded overlap record once (the recorded windows overlap).
const RECORDS = (function () {
  const seen = new Set();
  const out = [];
  ['10709-11509', '10793-11593', '10883-11685', '11102-11903', '13600-13820'].forEach(function (w) {
    JSON.parse(readText('variation/overlap_1_' + w + '.json')).forEach(function (r) {
      const k = JSON.stringify([r.id, r.start, r.end, r.alleles]);
      if (!seen.has(k)) {
        seen.add(k);
        out.push(r);
      }
    });
  });
  return out;
})();

function respond(body, status) {
  return new Response(typeof body === 'string' ? body : JSON.stringify(body), { status: status || 200 });
}

function refused() {
  return Object.assign(new TypeError('fetch failed'), { cause: Object.assign(new Error('connect ECONNREFUSED'), { code: 'ECONNREFUSED' }) });
}

// The recorded Ensembl service. opts: {lookups {id: body}, extraRecords [overlap records]}
function recordedEnsembl(opts) {
  opts = opts || {};
  const records = RECORDS.concat(opts.extraRecords || []);
  return function (url) {
    const u = new URL(url);
    let m = /^\/pansite-ensembl-115\/overlap\/region\/sorghum_bicolor\/([^/]+):(\d+)-(\d+)$/.exec(u.pathname);
    if (m) {
      const region = decodeURIComponent(m[1]);
      const s = Number(m[2]);
      const e = Number(m[3]);
      return respond(records.filter(function (r) {
        return String(r.seq_region_name) === region && (!Number.isSafeInteger(r.start) || (Math.min(r.start, r.end) <= e && Math.max(r.start, r.end) >= s));
      }));
    }
    m = /^\/pansite-ensembl-115\/variation\/sorghum_bicolor\/([^/]+)$/.exec(u.pathname);
    if (m) {
      const id = decodeURIComponent(m[1]);
      if (opts.lookups && opts.lookups[id]) return respond(opts.lookups[id]);
      const file = path.join(FIX, 'variation', 'variation_' + id + '.json');
      if (fs.existsSync(file)) return respond(fs.readFileSync(file, 'utf8'));
      const r = RECORDS.find(function (x) { return x.id === id; });
      if (r) {
        return respond({ name: id, synonyms: [], most_severe_consequence: r.consequence_type,
          mappings: [{ seq_region_name: r.seq_region_name, start: r.start, end: r.end, allele_string: r.alleles.join('/') }] });
      }
      return respond(readText('variation/variation_not_found.json'), 400);
    }
    return respond('<html><body>404 Not Found</body></html>', 404);
  };
}

// opts: {overrides (primers config), fetch (handler), lookups, extraRecords}
function world(opts) {
  opts = opts || {};
  const cfg = config._build({ env: {}, fileConfig: {}, overrides: opts.overrides || {} }).config;
  const calls = [];
  const reads = [];
  const resolves = [];
  const handler = opts.fetch || recordedEnsembl(opts);
  const client = createVariationClient({ cfg: cfg, log: SILENT, fetch: async function (url) { calls.push(url); return handler(url); } });
  const genomes = {
    sorghum_bicolor: { system_name: 'sorghum_bicolor', fasta: { dna: '/fake/sorghum_bicolor.fa.gz' } },
    sorghum_rio: { system_name: 'sorghum_rio', fasta: { dna: '/fake/sorghum_rio.fa.gz' } },
    sorghum_nofasta: { system_name: 'sorghum_nofasta', fasta: { dna: null } }
  };
  const deps = {
    cfg: cfg,
    client: client,
    log: SILENT,
    resolve: async function (name) {
      resolves.push(name);
      if (!genomes[name]) throw new PrimerHttpError(404, 'UNKNOWN_GENOME', 'unknown genome: ' + name, { system_name: name });
      return genomes[name];
    },
    sequence: {
      regionLength: async function (fasta, region) { return region === '1' ? CHR1_LENGTH : region === '2' ? REGION2_LENGTH : undefined; },
      fetch: async function (fasta, region, start, end) {
        reads.push([region, start, end]);
        let s = '';
        for (let p = start; p <= end; p++) s += region === '1' ? baseAt(p) : 'A';
        return s;
      }
    }
  };
  return { cfg: cfg, deps: deps, calls: calls, reads: reads, resolves: resolves };
}

async function rejection(promise) {
  try {
    await promise;
  } catch (e) {
    return e;
  }
  throw new Error('expected a rejection');
}

async function rejectsWith(promise, status, code, details) {
  const e = await rejection(promise);
  e.should.match({ status: status, code: code });
  if (details !== undefined) e.details.should.eql(details);
  return e;
}

function site(position, refBase, altBase, altMapsTo) {
  return { position: position, ref_base: refBase, alt_base: altBase, alt_maps_to: altMapsTo };
}

// ---- the documented bodies ------------------------------------------------------------------------------------

const LIST_EXAMPLE = {
  system_name: 'sorghum_bicolor', region: '1', start: 11180, end: 11290,
  source: { name: 'ensembl', release: '115' },
  total: 4, returned: 4, truncated: false,
  variants: [
    { key: '1:11182:A:G', ids: ['rs873026643'], synonyms: [], label: '1:11182 A/G', kind: 'snv', region: '1',
      vcf: { position: 11182, ref: 'A', alt: 'G' }, minimal: { start: 11182, end: 11182, ref: 'A', alt: 'G' },
      alleles: ['A', 'G'], multiallelic: null, shift: 0, zone: { start: 11182, end: 11182 },
      discriminating: { forward: site(11182, 'A', 'G', 11182), reverse: site(11182, 'A', 'G', 11182) },
      records: [{ id: 'rs873026643', source: 'EVA', ems: false }], ems: false,
      consequence: '3_prime_UTR_variant', ref_verified: true, designable: true, issues: [] },
    { key: '1:11193:C:T', ids: ['tmp_1_11193_C_T'], synonyms: [], label: '1:11193 C/T', kind: 'snv', region: '1',
      vcf: { position: 11193, ref: 'C', alt: 'T' }, minimal: { start: 11193, end: 11193, ref: 'C', alt: 'T' },
      alleles: ['C', 'T'], multiallelic: null, shift: 0, zone: { start: 11193, end: 11193 },
      discriminating: { forward: site(11193, 'C', 'T', 11193), reverse: site(11193, 'C', 'T', 11193) },
      records: [{ id: 'tmp_1_11193_C_T', source: 'EMS_PMID38100514_Jiao', ems: true }], ems: true,
      consequence: '3_prime_UTR_variant', ref_verified: true, designable: true, issues: [] },
    { key: '1:11203:C:T', ids: ['tmp_1_11203_C_T'], synonyms: [], label: '1:11203 C/T', kind: 'snv', region: '1',
      vcf: { position: 11203, ref: 'C', alt: 'T' }, minimal: { start: 11203, end: 11203, ref: 'C', alt: 'T' },
      alleles: ['C', 'T'], multiallelic: null, shift: 0, zone: { start: 11203, end: 11203 },
      discriminating: { forward: site(11203, 'C', 'T', 11203), reverse: site(11203, 'C', 'T', 11203) },
      records: [{ id: 'tmp_1_11203_C_T', source: 'EMS_PMID29378822_Addo-Qu', ems: true }], ems: true,
      consequence: '3_prime_UTR_variant', ref_verified: true, designable: true, issues: [] },
    { key: '1:11282:CA:C', ids: ['rs5413864115'], synonyms: [], label: '1:11283-11283 A/-', kind: 'deletion', region: '1',
      vcf: { position: 11282, ref: 'CA', alt: 'C' }, minimal: { start: 11283, end: 11283, ref: 'A', alt: '-' },
      alleles: ['A', '-'], multiallelic: null, shift: 2, zone: { start: 11282, end: 11286 },
      discriminating: { forward: site(11285, 'A', 'G', 11286), reverse: site(11283, 'A', 'C', 11282) },
      records: [{ id: 'rs5413864115', source: 'EVA', ems: false }], ems: false,
      consequence: '3_prime_UTR_variant', ref_verified: true, designable: true, issues: [] }
  ],
  warnings: []
};

const LOOKUP_EXAMPLE = {
  requested_id: 'tmp_1_11502_C_CGT',
  system_name: 'sorghum_bicolor',
  source: { name: 'ensembl', release: '115' },
  variants: [
    { key: '1:11502:C:CGT', ids: ['tmp_1_11502_C_CGT', 'rs5413863549'], synonyms: [], label: '1:11502^11503 -/GT',
      kind: 'insertion', region: '1',
      vcf: { position: 11502, ref: 'C', alt: 'CGT' }, minimal: { start: 11503, end: 11502, ref: '-', alt: 'GT' },
      alleles: ['-', 'GT'], multiallelic: null, shift: 0, zone: { start: 11502, end: 11503 },
      discriminating: { forward: site(11503, 'A', 'G', null), reverse: site(11502, 'C', 'T', null) },
      records: [{ id: 'tmp_1_11502_C_CGT', source: 'SAP_PMID35653240_Boatwri', ems: false },
        { id: 'rs5413863549', source: 'EVA', ems: false }],
      ems: false, consequence: '3_prime_UTR_variant', ref_verified: true, designable: true, issues: [] }
  ],
  warnings: [
    { code: 'DUPLICATE_VARIANT_IDS', message: '2 Ensembl ids describe the same event 1:11502:C:CGT; they were merged',
      details: { key: '1:11502:C:CGT', ids: ['tmp_1_11502_C_CGT', 'rs5413863549'] } }
  ]
};

// §2.9 variant block, without the design-only submission_sequence (M6 adds it).
const DESIGN_VARIANT_RS871475760 = {
  key: '1:11109:C:A', requested_id: 'rs871475760', ids: ['rs871475760'], synonyms: ['tmp_1_11109_C_A'], label: '1:11109 C/A',
  kind: 'snv', region: '1', vcf: { position: 11109, ref: 'C', alt: 'A' }, minimal: { start: 11109, end: 11109, ref: 'C', alt: 'A' },
  alleles: ['C', 'A'], multiallelic: null, shift: 0, zone: { start: 11109, end: 11109 },
  discriminating: { forward: site(11109, 'C', 'A', 11109), reverse: site(11109, 'C', 'A', 11109) },
  records: [{ id: 'rs871475760', source: 'EVA', ems: false }], ems: false, consequence: 'downstream_gene_variant',
  ref_verified: true, designable: true, issues: []
};

// ---- GET /primers/variants -------------------------------------------------------------------------------------

test('§2.3 example: 1:11180-11290 deep-equals the documented body, from one chunk fetch and one sequence read', async function () {
  const w = world();
  const body = await variation.listVariants({ system_name: 'sorghum_bicolor', region: '1', start: 11180, end: 11290 }, w.deps);
  body.should.eql(LIST_EXAMPLE);
  w.calls.should.eql(['https://data.gramene.org/pansite-ensembl-115/overlap/region/sorghum_bicolor/1:10001-20000?feature=variation;content-type=application/json']);
  w.reads.should.eql([['1', 11180 - 1250, 11290 + 1250]]);
});

test('list filters: types, include_ems, limit with truncated and VARIANTS_TRUNCATED; query strings are accepted; the chunk is reused', async function () {
  const w = world();
  const q = function (extra) { return Object.assign({ system_name: 'sorghum_bicolor', region: '1', start: '11180', end: '11290' }, extra); };
  const keys = function (body) { return body.variants.map(function (v) { return v.key; }); };
  keys(await variation.listVariants(q({ types: ['snv'] }), w.deps)).should.eql(['1:11182:A:G', '1:11193:C:T', '1:11203:C:T']);
  keys(await variation.listVariants(q({ types: 'deletion,insertion' }), w.deps)).should.eql(['1:11282:CA:C']);
  keys(await variation.listVariants(q({ include_ems: false }), w.deps)).should.eql(['1:11182:A:G', '1:11282:CA:C']);
  keys(await variation.listVariants(q({ include_ems: 'false' }), w.deps)).should.eql(['1:11182:A:G', '1:11282:CA:C']);
  const cut = await variation.listVariants(q({ limit: 2 }), w.deps);
  cut.should.match({ start: 11180, end: 11290, total: 4, returned: 2, truncated: true });
  keys(cut).should.eql(['1:11182:A:G', '1:11193:C:T']);
  cut.warnings.should.eql([{ code: 'VARIANTS_TRUNCATED', message: 'showing 2 of 4 variants; narrow the window', details: { returned: 2, total: 4, limit: 2 } }]);
  w.calls.should.have.length(1);
});

test('list errors: every rule before I/O is checked first, and nothing reaches Ensembl or the FASTA', async function () {
  const w = world({ overrides: { variation: { species: { sorghum_nofasta: 'sorghum_nofasta' } } } });
  const q = function (extra) { return Object.assign({ system_name: 'sorghum_bicolor', region: '1', start: 11180, end: 11290 }, extra); };
  await rejectsWith(variation.listVariants(q({ start: 11291 }), w.deps), 400, 'REGION_OUT_OF_BOUNDS', { region: '1', start: 11291, end: 11290 });
  await rejectsWith(variation.listVariants(q({ start: 0 }), w.deps), 400, 'REGION_OUT_OF_BOUNDS');
  await rejectsWith(variation.listVariants(q({ start: 1, end: 60000 }), w.deps), 400, 'VARIANT_WINDOW_TOO_LONG', { length: 60000, max: 50000 });
  await rejectsWith(variation.listVariants(q({ limit: 5001 }), w.deps), 400, 'INVALID_REQUEST', { field: 'limit' });
  await rejectsWith(variation.listVariants(q({ types: ['snp'] }), w.deps), 400, 'INVALID_REQUEST', { field: 'types' });
  await rejectsWith(variation.listVariants(q({ include_ems: 'maybe' }), w.deps), 400, 'INVALID_REQUEST', { field: 'include_ems' });
  await rejectsWith(variation.listVariants(q({ region: '' }), w.deps), 400, 'INVALID_REQUEST', { field: 'region' });
  w.resolves.should.eql([]);

  await rejectsWith(variation.listVariants(q({ system_name: 'sorghum_nope' }), w.deps), 404, 'UNKNOWN_GENOME');
  const noData = await rejectsWith(variation.listVariants(q({ system_name: 'sorghum_rio' }), w.deps), 422, 'NO_VARIATION_DATA', { system_name: 'sorghum_rio' });
  noData.message.should.equal('sorghum_rio has no known-variant data; enter the variant as region, position, ref and alt');
  await rejectsWith(variation.listVariants(q({ system_name: 'sorghum_nofasta' }), w.deps), 422, 'NO_SEQUENCE');
  await rejectsWith(variation.listVariants(q({ region: '99' }), w.deps), 404, 'UNKNOWN_REGION', { region: '99' });
  await rejectsWith(variation.listVariants(q({ region: '2', start: 4000, end: 5001 }), w.deps), 400, 'REGION_OUT_OF_BOUNDS',
    { region: '2', start: 4000, end: 5001, length: REGION2_LENGTH });
  w.calls.should.eql([]);
  w.reads.should.eql([]);

  const off = world({ overrides: { variation: { enabled: false } } });
  const e = await rejectsWith(variation.listVariants(q(), off.deps), 503, 'FEATURE_DISABLED');
  e.details.retry_after_s.should.equal(300);
  await rejectsWith(variation.lookupVariant({ system_name: 'sorghum_bicolor', variant_id: 'rs871475760' }, off.deps), 503, 'FEATURE_DISABLED');
  off.resolves.should.eql([]);
  off.calls.should.eql([]);
});

test('§3.8: Ensembl down answers 503 VARIATION_SOURCE_UNAVAILABLE with reason and retry_after_s; the repeat comes from the negative cache', async function () {
  const w = world({ fetch: function () { throw refused(); } });
  const q = { system_name: 'sorghum_bicolor', region: '1', start: 11180, end: 11290 };
  const e = await rejectsWith(variation.listVariants(q, w.deps), 503, 'VARIATION_SOURCE_UNAVAILABLE', { retry_after_s: 30, reason: 'transport' });
  e.message.should.equal('the Ensembl variation service is temporarily unavailable');
  await rejectsWith(variation.listVariants(q, w.deps), 503, 'VARIATION_SOURCE_UNAVAILABLE');
  await rejectsWith(variation.lookupVariant({ system_name: 'sorghum_bicolor', variant_id: 'rs871475760' }, w.deps), 503, 'VARIATION_SOURCE_UNAVAILABLE');
  w.calls.should.have.length(2); // the overlap chunk once, the lookup once
  w.reads.should.eql([]);
});

test('skipped records inside the window become VARIATION_RECORDS_SKIPPED; those outside it do not; REF_MISMATCHES counts unverified entries', async function () {
  const genomeBase = baseAt(11250);
  const wrong = genomeBase === 'A' ? 'C' : 'A';
  const alt = wrong === 'G' ? 'T' : 'G';
  const extra = [
    { feature_type: 'variation', seq_region_name: '1', start: 11200, end: 11200, id: 'bad_alleles', source: 'EVA', alleles: ['C', 'Q'] },
    { feature_type: 'variation', seq_region_name: '1', start: 15000, end: 15000, id: 'bad_far_away', source: 'EVA', alleles: ['X', 'Y'] },
    { feature_type: 'variation', seq_region_name: '1', start: 'x', end: 11250, id: 'no_position', source: 'EVA', alleles: ['A', 'C'] },
    { feature_type: 'variation', seq_region_name: '1', start: 11250, end: 11250, id: 'ref_differs', source: 'EVA', alleles: [wrong, alt] }
  ];
  const w = world({ extraRecords: extra });
  const body = await variation.listVariants({ system_name: 'sorghum_bicolor', region: '1', start: 11180, end: 11290 }, w.deps);
  body.total.should.equal(5);
  body.warnings.should.eql([
    { code: 'VARIATION_RECORDS_SKIPPED', message: '2 malformed Ensembl variation records were skipped',
      details: { count: 2, reasons: { invalid_alleles: 1, invalid_coordinates: 1 } } },
    { code: 'REF_MISMATCHES', message: '1 variant has a reference allele that differs from this genome', details: { count: 1 } }
  ]);
  const bad = body.variants.find(function (v) { return v.ids[0] === 'ref_differs'; });
  bad.should.match({ key: '1:11250:' + wrong + ':' + alt, ref_verified: false, designable: false });
  bad.issues.map(function (i) { return [i.code, i.details]; }).should.eql([['REF_MISMATCH', { region: '1', position: 11250, given: wrong, genome: genomeBase }]]);
});

// ---- GET /primers/variants/{variant_id} ----------------------------------------------------------------------------

test('§2.4 example: tmp_1_11502_C_CGT deep-equals the documented body; rs871475760 keeps synonym tmp_1_11109_C_A and drops "."', async function () {
  const w = world();
  (await variation.lookupVariant({ system_name: 'sorghum_bicolor', variant_id: 'tmp_1_11502_C_CGT' }, w.deps)).should.eql(LOOKUP_EXAMPLE);
  w.calls.should.eql([
    'https://data.gramene.org/pansite-ensembl-115/variation/sorghum_bicolor/tmp_1_11502_C_CGT?content-type=application/json',
    'https://data.gramene.org/pansite-ensembl-115/overlap/region/sorghum_bicolor/1:10001-20000?feature=variation;content-type=application/json'
  ]);
  w.reads.should.eql([['1', 11502 - 1250, 11503 + 1250]]);

  const rs = await variation.lookupVariant({ system_name: 'sorghum_bicolor', variant_id: 'rs871475760' }, w.deps);
  rs.variants.map(function (v) { return [v.key, v.ids, v.synonyms]; }).should.eql([['1:11109:C:A', ['rs871475760'], ['tmp_1_11109_C_A']]]);
  rs.warnings.should.eql([]);
  w.calls.should.have.length(3); // the chunk is shared

  // a deletion described once, with its synonym; a multi-allelic id returns one entry per alt
  const del = await variation.lookupVariant({ system_name: 'sorghum_bicolor', variant_id: 'rs5413864115' }, w.deps);
  del.variants.map(function (v) { return [v.key, v.synonyms]; }).should.eql([['1:11282:CA:C', ['tmp_1_11282_CA_C']]]);
  const multi = await variation.lookupVariant({ system_name: 'sorghum_bicolor', variant_id: 'rs5413863494' }, w.deps);
  multi.variants.map(function (v) { return [v.key, v.multiallelic.other_alts]; }).should.eql([['1:10718:C:G', ['T']], ['1:10718:C:T', ['G']]]);
});

test('lookup errors: unknown id, no mapping on this assembly, several mappings, bad ids; a lookup-only id falls back to its mapping', async function () {
  const lookups = {
    elsewhere: { name: 'elsewhere', synonyms: [], mappings: [{ seq_region_name: 'Pt', start: 10, end: 10, allele_string: 'A/G' }] },
    beyond_end: { name: 'beyond_end', synonyms: [], mappings: [{ seq_region_name: '1', start: CHR1_LENGTH + 5, end: CHR1_LENGTH + 5, allele_string: 'A/G' }] },
    two_places: { name: 'two_places', synonyms: [], mappings: [
      { seq_region_name: '1', start: 11109, end: 11109, allele_string: 'C/A' },
      { seq_region_name: 'Pt', start: 5, end: 5, allele_string: 'A/G' },
      { seq_region_name: '1', start: 20000, end: 20000, allele_string: 'A/G' }] },
    listed_twice: { name: 'listed_twice', synonyms: ['.'], mappings: [
      { seq_region_name: '1', start: 11193, end: 11193, allele_string: 'C/T' },
      { seq_region_name: '1', start: 11193, end: 11193, allele_string: 'C/T' }] },
    lookup_only: { name: 'lookup_only', synonyms: [], most_severe_consequence: 'intergenic_variant',
      mappings: [{ seq_region_name: '1', start: 11109, end: 11109, allele_string: 'C/A' }] }
  };
  const w = world({ lookups: lookups });
  const q = function (id) { return { system_name: 'sorghum_bicolor', variant_id: id }; };
  const unknown = await rejectsWith(variation.lookupVariant(q('rs0000000001'), w.deps), 404, 'UNKNOWN_VARIANT', { id: 'rs0000000001', system_name: 'sorghum_bicolor' });
  unknown.message.should.equal('unknown variant rs0000000001 for sorghum_bicolor');
  await rejectsWith(variation.lookupVariant(q('elsewhere'), w.deps), 422, 'VARIANT_NOT_ON_ASSEMBLY', { id: 'elsewhere' });
  await rejectsWith(variation.lookupVariant(q('beyond_end'), w.deps), 422, 'VARIANT_NOT_ON_ASSEMBLY', { id: 'beyond_end' });
  await rejectsWith(variation.lookupVariant(q('two_places'), w.deps), 422, 'AMBIGUOUS_VARIANT_MAPPING', {
    id: 'two_places',
    mappings: [{ region: '1', start: 11109, end: 11109, allele_string: 'C/A' }, { region: '1', start: 20000, end: 20000, allele_string: 'A/G' }]
  });
  (await variation.lookupVariant(q('listed_twice'), w.deps)).variants.map(function (v) { return v.ids; })
    .should.eql([['listed_twice', 'tmp_1_11193_C_T']]);

  const only = await variation.lookupVariant(q('lookup_only'), w.deps);
  only.variants.should.have.length(1);
  only.variants[0].should.match({ key: '1:11109:C:A', ids: ['lookup_only', 'rs871475760'], consequence: 'intergenic_variant', ems: false });
  only.variants[0].records.should.eql([{ id: 'lookup_only', source: null, ems: false }, { id: 'rs871475760', source: 'EVA', ems: false }]);
  only.warnings.map(function (x) { return x.code; }).should.eql(['DUPLICATE_VARIANT_IDS']);

  const before = w.calls.length;
  for (const bad of ['../etc', 'a/b', '', 'x'.repeat(256), null]) {
    await rejectsWith(variation.lookupVariant(q(bad), w.deps), 400, 'INVALID_REQUEST', { field: 'variant_id' });
  }
  await rejectsWith(variation.lookupVariant({ system_name: 'sorghum_rio', variant_id: 'rs871475760' }, w.deps), 422, 'NO_VARIATION_DATA');
  w.calls.should.have.length(before);
});

// ---- genotyping design: resolveDesignVariant --------------------------------------------------------------------

test('resolveDesignVariant by id: the §2.9 variant block with requested_id after key, from one sequence read covering the template flank', async function () {
  const w = world();
  const r = await variation.resolveDesignVariant({ system_name: 'sorghum_bicolor', variant: { id: 'rs871475760', alt: 'A' } }, w.deps);
  r.variant.should.eql(DESIGN_VARIANT_RS871475760);
  Object.keys(r.variant).slice(0, 3).should.eql(['key', 'requested_id', 'ids']);
  r.should.match({ region: '1', region_length: CHR1_LENGTH, species: 'sorghum_bicolor', has_variation: true, fasta: '/fake/sorghum_bicolor.fa.gz', warnings: [] });
  r.resolved.system_name.should.equal('sorghum_bicolor');
  w.reads.should.eql([['1', 11109 - 1650, 11109 + 1650]]);
  r.genome.start.should.equal(11109 - 1650);
  r.genome.bases(11109, 11109).should.equal('C');

  // a larger template flank widens the same single read; aliases give DUPLICATE_VARIANT_IDS
  const ins = await variation.resolveDesignVariant({ system_name: 'sorghum_bicolor', variant: { id: 'tmp_1_11502_C_CGT' }, flank: 1040 }, w.deps);
  ins.variant.should.match({ key: '1:11502:C:CGT', requested_id: 'tmp_1_11502_C_CGT', ids: ['tmp_1_11502_C_CGT', 'rs5413863549'] });
  ins.warnings.should.eql([LOOKUP_EXAMPLE.warnings[0]]);
  w.reads[1].should.eql(['1', 11502 - 2290, 11503 + 2290]);
});

test('resolveDesignVariant: ALT_REQUIRED in site order, ALT_NOT_AT_SITE, UNSUPPORTED_ALLELE, NO_VARIATION_DATA, FEATURE_DISABLED, and a 503 when Ensembl is down', async function () {
  const w = world();
  const byId = function (variant, system) { return variation.resolveDesignVariant({ system_name: system || 'sorghum_bicolor', variant: variant }, w.deps); };
  const required = await rejectsWith(byId({ id: 'rs5413863494' }), 400, 'ALT_REQUIRED', { id: 'rs5413863494', alts: ['T', 'G'] });
  required.message.should.equal('rs5413863494 has more than one alternative allele; choose one with variant.alt');
  (await byId({ id: 'rs5413863494', alt: 'g' })).variant.should.match({ key: '1:10718:C:G', multiallelic: { alleles: ['C', 'T', 'G'], other_alts: ['T'] } });
  await rejectsWith(byId({ id: 'rs5413863494', alt: 'A' }), 400, 'ALT_NOT_AT_SITE', { id: 'rs5413863494', alt: 'A', alleles: ['C', 'T', 'G'] });
  await rejectsWith(byId({ id: 'rs5413863494', alt: 'C' }), 400, 'ALT_NOT_AT_SITE');
  await rejectsWith(byId({ id: 'rs5413863901', alt: '*' }), 400, 'UNSUPPORTED_ALLELE', { allele: '*' });
  (await byId({ id: 'rs5413863901' })).variant.issues.map(function (i) { return i.code; }).should.eql(['STAR_ALLELE']);
  (await byId({ id: 'rs5413864115', alt: '-' })).variant.key.should.equal('1:11282:CA:C');
  await rejectsWith(byId({ id: 'rs0000000001' }), 404, 'UNKNOWN_VARIANT', { id: 'rs0000000001', system_name: 'sorghum_bicolor' });

  const calls = w.calls.length;
  await rejectsWith(byId({ id: 'rs871475760', region: '1' }), 400, 'INVALID_VARIANT', { reason: 'id_or_manual' });
  await rejectsWith(byId({ id: '../x' }), 400, 'INVALID_VARIANT', { reason: 'id_or_manual' });
  await rejectsWith(byId({ id: 'rs871475760' }, 'sorghum_rio'), 422, 'NO_VARIATION_DATA', { system_name: 'sorghum_rio' });
  w.calls.should.have.length(calls);

  const off = world({ overrides: { variation: { enabled: false } } });
  await rejectsWith(variation.resolveDesignVariant({ system_name: 'sorghum_bicolor', variant: { id: 'rs871475760' } }, off.deps), 503, 'FEATURE_DISABLED');
  (await variation.resolveDesignVariant({ system_name: 'sorghum_bicolor', variant: { region: '1', position: 11109, ref: 'C', alt: 'A' } }, off.deps))
    .should.match({ species: null, has_variation: false });
  off.calls.should.eql([]);

  const down = world({ fetch: function () { throw refused(); } });
  await rejectsWith(variation.resolveDesignVariant({ system_name: 'sorghum_bicolor', variant: { id: 'rs871475760' } }, down.deps), 503,
    'VARIATION_SOURCE_UNAVAILABLE', { retry_after_s: 30, reason: 'transport' });
});

test('resolveDesignVariant manual: the §2.10(e) REF_MISMATCH body; both input styles of one deletion with no Ensembl call; bounds and repeat limits', async function () {
  const w = world();
  const manual = function (variant, extra) {
    return variation.resolveDesignVariant(Object.assign({ system_name: 'sorghum_bicolor', variant: variant }, extra), w.deps);
  };
  const e = await rejection(manual({ region: '1', position: 11109, ref: 'A', alt: 'C' }, { template_only: true }));
  e.status.should.equal(400);
  e.toJSON().should.eql({ message: 'the reference allele A does not match the genome base C at 1:11109', code: 'REF_MISMATCH',
    details: { region: '1', position: 11109, given: 'A', genome: 'C' } });

  const ensemblStyle = await manual({ region: '1', position: 11283, ref: 'A', alt: '-' });
  const vcfStyle = await manual({ region: '1', position: 11282, ref: 'ca', alt: 'c' });
  ensemblStyle.variant.should.eql(vcfStyle.variant);
  ensemblStyle.variant.should.match({ key: '1:11282:CA:C', requested_id: null, ids: [], records: [], shift: 2, zone: { start: 11282, end: 11286 },
    alleles: ['A', '-'], ems: false, consequence: null, ref_verified: true, designable: true });
  ensemblStyle.should.match({ species: 'sorghum_bicolor', has_variation: true, region: '1', warnings: [] });
  w.calls.should.eql([]);
  // one read per call, the REF_MISMATCH one included (REF is verified against that read)
  w.reads.should.eql([['1', 11109 - 1650, 11109 + 1650], ['1', 11283 - 1650, 11283 + 1650], ['1', 11282 - 1650, 11283 + 1650]]);
  (await variation.resolveDesignVariant({ system_name: 'sorghum_rio', variant: { region: '1', position: 11109, ref: 'C', alt: 'A' } }, w.deps))
    .should.match({ species: null, has_variation: false, variant: { key: '1:11109:C:A' } });

  await rejectsWith(manual({ region: '99', position: 10, ref: 'A', alt: 'C' }), 404, 'UNKNOWN_REGION', { region: '99' });
  await rejectsWith(manual({ region: '2', position: 5001, ref: 'A', alt: 'C' }), 400, 'REGION_OUT_OF_BOUNDS', { region: '2', position: 5001, length: REGION2_LENGTH });
  await rejectsWith(manual({ region: '2', position: 1, ref: '-', alt: 'A' }), 400, 'REGION_OUT_OF_BOUNDS', { region: '2', position: 1, length: REGION2_LENGTH });
  await rejectsWith(manual({ region: '1', position: 11109, ref: '-', alt: '-' }), 400, 'INVALID_VARIANT', { reason: 'alleles' });
  await rejectsWith(manual({ region: '1', position: 11109, ref: 'C' }), 400, 'INVALID_VARIANT', { reason: 'id_or_manual' });
  await rejectsWith(manual({ region: '1', position: 11109, ref: 'C', alt: 'N' }), 400, 'UNSUPPORTED_ALLELE', { allele: 'N' });

  const tight = world({ overrides: { variation: { max_shift: 1 } } });
  const rep = await rejectsWith(variation.resolveDesignVariant({ system_name: 'sorghum_bicolor', variant: { region: '1', position: 11283, ref: 'A', alt: '-' } }, tight.deps),
    400, 'VARIANT_TOO_REPETITIVE', { region: '1', position: 11282, shift: 2, max: 1 });
  rep.message.should.equal('the event can slide more than 1 bp, so no stable discriminating base exists');
  await rejectsWith(variation.resolveDesignVariant({ system_name: 'sorghum_nope', variant: { region: '1', position: 1, ref: 'A', alt: 'C' } }, w.deps), 404, 'UNKNOWN_GENOME');
});

// ---- genotyping design: neighboursFor ----------------------------------------------------------------------------

test('neighboursFor: the §2.9-§2.10 summaries and DENSE_NEIGHBOURS from the cached chunk, over the design\'s own sequence read', async function () {
  const w = world();
  const cases = [
    [{ id: 'rs871475760', alt: 'A' }, { start: 10709, end: 11509 }, { variants: 57, non_ems: 45, ems: 12, dense_non_ems: 0 }],
    [{ id: 'tmp_1_11193_C_T' }, { start: 10793, end: 11593 }, { variants: 54, non_ems: 45, ems: 9, dense_non_ems: 1 }],
    [{ region: '1', position: 11283, ref: 'A', alt: '-' }, { start: 10883, end: 11685 }, { variants: 53, non_ems: 45, ems: 8, dense_non_ems: 4 }],
    [{ id: 'tmp_1_11502_C_CGT' }, { start: 11102, end: 11903 }, { variants: 59, non_ems: 49, ems: 10, dense_non_ems: 5 }]
  ];
  for (const [input, window, counts] of cases) {
    const reads = w.reads.length;
    const d = await variation.resolveDesignVariant({ system_name: 'sorghum_bicolor', variant: input }, w.deps);
    const n = await variation.neighboursFor({ system_name: 'sorghum_bicolor', variant: d.variant, window: window, resolved: d.resolved,
      region_length: d.region_length, genome: d.genome }, w.deps);
    n.summary.should.eql(Object.assign({ data: 'ensembl', window: window }, counts));
    n.data.should.equal('ensembl');
    n.entries.should.have.length(counts.variants);
    n.entries.some(function (e) { return e.key === d.variant.key; }).should.be.false();
    n.target.key.should.equal(d.variant.key);
    n.warnings.should.eql(counts.dense_non_ems > 2
      ? [{ code: 'DENSE_NEIGHBOURS', message: counts.dense_non_ems + ' known non-EMS variants lie within 30 bp of the variant', details: { count: counts.dense_non_ems, window: 30 } }]
      : []);
    w.reads.length.should.equal(reads + 1); // one read per design, neighbours included
  }
  // one overlap fetch for all four designs; the ids needed lookups
  w.calls.filter(function (u) { return /overlap/.test(u); }).should.have.length(1);
});

test('neighboursFor manual: ids and records come from Ensembl (§3.9); "none" without variation data; "unavailable" when Ensembl is down, a 503 for a design by id', async function () {
  const w = world();
  const d = await variation.resolveDesignVariant({ system_name: 'sorghum_bicolor', variant: { region: '1', position: 11502, ref: 'C', alt: 'CGT' } }, w.deps);
  d.variant.ids.should.eql([]);
  const n = await variation.neighboursFor({ system_name: 'sorghum_bicolor', variant: d.variant, window: { start: 11102, end: 11903 }, genome: d.genome }, w.deps);
  n.variant.should.match({ key: '1:11502:C:CGT', requested_id: null, ids: ['rs5413863549', 'tmp_1_11502_C_CGT'], ems: false, consequence: '3_prime_UTR_variant' });
  n.variant.records.map(function (r) { return r.source; }).should.eql(['EVA', 'SAP_PMID35653240_Boatwri']);
  n.warnings.map(function (x) { return x.code; }).should.eql(['DUPLICATE_VARIANT_IDS', 'DENSE_NEIGHBOURS']);
  d.variant.ids.should.eql([]); // the input is not mutated

  const rio = await variation.resolveDesignVariant({ system_name: 'sorghum_rio', variant: { region: '1', position: 11109, ref: 'C', alt: 'A' } }, w.deps);
  const none = await variation.neighboursFor({ system_name: 'sorghum_rio', variant: rio.variant, window: { start: 10709, end: 11509 } }, w.deps);
  none.should.match({ data: 'none', entries: [], target: null, variant: rio.variant });
  none.summary.should.eql({ data: 'none', window: { start: 10709, end: 11509 }, variants: 0, non_ems: 0, ems: 0, dense_non_ems: 0 });
  none.warnings.should.eql([{ code: 'NO_VARIATION_DATA', message: 'sorghum_rio has no known-variant data; neighbouring variants were not screened', details: { system_name: 'sorghum_rio' } }]);

  // Switched off: same code and data 'none', but the message must not claim that sorghum_bicolor lacks data.
  const off = world({ overrides: { variation: { enabled: false } } });
  const disabled = await variation.neighboursFor({ system_name: 'sorghum_bicolor', variant: d.variant, window: { start: 11102, end: 11903 }, genome: d.genome }, off.deps);
  disabled.should.match({ data: 'none', entries: [], target: null });
  disabled.warnings.should.eql([{ code: 'NO_VARIATION_DATA', message: 'known-variant lookups are disabled on this server; neighbouring variants were not screened', details: { system_name: 'sorghum_bicolor' } }]);
  off.calls.should.have.length(0);

  const down = world({ fetch: function () { throw refused(); } });
  const m = await variation.resolveDesignVariant({ system_name: 'sorghum_bicolor', variant: { region: '1', position: 11109, ref: 'C', alt: 'A' } }, down.deps);
  const req = { system_name: 'sorghum_bicolor', variant: m.variant, window: { start: 10709, end: 11509 }, genome: m.genome };
  const unavailable = await variation.neighboursFor(req, down.deps);
  unavailable.should.match({ data: 'unavailable', entries: [], variant: { ids: [] } });
  unavailable.summary.data.should.equal('unavailable');
  unavailable.warnings.should.eql([{ code: 'NEIGHBOURS_UNAVAILABLE', message: 'known variants near the primers could not be fetched from Ensembl and were not screened', details: { reason: 'transport' } }]);
  const byId = Object.assign({}, m.variant, { requested_id: 'rs871475760' });
  await rejectsWith(variation.neighboursFor(Object.assign({}, req, { variant: byId }), down.deps), 503, 'VARIATION_SOURCE_UNAVAILABLE');
  (await variation.neighboursFor(Object.assign({}, req, { variant: byId, strict: false }), down.deps)).data.should.equal('unavailable');
  (function () { variation.neighboursFor({ system_name: 'sorghum_bicolor', variant: m.variant, window: { start: 5, end: 4 } }, down.deps).catch(function () {}); }).should.not.throw();
  await rejection(variation.neighboursFor({ system_name: 'sorghum_bicolor', variant: m.variant, window: { start: 5, end: 4 } }, down.deps));
});

// ---- flags and helpers ------------------------------------------------------------------------------------------

test('hasVariation and variationInfo (§3.2): enabled, species listed and a FASTA; no resolver call for other genomes', async function () {
  const w = world({ overrides: { variation: { species: { sorghum_nofasta: 'sorghum_nofasta', sorghum_gone: 'sorghum_gone' } } } });
  (await variation.hasVariation('sorghum_bicolor', w.deps)).should.be.true();
  (await variation.hasVariation('sorghum_nofasta', w.deps)).should.be.false();
  (await variation.hasVariation('sorghum_gone', w.deps)).should.be.false();
  w.resolves.should.eql(['sorghum_bicolor', 'sorghum_nofasta', 'sorghum_gone']);
  (await variation.hasVariation('sorghum_rio', w.deps)).should.be.false();
  w.resolves.should.have.length(3);
  const off = world({ overrides: { variation: { enabled: false } } });
  (await variation.hasVariation('sorghum_bicolor', off.deps)).should.be.false();
  off.resolves.should.eql([]);

  const resolved = { fasta: { dna: '/x.fa.gz' } };
  variation.variationInfo('sorghum_bicolor', resolved, w.cfg).should.eql({ available: true, source: 'ensembl', release: '115' });
  variation.variationInfo('sorghum_bicolor', { fasta: { dna: null } }, w.cfg).should.eql({ available: false, source: null, release: null });
  variation.variationInfo('sorghum_bicolor', null, w.cfg).should.eql({ available: false, source: null, release: null });
  variation.variationInfo('sorghum_rio', resolved, w.cfg).should.eql({ available: false, source: null, release: null });
  variation.variationInfo('sorghum_bicolor', resolved, off.cfg).should.eql({ available: false, source: null, release: null });
  variation.variationInfo('sorghum_bicolor', resolved, { fasta_root: '/x' }).should.eql({ available: false, source: null, release: null });
});

test('withGenome re-reads a wider window after a SequenceWindowError, at most three times; chunkStarts covers [lo, hi]', async function () {
  const reads = [];
  const ctx = {
    cfg: {}, v: { max_shift: 10 }, fasta: '/x',
    sequence: { fetch: async function (fasta, region, s, e) { reads.push([s, e]); return 'A'.repeat(e - s + 1); } }
  };
  let failures = 1;
  const out = await variation._internal.withGenome(ctx, '1', 100000, { start: 1000, end: 2000 }, function (g) {
    if (failures-- > 0) throw new normalize.SequenceWindowError(2500, 2600, g);
    return g.bases(2500, 2600).length;
  });
  out.value.should.equal(101);
  reads.should.eql([[1000, 2000], [1000, 2860]]);
  const stuck = await rejection(variation._internal.withGenome(ctx, '1', 100000, { start: 1000, end: 2000 }, function (g) {
    throw new normalize.SequenceWindowError(g.end + 1, g.end + 1, g);
  }));
  stuck.should.be.instanceOf(normalize.SequenceWindowError);
  reads.length.should.equal(2 + 4);

  variation._internal.chunkStarts(11179, 11291, 10000).should.eql([10001]);
  variation._internal.chunkStarts(1, 25002, 10000).should.eql([1, 10001, 20001]);
  variation._internal.chunkStarts(10000, 10001, 10000).should.eql([1, 10001]);
  variation.KINDS.should.eql(['snv', 'mnv', 'insertion', 'deletion', 'complex']);
});
