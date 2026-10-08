'use strict';

// Integration: the genotyping endpoints against a running dev server with live Ensembl (genotyping spec §7.6):
// GET /primers/variants, GET /primers/variants/{variant_id}, the GET /primers/genomes variation fields and
// POST /primers/genotyping/design. Skipped unless PRIMERS_IT_BASE is set:
//
//   HOST=127.0.0.1 PORT=50112 SWAGGER_HOST=localhost:50112 SWAGGER_SCHEMES=http PRIMERS_SITE_KEY=sorghum_v11_geno \
//     PRIMERS_GLOBAL_MAX_JOBS=1 node app.js
//   PRIMERS_IT_BASE=http://127.0.0.1:50112/sorghum_v11 node --test --test-concurrency=1 test/primers/integration/genotyping_endpoints.test.js
//
// Read-only. No check job is submitted (POST /primers/check is never called), so no worker is needed. The captures in
// test/primers/fixtures/docs/capture-*.json were recorded over HTTP from this dev API with live Ensembl REST release 115,
// which is static, so a live response must equal its capture. Primer footprints are fetched from fastaIdx
// (PRIMERS_FASTAIDX, default http://localhost:8888) and must spell the primers. The fake-Ensembl cases (Ensembl down or
// slow, the breaker, the slot-starvation test) are in variation_down.test.js, which starts its own API.

const { describe, it } = require('node:test');
const should = require('should');
const fs = require('fs');
const path = require('path');

const BASE = (process.env.PRIMERS_IT_BASE || '').replace(/\/+$/, '');
const SKIP = BASE ? false : 'set PRIMERS_IT_BASE, e.g. http://127.0.0.1:50112/sorghum_v11';
const FASTAIDX = (process.env.PRIMERS_FASTAIDX || 'http://localhost:8888').replace(/\/+$/, '');
const CAPTURES = path.join(__dirname, '..', 'fixtures', 'docs');
const ORIGIN = 'http://example.org';

const LIST_KEYS = ['1:11182:A:G', '1:11193:C:T', '1:11203:C:T', '1:11282:CA:C'];
const COMPLEMENT = { A: 'T', C: 'G', G: 'C', T: 'A' };

// §2.9 check.request
const RS871475760_CHECK_REQUEST = {
  system_name: 'sorghum_bicolor',
  mode: 'region',
  checks: ['specificity', 'pangenome'],
  pairs: [
    { id: 'S1_REF', left: 'AGCTTCTCTAAGTGGTTATCCGA', right: 'ATCTTTGACTAGCGAGAAATTCAG', expected: { region: '1', start: 11068, end: 11132 } },
    { id: 'S1_ALT', left: 'AGCTTCTCTAAGTGGTTATCCGA', right: 'ATCTTTGACTAGCGAGAAATTCAT', expected: { region: '1', start: 11068, end: 11132 } },
    { id: 'S2_REF', left: 'GGTTATCCGAATATAGTCATACTCTATTC', right: 'TCTTTGTCTACTGAGAAATCCAGA', expected: { region: '1', start: 11081, end: 11172 } },
    { id: 'S2_ALT', left: 'GGTTATCCGAATATAGTCATACTCTATTA', right: 'TCTTTGTCTACTGAGAAATCCAGA', expected: { region: '1', start: 11081, end: 11172 } }
  ],
  genotyping: {
    variant: { region: '1', position: 11109, ref: 'C', alt: 'A' },
    sets: [
      { id: 'S1', ref_pair: 'S1_REF', alt_pair: 'S1_ALT' },
      { id: 'S2', ref_pair: 'S2_REF', alt_pair: 'S2_ALT' }
    ]
  }
};

async function http(method, p, opts) {
  opts = opts || {};
  const headers = Object.assign({ Origin: ORIGIN }, opts.headers || {});
  let body = opts.raw;
  if (opts.json !== undefined) {
    body = JSON.stringify(opts.json);
    if (!Object.keys(headers).some(function (h) { return h.toLowerCase() === 'content-type'; })) headers['Content-Type'] = 'application/json';
  }
  const res = await fetch(BASE + p, { method: method, headers: headers, body: body, redirect: 'manual' });
  const text = await res.text();
  let json = null;
  try { json = JSON.parse(text); } catch (e) { /* not JSON */ }
  return { status: res.status, headers: res.headers, json: json, text: text };
}

function get(p) {
  return http('GET', p);
}

function design(body) {
  return http('POST', '/primers/genotyping/design', { json: body });
}

// Cache-Control: no-store and CORS on every response of the new endpoints.
function headersOk(r) {
  should(r.headers.get('cache-control')).equal('no-store');
  should(r.headers.get('access-control-allow-origin')).equal('*');
}

function expectStatus(r, status) {
  r.status.should.equal(status, JSON.stringify(r.json || r.text).slice(0, 500));
  headersOk(r);
  return r;
}

function expectError(r, status, code) {
  expectStatus(r, status);
  should.exist(r.json, 'error bodies are JSON');
  r.json.code.should.equal(code, JSON.stringify(r.json).slice(0, 500));
  return r.json;
}

function validationCodes(json) {
  const out = [];
  const walk = function (errs) {
    (errs || []).forEach(function (e) {
      if (e.code) out.push(e.code);
      walk(e.errors);
    });
  };
  walk(json && json.errors);
  return out;
}

function expectValidation(r, code) {
  expectStatus(r, 400);
  r.json.message.should.equal('Validation errors');
  validationCodes(r.json).should.containEql(code);
}

function capture(name) {
  return JSON.parse(fs.readFileSync(path.join(CAPTURES, 'capture-' + name + '.json'), 'utf8'));
}

// Sends a capture's request and checks the status; the caller compares the body after its own, more telling assertions.
async function replay(name) {
  const c = capture(name);
  const basePath = new URL(BASE).pathname;
  c.request.path.indexOf(basePath + '/').should.equal(0, 'capture ' + name + ' was recorded under another basePath');
  const r = await http(c.request.method, c.request.path.slice(basePath.length), c.request.body !== undefined ? { json: c.request.body } : {});
  expectStatus(r, c.status);
  return { r: r, c: c };
}

const regionCache = new Map();
async function fastaIdx(systemName, region, start, end, strand) {
  const key = [systemName, region, start, end, strand].join(':');
  if (!regionCache.has(key)) {
    const url = FASTAIDX + '/sequence/region/' + systemName + '/' + region + ':' + start + '..' + end + ':' + strand;
    const res = await fetch(url);
    if (!res.ok) throw new Error('fastaIdx HTTP ' + res.status + ' for ' + url);
    regionCache.set(key, String((await res.json()).seq).toUpperCase());
  }
  return regionCache.get(key);
}

// The genome under an oligo's footprint, 5'->3' on its strand.
async function footprint(systemName, g) {
  const ordered = g.strand === 1 ? g.blocks : g.blocks.slice().reverse();
  let seq = '';
  for (const b of ordered) seq += await fastaIdx(systemName, g.region, b.start, b.end, g.strand);
  return seq;
}

// as_ref and common spell the reference exactly. For an SNV, as_alt spells it up to its 3' base, which is the ALT allele
// on the primer's strand.
async function verifySetPrimers(systemName, set, variant) {
  for (const role of ['as_ref', 'common']) {
    const o = set.primers[role];
    o.target_seq.should.equal(o.matched_seq);
    (await footprint(systemName, o.genomic)).should.equal(o.matched_seq, set.id + ' ' + role + ' ' + JSON.stringify(o.genomic));
  }
  if (variant.kind !== 'snv') return;
  const alt = set.primers.as_alt;
  const genome = await footprint(systemName, alt.genomic);
  alt.matched_seq.slice(0, -1).should.equal(genome.slice(0, -1), set.id + ' as_alt');
  alt.matched_seq.slice(-1).should.equal(alt.genomic.strand === 1 ? variant.vcf.alt : COMPLEMENT[variant.vcf.alt]);
  alt.order_seq.should.equal((alt.tail_seq || '') + alt.target_seq);
}

function codes(warnings) {
  return (warnings || []).map(function (w) { return w.code; });
}

describe('GET /primers/variants (dev server, live Ensembl)', { skip: SKIP }, function () {
  it('§2.3 window 1:11180-11290: keys 1:11182:A:G, 1:11193:C:T, 1:11203:C:T, 1:11282:CA:C; equal to the capture; no-store and CORS', async function () {
    const { r, c } = await replay('variants-list-1_11180-11290');
    r.json.variants.map(function (v) { return v.key; }).should.eql(LIST_KEYS);
    r.json.should.match({
      system_name: 'sorghum_bicolor', region: '1', start: 11180, end: 11290, source: { name: 'ensembl', release: '115' },
      total: 4, returned: 4, truncated: false, warnings: []
    });
    r.json.variants.forEach(function (v) {
      v.synonyms.should.eql([]);
      v.ref_verified.should.be.true();
      v.designable.should.be.true();
    });
    r.json.variants.filter(function (v) { return v.ems; }).map(function (v) { return v.key; }).should.eql(['1:11193:C:T', '1:11203:C:T']);
    r.json.variants[3].should.match({ kind: 'deletion', shift: 2, zone: { start: 11282, end: 11286 }, minimal: { start: 11283, end: 11283, ref: 'A', alt: '-' } });
    r.json.should.eql(c.response);
  });

  it('types, include_ems and limit: deletion only; snv only; EMS entries dropped; limit 2 truncates with VARIANTS_TRUNCATED', async function () {
    const q = '/primers/variants?system_name=sorghum_bicolor&region=1&start=11180&end=11290';
    const keys = function (r) { return r.json.variants.map(function (v) { return v.key; }); };
    let r = expectStatus(await get(q + '&types=deletion'), 200);
    keys(r).should.eql(['1:11282:CA:C']);
    r.json.should.match({ total: 1, returned: 1, truncated: false });
    r = expectStatus(await get(q + '&types=snv'), 200);
    keys(r).should.eql(LIST_KEYS.slice(0, 3));
    r = expectStatus(await get(q + '&include_ems=false'), 200);
    keys(r).should.eql(['1:11182:A:G', '1:11282:CA:C']);
    r.json.total.should.equal(2);
    r = expectStatus(await get(q + '&limit=2'), 200);
    keys(r).should.eql(LIST_KEYS.slice(0, 2));
    r.json.should.match({ total: 4, returned: 2, truncated: true });
    r.json.warnings.should.have.length(1);
    r.json.warnings[0].should.match({ code: 'VARIANTS_TRUNCATED', details: { returned: 2, total: 4, limit: 2 } });
  });

  it('multi-allelic and * sites: rs5413863494 C/T/G gives 1:10718:C:G and 1:10718:C:T; rs5413863901 C/T/* merges tmp_1_11318_C_T with STAR_ALLELE', async function () {
    let r = expectStatus(await get('/primers/variants?system_name=sorghum_bicolor&region=1&start=10718&end=10718'), 200);
    r.json.variants.map(function (v) { return v.key; }).should.eql(['1:10718:C:G', '1:10718:C:T']);
    const byKey = {};
    r.json.variants.forEach(function (v) { byKey[v.key] = v; });
    byKey['1:10718:C:T'].should.match({ ids: ['rs5413863494'], multiallelic: { alleles: ['C', 'T', 'G'], other_alts: ['G'] }, designable: true });
    byKey['1:10718:C:G'].should.match({ ids: ['rs5413863494'], multiallelic: { alleles: ['C', 'T', 'G'], other_alts: ['T'] }, designable: true });

    r = expectStatus(await get('/primers/variants?system_name=sorghum_bicolor&region=1&start=11318&end=11318'), 200);
    r.json.variants.should.have.length(1);
    const v = r.json.variants[0];
    v.should.match({ key: '1:11318:C:T', ids: ['rs5413863901', 'tmp_1_11318_C_T'], multiallelic: { alleles: ['C', 'T', '*'], other_alts: ['*'] }, designable: true });
    codes(v.issues).should.eql(['STAR_ALLELE']);
  });

  it('errors: 422 NO_VARIATION_DATA (capture), 400 VARIANT_WINDOW_TOO_LONG, 400 REGION_OUT_OF_BOUNDS, 404 UNKNOWN_REGION and UNKNOWN_GENOME, validator 400s', async function () {
    const { r, c } = await replay('variants-list-no-variation-data');
    r.json.should.match({ code: 'NO_VARIATION_DATA', details: { system_name: 'sorghum_rio' } });
    r.json.should.eql(c.response);

    const q = '/primers/variants?system_name=sorghum_bicolor&region=1';
    expectError(await get(q + '&start=1&end=60000'), 400, 'VARIANT_WINDOW_TOO_LONG').details.should.eql({ length: 60000, max: 50000 });
    expectError(await get(q + '&start=80884390&end=80884400'), 400, 'REGION_OUT_OF_BOUNDS').details.should.match({ region: '1', length: 80884392 });
    expectError(await get('/primers/variants?system_name=sorghum_bicolor&region=nope&start=1&end=100'), 404, 'UNKNOWN_REGION');
    expectError(await get('/primers/variants?system_name=no_such_genome&region=1&start=1&end=100'), 404, 'UNKNOWN_GENOME');

    expectValidation(await get(q + '&end=100'), 'REQUIRED');
    expectValidation(await get(q + '&start=0&end=100'), 'MINIMUM');
    expectValidation(await get('/primers/variants?system_name=../etc&region=1&start=1&end=100'), 'PATTERN');
    expectValidation(await get(q + '&start=1&end=100&types=bogus'), 'ENUM_MISMATCH');
    expectValidation(await get(q + '&start=1&end=100&limit=5001'), 'MAXIMUM');
  });
});

describe('GET /primers/variants/{variant_id} (dev server, live Ensembl)', { skip: SKIP }, function () {
  it('tmp_1_11502_C_CGT: one entry 1:11502:C:CGT with both aliases and DUPLICATE_VARIANT_IDS; equal to the capture', async function () {
    const { r, c } = await replay('variants-lookup-tmp_1_11502_C_CGT');
    r.json.should.match({ requested_id: 'tmp_1_11502_C_CGT', system_name: 'sorghum_bicolor', source: { name: 'ensembl', release: '115' } });
    r.json.variants.should.have.length(1);
    r.json.variants[0].should.match({
      key: '1:11502:C:CGT', ids: ['tmp_1_11502_C_CGT', 'rs5413863549'], kind: 'insertion', zone: { start: 11502, end: 11503 },
      minimal: { start: 11503, end: 11502, ref: '-', alt: 'GT' }
    });
    r.json.variants[0].records.map(function (x) { return x.source; }).should.eql(['SAP_PMID35653240_Boatwri', 'EVA']);
    r.json.warnings.should.have.length(1);
    r.json.warnings[0].should.match({ code: 'DUPLICATE_VARIANT_IDS', details: { key: '1:11502:C:CGT', ids: ['tmp_1_11502_C_CGT', 'rs5413863549'] } });
    r.json.should.eql(c.response);
  });

  it('rs871475760: one entry 1:11109:C:A with synonyms [tmp_1_11109_C_A] (Ensembl\'s "." dropped)', async function () {
    const r = expectStatus(await get('/primers/variants/rs871475760?system_name=sorghum_bicolor'), 200);
    r.json.requested_id.should.equal('rs871475760');
    r.json.variants.should.have.length(1);
    r.json.variants[0].should.match({ key: '1:11109:C:A', ids: ['rs871475760'], synonyms: ['tmp_1_11109_C_A'], kind: 'snv', designable: true });
    r.json.variants[0].synonyms.should.eql(['tmp_1_11109_C_A']);
    r.json.warnings.should.eql([]);
  });

  it('id syntax: tmp_1_13549_TTA_T,* (percent-encoded), 254- and 255-character ids resolve or answer 404 UNKNOWN_VARIANT, never a validator 400', async function () {
    let r = await get('/primers/variants/tmp_1_13549_TTA_T%2C%2A?system_name=sorghum_bicolor');
    [200, 404].should.containEql(r.status, r.text.slice(0, 300));
    headersOk(r);
    if (r.status === 200) {
      r.json.variants.length.should.be.above(0);
      r.json.variants.forEach(function (v) { v.ids.should.containEql('tmp_1_13549_TTA_T,*'); });
    } else {
      r.json.code.should.equal('UNKNOWN_VARIANT');
    }
    for (const id of ['rs' + '1'.repeat(252), 'A'.repeat(255)]) {
      r = await get('/primers/variants/' + id + '?system_name=sorghum_bicolor');
      [200, 404].should.containEql(r.status, id.length + ' characters: ' + r.text.slice(0, 200));
      if (r.status === 404) r.json.should.match({ code: 'UNKNOWN_VARIANT', details: { id: id, system_name: 'sorghum_bicolor' } });
      headersOk(r);
    }
  });

  it('errors: 404 UNKNOWN_VARIANT (capture), 422 NO_VARIATION_DATA, validator 400 for a%2Fb, a 256-character id and a missing system_name', async function () {
    const { r, c } = await replay('variants-lookup-unknown-variant');
    r.json.should.eql({ message: 'unknown variant rs0000000001 for sorghum_bicolor', code: 'UNKNOWN_VARIANT', details: { id: 'rs0000000001', system_name: 'sorghum_bicolor' } });
    r.json.should.eql(c.response);

    expectError(await get('/primers/variants/rs871475760?system_name=sorghum_rio'), 422, 'NO_VARIATION_DATA').details.should.eql({ system_name: 'sorghum_rio' });
    expectError(await get('/primers/variants/rs871475760?system_name=no_such_genome'), 404, 'UNKNOWN_GENOME');
    expectValidation(await get('/primers/variants/a%2Fb?system_name=sorghum_bicolor'), 'PATTERN');
    expectValidation(await get('/primers/variants/' + 'A'.repeat(256) + '?system_name=sorghum_bicolor'), 'PATTERN');
    expectValidation(await get('/primers/variants/rs871475760'), 'REQUIRED');
  });
});

describe('GET /primers/genomes: variation fields (dev server)', { skip: SKIP }, function () {
  it('sorghum_bicolor: variation {available: true, source: ensembl, release: 115}; has_variation only on sorghum_bicolor; equal to the capture', async function () {
    const { r, c } = await replay('genomes-sorghum_bicolor');
    r.json.variation.should.eql({ available: true, source: 'ensembl', release: '115' });
    r.json.genomes.length.should.equal(120);
    r.json.genomes[0].should.match({ system_name: 'sorghum_bicolor', is_query: true, has_variation: true });
    r.json.genomes.filter(function (g) { return g.has_variation; }).map(function (g) { return g.system_name; }).should.eql(['sorghum_bicolor']);
    r.json.genomes.find(function (g) { return g.system_name === 'sorghum_rio'; }).has_variation.should.be.false();
    r.json.should.eql(c.response);
  });

  it('zea_maysb73 (no variation data): variation {available: false, source: null, release: null}; has_variation false', async function () {
    const r = expectStatus(await get('/primers/genomes?system_name=zea_maysb73'), 200);
    r.json.variation.should.eql({ available: false, source: null, release: null });
    r.json.genomes.length.should.be.above(0);
    r.json.genomes.forEach(function (g) { g.has_variation.should.be.false(); });
  });
});

describe('POST /primers/genotyping/design (dev server, live Ensembl)', { skip: SKIP }, function () {
  it('§2.9 rs871475760 KASP num_sets 2: set keys f9df650ad116 and 1accc54c262d; check.request as documented; equal to the capture; primers spell the genome', async function (t) {
    const started = Date.now();
    const { r, c } = await replay('genotyping-design-rs871475760-kasp');
    t.diagnostic('rs871475760 KASP design in ' + (Date.now() - started) + ' ms');
    const res = r.json;
    res.sets.map(function (s) { return s.key; }).should.eql(['f9df650ad116', '1accc54c262d']);
    res.sets.map(function (s) { return [s.id, s.orientation, s.quality, s.score]; }).should.eql([['S1', 'reverse', 'usable', 9.18], ['S2', 'forward', 'poor', 19.63]]);
    res.check.request.should.eql(RS871475760_CHECK_REQUEST);
    res.check.should.match({ set_ids: ['S1', 'S2'], unique_primers: 6, omitted_set_ids: [] });
    res.variant.should.match({ key: '1:11109:C:A', requested_id: 'rs871475760', ids: ['rs871475760'], synonyms: ['tmp_1_11109_C_A'] });
    res.variant.submission_sequence.should.equal('YCCTCAAAAAGCTTCTCTAAGTGGTTATCCGAATATAGTCATACTCTATT[C/A]TGAATTTCTCGCTAGTCAAAGATAACAAAAATAGCATATTCTGGATTTCT');
    res.template.should.match({ region: '1', start: 10709, end: 11509, strand: 1, length: 801, alt_length: 801 });
    res.template.seq.charAt(400).should.equal('C');
    res.template.alt_seq.charAt(400).should.equal('A');
    res.neighbours.should.eql({ data: 'ensembl', window: { start: 10709, end: 11509 }, variants: 57, non_ems: 45, ems: 12, dense_non_ems: 0 });
    [res.orientations.forward.status, res.orientations.reverse.status].should.eql(['ok', 'ok']);
    res.engine.should.eql({ primer3: '2.6.1', thermo: 'ntthal 2.6.1', genotyping_design: '1', variation_source: 'ensembl 115' });
    res.sets[0].order.map(function (o) { return o.name; }).should.eql(['rs871475760_S1_REF_FAM', 'rs871475760_S1_ALT_HEX', 'rs871475760_S1_COM']);
    res.warnings.should.eql([]);
    for (const s of res.sets) await verifySetPrimers('sorghum_bicolor', s, res.variant);
    res.should.eql(c.response);
  });

  it('§2.10(c) tmp_1_11502_C_CGT num_sets 1: forward blocked by rs5413863234 4 nt from the 3′ end; reverse set 53942cb55348; equal to the capture', async function () {
    const { r, c } = await replay('genotyping-design-tmp_1_11502_C_CGT');
    const res = r.json;
    res.variant.should.match({ key: '1:11502:C:CGT', requested_id: 'tmp_1_11502_C_CGT', ids: ['tmp_1_11502_C_CGT', 'rs5413863549'], kind: 'insertion' });
    res.orientations.forward.should.match({ status: 'blocked', reason: 'neighbour_at_3p', discriminating_position: 11503, relaxation_level: null, sets_found: 0 });
    res.orientations.forward.attempts.should.eql([]);
    res.orientations.forward.blockers.should.eql([
      { key: '1:11500:G:A', ids: ['rs5413863234'], label: '1:11500 G/A', start: 11500, end: 11500, alleles: 'G/A', ems: false, distance_from_3p: 4 }
    ]);
    res.orientations.reverse.status.should.equal('ok');
    res.sets.map(function (s) { return [s.key, s.orientation]; }).should.eql([['53942cb55348', 'reverse']]);
    codes(res.warnings).should.eql(['DUPLICATE_VARIANT_IDS', 'ORIENTATION_BLOCKED', 'DENSE_NEIGHBOURS']);
    res.warnings[1].details.should.eql({ orientation: 'forward', ids: ['rs5413863234'], distances: [4] });
    await verifySetPrimers('sorghum_bicolor', res.sets[0], res.variant);
    res.should.eql(c.response);
  });

  it('§2.10(b) manual Ensembl-style deletion 1:11283 A/-: key 1:11282:CA:C, ids filled from Ensembl, sets a5277232d8ab and cb3ef66afd37; equal to the capture', async function () {
    const { r, c } = await replay('genotyping-design-manual-deletion');
    const res = r.json;
    res.variant.should.match({ key: '1:11282:CA:C', requested_id: null, ids: ['rs5413864115'], kind: 'deletion', shift: 2, zone: { start: 11282, end: 11286 } });
    res.neighbours.data.should.equal('ensembl');
    res.sets.map(function (s) { return s.key; }).should.eql(['a5277232d8ab', 'cb3ef66afd37']);
    codes(res.warnings).should.eql(['SHIFTABLE_INDEL', 'DENSE_NEIGHBOURS', 'RELAXED_CONSTRAINTS', 'RELAXED_CONSTRAINTS']);
    res.check.request.genotyping.variant.should.eql({ region: '1', position: 11282, ref: 'CA', alt: 'C' });
    for (const s of res.sets) await verifySetPrimers('sorghum_bicolor', s, res.variant);
    res.should.eql(c.response);
  });

  it('manual SNV on sorghum_rio (no variation data): REF_MISMATCH names the genome base; with it, 200 with neighbours.data none and NO_VARIATION_DATA', async function (t) {
    const variant = { region: '1', position: 20000, ref: 'A', alt: 'C' };
    let r = await design({ system_name: 'sorghum_rio', variant: variant, assay: { num_sets: 1 } });
    let base = 'A';
    if (r.status === 400) {
      const e = expectError(r, 400, 'REF_MISMATCH');
      e.details.should.match({ region: '1', position: 20000, given: 'A' });
      base = e.details.genome;
    }
    base.should.match(/^[ACGT]$/);
    const alt = base === 'C' ? 'T' : 'C';
    r = expectStatus(await design({ system_name: 'sorghum_rio', variant: { region: '1', position: 20000, ref: base, alt: alt }, assay: { num_sets: 1 } }), 200);
    const res = r.json;
    res.variant.should.match({ key: '1:20000:' + base + ':' + alt, requested_id: null, ids: [], records: [] });
    res.neighbours.should.match({ data: 'none', variants: 0, non_ems: 0, ems: 0 });
    const w = res.warnings.find(function (x) { return x.code === 'NO_VARIATION_DATA'; });
    should.exist(w, JSON.stringify(res.warnings));
    w.details.should.eql({ system_name: 'sorghum_rio' });
    should(res.engine.variation_source).be.null();
    should.exist(res.orientations);
    t.diagnostic('sorghum_rio 1:20000 ' + base + '/' + alt + ': ' + res.sets.length + ' set(s), warnings ' + JSON.stringify(codes(res.warnings)));
    if (res.sets.length) {
      res.check.request.system_name.should.equal('sorghum_rio');
      res.check.request.genotyping.variant.should.eql({ region: '1', position: 20000, ref: base, alt: alt });
      await verifySetPrimers('sorghum_rio', res.sets[0], res.variant);
    }
  });

  it('design errors: 400 REF_MISMATCH (capture), 404 UNKNOWN_VARIANT, 400 ALT_REQUIRED and ALT_NOT_AT_SITE, 422 NO_VARIATION_DATA, 400 INVALID_VARIANT and INVALID_PARAMS, 404 UNKNOWN_GENOME and UNKNOWN_REGION, 400 REGION_OUT_OF_BOUNDS', async function () {
    const { r, c } = await replay('genotyping-design-ref-mismatch');
    r.json.should.eql({
      message: 'the reference allele A does not match the genome base C at 1:11109', code: 'REF_MISMATCH',
      details: { region: '1', position: 11109, given: 'A', genome: 'C' }
    });
    r.json.should.eql(c.response);

    const sb = 'sorghum_bicolor';
    expectError(await design({ system_name: sb, variant: { id: 'rs0000000001' } }), 404, 'UNKNOWN_VARIANT');
    (await design({ system_name: sb, variant: { id: 'rs0000000001' } })).json.should.eql({
      message: 'unknown variant rs0000000001 for sorghum_bicolor', code: 'UNKNOWN_VARIANT', details: { id: 'rs0000000001', system_name: sb }
    });
    expectError(await design({ system_name: sb, variant: { id: 'rs5413863494' } }), 400, 'ALT_REQUIRED').details.should.eql({ id: 'rs5413863494', alts: ['T', 'G'] });
    expectError(await design({ system_name: sb, variant: { id: 'rs5413863494', alt: 'A' } }), 400, 'ALT_NOT_AT_SITE')
      .details.should.eql({ id: 'rs5413863494', alt: 'A', alleles: ['C', 'T', 'G'] });
    expectError(await design({ system_name: 'sorghum_rio', variant: { id: 'rs871475760' } }), 422, 'NO_VARIATION_DATA').details.should.eql({ system_name: 'sorghum_rio' });
    expectError(await design({ system_name: sb, variant: { id: 'rs871475760', region: '1' } }), 400, 'INVALID_VARIANT').details.should.eql({ reason: 'id_or_manual' });
    expectError(await design({ system_name: sb, variant: { region: '1', position: 11109, ref: 'C' } }), 400, 'INVALID_VARIANT').details.should.eql({ reason: 'id_or_manual' });
    expectError(await design({ system_name: sb, variant: { region: '1', position: 11109, ref: '-', alt: '-' } }), 400, 'INVALID_VARIANT').details.should.eql({ reason: 'alleles' });
    expectError(await design({ system_name: sb, variant: { id: 'rs871475760' }, params: { min_size: 25, opt_size: 22, max_size: 20 } }), 400, 'INVALID_PARAMS');
    expectError(await design({ system_name: 'no_such_genome', variant: { region: '1', position: 11109, ref: 'C', alt: 'A' } }), 404, 'UNKNOWN_GENOME');
    expectError(await design({ system_name: sb, variant: { region: 'nope', position: 11109, ref: 'C', alt: 'A' } }), 404, 'UNKNOWN_REGION');
    expectError(await design({ system_name: sb, variant: { region: '1', position: 90000000, ref: 'C', alt: 'A' } }), 400, 'REGION_OUT_OF_BOUNDS');
  });

  it('validator 400s: unknown top-level and nested fields, variant.id and ref patterns, num_sets 11, assay.type enum, missing variant, no JSON Content-Type', async function () {
    const sb = 'sorghum_bicolor';
    expectValidation(await design({ system_name: sb, variant: { id: 'rs871475760' }, bogus: 1 }), 'OBJECT_ADDITIONAL_PROPERTIES');
    expectValidation(await design({ system_name: sb, variant: { id: 'rs871475760', extra: 1 } }), 'OBJECT_ADDITIONAL_PROPERTIES');
    expectValidation(await design({ system_name: sb, variant: { id: '../etc' } }), 'PATTERN');
    expectValidation(await design({ system_name: sb, variant: { region: '1', position: 11109, ref: 'N', alt: 'A' } }), 'PATTERN');
    expectValidation(await design({ system_name: sb, variant: { id: 'rs871475760' }, assay: { num_sets: 11 } }), 'MAXIMUM');
    expectValidation(await design({ system_name: sb, variant: { id: 'rs871475760' }, assay: { type: 'caps' } }), 'ENUM_MISMATCH');
    expectValidation(await design({ system_name: sb }), 'OBJECT_MISSING_REQUIRED_PROPERTY');
    const r = await http('POST', '/primers/genotyping/design', { raw: JSON.stringify({ system_name: sb, variant: { id: 'rs871475760' } }) });
    expectValidation(r, 'INVALID_CONTENT_TYPE');
  });
});
