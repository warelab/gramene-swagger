'use strict';

// Integration: POST /primers/design and the /primers middleware against a running dev server
// (spec §10.3, plan §V3 rows 1-8 and 10). Skipped unless PRIMERS_IT_BASE is set:
//
//   PRIMERS_IT_BASE=http://127.0.0.1:50111/sorghum_v11 node --test --test-concurrency=1 "test/primers/integration/**/*.test.js"
//
// Read-only. Every returned primer's genomic.blocks are fetched from fastaIdx (PRIMERS_FASTAIDX, default
// http://localhost:8888) on the primer's strand and must spell the primer. Two requests run megablast in
// the server (BLAST depth mask); --test-concurrency=1 keeps them from overlapping check_jobs' BLAST runs.

const { describe, it } = require('node:test');
const should = require('should');

const coords = require('../../../api/helpers/primers/coords');
const designLib = require('../../../api/helpers/primers/design');

const BASE = (process.env.PRIMERS_IT_BASE || '').replace(/\/+$/, '');
const SKIP = BASE ? false : 'set PRIMERS_IT_BASE, e.g. http://127.0.0.1:50111/sorghum_v11';
const FASTAIDX = (process.env.PRIMERS_FASTAIDX || 'http://localhost:8888').replace(/\/+$/, '');

const G200 = 'SORBI_3001G000200';
const G700 = 'SORBI_3001G000700';
const G087700 = 'SORBI_3004G087700';

async function http(method, path, opts) {
  opts = opts || {};
  const headers = Object.assign({}, opts.headers || {});
  let body = opts.raw;
  if (opts.json !== undefined) {
    body = JSON.stringify(opts.json);
    if (!Object.keys(headers).some(function (h) { return h.toLowerCase() === 'content-type'; })) headers['Content-Type'] = 'application/json';
  }
  const res = await fetch(BASE + path, { method: method, headers: headers, body: body, redirect: 'manual' });
  const text = await res.text();
  let json = null;
  try { json = JSON.parse(text); } catch (e) { /* not JSON */ }
  return { status: res.status, headers: res.headers, json: json, text: text };
}

async function design(body, expectStatus) {
  const r = await http('POST', '/primers/design', { json: body });
  if (expectStatus !== undefined) r.status.should.equal(expectStatus, JSON.stringify(r.json || r.text).slice(0, 500));
  r.headers.get('cache-control').should.equal('no-store');
  return r;
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

// Blocks fetched on the primer's genomic strand, joined 5'->3', must equal the primer.
async function verifyPrimer(systemName, primer) {
  const g = primer.genomic;
  should.exist(g, 'primer has no genomic location');
  g.blocks.length.should.be.aboveOrEqual(1);
  (g.strand === 1 || g.strand === -1).should.be.true();
  const ordered = g.strand === 1 ? g.blocks : g.blocks.slice().reverse();
  let seq = '';
  for (const b of ordered) seq += await fastaIdx(systemName, g.region, b.start, b.end, g.strand);
  seq.should.equal(primer.seq, systemName + ' ' + JSON.stringify(g));
}

async function verifyAllPrimers(systemName, r) {
  r.pairs.length.should.be.above(0);
  for (const p of r.pairs) {
    await verifyPrimer(systemName, p.left);
    await verifyPrimer(systemName, p.right);
  }
}

// Each pair has a primer spanning a junction by the §A.3 rule with the effective overlaps.
function checkJunctionPairs(r) {
  const min5 = r.settings.params.min_5_prime_overlap_of_junction;
  const min3 = r.settings.params.min_3_prime_overlap_of_junction;
  r.pairs.forEach(function (p) {
    (p.left.junction !== null || p.right.junction !== null).should.be.true('pair ' + p.rank + ' spans no junction');
    ['left', 'right'].forEach(function (side) {
      const o = p[side];
      if (!o.junction) return;
      r.template.features.junctions.should.containEql(o.junction.position);
      coords.spansJunction(side, o.start, o.end, o.junction.position, min5, min3).should.be.true();
      o.genomic.blocks.length.should.equal(2);
    });
    p.product.genomic_size.should.equal(p.product.genomic.end - p.product.genomic.start + 1);
  });
}

describe('POST /primers/design (dev server)', { skip: SKIP }, function () {
  it('V3.1 transcript SORBI_3001G000200: 1982 nt, junctions [397,493,597,...], junction-spanning pairs of 70-150 bp, blocks verify', async function () {
    const r = (await design({ mode: 'transcript', gene_id: G200 }, 200)).json;
    r.template.mode.should.equal('transcript');
    r.template.length.should.equal(1982);
    r.template.transcript_id.should.equal('SORBI_3001G000200.1');
    r.template.features.junctions.slice(0, 3).should.eql([397, 493, 597]);
    r.settings.preset.should.equal('qpcr');
    r.settings.junction_spanning.should.be.true();
    r.pairs.length.should.be.within(1, 5);
    r.pairs.forEach(function (p) {
      p.product_size.should.be.within(70, 150);
      p.product_tm.should.be.a.Number();
      p.product_size.should.equal(p.right.end - p.left.start + 1);
      p.left.seq.should.match(/^[ACGT]+$/);
      p.right.seq.should.match(/^[ACGT]+$/);
    });
    checkJunctionPairs(r);
    r.engine.primer3.should.equal('2.6.1');
    await verifyAllPrimers('sorghum_bicolor', r);
  });

  it('V3.2 gene SORBI_3001G000200 with flanks 200/100: 1:11080-15099(-), 4020 bp, exon1 201-597, CDS start 499, blocks verify', async function () {
    const r = (await design({ mode: 'gene', gene_id: G200, flank_up: 200, flank_down: 100 }, 200)).json;
    const t = r.template;
    [t.region, t.start, t.end, t.strand, t.length].should.eql(['1', 11080, 15099, -1, 4020]);
    [t.features.exons[0].start, t.features.exons[0].end].should.eql([201, 597]);
    t.features.cds.start.should.equal(499);
    t.seq.length.should.equal(4020);
    r.pairs.forEach(function (p) { p.product.should.not.have.property('genomic_size'); });
    await verifyAllPrimers('sorghum_bicolor', r);
  });

  it('V3.4 + strand SORBI_3001G000700 in gene and transcript modes: blocks verify', async function () {
    const g = (await design({ mode: 'gene', gene_id: G700 }, 200)).json;
    g.template.strand.should.equal(1);
    await verifyAllPrimers('sorghum_bicolor', g);
    const t = (await design({ mode: 'transcript', gene_id: G700 }, 200)).json;
    t.template.strand.should.equal(1);
    t.template.features.junctions.length.should.be.above(0);
    checkJunctionPairs(t);
    await verifyAllPrimers('sorghum_bicolor', t);
  });

  it('V3.5 region 1:11080-15099(-) with target, included, excluded and params: every pair obeys every constraint; settings echo', async function () {
    const params = { min_size: 19, max_size: 22, min_tm: 58, max_tm: 61, min_gc: 40, max_gc: 60, max_tm_diff: 2, product_size_ranges: [[200, 300]] };
    const r = (await design({
      mode: 'region', system_name: 'sorghum_bicolor', region: { region: '1', start: 11080, end: 15099, strand: -1 },
      target: [499, 50], included: [100, 3000], excluded: [[1000, 40]], params: params
    }, 200)).json;
    r.template.length.should.equal(4020);
    r.settings.params.should.have.properties(params);
    r.pairs.length.should.be.above(0);
    const eps = 1e-3;
    r.pairs.forEach(function (p) {
      p.left.end.should.be.below(499);
      p.right.start.should.be.above(548);
      p.left.start.should.be.aboveOrEqual(100);
      p.right.end.should.be.belowOrEqual(3099);
      coords.overlapsRuns(p.left.start, p.left.end, [[1000, 40]]).should.be.false();
      coords.overlapsRuns(p.right.start, p.right.end, [[1000, 40]]).should.be.false();
      p.product_size.should.be.within(200, 300);
      ['left', 'right'].forEach(function (side) {
        p[side].len.should.be.within(19, 22);
        p[side].tm.should.be.within(58 - eps, 61 + eps);
        p[side].gc.should.be.within(40 - eps, 60 + eps);
      });
      Math.abs(p.left.tm - p.right.tm).should.be.belowOrEqual(2 + eps);
    });
    await verifyAllPrimers('sorghum_bicolor', r);
  });

  it('V3.6 sequence with IUPAC codes: 200 with IUPAC_CONVERTED (not PRIMER3_INPUT_ERROR); genomic fields null', async function () {
    // 600 nt of sorghum_bicolor 1:14300-14899 on the - strand (exon 1 of SORBI_3001G000200 and intron), from
    // fastaIdx, with R, Y and K planted
    const cdna = (await fastaIdx('sorghum_bicolor', '1', 14300, 14899, -1));
    const seq = cdna.slice(0, 49) + 'R' + cdna.slice(50, 299) + 'Y' + cdna.slice(300, 449) + 'K' + cdna.slice(450, 600);
    const r = (await design({ mode: 'sequence', sequence: '>iupac_test\n' + seq.match(/.{1,60}/g).join('\n') }, 200)).json;
    r.warnings.map(function (w) { return w.code; }).should.containEql('IUPAC_CONVERTED');
    r.template.length.should.equal(600);
    [r.template.seq.charAt(49), r.template.seq.charAt(299), r.template.seq.charAt(449)].should.eql(['N', 'N', 'N']);
    should(r.template.region).be.null();
    r.pairs.length.should.be.above(0);
    r.pairs.forEach(function (p) {
      should(p.left.genomic).be.null();
      should(p.right.genomic).be.null();
      should(p.product.genomic).be.null();
      p.left.seq.should.not.match(/N/);
      // primers are template substrings: left as is, right reverse-complemented
      r.template.seq.slice(p.left.start - 1, p.left.end).should.equal(p.left.seq);
    });
  });

  it('sequence mode with two FASTA records: 200, one joined template, warning MULTIPLE_RECORDS; one record (plus an empty header) has none', async function () {
    const a = await fastaIdx('sorghum_bicolor', '1', 14300, 14599, -1);
    const b = await fastaIdx('sorghum_bicolor', '1', 13300, 13599, -1);
    const r = (await design({ mode: 'sequence', sequence: '>rec_a\n' + a.match(/.{1,60}/g).join('\n') + '\n>rec_b\n' + b.match(/.{1,60}/g).join('\n') + '\n' }, 200)).json;
    const w = r.warnings.filter(function (x) { return x.code === 'MULTIPLE_RECORDS'; });
    w.should.have.length(1);
    w[0].message.should.equal('2 FASTA records were joined into one template; primers may span the joins');
    r.template.length.should.equal(600);
    r.template.seq.should.equal(a + b);
    const one = (await design({ mode: 'sequence', sequence: '>empty\n>rec_a\n' + a }, 200)).json;
    one.warnings.map(function (x) { return x.code; }).should.not.containEql('MULTIPLE_RECORDS');
    one.template.length.should.equal(300);
  });

  it('template_only with a sequence shorter than every product size range: 200 preview; without template_only 400 INVALID_PARAMS; intervals still checked', async function () {
    const seq = await fastaIdx('sorghum_bicolor', '1', 14300, 14379, -1); // 80 nt; pcr preset product range [100, 1000]
    const r = (await design({ mode: 'sequence', sequence: seq, template_only: true }, 200)).json;
    r.template.length.should.equal(80);
    r.pairs.should.eql([]);
    should(r.explain).be.null();
    const full = await design({ mode: 'sequence', sequence: seq }, 400);
    full.json.should.match({ code: 'INVALID_PARAMS', details: { param: 'product_size_ranges', template_length: 80 } });
    const outside = await design({ mode: 'sequence', sequence: seq, template_only: true, target: [70, 20] }, 400);
    outside.json.code.should.equal('INTERVAL_OUT_OF_BOUNDS');
  });

  it('n_mask on sorghum_tx436pac 4:7547610-7564601 with target [2483,150]: product_tm is long_seq_tm of the UNMASKED product, also for products covering the mask', async function () {
    const r = (await design({
      mode: 'region', system_name: 'sorghum_tx436pac', region: { region: '4', start: 7547610, end: 7564601 },
      target: [2483, 150], avoid_repeats: true, repeat_mask_mode: 'n_mask'
    }, 200)).json;
    r.template.mask_source.should.equal('softmask');
    r.settings.repeat_mask_mode.should.equal('n_mask');
    r.pairs.length.should.be.above(0);
    let covering = 0;
    r.pairs.forEach(function (p) {
      p.product_tm.should.be.a.Number();
      const unmasked = r.template.seq.slice(p.product.start - 1, p.product.end);
      unmasked.length.should.equal(p.product_size);
      const expected = designLib.longSeqTm(unmasked, r.settings.params);
      Math.abs(p.product_tm - expected).should.be.belowOrEqual(1e-3, 'pair ' + p.rank + ' product_tm ' + p.product_tm + ' vs unmasked ' + expected);
      if (!coords.overlapsRuns(p.product.start, p.product.end, r.template.mask)) return;
      covering++;
      // Primer3 itself sees the N-masked product, whose long_seq_tm is several degrees lower
      const masked = unmasked.split('');
      r.template.mask.forEach(function (run) {
        for (let i = Math.max(run[0], p.product.start); i <= Math.min(run[0] + run[1] - 1, p.product.end); i++) masked[i - p.product.start] = 'N';
      });
      (expected - designLib.longSeqTm(masked.join(''), r.settings.params)).should.be.above(1);
    });
    covering.should.be.above(0, 'no product covers the mask; the n_mask correction is not exercised');
  });

  it('V3.8 avoid_repeats on sorghum_tx436pac 4:7547610-7564601 uses the real soft-mask; no primer overlaps it; blocks verify', async function () {
    const r = (await design({ mode: 'region', system_name: 'sorghum_tx436pac', region: { region: '4', start: 7547610, end: 7564601 }, avoid_repeats: true }, 200)).json;
    r.template.mask_source.should.equal('softmask');
    r.template.masked.should.be.true();
    r.template.masked_fraction.should.be.above(0);
    r.warnings.map(function (w) { return w.code; }).should.not.containEql('BLAST_DEPTH_MASK');
    r.pairs.length.should.be.above(0);
    r.pairs.forEach(function (p) {
      coords.overlapsRuns(p.left.start, p.left.end, r.template.mask).should.be.false();
      coords.overlapsRuns(p.right.start, p.right.end, r.template.mask).should.be.false();
    });
    await verifyAllPrimers('sorghum_tx436pac', r);
  });

  it('V3.8 SORBI_3004G087700 avoid_repeats + template_only: blast_depth mask, BLAST_DEPTH_MASK, pairs []', async function () {
    const r = (await design({ mode: 'gene', gene_id: G087700, avoid_repeats: true, template_only: true }, 200)).json;
    r.template.mask_source.should.equal('blast_depth');
    r.warnings.map(function (w) { return w.code; }).should.containEql('BLAST_DEPTH_MASK');
    r.pairs.should.eql([]);
    should(r.explain).be.null();
    r.template.mask.forEach(function (run) {
      run[0].should.be.aboveOrEqual(1);
      (run[0] + run[1] - 1).should.be.belowOrEqual(r.template.length);
    });
  });

  it('V3.7 error cases: missing Content-Type, unknown field, unknown gene, ../etc, junction overlap, bad JSON, limits', async function () {
    const body = { mode: 'transcript', gene_id: G200 };
    // no JSON Content-Type: a string body is sent as text/plain, a byte body with no Content-Type at all
    for (const opts of [{ raw: JSON.stringify(body) }, { raw: new TextEncoder().encode(JSON.stringify(body)) },
      { raw: 'mode=gene', headers: { 'Content-Type': 'application/x-www-form-urlencoded' } }]) {
      const r = await http('POST', '/primers/design', opts);
      r.status.should.equal(400);
      r.json.message.should.equal('Validation errors');
      validationCodes(r.json).should.containEql('INVALID_CONTENT_TYPE');
      r.headers.get('cache-control').should.equal('no-store');
    }

    let r = await design({ mode: 'gene', gene_id: 'X', bogus: 1 }, 400);
    r.json.message.should.equal('Validation errors');
    validationCodes(r.json).should.containEql('OBJECT_ADDITIONAL_PROPERTIES');

    r = await design({ mode: 'gene', gene_id: 'NOPE' }, 404);
    r.json.should.match({ code: 'UNKNOWN_GENE' });
    r.json.should.have.property('details');

    r = await design({ mode: 'region', system_name: '../etc', region: { region: '1', start: 1, end: 100 } }, 400);
    r.json.message.should.equal('Validation errors');
    validationCodes(r.json).should.containEql('PATTERN');

    r = await design({ mode: 'transcript', gene_id: G200, params: { max_size: 24, min_5_prime_overlap_of_junction: 13 } }, 400);
    r.json.code.should.equal('INVALID_PARAMS');
    r.json.details.param.should.equal('min_5_prime_overlap_of_junction');

    r = await http('POST', '/primers/design', { raw: '{"mode":', headers: { 'Content-Type': 'application/json' } });
    r.status.should.equal(400);
    should.exist(r.json, 'malformed JSON must still get a JSON error body');
    r.headers.get('cache-control').should.equal('no-store');

    r = await design({ mode: 'region', system_name: 'sorghum_bicolor', region: { region: '1', start: 1, end: 60000 } }, 400);
    r.json.code.should.equal('TEMPLATE_TOO_LONG');

    r = await design({ mode: 'region', system_name: 'no_such_genome', region: { region: '1', start: 1, end: 100 } }, 404);
    r.json.code.should.equal('UNKNOWN_GENOME');

    r = await design({ mode: 'gene', gene_id: G200, target: [5000, 10] }, 400);
    r.json.code.should.equal('INTERVAL_OUT_OF_BOUNDS');

    r = await design({ mode: 'transcript', gene_id: G200, transcript_id: 'SORBI_3001G000200.9' }, 404);
    r.json.code.should.equal('UNKNOWN_TRANSCRIPT');

    // body over 100 kb: body-parser 413, still JSON and no-store
    r = await http('POST', '/primers/design', { json: { mode: 'sequence', sequence: 'A'.repeat(60000), excluded: Array(50).fill([1, 1]), gene_id: 'x'.repeat(255), pad: 'x'.repeat(50000) } });
    r.status.should.equal(413);
    r.headers.get('cache-control').should.equal('no-store');
    should.exist(r.json);
  });

  it('V3.10 headers and routing: no-store on 200 and 400, JSON 404 for bare/unknown /primers paths, JSON 405 for wrong methods, Retry-After exposed', async function () {
    let r = await http('POST', '/primers/design', { json: { mode: 'transcript', gene_id: G200 }, headers: { Origin: 'http://example.org' } });
    r.status.should.equal(200);
    r.headers.get('cache-control').should.equal('no-store');
    r.headers.get('access-control-allow-origin').should.equal('*');
    r.headers.get('access-control-expose-headers').should.match(/Retry-After/i);

    r = await http('POST', '/primers/design', { json: { mode: 'gene', gene_id: 'X', bogus: 1 } });
    r.status.should.equal(400);
    r.headers.get('cache-control').should.equal('no-store');

    for (const [method, path] of [['GET', '/primers'], ['POST', '/primers'], ['GET', '/primers/'], ['GET', '/primers/nope'], ['POST', '/primers/design/extra']]) {
      r = await http(method, path);
      r.status.should.equal(404, method + ' ' + path + ': ' + r.text.slice(0, 200));
      r.json.code.should.equal('NOT_FOUND');
      r.headers.get('content-type').should.match(/application\/json/);
      r.headers.get('cache-control').should.equal('no-store');
    }

    r = await http('PUT', '/primers/check/0123456789abcdef0123456789abcdef', { json: {} });
    r.status.should.equal(405);
    r.json.code.should.equal('METHOD_NOT_ALLOWED');
    r.json.details.allowed_methods.should.eql(['GET']);
    r.headers.get('allow').should.equal('GET');
    r.headers.get('cache-control').should.equal('no-store');
    r.text.should.not.match(/<html|<pre>|at .*\.js:\d+/i);

    r = await http('GET', '/primers/design');
    r.status.should.equal(405);
    r.json.code.should.equal('METHOD_NOT_ALLOWED');
    r.headers.get('allow').should.equal('POST');

    r = await http('DELETE', '/primers/genomes?system_name=sorghum_bicolor');
    r.status.should.equal(405);
    r.json.code.should.equal('METHOD_NOT_ALLOWED');

    // CORS preflight for a JSON POST is answered by the global cors middleware
    r = await http('OPTIONS', '/primers/design', { headers: { Origin: 'http://example.org', 'Access-Control-Request-Method': 'POST', 'Access-Control-Request-Headers': 'content-type' } });
    r.status.should.be.within(200, 204);
    r.headers.get('access-control-allow-origin').should.equal('*');
  });
});
