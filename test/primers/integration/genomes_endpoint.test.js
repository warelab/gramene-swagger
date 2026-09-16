'use strict';

// Integration: GET /primers/genomes against a running dev server (spec §10.3 row 13, plan §V3 row 9), plus a
// smoke test that the pre-existing endpoints still answer on the same instance. Skipped unless PRIMERS_IT_BASE:
//
//   PRIMERS_IT_BASE=http://127.0.0.1:50111/sorghum_v11 node --test --test-concurrency=1 "test/primers/integration/**/*.test.js"
//
// Optional: PRIMERS_IT_SWAGGER_HOST=localhost:50111 asserts the SWAGGER_HOST override in GET /swagger.

const { describe, it } = require('node:test');
const should = require('should');

const BASE = (process.env.PRIMERS_IT_BASE || '').replace(/\/+$/, '');
const SKIP = BASE ? false : 'set PRIMERS_IT_BASE, e.g. http://127.0.0.1:50111/sorghum_v11';

// The four sorghum assemblies with a real soft-masked dna_sm (checked by sampling for lowercase).
const SOFT_MASKED_SORGHUM = ['sorghum_rio', 'sorghum_tx2783pac', 'sorghum_tx430nano', 'sorghum_tx436pac'];
// has_variation: the additive genotyping field (genotyping spec §2.2); genotyping_endpoints.test.js checks its values.
const GENOME_KEYS = ['system_name', 'display_name', 'taxon_id', 'map_id', 'is_query', 'has_sequence', 'has_blastdb',
  'has_cdna_blastdb', 'has_variation', 'repeat_masking', 'total_bases', 'warnings'];

async function http(method, path, opts) {
  opts = opts || {};
  const res = await fetch(BASE + path, { method: method, headers: opts.headers || {}, redirect: 'manual' });
  const text = await res.text();
  let json = null;
  try { json = JSON.parse(text); } catch (e) { /* not JSON */ }
  return { status: res.status, headers: res.headers, json: json, text: text };
}

function compareGenomes(a, b) {
  const c = String(a.display_name).localeCompare(String(b.display_name), 'en', { numeric: true, sensitivity: 'base' });
  return c !== 0 ? c : (a.system_name < b.system_name ? -1 : a.system_name > b.system_name ? 1 : 0);
}

describe('GET /primers/genomes (dev server)', { skip: SKIP }, function () {
  it('V3.9 sorghum_bicolor: 120 genomes, query first then by display_name, sorghum_rio soft_masked, no filesystem paths', async function () {
    const r = await http('GET', '/primers/genomes?system_name=sorghum_bicolor');
    r.status.should.equal(200, r.text.slice(0, 300));
    r.headers.get('cache-control').should.equal('no-store');
    const g = r.json;
    g.system_name.should.equal('sorghum_bicolor');
    g.species.should.eql({ taxon_id: 4558, name: 'Sorghum bicolor' });
    g.counts.total.should.equal(120);
    g.genomes.length.should.equal(120);
    g.genomes[0].system_name.should.equal('sorghum_bicolor');
    g.genomes[0].is_query.should.be.true();
    g.genomes.slice(1).forEach(function (x) { x.is_query.should.be.false(); });
    const rest = g.genomes.slice(1);
    rest.slice().sort(compareGenomes).map(function (x) { return x.system_name; }).should.eql(rest.map(function (x) { return x.system_name; }));
    g.genomes.forEach(function (x) {
      Object.keys(x).sort().should.eql(GENOME_KEYS.slice().sort());
      ['soft_masked', 'unmasked_copy', 'absent'].should.containEql(x.repeat_masking);
      x.system_name.should.match(/^sorghum_/);
    });
    g.counts.with_blastdb.should.equal(g.genomes.filter(function (x) { return x.has_blastdb; }).length);
    g.counts.with_cdna_blastdb.should.equal(g.genomes.filter(function (x) { return x.has_cdna_blastdb; }).length);
    g.counts.with_blastdb.should.equal(120);
    const rio = g.genomes.find(function (x) { return x.system_name === 'sorghum_rio'; });
    rio.repeat_masking.should.equal('soft_masked');
    g.genomes.filter(function (x) { return x.repeat_masking === 'soft_masked'; }).map(function (x) { return x.system_name; }).sort()
      .should.eql(SOFT_MASKED_SORGHUM);
    const bicolor = g.genomes[0];
    bicolor.should.match({ map_id: 'GCA_000003195.3', total_bases: 708735318, has_sequence: true, has_cdna_blastdb: true });
    // synthetic bins (UNANCHORED, UNPLACED, ...) are not real map regions: no sorghum assembly is ASSEMBLY_MISMATCH,
    // including sorghum_tx430nano whose map is a single UNANCHORED bin
    g.genomes.forEach(function (x) {
      x.warnings.map(function (w) { return w.code; }).should.not.containEql('ASSEMBLY_MISMATCH', x.system_name);
    });
    g.genomes.find(function (x) { return x.system_name === 'sorghum_tx430nano'; }).warnings.should.eql([]);
    // no filesystem paths or file names anywhere in the body
    r.text.should.not.match(/\/scratch|\/home|\/usr\/|\.fa\.gz|\.fai|\.nal\b|\.nin\b|toplevel/);
  });

  it('other species and error cases: zea_maysb73 alone, 404 UNKNOWN_GENOME, 400 validation for a bad or missing system_name', async function () {
    let r = await http('GET', '/primers/genomes?system_name=zea_maysb73');
    r.status.should.equal(200);
    r.json.genomes[0].system_name.should.equal('zea_maysb73');
    r.json.genomes.every(function (x) { return !/^sorghum_/.test(x.system_name); }).should.be.true();

    // a map whose only region is a synthetic bin gets no ASSEMBLY_MISMATCH
    r = await http('GET', '/primers/genomes?system_name=selaginella_moellendorffii');
    r.status.should.equal(200);
    r.json.genomes[0].system_name.should.equal('selaginella_moellendorffii');
    r.json.genomes[0].warnings.map(function (w) { return w.code; }).should.not.containEql('ASSEMBLY_MISMATCH');

    r = await http('GET', '/primers/genomes?system_name=no_such_genome');
    r.status.should.equal(404);
    r.json.code.should.equal('UNKNOWN_GENOME');
    r.headers.get('cache-control').should.equal('no-store');

    for (const q of ['?system_name=../etc', '?system_name=Sorghum_Bicolor', '']) {
      r = await http('GET', '/primers/genomes' + q);
      r.status.should.equal(400, q);
      r.json.message.should.equal('Validation errors');
      r.headers.get('cache-control').should.equal('no-store');
    }
  });

  it('pre-existing endpoints still work on the same instance; /swagger documents the primers paths', async function () {
    let r = await http('GET', '/genes?idList=SORBI_3001G000200');
    r.status.should.equal(200);
    const genes = r.json;
    should(Array.isArray(genes)).be.true();
    genes.map(function (x) { return x._id; }).should.containEql('SORBI_3001G000200');
    should(r.headers.get('cache-control')).not.equal('no-store');

    r = await http('GET', '/maps?rows=1');
    r.status.should.equal(200);
    should(Array.isArray(r.json)).be.true();
    r.json.length.should.equal(1);

    r = await http('GET', '/docs');
    [200, 301, 302].should.containEql(r.status);
    r = await http('GET', '/docs/');
    r.status.should.equal(200);
    r.text.should.match(/swagger/i);

    r = await http('GET', '/swagger');
    r.status.should.equal(200);
    const doc = r.json;
    ['/primers/design', '/primers/genomes', '/primers/check', '/primers/check/{job_id}'].forEach(function (p) {
      doc.paths.should.have.property(p);
    });
    doc.paths['/primers/design'].post.tags.should.eql(['Primer design']);
    if (process.env.PRIMERS_IT_SWAGGER_HOST) {
      doc.host.should.equal(process.env.PRIMERS_IT_SWAGGER_HOST);
      doc.schemes.should.eql(['http']);
    }
  });
});
