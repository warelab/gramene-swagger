'use strict';

// Integration: POST /primers/check + GET /primers/check/{job_id} against a running dev server AND its check
// worker (spec §10.4, plan §V4). Skipped unless PRIMERS_IT_BASE is set:
//
//   HOST=127.0.0.1 PORT=50111 SWAGGER_HOST=localhost:50111 SWAGGER_SCHEMES=http PRIMERS_SITE_KEY=sorghum_v11_dev node app.js
//   PRIMERS_SITE_KEY=sorghum_v11_dev node api/helpers/primers/jobs/worker_main.js
//   PRIMERS_IT_BASE=http://127.0.0.1:50111/sorghum_v11 node --test --test-concurrency=1 "test/primers/integration/**/*.test.js"
//
// The checks run real BLAST in the worker (sorghum_bicolor reference at word size 5; three or four pan-genome
// assemblies at word size 6). Tests in this file run one after another and each waits for its job, so one check
// runs at a time. Job ids are deterministic, so a re-run within 24 h gets the finished jobs back (POST 200).
//
// Expectations are those of check algorithm version 2 (plan stringency decision): a product amplifies only if each
// primer site has <= max_amplifying_mismatches (default 3) edits and passes the Primer-BLAST 3' rule; 4-5 edit
// products are 'unlikely'; products whose two primer footprints overlap are discarded.
//
// Restart recovery and Redis-down behaviour are manual checks (plan Overrides), not tests here.
// JOB_TOO_LARGE cannot be reached with the default limit (6000 CPU-s; the largest sorghum_v11 job, 20 primers over
// all 119 other assemblies, is estimated at ~5200), so that test runs only with PRIMERS_IT_EXPECT_JOB_TOO_LARGE=1
// against a server started with a lower limit, e.g.
//   NODE_CONFIG='{"primers":{"check":{"max_job_cpu_s":1000}}}'  (never set PRIMERS_IT_EXPECT_JOB_TOO_LARGE against a
//   default server: the request would queue a ~5200 CPU-s job).

const { describe, it } = require('node:test');
const should = require('should');

const BASE = (process.env.PRIMERS_IT_BASE || '').replace(/\/+$/, '');
const SKIP = BASE ? false : 'set PRIMERS_IT_BASE, e.g. http://127.0.0.1:50111/sorghum_v11';
const JOB_TIMEOUT_MS = Number(process.env.PRIMERS_IT_JOB_TIMEOUT_MS) || 20 * 60 * 1000;
const POLL_MS = 2000;
const CAP = 3; // default params.max_amplifying_mismatches
const PARAM_KEYS = ['ignore_mismatches', 'include_unlikely', 'max_amplifying_mismatches', 'max_product_size', 'min_3p_mismatches',
  'min_total_mismatches', 'repeat_site_threshold', 'three_prime_window'];

// §10.4 primers
const P = Object.freeze({
  P1_L: 'GATCGACAATCCGACGATAGAAG',
  P1_R: 'GTGAACATCATGCTGCCCGATG',
  P2_L: 'GGACAGCTCCACAACATATCAG',
  P2_R: 'GGACATTTGAAGCCCATGGCC',
  P3_L: 'GATATCAGTGGAATCATAAGACCG',
  P3_R: 'CATCGATATCAGGATCTGGCTT',
  P5_L: 'CGACGAGGAGGCGCTGCGATG',
  m20_7_14: 'GGACAGATCCACATCATATC',
  del18: 'GGACAGCTCCACAACATTCAG',
  J_L: 'CCAACAAAGTCATGGATGCACT',
  PA_L: 'TTACCTCTCCTCCTTCCCGATC',
  PA_R: 'GCATAGTCTCATTAGCTAGCAC'
});

async function http(method, path, opts) {
  opts = opts || {};
  const headers = Object.assign({}, opts.headers || {});
  let body;
  if (opts.json !== undefined) {
    body = JSON.stringify(opts.json);
    headers['Content-Type'] = 'application/json';
  }
  const res = await fetch(BASE + path, { method: method, headers: headers, body: body, redirect: 'manual' });
  const text = await res.text();
  let json = null;
  try { json = JSON.parse(text); } catch (e) { /* not JSON */ }
  return { status: res.status, headers: res.headers, json: json, text: text };
}

const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));

async function submit(body) {
  const r = await http('POST', '/primers/check', { json: body });
  r.headers.get('cache-control').should.equal('no-store');
  return r;
}

// Poll until done or error; returns the final status document and the statuses seen.
async function poll(jobId, t) {
  const started = Date.now();
  const seen = [];
  let last = null;
  for (;;) {
    const r = await http('GET', '/primers/check/' + jobId);
    r.status.should.equal(200, r.text.slice(0, 300));
    r.headers.get('cache-control').should.equal('no-store');
    last = r.json;
    if (seen[seen.length - 1] !== last.status) seen.push(last.status);
    if (last.status === 'done' || last.status === 'error') break;
    if (Date.now() - started > JOB_TIMEOUT_MS) throw new Error('job ' + jobId + ' still ' + last.status + ' after ' + JOB_TIMEOUT_MS + ' ms');
    await sleep(POLL_MS);
  }
  if (t) t.diagnostic(jobId + ' ' + seen.join(' -> ') + ' in ' + Math.round((Date.now() - started) / 1000) + ' s; attempts ' + last.attempts + '; timings_ms ' + JSON.stringify(last.results && last.results.timings_ms));
  return { job: last, seen: seen };
}

// Every job: algorithm version 2 and all 8 params filled in both the normalized request and the results.
function checkEngineAndParams(job) {
  job.results.engine.should.match({ algorithm_version: '2', blast: '2.13.0' });
  Object.keys(job.results.params).sort().should.eql(PARAM_KEYS);
  Object.keys(job.request.params).sort().should.eql(PARAM_KEYS);
  job.results.params.should.eql(job.request.params);
}

async function runCheck(body, t) {
  const s = await submit(body);
  [200, 202].should.containEql(s.status, JSON.stringify(s.json).slice(0, 400));
  s.json.job_id.should.match(/^[0-9a-f]{32}$/);
  s.json.should.not.have.property('created');
  s.json.should.not.have.property('resolved');
  const p = await poll(s.json.job_id, t);
  p.job.status.should.equal('done', JSON.stringify(p.job.error));
  p.job.partial.should.be.false();
  should.exist(p.job.results);
  checkEngineAndParams(p.job);
  return { submitted: s, job: p.job, results: p.job.results };
}

const find = (list, region, start, end, orientation) => (list || []).find((a) => a.region === region && a.start === start && a.end === end && (!orientation || a.orientation === orientation));
const geneIds = (a) => (a && Array.isArray(a.genes) ? a.genes.map((g) => g.id) : []);
const amplifying = (a) => a.likelihood === 'likely' || a.likelihood === 'likely_weak';

// Stringency: every listed amplifying product has each primer within the cap; unlikely[] holds only unlikely ones.
function checkStringency(pair, cap) {
  const listed = (pair.on_target ? [pair.on_target] : []).concat(pair.off_targets || [], pair.gdna_products || []);
  listed.forEach(function (a) {
    amplifying(a).should.be.true(JSON.stringify(a).slice(0, 200));
    a.left_mm.should.be.belowOrEqual(cap);
    a.right_mm.should.be.belowOrEqual(cap);
  });
  (pair.unlikely || []).forEach(function (a) { a.likelihood.should.equal('unlikely'); });
}

function checkSummaries(results) {
  (results.pangenome ? results.pangenome.pairs : []).forEach(function (pair) {
    const s = pair.summary;
    (s.single_perfect + s.single_mismatch + s.multiple + s.no_amplicon + s.db_unavailable + s.error).should.equal(s.genomes_total);
    s.amplifies.should.equal(s.single_perfect + s.single_mismatch + s.multiple);
    s.truncated.should.equal(pair.genomes.filter((g) => g.truncated === true).length);
    s.genomes_total.should.equal(pair.genomes.length);
  });
}

function presenceStatuses(pair) {
  const by = {};
  pair.genomes.forEach(function (x) { by[x.system_name] = x; });
  // §10.4: 353 and is12661 no_amplicon (their best products need 4+ edits in a primer: unlikely, shown as nearest)
  ['sorghum_353', 'sorghum_is12661'].forEach(function (sys) {
    by[sys].should.match({ status: 'no_amplicon', ortholog_annotated: false });
    should(by[sys].primary).be.null();
    should.exist(by[sys].nearest, sys + ' nearest');
    by[sys].nearest.likelihood.should.equal('unlikely');
    Math.max(by[sys].nearest.left_mm, by[sys].nearest.right_mm).should.be.above(CAP);
  });
  by.sorghum_ji2731.should.match({ status: 'single_perfect', primary: { size_delta: -177, left_mm: 0, right_mm: 0 } });
  by.sorghum_leoti.should.match({ status: 'single_perfect', primary: { left_mm: 0, right_mm: 0, size_delta: 0 } });
  pair.summary.should.eql({ genomes_total: 4, single_perfect: 2, single_mismatch: 0, multiple: 0, no_amplicon: 2, db_unavailable: 0, error: 0, amplifies: 2, truncated: 0 });
  return by;
}

describe('POST /primers/check and GET /primers/check/{job_id} (dev server + worker)', { skip: SKIP }, function () {
  let presenceOnTarget = null;

  it('P2: POST 202 (or 200 if a previous run left the job), re-POST 200 with the same job_id; the two §10.4 off_targets and nothing else', { timeout: JOB_TIMEOUT_MS + 60000 }, async function (t) {
    const body = {
      system_name: 'sorghum_bicolor', mode: 'gene', gene_id: 'SORBI_3004G087700',
      pairs: [{ id: 'P2', left: P.P2_L, right: P.P2_R, expected: { region: '4', start: 7423537, end: 7423746 } }]
    };
    const first = await submit(body);
    [200, 202].should.containEql(first.status);
    t.diagnostic('first POST ' + first.status + ', estimate ' + JSON.stringify(first.json.estimate));
    first.json.should.match({ kind: 'specificity' });
    first.json.estimate.cpu_s.should.be.above(0);
    Object.keys(first.json).sort().should.eql(['created_at', 'estimate', 'job_id', 'kind', 'progress', 'queue_position', 'status', 'warnings']);
    // key order, case and an explicit default do not change the job id
    const second = await submit(Object.assign({}, body, {
      params: { max_amplifying_mismatches: CAP },
      pairs: [{ id: 'P2', left: P.P2_L.toLowerCase(), right: P.P2_R, expected: body.pairs[0].expected }]
    }));
    second.status.should.equal(200);
    second.json.job_id.should.equal(first.json.job_id);

    const { job } = await poll(first.json.job_id, t);
    job.status.should.equal('done', JSON.stringify(job.error));
    job.attempts.should.be.aboveOrEqual(1);
    should.exist(job.finished_at);
    job.request.pairs[0].left.should.equal(P.P2_L);
    job.request.checks.should.eql(['specificity']);
    job.request.params.max_amplifying_mismatches.should.equal(CAP);
    job.should.not.have.property('resolved');
    job.should.not.have.property('worker');
    checkEngineAndParams(job);
    const res = job.results;
    res.reference.should.eql({ system_name: 'sorghum_bicolor', map_id: 'GCA_000003195.3', total_bases: 708735318 });
    res.sensitivity_note.should.match(/at most 3 mismatches \(max_amplifying_mismatches\)/);
    res.sensitivity_note.should.match(/reference search \(word size 5\) finds every site with up to 3 mismatches for these 21-22 nt primers, which covers the cap/);
    res.sensitivity_note.should.not.match(/pan-genome search/);
    const pair = res.specificity.pairs[0];
    checkStringency(pair, CAP);
    pair.should.match({ verdict: 'off_targets', off_target_count: 2 });
    pair.on_target.should.match({ region: '4', start: 7423537, end: 7423746, orientation: 'LR', left_mm: 0, right_mm: 0, likelihood: 'likely' });
    geneIds(pair.on_target).should.containEql('SORBI_3004G087700');
    const para = find(pair.off_targets, '4', 7437317, 7437526, 'LR');
    should.exist(para, 'off-target 4:7437317-7437526');
    para.should.match({ left_mm: 0, right_mm: 0, likelihood: 'likely' });
    geneIds(para).should.containEql('SORBI_3004G087800');
    const chr5 = find(pair.off_targets, '5', 66890615, 66890824, 'RL');
    should.exist(chr5, 'off-target 5:66890615-66890824');
    chr5.should.match({ left_mm_pos: [10, 7], right_mm_pos: [15, 10], likelihood: 'likely' });
    geneIds(chr5).should.containEql('SORBI_3005G183900');
    res.primers[P.P2_L].should.match({ near_perfect_sites: 2, repetitive: false });
    res.primers[P.P2_L].sensitivity.reference.should.eql({ word_size: 5, guaranteed_max_mismatches: 3 });
  });

  it('P1, P3, P5_L, m20_7_14, del18 in one gene-mode job: §10.4 off-targets (4-edit sites now unlikely), repetitive, word-size-5 and gap sites', { timeout: JOB_TIMEOUT_MS + 60000 }, async function (t) {
    const { results } = await runCheck({
      system_name: 'sorghum_bicolor', mode: 'gene', gene_id: 'SORBI_3004G087700', params: { include_unlikely: true },
      pairs: [
        { id: 'P1', left: P.P1_L, right: P.P1_R, expected: { region: '4', start: 7422190, end: 7422843 } },
        { id: 'P3', left: P.P3_L, right: P.P3_R, expected: { region: '4', start: 7422482, end: 7423061 } },
        { id: 'P5_L', left: P.P5_L, right: P.P2_R },
        { id: 'm20_7_14', left: P.m20_7_14, right: P.P2_R },
        { id: 'del18', left: P.del18, right: P.P2_R }
      ]
    }, t);
    const pair = {};
    results.specificity.pairs.forEach(function (p) { pair[p.id] = p; });
    t.diagnostic('verdicts [id, verdict, off_target_count, unlikely_count] ' + JSON.stringify(results.specificity.pairs.map((p) => [p.id, p.verdict, p.off_target_count, p.unlikely_count])));
    results.specificity.pairs.forEach(function (p) { checkStringency(p, CAP); });

    // P1: 4:7435938-7436628 (691 bp) still amplifies (L 3 / R 2 edits); 5:66891530-66892235 (706 bp) has a 4-edit
    // L site (mm_pos [20,17,14,11]) and is unlikely (not re-aligned: its lower bound is over the cap)
    pair.P1.on_target.should.match({ region: '4', start: 7422190, end: 7422843, size: 654 });
    const p1a = find(pair.P1.off_targets, '4', 7435938, 7436628);
    should.exist(p1a, 'P1 off-target 4:7435938-7436628');
    p1a.should.match({ size: 691, likelihood: 'likely', left_mm: 3, right_mm: 2, approx: false });
    should(find(pair.P1.off_targets, '5', 66891530, 66892235)).be.undefined();
    const p1b = find(pair.P1.unlikely, '5', 66891530, 66892235);
    should.exist(p1b, 'P1 unlikely 5:66891530-66892235');
    p1b.should.match({ size: 706, likelihood: 'unlikely', left_mm: 4, right_mm: 0, approx: true });
    pair.P1.should.match({ verdict: 'off_targets', off_target_count: 2 });
    const nearP1L = (a) => a.region === '4' && (Math.abs(a.start - 7432183) <= 3 || Math.abs(a.end - 7432183) <= 3);
    pair.P1.off_targets.concat(pair.P1.unlikely || []).some(nearP1L).should.be.false('the P1_L site at 4:7432183 must be ignored');

    // P3: the §10.4 off-target 4:7436221-7436845 has a 4-edit L site (mm_pos [22,12,6,1]): unlikely, so P3 is specific
    pair.P3.on_target.should.match({ region: '4', start: 7422482, end: 7423061, size: 580 });
    should(find(pair.P3.off_targets, '4', 7436221, 7436845)).be.undefined();
    find(pair.P3.unlikely, '4', 7436221, 7436845).should.match({ likelihood: 'unlikely', left_mm: 4, right_mm: 1, right_mm_pos: [1], terminal_mismatch: true, approx: true });
    pair.P3.should.match({ verdict: 'specific', off_target_count: 0 });

    results.primers[P.P5_L].repetitive.should.be.true();
    results.primers[P.P5_L].near_perfect_sites.should.be.above(results.params.repeat_site_threshold);

    const m20 = pair.m20_7_14.off_targets.concat(pair.m20_7_14.on_target ? [pair.m20_7_14.on_target] : []);
    find(m20, '4', 7423537, 7423746, 'LR').should.match({ left_mm: 2, left_mm_pos: [14, 7], right_mm: 0 });
    const del = pair.del18.off_targets.concat(pair.del18.on_target ? [pair.del18.on_target] : []);
    find(del, '4', 7423537, 7423746, 'LR').should.match({ left_mm: 1, likelihood: 'likely' });
  });

  it('qPCR J_L + P1_R on SORBI_3004G087700 (transcript mode): genome verdict not on_target_missing; cDNA on-target .1/.2/.3 at 278 bp; 2 off-target gene groups', { timeout: JOB_TIMEOUT_MS + 60000 }, async function (t) {
    const { results, job } = await runCheck({
      system_name: 'sorghum_bicolor', mode: 'transcript', gene_id: 'SORBI_3004G087700',
      pairs: [{ id: 'J', left: P.J_L, right: P.P1_R }]
    }, t);
    job.request.transcript_id.should.match(/^SORBI_3004G087700\.\d+$/);
    const g = results.specificity.pairs[0];
    g.verdict.should.not.equal('on_target_missing');
    checkStringency(g, CAP);
    const c = results.transcriptome.pairs[0];
    c.on_target.gene_id.should.equal('SORBI_3004G087700');
    c.on_target.isoforms.map((i) => [i.transcript_id, i.size]).sort().should.eql([
      ['SORBI_3004G087700.1', 278], ['SORBI_3004G087700.2', 278], ['SORBI_3004G087700.3', 278]
    ]);
    c.off_target_count.should.equal(2);
    c.off_targets.map((x) => x.gene_id).sort().should.eql(['SORBI_3004G087800', 'SORBI_3005G183900']);
    c.off_targets.forEach(function (x) {
      x.left_mm.should.be.belowOrEqual(CAP);
      x.right_mm.should.be.belowOrEqual(CAP);
    });
    c.off_targets.find((x) => x.gene_id === 'SORBI_3005G183900').isoforms.map((i) => i.transcript_id).sort()
      .should.eql(['SORBI_3005G183900.1', 'SORBI_3005G183900.2']);
    t.diagnostic('genome verdict ' + g.verdict + ', off_target_count ' + g.off_target_count + ', gdna_products ' + g.gdna_products.length);
  });

  it('pan-genome P3 over sorghum_353, sorghum_grassl, sorghum_leoti: 580/580/581 bp primaries with terminal-mismatch flags, single_mismatch in all three', { timeout: JOB_TIMEOUT_MS + 60000 }, async function (t) {
    const { results, submitted, job } = await runCheck({
      system_name: 'sorghum_bicolor', mode: 'gene', gene_id: 'SORBI_3004G087700',
      checks: ['specificity', 'pangenome'], genomes: ['sorghum_353', 'sorghum_grassl', 'sorghum_leoti'],
      pairs: [{ id: 'P3', left: P.P3_L, right: P.P3_R, expected: { region: '4', start: 7422482, end: 7423061 } }]
    }, t);
    submitted.json.kind.should.equal('pangenome');
    job.progress.total.should.equal(4);
    const pair = results.pangenome.pairs[0];
    pair.reference_size.should.equal(580);
    const by = {};
    pair.genomes.forEach(function (x) { by[x.system_name] = x; });
    t.diagnostic('statuses ' + JSON.stringify(pair.genomes.map((x) => [x.system_name, x.status, x.other_amplicons])) + ' summary ' + JSON.stringify(pair.summary));
    by.sorghum_353.primary.should.match({ region: '4', start: 7499931, end: 7500510, size: 580, size_delta: 0, terminal_mismatch: true, likelihood: 'likely_weak', left_mm_pos: [1] });
    by.sorghum_grassl.primary.should.match({ size: 580, terminal_mismatch: true, likelihood: 'likely_weak' });
    by.sorghum_leoti.primary.should.match({ size: 581, size_delta: 1, right_mm: 1 });
    // the 626 bp (4-edit L site) and ~3.96 kb (5/4 edits) products no longer amplify, so no genome is 'multiple'
    pair.genomes.map((x) => [x.system_name, x.status]).sort().should.eql([
      ['sorghum_353', 'single_mismatch'], ['sorghum_grassl', 'single_mismatch'], ['sorghum_leoti', 'single_mismatch']
    ]);
    pair.summary.should.eql({ genomes_total: 3, single_perfect: 0, single_mismatch: 3, multiple: 0, no_amplicon: 0, db_unavailable: 0, error: 0, amplifies: 3, truncated: 0 });
    checkStringency(results.specificity.pairs[0], CAP);
    checkSummaries(results);
  });

  it('presence/absence SORBI_3001G046200 (default params): 353 and is12661 no_amplicon with ortholog_annotated false; ji2731 -177; leoti 0 mm', { timeout: JOB_TIMEOUT_MS + 60000 }, async function (t) {
    const { results } = await runCheck({
      system_name: 'sorghum_bicolor', mode: 'gene', gene_id: 'SORBI_3001G046200',
      checks: ['pangenome'], genomes: ['sorghum_353', 'sorghum_is12661', 'sorghum_ji2731', 'sorghum_leoti'],
      pairs: [{ id: 'PA', left: P.PA_L, right: P.PA_R }]
    }, t);
    const spec = results.specificity.pairs[0];
    spec.should.match({ on_target_inferred: true, on_target: { region: '1', size: 3621 } });
    checkStringency(spec, CAP);
    presenceOnTarget = spec.on_target;
    const pair = results.pangenome.pairs[0];
    pair.should.match({ reference_size: 3621, max_size: 5932 });
    t.diagnostic('statuses ' + JSON.stringify(pair.genomes.map((x) => [x.system_name, x.status, x.nearest && [x.nearest.left_mm, x.nearest.right_mm, x.nearest.approx]])));
    presenceStatuses(pair);
    checkSummaries(results);
  });

  it('expected product > 4 kb (default params): MAX_PRODUCT_SIZE_RAISED to 4346; on-target found; 353/is12661 no_amplicon, ji2731/leoti amplify', { timeout: JOB_TIMEOUT_MS + 60000 }, async function (t) {
    should.exist(presenceOnTarget, 'needs the on-target from the previous test');
    const expected = { region: presenceOnTarget.region, start: presenceOnTarget.start, end: presenceOnTarget.end };
    const body = {
      system_name: 'sorghum_bicolor', mode: 'gene', gene_id: 'SORBI_3001G046200', checks: ['pangenome'],
      genomes: ['sorghum_353', 'sorghum_is12661', 'sorghum_ji2731', 'sorghum_leoti'],
      pairs: [{ id: 'PA', left: P.PA_L, right: P.PA_R, expected: expected }]
    };
    const s = await submit(body);
    [200, 202].should.containEql(s.status, JSON.stringify(s.json).slice(0, 300));
    s.json.warnings.map((w) => w.code).should.containEql('MAX_PRODUCT_SIZE_RAISED');
    const { job } = await poll(s.json.job_id, t);
    job.status.should.equal('done', JSON.stringify(job.error));
    checkEngineAndParams(job);
    job.request.params.max_product_size.should.equal(Math.ceil(1.2 * 3621));
    const res = job.results;
    const spec = res.specificity.pairs[0];
    spec.verdict.should.not.equal('on_target_missing');
    spec.on_target.should.match(expected);
    const pair = res.pangenome.pairs[0];
    t.diagnostic('statuses ' + JSON.stringify(pair.genomes.map((x) => [x.system_name, x.status])));
    presenceStatuses(pair);
    pair.summary.amplifies.should.be.above(0);
    checkSummaries(res);
  });

  it('requests refused before queueing: validation 400s, INVALID_PARAMS, GENOME_NOT_CHECKABLE, UNKNOWN_*, DUPLICATE_PAIR_ID, SYSTEM_NAME_MISMATCH, PRODUCT_TOO_LONG_TO_CHECK; 404/400/405 on status', async function () {
    const pair = { id: 'P2', left: P.P2_L, right: P.P2_R };
    let r = await submit({ system_name: 'sorghum_bicolor', max_mismatches: 3, pairs: [pair] });
    r.status.should.equal(400);
    r.json.message.should.equal('Validation errors');
    r = await submit({ system_name: 'sorghum_bicolor', checks: ['transcriptome'], pairs: [pair] });
    r.status.should.equal(400);
    r.json.message.should.equal('Validation errors');
    r = await http('POST', '/primers/check', { json: undefined, headers: { 'Content-Type': 'text/plain' } });
    r.status.should.equal(400);

    // stringency cap: 0-5 in the contract, below ignore_mismatches in the handler
    r = await submit({ system_name: 'sorghum_bicolor', params: { max_amplifying_mismatches: 6 }, pairs: [pair] });
    r.status.should.equal(400);
    r.json.message.should.equal('Validation errors');
    r = await submit({ system_name: 'sorghum_bicolor', params: { ignore_mismatches: 3, max_amplifying_mismatches: 3 }, pairs: [pair] });
    r.status.should.equal(400);
    r.json.should.match({ code: 'INVALID_PARAMS', details: { field: 'params.max_amplifying_mismatches', max_amplifying_mismatches: 3, ignore_mismatches: 3 } });

    r = await submit({ system_name: 'sorghum_bicolor', checks: ['pangenome'], genomes: ['zea_maysb73'], pairs: [pair] });
    r.status.should.equal(400);
    r.json.code.should.equal('GENOME_NOT_CHECKABLE');
    r.json.details.genomes.should.containEql('zea_maysb73');

    r = await submit({ system_name: 'sorghum_bicolor', checks: ['pangenome'], genomes: ['sorghum_no_such'], pairs: [pair] });
    r.status.should.equal(404);
    r.json.code.should.equal('UNKNOWN_GENOME');

    r = await submit({ system_name: 'no_such_genome', pairs: [pair] });
    r.status.should.equal(404);
    r.json.code.should.equal('UNKNOWN_GENOME');

    r = await submit({ system_name: 'sorghum_bicolor', mode: 'gene', gene_id: 'NOPE', pairs: [pair] });
    r.status.should.equal(404);
    r.json.code.should.equal('UNKNOWN_GENE');

    r = await submit({ system_name: 'sorghum_rio', mode: 'gene', gene_id: 'SORBI_3004G087700', pairs: [pair] });
    r.status.should.equal(400);
    r.json.code.should.equal('SYSTEM_NAME_MISMATCH');

    r = await submit({ system_name: 'sorghum_bicolor', pairs: [pair, Object.assign({}, pair, { left: P.P1_L })] });
    r.status.should.equal(400);
    r.json.code.should.equal('DUPLICATE_PAIR_ID');

    r = await submit({ system_name: 'sorghum_bicolor', mode: 'region', pairs: [Object.assign({}, pair, { expected: { region: '4', start: 1, end: 10002 } })] });
    r.status.should.equal(400);
    r.json.code.should.equal('PRODUCT_TOO_LONG_TO_CHECK');

    r = await http('GET', '/primers/check/ffffffffffffffffffffffffffffffff');
    r.status.should.equal(404);
    r.json.code.should.equal('UNKNOWN_JOB');
    r.headers.get('cache-control').should.equal('no-store');
    r = await http('GET', '/primers/check/not-a-job-id');
    r.status.should.equal(400);
    r.json.message.should.equal('Validation errors');
    r = await http('PUT', '/primers/check/ffffffffffffffffffffffffffffffff', { json: {} });
    r.status.should.equal(405);
    r.json.code.should.equal('METHOD_NOT_ALLOWED');
  });

  it('cost guard: 10 pairs x 20 primers over every same-species genome -> 422 JOB_TOO_LARGE with the estimate',
    { skip: process.env.PRIMERS_IT_EXPECT_JOB_TOO_LARGE === '1' ? false : 'set PRIMERS_IT_EXPECT_JOB_TOO_LARGE=1 against a server started with a lower primers.check.max_job_cpu_s (see the file header)' },
    async function (t) {
      const bases = 'ACGT';
      // 20 distinct 22-mers: LCG high bits (the low bits of a power-of-two LCG cycle with period 4)
      const primer = function (i, j) {
        let s = '';
        let x = i * 7919 + j * 104729 + 13;
        for (let k = 0; k < 22; k++) { x = (x * 1103515245 + 12345) % 2147483648; s += bases[Math.floor(x / 65536) % 4]; }
        return s;
      };
      const pairs = Array.from({ length: 10 }, function (_, i) { return { id: 'C' + i, left: primer(i, 1), right: primer(i, 2) }; });
      new Set(pairs.map((p) => p.left).concat(pairs.map((p) => p.right))).size.should.equal(20);
      const r = await submit({ system_name: 'sorghum_bicolor', checks: ['specificity', 'pangenome'], pairs: pairs });
      r.status.should.equal(422, JSON.stringify(r.json).slice(0, 300));
      r.json.code.should.equal('JOB_TOO_LARGE');
      r.json.details.estimate_cpu_s.should.be.above(r.json.details.limit);
      // 20 primers x (reference + 119 assemblies at word size 6 + 120 genome re-alignment tasks x 0.6 CPU-s)
      r.json.details.estimate_cpu_s.should.be.within(8000, 10000); // 20 primers × all 119 genomes with pangenome_cpu_factor 2
      t.diagnostic('JOB_TOO_LARGE details ' + JSON.stringify(r.json.details));
    });
});
