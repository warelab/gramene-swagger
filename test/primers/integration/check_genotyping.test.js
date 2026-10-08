'use strict';

// Integration: a genotyping check job end to end against a running dev API AND its worker (genotyping spec §7.6 genotyping check
// rows, §2.11-§2.13, §5). The check.request of a live rs871475760 KASP design (num_sets 2), narrowed to the 11 research
// assemblies, is submitted, polled through its partial results, and its results.genotyping is compared with the real-data
// expectations of M8/M8b (test/primers/unit/genotype_realdata.test.js) and the verified §2.13 entries. Skipped unless
// PRIMERS_IT_BASE is set and PRIMERS_IT_WORKER=1:
//
//   HOST=127.0.0.1 PORT=50112 SWAGGER_HOST=localhost:50112 SWAGGER_SCHEMES=http PRIMERS_SITE_KEY=sorghum_v11_geno \
//     PRIMERS_GLOBAL_MAX_JOBS=1 node app.js
//   PRIMERS_SITE_KEY=sorghum_v11_geno PRIMERS_GLOBAL_MAX_JOBS=1 node api/helpers/primers/jobs/worker_main.js
//   PRIMERS_IT_BASE=http://127.0.0.1:50112/sorghum_v11 PRIMERS_IT_WORKER=1 node --test --test-concurrency=1 \
//     test/primers/integration/check_genotyping.test.js
//
// The job must be new: the first POST has to answer 202 and partial results have to be seen. Job ids are deterministic and
// done jobs are kept for ttl_done_s (24 h), so a re-run first needs the site's keys deleted:
//   redis-cli -p 6380 -n 1 --scan --pattern 'primers:sorghum_v11_geno:*' | xargs -r redis-cli -p 6380 -n 1 DEL
// The worker runs real BLAST (the reference at word size 5, the 11 assemblies at word size 6) plus the allele caller; the
// design call goes to live Ensembl REST (release 115). The refusal bodies ask for specificity only, so a regression that let one
// through would queue only a small job.
//
// Optional:
//   PRIMERS_IT_REDIS_URL (e.g. redis://127.0.0.1:6380/1) with PRIMERS_IT_SITE_KEY (default sorghum_v11_geno, the dev API's key):
//     the same pairs without genotyping are submitted only to compare job ids. With these set, that twin job is removed while
//     it is still queued behind the genotyping job (its job doc and its queue entry, and only while its status is queued), so
//     the worker never runs it. Without them it stays queued, and the worker runs it after the genotyping job unless the site's
//     keys are deleted first.
//   PRIMERS_IT_CAPTURE_DIR: write capture-check-genotyping-submit.json (the 202) and capture-check-genotyping-result.json (the
//     final GET) there, in the test/primers/fixtures/docs capture format; PRIMERS_IT_CAPTURE_SOURCE replaces the source note.
//   PRIMERS_IT_JOB_TIMEOUT_MS (default 20 min).

const { describe, it } = require('node:test');
const should = require('should');
const fs = require('fs');
const path = require('path');

const BASE = (process.env.PRIMERS_IT_BASE || '').replace(/\/+$/, '');
const SKIP = !BASE ? 'set PRIMERS_IT_BASE, e.g. http://127.0.0.1:50112/sorghum_v11'
  : process.env.PRIMERS_IT_WORKER !== '1' ? 'set PRIMERS_IT_WORKER=1 when the dev API has its own worker (the job runs real BLAST)' : false;
const JOB_TIMEOUT_MS = Number(process.env.PRIMERS_IT_JOB_TIMEOUT_MS) || 20 * 60 * 1000;
const POLL_MS = 1000;
const REDIS_URL = process.env.PRIMERS_IT_REDIS_URL || '';
const SITE_KEY = process.env.PRIMERS_IT_SITE_KEY || 'sorghum_v11_geno';
const CAPTURE_DIR = process.env.PRIMERS_IT_CAPTURE_DIR || '';

// §2.13 genotyping block: copies, calls and primer statuses verified against the real assembly FASTAs (sets S1 and S2).
const EXAMPLE = require('../fixtures/check_core/genotype/results_2_13.json');

const REFERENCE = 'sorghum_bicolor';
const DESIGN_BODY = Object.freeze({ system_name: REFERENCE, variant: { id: 'rs871475760', alt: 'A' }, assay: { type: 'kasp', num_sets: 2 } });

// §7.5 / M8 real data: the 11 research assemblies, [allele, orthologous copies]
const TABLE = Object.freeze({
  sorghum_bicolorv5: ['ref', 1], sorghum_austrcf317961: ['ref', 1], sorghum_pi565121: ['ref', 1], sorghum_pi655972: ['ref', 1],
  sorghum_pi180348: ['alt', 2], sorghum_pi276837: ['alt', 2], sorghum_pi656027: ['alt', 2], sorghum_pi510757: ['alt', 1],
  sorghum_rio: ['alt', 2], sorghum_pi329250: ['ref', 1], sorghum_s3691: ['ref', 1]
});
const GENOMES = Object.keys(TABLE);

// §2.9 check.request of the design above
const CHECK_REQUEST_2_9 = Object.freeze({
  system_name: REFERENCE,
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
});

const VARIANT_INFO = Object.freeze({
  key: '1:11109:C:A', region: '1', position: 11109, ref: 'C', alt: 'A', shift: 0, zone: { start: 11109, end: 11109 }, flank: 15,
  core: { ref: 'TCT', alt: 'TAT' },
  haplotypes: { ref: 'ATAGTCATACTCTATTCTGAATTTCTCGCTAGT', alt: 'ATAGTCATACTCTATTATGAATTTCTCGCTAGT' }
});

const ALLELES = ['ref', 'alt', 'other', 'ambiguous', 'missing', 'unavailable'];
const PREDICTED = Object.freeze({ ref: 'predicted_ref', alt: 'predicted_alt', both: 'both', none: 'none', no_call: 'no_call', unknown: 'unknown' });

async function http(method, p, opts) {
  opts = opts || {};
  const headers = Object.assign({}, opts.headers || {});
  let body;
  if (opts.json !== undefined) {
    body = JSON.stringify(opts.json);
    headers['Content-Type'] = 'application/json';
  }
  const res = await fetch(BASE + p, { method: method, headers: headers, body: body, redirect: 'manual' });
  const text = await res.text();
  let json = null;
  try { json = JSON.parse(text); } catch (e) { /* not JSON */ }
  return { status: res.status, headers: res.headers, json: json, text: text };
}

const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));
const clone = (x) => JSON.parse(JSON.stringify(x));
const codes = (warnings) => (warnings || []).map((w) => w.code);

function byName(list) {
  const m = {};
  (list || []).forEach(function (x) { m[x.system_name] = x; });
  return m;
}

async function submit(body) {
  const r = await http('POST', '/primers/check', { json: body });
  should(r.headers.get('cache-control')).equal('no-store');
  return r;
}

function expectError(r, status, code) {
  r.status.should.equal(status, r.text.slice(0, 500));
  should(r.headers.get('cache-control')).equal('no-store');
  should.exist(r.json, 'error bodies are JSON');
  r.json.code.should.equal(code, JSON.stringify(r.json).slice(0, 500));
  return r.json;
}

function basePath() {
  return new URL(BASE).pathname.replace(/\/+$/, '');
}

// The test/primers/fixtures/docs capture format: {source, request {method, path, body?}, status, response}.
function writeCapture(name, request, status, response) {
  if (!CAPTURE_DIR) return null;
  const source = process.env.PRIMERS_IT_CAPTURE_SOURCE ||
    ('recorded ' + new Date().toISOString().slice(0, 10) + ' over HTTP from the dev API and its worker on ' + new URL(BASE).host +
      ' by test/primers/integration/check_genotyping.test.js');
  const text = '{\n  "source": ' + JSON.stringify(source) + ',\n  "request": ' + JSON.stringify(request) + ',\n  "status": ' + status +
    ',\n  "response": ' + JSON.stringify(response) + '\n}\n';
  const file = path.join(CAPTURE_DIR, 'capture-' + name + '.json');
  fs.writeFileSync(file, text);
  return { file: file, bytes: Buffer.byteLength(text) };
}

// Removes a job that is still queued: its job doc and its queue entry, atomically and only while the doc says queued.
// -> 0 no such doc under this site key, -1 not queued (left alone), 1 + the queue entries removed otherwise.
const DROP_QUEUED_LUA = [
  "local raw = redis.call('GET', KEYS[1])",
  'if not raw then return 0 end',
  'local ok, doc = pcall(cjson.decode, raw)',
  "if not ok or type(doc) ~= 'table' or doc.status ~= 'queued' then return -1 end",
  "local removed = redis.call('LREM', KEYS[2], 0, ARGV[1])",
  "redis.call('DEL', KEYS[1])",
  'return 1 + removed'
].join('\n');

async function dropQueuedJob(id, kind) {
  const Redis = require('ioredis');
  const client = new Redis(REDIS_URL, { lazyConnect: true, enableOfflineQueue: false, maxRetriesPerRequest: 1, connectTimeout: 2000 });
  await client.connect();
  try {
    const P = 'primers:' + SITE_KEY + ':';
    return await client.eval(DROP_QUEUED_LUA, 2, P + 'job:' + id, P + 'queue:' + (kind === 'pangenome' ? 'pan' : 'spec'), id);
  } finally {
    client.disconnect();
  }
}

// §5.8: the genome summary counts the pan-genome genomes only and adds up; each set lists the same genomes in the same order,
// and its summary adds up both ways and matches its rows. Holds for partial results too.
function checkGenotypingSummaries(g, label) {
  g.genomes.length.should.be.above(0, label);
  g.genomes[0].should.match({ system_name: REFERENCE, is_reference: true }, label);
  g.genomes.filter((x) => x.is_reference).length.should.equal(1, label);
  const pan = g.genomes.slice(1);
  const s = g.summary;
  s.genomes_total.should.equal(pan.length, label + ': summary.genomes_total');
  ALLELES.reduce((n, a) => n + s[a], 0).should.equal(s.genomes_total, label + ': summary parts');
  ALLELES.forEach((a) => s[a].should.equal(pan.filter((x) => x.allele === a).length, label + ': summary.' + a));
  g.sets.forEach(function (set) {
    const ss = set.summary;
    const where = label + ' ' + set.id;
    set.reference.system_name.should.equal(REFERENCE, where);
    set.genomes.map((x) => x.system_name).should.eql(pan.map((x) => x.system_name), where + ': genomes');
    ss.genomes_total.should.equal(s.genomes_total, where + ': genomes_total');
    (ss.predicted_ref + ss.predicted_alt + ss.both + ss.none + ss.no_call + ss.unknown).should.equal(ss.genomes_total, where + ': predictions');
    (ss.agree + ss.disagree + ss.not_comparable).should.equal(ss.genomes_total, where + ': agreement');
    Object.keys(PREDICTED).forEach((p) => ss[PREDICTED[p]].should.equal(set.genomes.filter((x) => x.predicted === p).length, where + ': ' + p));
    ss.agree.should.equal(set.genomes.filter((x) => x.agrees === true).length, where + ': agree');
    ss.disagree.should.equal(set.genomes.filter((x) => x.agrees === false).length, where + ': disagree');
    ss.not_comparable.should.equal(set.genomes.filter((x) => x.agrees === null).length, where + ': not_comparable');
    ss.weak.should.equal(set.genomes.filter((x) => x.strength === 'weak').length, where + ': weak');
  });
}

function checkPangenomeSummaries(results) {
  results.pangenome.pairs.forEach(function (pair) {
    const s = pair.summary;
    (s.single_perfect + s.single_mismatch + s.multiple + s.no_amplicon + s.db_unavailable + s.error).should.equal(s.genomes_total, pair.id);
    s.amplifies.should.equal(s.single_perfect + s.single_mismatch + s.multiple, pair.id);
    s.genomes_total.should.equal(pair.genomes.length, pair.id);
  });
}

// Polls every second until done or error. Partial results.genotyping blocks seen while running are checked for consistency.
async function poll(jobId) {
  const started = Date.now();
  const seen = [];
  const partials = [];
  const problems = [];
  const done = [];
  let last = null;
  for (;;) {
    const r = await http('GET', '/primers/check/' + jobId);
    r.status.should.equal(200, r.text.slice(0, 300));
    should(r.headers.get('cache-control')).equal('no-store');
    last = r;
    const job = r.json;
    if (seen[seen.length - 1] !== job.status) seen.push(job.status);
    if (done[done.length - 1] !== job.progress.done) done.push(job.progress.done);
    if (job.status === 'running' && job.partial === true && job.results && job.results.genotyping) {
      const g = job.results.genotyping;
      const at = Date.now() - started;
      partials.push({ at_ms: at, done: job.progress.done, total: job.progress.total, genomes: g.genomes.map((x) => x.system_name) });
      try {
        checkGenotypingSummaries(g, 'partial at ' + at + ' ms');
      } catch (e) {
        problems.push(e.message);
      }
    }
    if (job.status === 'done' || job.status === 'error') break;
    if (Date.now() - started > JOB_TIMEOUT_MS) throw new Error('job ' + jobId + ' still ' + job.status + ' after ' + JOB_TIMEOUT_MS + ' ms');
    await sleep(POLL_MS);
  }
  return { response: last, job: last.json, seen: seen, partials: partials, problems: problems, done: done, wall_ms: Date.now() - started };
}

describe('genotyping check job: POST /primers/check with genotyping and GET /primers/check/{job_id} (dev server + worker)', { skip: SKIP }, function () {
  const state = { body: null, jobId: null, estimate: null, twinId: null, twinDropped: false, job: null };

  it('rs871475760 KASP design → check.request over the 11 research assemblies: POST 202, re-POST 200 with the same id; the same pairs without genotyping get another id', async function (t) {
    const d = await http('POST', '/primers/genotyping/design', { json: DESIGN_BODY });
    d.status.should.equal(200, d.text.slice(0, 500));
    d.json.sets.map((s) => s.key).should.eql(['f9df650ad116', '1accc54c262d']);
    const cr = d.json.check.request;
    cr.should.eql(CHECK_REQUEST_2_9);

    const catalog = await http('GET', '/primers/genomes?system_name=' + REFERENCE);
    catalog.status.should.equal(200, catalog.text.slice(0, 300));
    const names = new Set(catalog.json.genomes.map((x) => x.system_name));
    GENOMES.forEach((sys) => names.has(sys).should.be.true(sys + ' is not a system_name of /primers/genomes'));

    const body = { system_name: cr.system_name, mode: cr.mode, checks: cr.checks, genomes: GENOMES.slice(), pairs: cr.pairs, genotyping: cr.genotyping };
    state.body = body;
    const first = await submit(body);
    first.status.should.equal(202, 'the first POST must create the job (delete the primers:' + SITE_KEY + ':* keys of an earlier run): ' + first.text.slice(0, 400));
    Object.keys(first.json).sort().should.eql(['created_at', 'estimate', 'job_id', 'kind', 'progress', 'queue_position', 'status', 'warnings']);
    first.json.job_id.should.match(/^[0-9a-f]{32}$/);
    first.json.should.match({ status: 'queued', kind: 'pangenome', progress: { done: 0, total: 1 + GENOMES.length, stage: 'queued', running: [] }, warnings: [] });
    first.json.estimate.cpu_s.should.be.above(0);
    state.jobId = first.json.job_id;
    state.estimate = first.json.estimate;
    const cap = writeCapture('check-genotyping-submit', { method: 'POST', path: basePath() + '/primers/check', body: clone(body) }, first.status, first.json);
    t.diagnostic('POST ' + first.status + ' ' + JSON.stringify(first.json) + (cap ? '; wrote ' + cap.file + ' (' + cap.bytes + ' bytes)' : ''));

    const again = await submit(clone(body));
    again.status.should.equal(200, again.text.slice(0, 300));
    again.json.job_id.should.equal(state.jobId);
    ['queued', 'running'].should.containEql(again.json.status);
    // §5.3: alleles are uppercased before hashing, so their case does not change the id
    const lower = clone(body);
    lower.genotyping.variant.ref = 'c';
    lower.genotyping.variant.alt = 'a';
    const cased = await submit(lower);
    cased.status.should.equal(200, cased.text.slice(0, 300));
    cased.json.job_id.should.equal(state.jobId);

    // §5.3: the same request without the genotyping block hashes with algorithm '2', not '2+g1'
    const twinBody = clone(body);
    delete twinBody.genotyping;
    const twin = await submit(twinBody);
    [200, 202].should.containEql(twin.status, twin.text.slice(0, 300));
    twin.json.kind.should.equal('pangenome');
    twin.json.job_id.should.match(/^[0-9a-f]{32}$/);
    twin.json.job_id.should.not.equal(state.jobId);
    // §5.4: the genotyping term only adds to the estimate
    twin.json.estimate.cpu_s.should.be.belowOrEqual(state.estimate.cpu_s);
    state.twinId = twin.json.job_id;
    t.diagnostic('without genotyping: POST ' + twin.status + ', job ' + twin.json.job_id + ' (status ' + twin.json.status + '), estimate ' +
      JSON.stringify(twin.json.estimate) + ' vs ' + JSON.stringify(state.estimate) + ' with genotyping');
    if (REDIS_URL) {
      const dropped = await dropQueuedJob(state.twinId, twin.json.kind);
      state.twinDropped = dropped >= 1;
      t.diagnostic('twin job ' + state.twinId + ' removed while queued: ' + dropped);
      if (twin.status === 202) dropped.should.be.aboveOrEqual(1, 'the twin job was expected to be queued behind the genotyping job under site key ' + SITE_KEY);
    } else {
      t.diagnostic('PRIMERS_IT_REDIS_URL is not set: the twin job ' + state.twinId + ' stays queued; delete the site\'s keys after the run');
    }
  });

  it('the worker runs it: running with partial results.genotyping that grow genome by genome and stay consistent, then done', { timeout: JOB_TIMEOUT_MS + 60000 }, async function (t) {
    should.exist(state.jobId, 'needs the job submitted by the previous test');
    const p = await poll(state.jobId);
    const job = p.job;
    state.job = job;
    t.diagnostic('statuses ' + p.seen.join(' -> ') + ' in ' + Math.round(p.wall_ms / 1000) + ' s; progress.done ' + p.done.join(',') + '; attempts ' + job.attempts);
    t.diagnostic('created ' + job.created_at + ', started ' + job.started_at + ', finished ' + job.finished_at + '; estimate ' + JSON.stringify(job.estimate) +
      '; timings_ms ' + JSON.stringify(job.results && job.results.timings_ms));
    const sizes = [];
    p.partials.forEach((x) => { if (sizes[sizes.length - 1] !== x.genomes.length) sizes.push(x.genomes.length); });
    t.diagnostic('partial results.genotyping.genomes lengths ' + sizes.join(',') + ' over ' + p.partials.length + ' polls');
    job.status.should.equal('done', JSON.stringify(job.error));
    job.partial.should.be.false();
    job.attempts.should.be.aboveOrEqual(1);
    p.seen.should.containEql('running');
    p.partials.length.should.be.above(0, 'no partial results.genotyping was seen while the job ran (statuses ' + p.seen.join(' -> ') + ')');
    p.partials.some((x) => x.genomes.length < 1 + GENOMES.length && x.done < x.total).should.be.true('no intermediate partial result: ' + sizes.join(','));
    p.partials.forEach(function (x) {
      x.genomes[0].should.equal(REFERENCE);
      x.genomes.slice(1).forEach((sys) => GENOMES.should.containEql(sys));
    });
    for (let i = 1; i < p.partials.length; i++) p.partials[i].genomes.length.should.be.aboveOrEqual(p.partials[i - 1].genomes.length, 'partial genomes never shrink');
    p.problems.should.eql([], 'partial results.genotyping inconsistencies');
    const cap = writeCapture('check-genotyping-result', { method: 'GET', path: basePath() + '/primers/check/' + state.jobId }, p.response.status, job);
    if (cap) t.diagnostic('wrote ' + cap.file + ' (' + cap.bytes + ' bytes)');
  });

  it('results.genotyping: 6 ref and 5 alt; S1 agrees 11/11, S2 10/11 with both on pi180348; reference controls pass; one WEAK_OFF_TARGETS; consistent summaries', function (t) {
    const job = state.job;
    should.exist(job, 'needs the finished job of the previous test');
    const results = job.results;
    const g = results.genotyping;
    t.diagnostic('calls ' + JSON.stringify(g.genomes.map((x) => [x.system_name, x.allele, x.observed, x.source, x.orthologous_copies, x.paralog_copies])));
    t.diagnostic('predictions ' + JSON.stringify(g.sets.map((s) => [s.id, s.control, s.specificity, s.summary, s.genomes.filter((x) => x.agrees !== true).map((x) => [x.system_name, x.predicted, x.reasons])])));
    t.diagnostic('results.warnings ' + JSON.stringify(results.warnings));

    // the stored request keeps only the client's keys; engine and caller versions
    job.warnings.should.eql([]);
    job.request.genotyping.should.eql(CHECK_REQUEST_2_9.genotyping);
    job.request.pairs.should.eql(CHECK_REQUEST_2_9.pairs);
    job.request.genomes.should.eql(GENOMES.slice().sort());
    job.progress.should.match({ done: 1 + GENOMES.length, total: 1 + GENOMES.length });
    results.engine.algorithm_version.should.equal('2');
    g.algorithm_version.should.equal('g1');
    g.variant.should.match(VARIANT_INFO);

    // allele calls: the reference, then the pan-genome genomes in request order
    g.genomes.map((x) => x.system_name).should.eql([REFERENCE].concat(job.request.genomes));
    const by = byName(g.genomes);
    by[REFERENCE].should.match({ is_reference: true, allele: 'ref', observed: 'TCT', source: 'amplicon', orthologous_copies: 1, reason: null });
    GENOMES.forEach(function (sys) {
      should([by[sys].is_reference, by[sys].allele, by[sys].orthologous_copies, by[sys].source, by[sys].observed, by[sys].reason])
        .eql([false, TABLE[sys][0], TABLE[sys][1], 'amplicon', TABLE[sys][0] === 'ref' ? 'TCT' : 'TAT', null], sys);
    });
    g.summary.should.eql({ genomes_total: 11, ref: 6, alt: 5, other: 0, ambiguous: 0, missing: 0, unavailable: 0 });
    // the verified §2.13 entries (reference, bicolorv5, pi180348, pi329250), except that pi180348 has one paralog of the locus at
    // 1:72.9 Mb, which the megablast loci behind the example did not include
    const PARALOGS = { sorghum_bicolor: 0, sorghum_bicolorv5: 0, sorghum_pi180348: 1, sorghum_pi329250: 0 };
    EXAMPLE.genomes.forEach(function (e) {
      should(by[e.system_name]).eql(Object.assign({}, e, { paralog_copies: PARALOGS[e.system_name] }), e.system_name);
    });
    by.sorghum_pi180348.copies.map((k) => [k.start, k.end, k.strand]).should.eql([[15028, 15132, -1], [38342, 38446, 1]]);

    // sets: derived orientation, no deliberate mismatch, §2.13 primer calls
    g.sets.map((s) => [s.id, s.ref_pair, s.alt_pair, s.orientation, s.deliberate_mismatch_positions]).should.eql([
      ['S1', 'S1_REF', 'S1_ALT', 'reverse', []], ['S2', 'S2_REF', 'S2_ALT', 'forward', []]
    ]);
    const primerCalls = (x) => [x.system_name, x.ref_primer, x.alt_primer, x.common_primer];
    const fullCalls = (x) => primerCalls(x).concat([x.predicted, x.strength, x.agrees, x.reasons, x.off_locus_products]);
    EXAMPLE.sets.forEach(function (s) {
      const got = g.sets.find((x) => x.id === s.id);
      should(fullCalls(got.reference)).eql(fullCalls(s.reference), s.id + ' reference');
      const mine = byName(got.genomes);
      s.genomes.forEach(function (e) {
        const calls = s.id === 'S2' && e.system_name === 'sorghum_pi180348' ? primerCalls : fullCalls;
        should(calls(mine[e.system_name])).eql(calls(e), s.id + ' ' + e.system_name);
      });
    });

    // predictions (§5.7 with the M8b off-locus threshold): every genome reads its own allele, except S2 on pi180348, whose REF pair
    // amplifies the 1:72.9 Mb paralog carrying C (1 and 2 mismatches): both
    g.sets.forEach(function (s) {
      const rows = byName(s.genomes);
      GENOMES.forEach(function (sys) {
        const paralog = s.id === 'S2' && sys === 'sorghum_pi180348';
        should([rows[sys].predicted, rows[sys].agrees, rows[sys].reasons, rows[sys].off_locus_products])
          .eql(paralog ? ['both', false, ['ref_signal_off_locus'], 1] : [TABLE[sys][0], true, [], 0], s.id + ' ' + sys);
      });
      // the reference control
      should([s.reference.predicted, s.control]).eql(['ref', { status: 'pass', allele: 'ref', reasons: [] }], s.id + ' control');
    });
    g.sets.map((s) => [s.id, s.summary]).should.eql([
      ['S1', { genomes_total: 11, predicted_ref: 6, predicted_alt: 5, both: 0, none: 0, no_call: 0, unknown: 0, weak: 0, agree: 11, disagree: 0, not_comparable: 0 }],
      ['S2', { genomes_total: 11, predicted_ref: 6, predicted_alt: 4, both: 1, none: 0, no_call: 0, unknown: 0, weak: 0, agree: 10, disagree: 1, not_comparable: 0 }]
    ]);

    // specificity is unchanged by genotyping: S1_REF keeps its 3 weak genome-wide off-targets, which the sets report as is
    const spec = {};
    results.specificity.pairs.forEach((p) => { spec[p.id] = p; });
    ['S1_REF', 'S1_ALT', 'S2_REF', 'S2_ALT'].map((id) => [id, spec[id].verdict, spec[id].off_target_count]).should.eql([
      ['S1_REF', 'off_targets', 3], ['S1_ALT', 'specific', 0], ['S2_REF', 'specific', 0], ['S2_ALT', 'specific', 0]
    ]);
    g.sets.forEach(function (s) {
      s.specificity.ref_pair.verdict.should.equal(spec[s.ref_pair].verdict, s.id);
      s.specificity.alt_pair.verdict.should.equal(spec[s.alt_pair].verdict, s.id);
    });
    g.sets.map((s) => [s.id, s.specificity.off_target_count]).should.eql([['S1', 3], ['S2', 0]]);

    // one job-level WEAK_OFF_TARGETS: S1_REF's off-targets carry 3 mismatches in a primer, over genotype_offlocus_max_mismatches
    const weak = results.warnings.filter((w) => w.code === 'WEAK_OFF_TARGETS');
    weak.should.have.length(1);
    weak[0].details.should.match({ max_mismatches: 2 });
    weak[0].details.count.should.be.above(0);
    weak[0].details.examples.length.should.be.within(1, 5);
    weak[0].details.examples.forEach(function (e) {
      e.pair_id.should.equal('S1_REF');
      Math.max(e.left_mm, e.right_mm).should.be.above(2);
    });
    weak[0].details.examples.slice(0, 3).map((e) => [e.system_name, e.pair_id, e.region + ':' + e.start + '-' + e.end, e.size, e.left_mm, e.right_mm]).should.eql([
      ['sorghum_bicolor', 'S1_REF', '1:372590-372654', 65, 3, 3],
      ['sorghum_bicolor', 'S1_REF', '4:46360759-46360823', 65, 3, 3],
      ['sorghum_bicolor', 'S1_REF', '9:14931544-14931612', 69, 3, 3]
    ]);
    codes(results.warnings).filter((c) => /^(GENOTYPE_|REFERENCE_CONTROL_FAILED)/.test(c)).should.eql([]);

    // summaries: genotyping (§5.8) and pan-genome, which count the same genomes in the same order
    checkGenotypingSummaries(g, 'final');
    checkPangenomeSummaries(results);
    results.pangenome.pairs.map((p) => p.id).should.eql(['S1_REF', 'S1_ALT', 'S2_REF', 'S2_ALT']);
    results.pangenome.pairs.forEach(function (p) {
      p.summary.genomes_total.should.equal(g.summary.genomes_total, p.id);
      p.genomes.map((x) => x.system_name).should.eql(g.genomes.slice(1).map((x) => x.system_name), p.id);
    });
  });

  it('refusals: GET of an unknown job 404 UNKNOWN_JOB; GENOTYPING_SET_INVALID alleles_swapped and common_in_zone (§5.2) answer 400', async function () {
    const unknown = 'ffffffffffffffffffffffffffffffff';
    expectError(await http('GET', '/primers/check/' + unknown), 404, 'UNKNOWN_JOB').details.should.eql({ job_id: unknown });
    if (state.twinDropped) expectError(await http('GET', '/primers/check/' + state.twinId), 404, 'UNKNOWN_JOB');

    // ref_pair and alt_pair exchanged: the REF-specific primer ends in the ALT base
    const swapped = clone(CHECK_REQUEST_2_9);
    swapped.checks = ['specificity'];
    swapped.genotyping.sets[0] = { id: 'S1', ref_pair: 'S1_ALT', alt_pair: 'S1_REF' };
    let e = expectError(await submit(swapped), 400, 'GENOTYPING_SET_INVALID');
    e.details.should.eql({ set_id: 'S1', reason: 'alleles_swapped', pair_ids: ['S1_ALT'], primer: 'ref_pair', sequence: 'ATCTTTGACTAGCGAGAAATTCAT' });
    e.message.should.match(/^set S1: .*are ref_pair and alt_pair swapped\?$/);

    // a forward set whose common primer (1:11108-11131) lies over the variant
    const expected = { region: '1', start: 11081, end: 11131 };
    const common = {
      system_name: REFERENCE,
      mode: 'region',
      checks: ['specificity'],
      pairs: [
        { id: 'C_REF', left: 'GGTTATCCGAATATAGTCATACTCTATTC', right: 'TCTTTGACTAGCGAGAAATTCAGA', expected: expected },
        { id: 'C_ALT', left: 'GGTTATCCGAATATAGTCATACTCTATTA', right: 'TCTTTGACTAGCGAGAAATTCAGA', expected: expected }
      ],
      genotyping: { variant: clone(CHECK_REQUEST_2_9.genotyping.variant), sets: [{ id: 'C', ref_pair: 'C_REF', alt_pair: 'C_ALT' }] }
    };
    e = expectError(await submit(common), 400, 'GENOTYPING_SET_INVALID');
    e.details.should.eql({ set_id: 'C', reason: 'common_in_zone', pair_ids: ['C_REF', 'C_ALT'], primer: 'common', sequence: 'TCTTTGACTAGCGAGAAATTCAGA' });
    e.message.should.match(/^set C: the common primer TCTTTGACTAGCGAGAAATTCAGA at 1:11108-11131 must lie after the variant zone 1:11109-11109/);
  });
});
