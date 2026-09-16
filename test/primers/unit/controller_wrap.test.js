'use strict';

// api/controllers/primers.js: every handler is promise-wrapped (a throw or rejection becomes a JSON error,
// never an unhandled rejection), every response is Cache-Control: no-store, and the /primers middleware
// used by app.js answers JSON 404/405. Driven through a real express app on an ephemeral 127.0.0.1 port
// with fake helper modules; app.js is never required.

require('../../../api/helpers/primers/node_compat');

const test = require('node:test');
const should = require('should');
const express = require('express');

const CONTROLLER = '../../../api/controllers/primers';
const { PrimerHttpError } = require('../../../api/helpers/primers/errors');

function quietLog() {
  const lines = [];
  const push = function (m) { lines.push(String(m)); };
  return { lines: lines, log: push, info: push, warn: push, error: push };
}

function fakeConfig(overrides) {
  const cfg = Object.assign({ enabled: true, check: { store: 'memory', global_prefix: 'primers:test:' } }, overrides || {});
  return {
    basePathSet: null,
    get: function () { return cfg; },
    setBasePath: function (bp) { this.basePathSet = bp; },
    basePath: function () { return this.basePathSet || '/sorghum_v11'; },
    siteKey: function () { return 'sorghum_v11:test'; }
  };
}

// An express app wired like app.js: noStore + bare-path 404 before the routes, 404 + error handler after.
// Routes stand in for the swagger router: they put req.swagger.params in place and call the handler the way
// swagger_router does (req, res, cb).
function buildApp(ctrl, extra) {
  const app = express();
  const base = '/sorghum_v11/primers';
  app.use(base, ctrl.noStore);
  app.all(base, ctrl.notFound);
  app.use(express.json());
  const swagger = function (params) {
    return function (req, res, next) {
      req.swagger = { params: {} };
      Object.keys(params).forEach(function (k) { req.swagger.params[k] = { value: params[k](req) }; });
      next();
    };
  };
  app.post(base + '/design', swagger({ body: (r) => r.body }), function (req, res) { ctrl.designPrimers(req, res, function () {}); });
  app.get(base + '/genomes', swagger({ system_name: (r) => r.query.system_name }), function (req, res) { ctrl.primerGenomes(req, res, function () {}); });
  app.post(base + '/check', swagger({ body: (r) => r.body }), function (req, res) { ctrl.submitPrimerCheck(req, res, function () {}); });
  app.get(base + '/check/:job_id', swagger({ job_id: (r) => r.params.job_id }), function (req, res) { ctrl.getPrimerCheck(req, res, function () {}); });
  // swagger-node-runner's 405 for a path that is defined but lacks the method
  app.all(base + '/check/:job_id', function (req, res, next) {
    const err = new Error('Path [/primers/check/{job_id}] defined in Swagger, but ' + req.method + ' operation is not.');
    err.statusCode = 405;
    err.allowedMethods = ['GET'];
    res.setHeader('Allow', 'GET');
    next(err);
  });
  if (extra) extra(app, base);
  app.use(base, ctrl.notFound);
  app.use(base, ctrl.errorHandler);
  return app;
}

async function withServer(app, fn) {
  const server = await new Promise(function (resolve) {
    const s = app.listen(0, '127.0.0.1', function () { resolve(s); });
  });
  const baseUrl = 'http://127.0.0.1:' + server.address().port + '/sorghum_v11';
  try {
    return await fn(baseUrl);
  } finally {
    await new Promise(function (resolve) { server.close(resolve); if (server.closeAllConnections) server.closeAllConnections(); });
  }
}

async function call(url, opts) {
  const res = await fetch(url, opts);
  const text = await res.text();
  let body = null;
  try { body = JSON.parse(text); } catch (e) { body = text; }
  return { status: res.status, headers: res.headers, body: body };
}

function postJson(url, body) {
  return call(url, { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body) });
}

// Records unhandled rejections for the duration of a test.
function watchUnhandled(t) {
  const seen = [];
  const onRej = function (reason) { seen.push(reason); };
  process.on('unhandledRejection', onRej);
  t.after(function () { process.removeListener('unhandledRejection', onRej); });
  return seen;
}

const tick = () => new Promise((resolve) => setImmediate(resolve));

test('loading the controller exports the four operations and the app.js middleware without loading heavy helpers', function () {
  const ctrl = require(CONTROLLER);
  ['designPrimers', 'primerGenomes', 'submitPrimerCheck', 'getPrimerCheck', 'noStore', 'notFound', 'errorHandler', 'startup', 'createController']
    .forEach(function (k) { ctrl[k].should.be.a.Function(); });
  const loaded = Object.keys(require.cache);
  ['design.js', 'jobs/index.js', 'genomes.js', 'assemblies.js', 'redis_store.js'].forEach(function (f) {
    loaded.some(function (p) { return p.endsWith('/api/helpers/primers/' + f); }).should.equal(false, f + ' loaded eagerly');
  });
});

test('a rejected design promise becomes a 500 JSON INTERNAL error with no-store, never an unhandled rejection', async function (t) {
  const unhandled = watchUnhandled(t);
  const log = quietLog();
  const ctrl = require(CONTROLLER).createController({
    log: log,
    config: fakeConfig(),
    design: { design: async function () { throw new Error('boom at /scratch/olson/fasta/secret.fa.gz'); } }
  });
  await withServer(buildApp(ctrl), async function (base) {
    const r = await postJson(base + '/primers/design', { mode: 'gene', gene_id: 'X' });
    r.status.should.equal(500);
    r.body.should.eql({ message: 'internal error', code: 'INTERNAL', details: {} });
    r.headers.get('cache-control').should.equal('no-store');
    r.headers.get('content-type').should.match(/application\/json/);
    JSON.stringify(r.body).should.not.match(/scratch|secret/);
  });
  await tick();
  unhandled.should.have.length(0);
  log.lines.join('\n').should.match(/boom/); // the stack goes to the server log only
});

test('a synchronous throw inside a handler (before any await) is also a JSON error', async function (t) {
  const unhandled = watchUnhandled(t);
  const ctrl = require(CONTROLLER).createController({
    log: quietLog(),
    config: { get: function () { throw new TypeError('config exploded'); } }
  });
  await withServer(buildApp(ctrl), async function (base) {
    const r = await call(base + '/primers/genomes?system_name=sorghum_bicolor');
    r.status.should.equal(500);
    r.body.code.should.equal('INTERNAL');
    r.headers.get('cache-control').should.equal('no-store');
  });
  // called directly with a bare response object: the returned promise resolves, the error is sent
  const res = {
    headers: {}, statusCode: 200, headersSent: false, body: null,
    setHeader: function (k, v) { this.headers[k.toLowerCase()] = v; },
    set: function (k, v) { this.headers[k.toLowerCase()] = v; return this; },
    status: function (c) { this.statusCode = c; return this; },
    json: function (b) { this.body = b; this.headersSent = true; return this; }
  };
  const p = ctrl.designPrimers({}, res);
  (typeof p.then).should.equal('function');
  await p;
  res.statusCode.should.equal(500);
  res.body.code.should.equal('INTERNAL');
  await tick();
  unhandled.should.have.length(0);
});

test('PrimerHttpError status, code and details pass through; 503 adds Retry-After', async function () {
  const errs = {
    design: new PrimerHttpError(404, 'UNKNOWN_GENE', 'unknown gene NOPE', { gene_id: 'NOPE' }),
    check: new PrimerHttpError(503, 'JOB_STORE_UNAVAILABLE', 'the job store is unavailable', { retry_after_s: 5 }),
    status: new PrimerHttpError(404, 'UNKNOWN_JOB', 'unknown job', { job_id: '0123456789abcdef0123456789abcdef' })
  };
  const ctrl = require(CONTROLLER).createController({
    log: quietLog(),
    config: fakeConfig(),
    design: { design: async function () { throw errs.design; } },
    jobs: {
      submit: async function () { throw errs.check; },
      status: async function () { throw errs.status; }
    }
  });
  await withServer(buildApp(ctrl), async function (base) {
    let r = await postJson(base + '/primers/design', { mode: 'gene', gene_id: 'NOPE' });
    r.status.should.equal(404);
    r.body.should.eql({ message: 'unknown gene NOPE', code: 'UNKNOWN_GENE', details: { gene_id: 'NOPE' } });
    r = await postJson(base + '/primers/check', { system_name: 'sorghum_bicolor', pairs: [] });
    r.status.should.equal(503);
    r.body.code.should.equal('JOB_STORE_UNAVAILABLE');
    r.body.details.retry_after_s.should.equal(5);
    r.headers.get('retry-after').should.equal('5');
    r.headers.get('cache-control').should.equal('no-store');
    r = await call(base + '/primers/check/0123456789abcdef0123456789abcdef');
    r.status.should.equal(404);
    r.body.code.should.equal('UNKNOWN_JOB');
  });
});

test('success paths: design 200, genomes 200, check 202 when created and 200 when existing (created dropped), status 200', async function () {
  const seen = {};
  let created = true;
  const ctrl = require(CONTROLLER).createController({
    log: quietLog(),
    config: fakeConfig(),
    design: { design: async function (body, deps) { seen.design = { body: body, signal: deps.signal }; return { pairs: [], warnings: [] }; } },
    genomes: { genomesResponse: async function (sys) { seen.genomes = sys; return { system_name: sys, genomes: [] }; } },
    jobs: {
      submit: async function (body) {
        seen.submit = body;
        return { job_id: 'a'.repeat(32), status: 'queued', kind: 'specificity', queue_position: 0, created: created };
      },
      status: async function (id) { seen.status = id; return { job_id: id, status: 'done', results: {} }; }
    }
  });
  await withServer(buildApp(ctrl), async function (base) {
    let r = await postJson(base + '/primers/design', { mode: 'transcript', gene_id: 'SORBI_3001G000200' });
    r.status.should.equal(200);
    r.body.should.eql({ pairs: [], warnings: [] });
    r.headers.get('cache-control').should.equal('no-store');
    seen.design.body.should.eql({ mode: 'transcript', gene_id: 'SORBI_3001G000200' });
    seen.design.signal.should.be.instanceOf(AbortSignal);
    seen.design.signal.aborted.should.equal(false);

    r = await call(base + '/primers/genomes?system_name=sorghum_bicolor');
    r.status.should.equal(200);
    seen.genomes.should.equal('sorghum_bicolor');

    r = await postJson(base + '/primers/check', { system_name: 'sorghum_bicolor' });
    r.status.should.equal(202);
    r.body.should.not.have.property('created');
    r.body.job_id.should.equal('a'.repeat(32));
    created = false;
    r = await postJson(base + '/primers/check', { system_name: 'sorghum_bicolor' });
    r.status.should.equal(200);
    r.body.should.not.have.property('created');
    r.headers.get('cache-control').should.equal('no-store');

    r = await call(base + '/primers/check/' + 'b'.repeat(32));
    r.status.should.equal(200);
    seen.status.should.equal('b'.repeat(32));
  });
});

test('FEATURE_DISABLED (503, retry_after_s 300) on every endpoint when primers are disabled', async function () {
  const ctrl = require(CONTROLLER).createController({
    log: quietLog(),
    config: fakeConfig({ enabled: false }),
    design: { design: async function () { throw new Error('must not be called'); } },
    genomes: { genomesResponse: async function () { throw new Error('must not be called'); } },
    jobs: { submit: async function () { throw new Error('must not be called'); }, status: async function () { throw new Error('must not be called'); } }
  });
  await withServer(buildApp(ctrl), async function (base) {
    const rs = [
      await postJson(base + '/primers/design', {}),
      await call(base + '/primers/genomes?system_name=x'),
      await postJson(base + '/primers/check', {}),
      await call(base + '/primers/check/' + 'c'.repeat(32))
    ];
    rs.forEach(function (r) {
      r.status.should.equal(503);
      r.body.code.should.equal('FEATURE_DISABLED');
      r.body.details.retry_after_s.should.equal(300);
      r.headers.get('retry-after').should.equal('300');
    });
  });
});

test('a client that disconnects aborts the design signal and nothing is sent or thrown', async function (t) {
  const unhandled = watchUnhandled(t);
  const log = quietLog();
  let resolveStarted;
  const started = new Promise((resolve) => { resolveStarted = resolve; });
  let abortedWith = null;
  const ctrl = require(CONTROLLER).createController({
    log: log,
    config: fakeConfig(),
    design: {
      design: function (body, deps) {
        resolveStarted();
        return new Promise(function (resolve, reject) {
          deps.signal.addEventListener('abort', function () { abortedWith = deps.signal.reason; reject(deps.signal.reason); });
        });
      }
    }
  });
  await withServer(buildApp(ctrl), async function (base) {
    const ac = new AbortController();
    const pending = fetch(base + '/primers/design', {
      method: 'POST', headers: { 'Content-Type': 'application/json' }, body: '{"mode":"gene","gene_id":"X"}', signal: ac.signal
    }).catch(function (e) { return e; });
    await started;
    ac.abort();
    (await pending).name.should.equal('AbortError');
    for (let i = 0; i < 50 && !abortedWith; i++) await new Promise((resolve) => setTimeout(resolve, 10));
  });
  should.exist(abortedWith);
  abortedWith.code.should.equal('CLIENT_CLOSED_REQUEST');
  await tick();
  unhandled.should.have.length(0);
});

test('app.js middleware: bare /primers and unknown /primers paths are JSON 404, wrong method is JSON 405, all no-store', async function () {
  const ctrl = require(CONTROLLER).createController({ log: quietLog(), config: fakeConfig() });
  await withServer(buildApp(ctrl), async function (base) {
    for (const [method, path] of [['GET', '/primers'], ['POST', '/primers'], ['GET', '/primers/'], ['GET', '/primers/nope'], ['DELETE', '/primers/design/x']]) {
      const r = await call(base + path, { method: method });
      r.status.should.equal(404, method + ' ' + path);
      r.body.code.should.equal('NOT_FOUND');
      r.body.details.method.should.equal(method);
      r.headers.get('cache-control').should.equal('no-store');
      r.headers.get('content-type').should.match(/application\/json/);
    }
    const r = await call(base + '/primers/check/0123456789abcdef0123456789abcdef', { method: 'PUT' });
    r.status.should.equal(405);
    r.body.code.should.equal('METHOD_NOT_ALLOWED');
    r.body.details.should.eql({ method: 'PUT', allowed_methods: ['GET'] });
    r.headers.get('allow').should.equal('GET');
    r.headers.get('cache-control').should.equal('no-store');
  });
});

test('errorHandler: generic errors become 500 INTERNAL, validator-shaped 4xx keep errors[], answered requests are left to express', async function () {
  const log = quietLog();
  const ctrl = require(CONTROLLER).createController({ log: log, config: fakeConfig() });
  const app = buildApp(ctrl, function (app, base) {
    app.get(base + '/crash', function (req, res, next) { next(new Error('kaboom /home/olson/x')); });
    app.get(base + '/invalid', function (req, res, next) {
      const e = new Error('Validation errors');
      e.statusCode = 400;
      e.errors = [{ code: 'PATTERN', message: 'bad', path: ['system_name'] }];
      next(e);
    });
    app.get(base + '/answered', function (req, res, next) {
      res.status(400).json({ message: 'Validation errors', errors: [] });
      next(); // what swagger-node-runner's _finish does after writing a validator error
    });
  });
  await withServer(app, async function (base) {
    let r = await call(base + '/primers/crash');
    r.status.should.equal(500);
    r.body.should.eql({ message: 'internal error', code: 'INTERNAL', details: {} });
    r = await call(base + '/primers/invalid');
    r.status.should.equal(400);
    r.body.should.eql({ message: 'Validation errors', code: 'INVALID_REQUEST', errors: [{ code: 'PATTERN', message: 'bad', path: ['system_name'] }] });
    r.headers.get('cache-control').should.equal('no-store');
    r = await call(base + '/primers/answered');
    r.status.should.equal(400);
    r.body.should.eql({ message: 'Validation errors', errors: [] });
  });
  log.lines.join('\n').should.match(/kaboom/);
});

test('startup(basePath) sets the basePath, logs "primers site_key=<key>" and never throws', function () {
  const log = quietLog();
  const config = fakeConfig();
  const ctrl = require(CONTROLLER).createController({
    log: log, config: config, jobs: require('../../../api/helpers/primers/jobs')
  });
  ctrl.startup('/sorghum_v11').should.equal('sorghum_v11:test');
  config.basePathSet.should.equal('/sorghum_v11');
  log.lines[0].should.equal('primers site_key=sorghum_v11:test');

  const bad = quietLog();
  const broken = require(CONTROLLER).createController({
    log: bad, config: { setBasePath: function () {}, get: function () { return {}; }, siteKey: function () { throw new Error('no mongo config'); } }
  });
  should(broken.startup('/sorghum_v11')).equal(null);
  bad.lines.join('\n').should.match(/startup check failed.*no mongo config/);
});

test('paramValue reads req.swagger.params, falling back to body, query and path', function () {
  const { paramValue } = require(CONTROLLER);
  should(paramValue({ swagger: { params: { body: { value: { a: 1 } } } } }, 'body', 'body')).eql({ a: 1 });
  should(paramValue({ swagger: { params: { job_id: { value: undefined } } }, params: { job_id: 'x' } }, 'job_id', 'path')).equal(undefined);
  should(paramValue({ body: { b: 2 } }, 'body', 'body')).eql({ b: 2 });
  should(paramValue({ query: { system_name: 's' } }, 'system_name', 'query')).equal('s');
  should(paramValue({ params: { job_id: 'y' } }, 'job_id', 'path')).equal('y');
  should(paramValue(undefined, 'body', 'body')).equal(undefined);
});
