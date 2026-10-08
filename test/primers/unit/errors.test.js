'use strict';

const test = require('node:test');
const should = require('should');

const errors = require('../../../api/helpers/primers/errors');
const { PrimerHttpError, sendError, wrap, redactPaths, isPrimerHttpError, toPublicError } = errors;

function fakeRes() {
  return {
    headers: {},
    statusCode: 200,
    body: undefined,
    headersSent: false,
    status: function (c) { this.statusCode = c; return this; },
    set: function (k, v) { this.headers[k.toLowerCase()] = v; return this; },
    json: function (b) { this.body = b; this.headersSent = true; return this; }
  };
}

function quietLog() {
  const lines = [];
  return { lines: lines, error: function (m) { lines.push(String(m)); } };
}

test('PrimerHttpError carries status, code, message and details', function () {
  const e = new PrimerHttpError(404, 'UNKNOWN_GENE', 'no such gene', { gene_id: 'X' });
  e.should.be.instanceOf(Error);
  e.name.should.equal('PrimerHttpError');
  e.status.should.equal(404);
  e.statusCode.should.equal(404);
  e.code.should.equal('UNKNOWN_GENE');
  e.message.should.equal('no such gene');
  e.details.should.eql({ gene_id: 'X' });
  JSON.parse(JSON.stringify(e)).should.eql({ message: 'no such gene', code: 'UNKNOWN_GENE', details: { gene_id: 'X' } });
  isPrimerHttpError(e).should.equal(true);
  isPrimerHttpError(new Error('x')).should.equal(false);
});

test('PrimerHttpError normalizes bad status and fills 503 retry_after_s', function () {
  new PrimerHttpError(200, 'X', 'm').status.should.equal(500);
  new PrimerHttpError('nope', 'X', 'm').status.should.equal(500);
  new PrimerHttpError(400, 'X', 'm').details.should.eql({});
  new PrimerHttpError(503, 'BUSY', 'm').details.retry_after_s.should.equal(5);
  new PrimerHttpError(503, 'QUEUE_FULL', 'm', { retry_after_s: 60 }).details.retry_after_s.should.equal(60);
});

test('sendError writes {message, code, details} with Cache-Control: no-store', function () {
  const res = fakeRes();
  sendError(res, new PrimerHttpError(404, 'UNKNOWN_REGION', 'region "9" is not in this assembly', { region: '9' }), { log: quietLog() }).should.equal(true);
  res.statusCode.should.equal(404);
  res.headers['cache-control'].should.equal('no-store');
  should(res.headers['retry-after']).be.undefined();
  res.body.should.eql({ message: 'region "9" is not in this assembly', code: 'UNKNOWN_REGION', details: { region: '9' } });
});

test('503 responses carry details.retry_after_s and a Retry-After header', function () {
  let res = fakeRes();
  sendError(res, new PrimerHttpError(503, 'BUSY', 'busy', { retry_after_s: 5 }), { log: quietLog() });
  res.statusCode.should.equal(503);
  res.headers['retry-after'].should.equal('5');
  res.body.details.retry_after_s.should.equal(5);

  res = fakeRes();
  sendError(res, new PrimerHttpError(503, 'QUEUE_FULL', 'full', { retry_after_s: 2.2 }), { log: quietLog() });
  res.headers['retry-after'].should.equal('3');

  res = fakeRes();
  sendError(res, new PrimerHttpError(503, 'PRIMER3_UNAVAILABLE', 'gone'), { log: quietLog() });
  res.headers['retry-after'].should.equal('5');
  res.body.details.retry_after_s.should.equal(5);
});

test('unknown errors become 500 INTERNAL without stack traces or paths', function () {
  const log = quietLog();
  const res = fakeRes();
  const boom = new Error("ENOENT: no such file or directory, open '/scratch/olson/fasta/sorghum_bicolor/dna/x.fa.gz.fai'");
  sendError(res, boom, { log: log });
  res.statusCode.should.equal(500);
  res.body.should.eql({ message: 'internal error', code: 'INTERNAL', details: {} });
  const text = JSON.stringify(res.body);
  text.should.not.match(/scratch|\.fai|\bat\s|ENOENT/);
  res.headers['cache-control'].should.equal('no-store');
  // ...but the server log keeps the full stack for operators
  log.lines.join('\n').should.match(/ENOENT/);
  log.lines.join('\n').should.match(/errors\.test\.js/);

  const res2 = fakeRes();
  sendError(res2, 'a string', { log: quietLog() });
  res2.body.code.should.equal('INTERNAL');
  const res3 = fakeRes();
  sendError(res3, undefined, { log: quietLog() });
  res3.statusCode.should.equal(500);
});

test('PrimerHttpError messages and details have filesystem paths reduced to basenames', function () {
  const res = fakeRes();
  sendError(res, new PrimerHttpError(422, 'NO_SEQUENCE', 'missing /scratch/olson/fasta/sorghum_x/dna/X.dna.toplevel.fa.gz',
    { file: '/home/olson/bin/primer3_core', nested: ['see /tmp/abc/def.txt'] }), { log: quietLog() });
  res.body.message.should.equal('missing X.dna.toplevel.fa.gz');
  res.body.details.should.eql({ file: 'primer3_core', nested: ['see def.txt'] });
});

test('redactPaths keeps API paths and plain text', function () {
  redactPaths('spawn /home/olson/bin/primer3_core ENOENT').should.equal('spawn primer3_core ENOENT');
  redactPaths('POST /sorghum_v11/primers/design failed').should.equal('POST /sorghum_v11/primers/design failed');
  redactPaths('region 4:7400001-7450000').should.equal('region 4:7400001-7450000');
  redactPaths("open '/usr/local/gramene/x/y.js'").should.equal("open 'y.js'");
  should(redactPaths(null)).equal(null);
});

test('http-errors style client errors keep their 4xx status', function () {
  const tooBig = Object.assign(new Error('request entity too large'), { status: 413, expose: true, type: 'entity.too.large' });
  const e = toPublicError(tooBig, quietLog());
  e.status.should.equal(413);
  e.code.should.equal('PAYLOAD_TOO_LARGE');
  const parse = Object.assign(new Error('Unexpected token'), { statusCode: 400, expose: true });
  toPublicError(parse, quietLog()).code.should.equal('INVALID_REQUEST');
  const internal = Object.assign(new Error('secret'), { status: 500, expose: false });
  toPublicError(internal, quietLog()).code.should.equal('INTERNAL');
});

test('sendError works with a bare http.ServerResponse-like object', function () {
  const res = {
    statusCode: 200,
    headers: {},
    headersSent: false,
    setHeader: function (k, v) { this.headers[k.toLowerCase()] = v; },
    end: function (b) { this.body = b; this.headersSent = true; }
  };
  sendError(res, new PrimerHttpError(400, 'INVALID_PARAMS', 'min_size > max_size'), { log: quietLog() });
  res.statusCode.should.equal(400);
  res.headers['content-type'].should.match(/application\/json/);
  res.headers['cache-control'].should.equal('no-store');
  JSON.parse(res.body).should.eql({ message: 'min_size > max_size', code: 'INVALID_PARAMS', details: {} });
});

test('sendError after headers were sent ends the response and returns false', function () {
  let ended = false;
  const res = fakeRes();
  res.headersSent = true;
  res.end = function () { ended = true; };
  sendError(res, new Error('late'), { log: quietLog() }).should.equal(false);
  ended.should.equal(true);
  should(res.body).be.undefined();
  sendError(null, new Error('x'), { log: quietLog() }).should.equal(false);
});

test('wrap turns sync throws and async rejections into JSON errors, never unhandled rejections', async function () {
  const unhandled = [];
  const onUnhandled = function (r) { unhandled.push(r); };
  process.on('unhandledRejection', onUnhandled);
  try {
    const log = quietLog();
    const res1 = fakeRes();
    await wrap(function () { throw new PrimerHttpError(400, 'INVALID_REQUEST', 'bad'); }, { log: log })({}, res1);
    res1.statusCode.should.equal(400);
    res1.body.code.should.equal('INVALID_REQUEST');

    const res2 = fakeRes();
    await wrap(async function () {
      await new Promise(function (r) { setImmediate(r); });
      throw new Error('kaboom /home/olson/secret');
    }, { log: log })({}, res2);
    res2.statusCode.should.equal(500);
    res2.body.should.eql({ message: 'internal error', code: 'INTERNAL', details: {} });

    const res3 = fakeRes();
    await wrap(async function (req, res) { res.status(200).json({ ok: true }); }, { log: log })({}, res3);
    res3.body.should.eql({ ok: true });

    await new Promise(function (r) { setImmediate(r); });
    unhandled.length.should.equal(0);
  } finally {
    process.removeListener('unhandledRejection', onUnhandled);
  }
});
