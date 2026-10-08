'use strict';

// A loopback fake of the Ensembl REST variation service for the genotyping integration tests (genotyping spec §3.3, §3.8,
// §6.4, §7.6). Point an API at it with PRIMERS_VARIATION_URL=<fake.url>; config.js accepts plain http only on 127.0.0.1.
//
// It serves the recorded Ensembl bodies of test/primers/fixtures/variation/, under the path prefix /pansite-ensembl-115:
//   GET .../overlap/region/{species}/{region}:{start}-{end}  every recorded overlap record (the union of the overlap_*.json
//       windows) whose span overlaps start..end
//   GET .../variation/{species}/{id}  variation_{id}.json when recorded, else a lookup built from the recorded overlap record
//       with that id, else Ensembl's HTTP 400 {"error": "<id> not found for <species>"}
// Any species but sorghum_bicolor gets Ensembl's HTTP 400, and any other path an Apache-style HTML 404.
//
// Modes, set with setMode(mode, options):
//   normal     the recorded answers
//   slow       the recorded answers after options.delay_ms (default 9000, above the client's 8 s timeout)
//   5xx        HTTP options.status (default 500) with an HTML body
//   malformed  HTTP 200 application/json with a truncated JSON body
//   html404    the Apache-style HTML 404 of a mistyped base URL, for every path
//   oversize   the recorded answer padded with whitespace to options.bytes (default 1,500,000). It stays valid JSON, so only
//              a byte cap can reject it. It is streamed in 64 KB pieces without Content-Length, unless
//              options.declare_length is set.
//   down       closes the listening socket and every open connection, so connections are refused. Any other mode
//              reopens the same port.
// The aliases 500 (for 5xx) and refused (for down) match the spec §6.4 names.
//
//   const fake = await createFakeEnsembl().start();
//   await fake.setMode('slow', {delay_ms: 7000});
//   fake.requests       [{kind: 'overlap' | 'variation' | 'other', key, path, mode, status, at, done_at, aborted}]
//                       key mirrors the client's cache keys: o|species|region|chunkStart and v|species|id
//   fake.count(filter)  requests matching {kind, key, mode}; fake.inflight(filter) those still open
//   await fake.waitFor(predicate, timeoutMs, what)
//   await fake.stop()
//
// CLI for manual runs against a dev API (spec §6.4). The control port takes POST /mode/<mode>?delay_ms=&status=&bytes=&declare_length=
//   node test/primers/integration/helpers/fake_ensembl.js --port 50199 --control-port 50198 --mode normal
//   curl -s -X POST 'http://127.0.0.1:50198/mode/slow?delay_ms=9000'

const http = require('http');
const fs = require('fs');
const path = require('path');

const FIXTURES = path.resolve(__dirname, '..', '..', 'fixtures', 'variation');
const HOST = '127.0.0.1';
const PREFIX = '/pansite-ensembl-115';
const SPECIES = 'sorghum_bicolor';
const MODES = Object.freeze(['normal', 'slow', '5xx', 'malformed', 'html404', 'oversize', 'down']);
const ALIASES = Object.freeze({ 500: '5xx', refused: 'down' });
const DEFAULTS = Object.freeze({ delay_ms: 9000, status: 500, bytes: 1500000, declare_length: false });
const STREAM_PIECE = 65536;
const JSON_TYPE = 'application/json';
const HTML_TYPE = 'text/html; charset=iso-8859-1';

const APACHE_404 = '<!DOCTYPE HTML PUBLIC "-//IETF//DTD HTML 2.0//EN">\n<html><head>\n<title>404 Not Found</title>\n</head><body>\n' +
  '<h1>Not Found</h1>\n<p>The requested URL was not found on this server.</p>\n</body></html>\n';
const HTML_5XX = '<html><head><title>500 Internal Server Error</title></head><body><h1>Internal Server Error</h1></body></html>\n';
const TRUNCATED_JSON = '[{"feature_type":"variation","seq_region_name":"1","start":11109,"end":11109,"id":"rs871475760","alleles":["C",';

// Every recorded overlap record once (the recorded windows overlap).
function loadRecords(dir) {
  const seen = new Set();
  const out = [];
  fs.readdirSync(dir).filter(function (f) { return /^overlap_.*\.json$/.test(f); }).sort().forEach(function (f) {
    JSON.parse(fs.readFileSync(path.join(dir, f), 'utf8')).forEach(function (r) {
      const k = JSON.stringify([r.id, r.seq_region_name, r.start, r.end, r.alleles]);
      if (seen.has(k)) return;
      seen.add(k);
      out.push(r);
    });
  });
  return out;
}

function loadLookups(dir) {
  const out = new Map();
  fs.readdirSync(dir).forEach(function (f) {
    const m = /^variation_(.+)\.json$/.exec(f);
    if (m && m[1] !== 'not_found') out.set(m[1], fs.readFileSync(path.join(dir, f), 'utf8'));
  });
  return out;
}

function decode(s) {
  try {
    return decodeURIComponent(s);
  } catch (e) {
    return null;
  }
}

function json(status, body, kind, key) {
  return { kind: kind, key: key, status: status, type: JSON_TYPE, body: typeof body === 'string' ? body : JSON.stringify(body) };
}

function other() {
  return { kind: 'other', key: null, status: 404, type: HTML_TYPE, body: APACHE_404 };
}

// Ensembl's answer to a species it does not know.
function unknownSpecies(species, kind, key) {
  return json(400, { error: 'Can not find internal name for species \'' + species + '\'' }, kind, key);
}

// The recorded answer for a path: {kind, key, status, type, body}.
function recordedAnswer(pathname, data) {
  if (pathname.indexOf(PREFIX + '/') !== 0) return other();
  const rest = pathname.slice(PREFIX.length);
  let m = /^\/overlap\/region\/([^/]+)\/([^/]+):(\d+)-(\d+)$/.exec(rest);
  if (m) {
    const species = decode(m[1]);
    const region = decode(m[2]);
    const start = Number(m[3]);
    const end = Number(m[4]);
    if (species === null || region === null) return other();
    const key = 'o|' + species + '|' + region + '|' + start;
    if (species !== SPECIES) return unknownSpecies(species, 'overlap', key);
    return json(200, data.records.filter(function (r) {
      return String(r.seq_region_name) === region && Math.min(r.start, r.end) <= end && Math.max(r.start, r.end) >= start;
    }), 'overlap', key);
  }
  m = /^\/variation\/([^/]+)\/([^/]+)$/.exec(rest);
  if (m) {
    const species = decode(m[1]);
    const id = decode(m[2]);
    if (species === null || id === null) return other();
    const key = 'v|' + species + '|' + id;
    if (species !== SPECIES) return unknownSpecies(species, 'variation', key);
    if (data.lookups.has(id)) return json(200, data.lookups.get(id), 'variation', key);
    const r = data.records.find(function (x) { return x.id === id; });
    if (r) {
      return json(200, {
        name: id, synonyms: [], most_severe_consequence: r.consequence_type,
        mappings: [{ seq_region_name: String(r.seq_region_name), start: r.start, end: r.end, allele_string: r.alleles.join('/') }]
      }, 'variation', key);
    }
    return json(400, { error: id + ' not found for ' + species }, 'variation', key);
  }
  return other();
}

// Whitespace after the opening bracket keeps a JSON body valid at any size.
function padded(body, bytes) {
  const length = Buffer.byteLength(body, 'utf8');
  if (length >= bytes || !/^[[{]/.test(body)) return body;
  return body.charAt(0) + ' '.repeat(bytes - length) + body.slice(1);
}

function send(res, entry, a) {
  if (res.destroyed || res.headersSent) return;
  const body = Buffer.from(a.body, 'utf8');
  entry.status = a.status;
  res.writeHead(a.status, { 'Content-Type': a.type, 'Content-Length': body.length });
  res.end(body);
}

function stream(res, entry, a, bytes, declareLength) {
  if (res.destroyed || res.headersSent) return;
  const body = Buffer.from(padded(a.body, bytes), 'utf8');
  const headers = { 'Content-Type': a.type };
  if (declareLength) headers['Content-Length'] = body.length;
  entry.status = a.status;
  res.writeHead(a.status, headers);
  let offset = 0;
  const pump = function () {
    while (offset < body.length) {
      if (res.destroyed) return;
      const piece = body.subarray(offset, offset + STREAM_PIECE);
      offset += piece.length;
      if (!res.write(piece)) {
        res.once('drain', pump);
        return;
      }
    }
    res.end();
  };
  pump();
}

function matches(r, filter) {
  filter = filter || {};
  return Object.keys(filter).every(function (k) { return r[k] === filter[k]; });
}

function createFakeEnsembl(opts) {
  opts = opts || {};
  const dir = opts.fixtures || FIXTURES;
  const data = { records: loadRecords(dir), lookups: loadLookups(dir) };
  const state = { mode: 'normal', options: Object.assign({}, DEFAULTS) };
  const requests = [];
  let port = opts.port || 0;
  let server = null;
  let listening = false;
  let started = false;

  function handle(req, res) {
    req.on('error', function () { /* the client went away */ });
    res.on('error', function () { /* the client went away */ });
    const pathname = new URL(req.url, 'http://' + HOST).pathname;
    const answer = recordedAnswer(pathname, data);
    const entry = { kind: answer.kind, key: answer.key, path: pathname, mode: state.mode, status: null, at: Date.now(), done_at: null, aborted: false };
    requests.push(entry);
    res.on('close', function () {
      entry.done_at = Date.now();
      entry.aborted = !res.writableFinished;
    });
    const o = state.options;
    switch (state.mode) {
      case 'slow': {
        const timer = setTimeout(function () { send(res, entry, answer); }, o.delay_ms);
        res.on('close', function () { clearTimeout(timer); });
        return;
      }
      case '5xx':
        return send(res, entry, { status: o.status, type: HTML_TYPE, body: HTML_5XX });
      case 'malformed':
        return send(res, entry, { status: 200, type: JSON_TYPE, body: TRUNCATED_JSON });
      case 'html404':
        return send(res, entry, other());
      case 'oversize':
        return stream(res, entry, answer, o.bytes, o.declare_length);
      default:
        return send(res, entry, answer);
    }
  }

  function listen() {
    return new Promise(function (resolve, reject) {
      server = http.createServer(handle);
      server.on('clientError', function (err, socket) { socket.destroy(); });
      server.once('error', reject);
      server.listen(port, HOST, function () {
        server.removeListener('error', reject);
        port = server.address().port;
        listening = true;
        resolve();
      });
    });
  }

  function close() {
    if (!listening) return Promise.resolve();
    listening = false;
    const s = server;
    return new Promise(function (resolve) {
      s.close(function () { resolve(); });
      s.closeAllConnections();
    });
  }

  const fake = {
    get url() { return 'http://' + HOST + ':' + port + PREFIX; },
    get port() { return port; },
    get mode() { return state.mode; },
    get options() { return Object.assign({}, state.options); },
    get listening() { return listening; },
    requests: requests,

    start: async function () {
      if (!started) {
        started = true;
        await listen();
      }
      return fake;
    },

    stop: async function () {
      started = false;
      await close();
    },

    setMode: async function (mode, options) {
      const m = ALIASES[mode] || mode;
      if (MODES.indexOf(m) < 0) throw new TypeError('unknown fake Ensembl mode ' + JSON.stringify(mode) + '; modes: ' + MODES.join(', '));
      state.mode = m;
      state.options = Object.assign({}, DEFAULTS, options || {});
      if (m === 'down') await close();
      else if (started && !listening) await listen();
      return fake;
    },

    reset: function () {
      requests.length = 0;
    },

    count: function (filter) {
      return requests.filter(function (r) { return matches(r, filter); }).length;
    },

    inflight: function (filter) {
      return requests.filter(function (r) { return r.done_at === null && matches(r, filter); }).length;
    },

    waitFor: async function (predicate, timeoutMs, what) {
      const until = Date.now() + (timeoutMs || 10000);
      while (!predicate()) {
        if (Date.now() > until) throw new Error('fake Ensembl: timed out waiting for ' + (what || 'a condition'));
        await new Promise(function (resolve) { setTimeout(resolve, 20); });
      }
    }
  };
  return fake;
}

// ---- CLI --------------------------------------------------------------------------------------------------------------

function cliOptions(o) {
  const out = {};
  const num = function (v) { return v === undefined || v === '' ? undefined : Number(v); };
  if (num(o.delay_ms) !== undefined) out.delay_ms = num(o.delay_ms);
  if (num(o.status) !== undefined) out.status = num(o.status);
  if (num(o.bytes) !== undefined) out.bytes = num(o.bytes);
  if (o.declare_length !== undefined) out.declare_length = /^(1|true|yes|)$/i.test(String(o.declare_length));
  return out;
}

if (require.main === module) {
  const args = {};
  const argv = process.argv.slice(2);
  for (let i = 0; i < argv.length; i++) {
    const m = /^--([a-z-]+)(?:=(.*))?$/.exec(argv[i]);
    if (!m) continue;
    const name = m[1].replace(/-/g, '_');
    if (m[2] !== undefined) args[name] = m[2];
    else if (i + 1 < argv.length && !/^--/.test(argv[i + 1])) args[name] = argv[++i];
    else args[name] = '';
  }
  (async function () {
    const fake = await createFakeEnsembl({ port: Number(args.port) || 50199 }).start();
    await fake.setMode(args.mode || 'normal', cliOptions(args));
    const controlPort = Number(args.control_port) || 50198;
    const control = http.createServer(function (req, res) {
      const u = new URL(req.url, 'http://' + HOST);
      const m = /^\/mode\/([A-Za-z0-9]+)$/.exec(u.pathname);
      if (req.method !== 'POST' || !m) {
        res.writeHead(404, { 'Content-Type': 'text/plain' });
        return res.end('POST /mode/<' + MODES.join('|') + '>\n');
      }
      fake.setMode(m[1], cliOptions(Object.fromEntries(u.searchParams))).then(function () {
        res.writeHead(200, { 'Content-Type': JSON_TYPE });
        res.end(JSON.stringify({ mode: fake.mode, options: fake.options, url: fake.url }) + '\n');
      }, function (err) {
        res.writeHead(400, { 'Content-Type': 'text/plain' });
        res.end(err.message + '\n');
      });
    });
    control.listen(controlPort, HOST, function () {
      console.log('fake Ensembl ' + fake.url + ' mode ' + fake.mode + '; control http://' + HOST + ':' + controlPort + '/mode/<mode>');
    });
    const shutdown = function () {
      control.close();
      fake.stop().then(function () { process.exit(0); });
    };
    process.on('SIGTERM', shutdown);
    process.on('SIGINT', shutdown);
  })().catch(function (err) {
    console.error(err && err.stack ? err.stack : err);
    process.exit(1);
  });
}

module.exports = { createFakeEnsembl, MODES, PREFIX, _internal: { recordedAnswer, loadRecords, loadLookups, padded } };
