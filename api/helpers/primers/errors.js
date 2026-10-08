'use strict';

// Error type and JSON error responses for the /primers endpoints.
// Body shape: {message, code, details}. Every error response is Cache-Control: no-store.
// 503 responses carry details.retry_after_s and a matching Retry-After header.
// Anything that is not a PrimerHttpError becomes 500 INTERNAL with a generic message:
// stack traces and filesystem paths never reach the client.

const path = require('path');

const DEFAULT_RETRY_AFTER_S = 5;

// Absolute paths under these top-level directories are reduced to their basename.
// API paths such as /sorghum_v11/primers/design are left alone.
const FS_ROOTS = ['home', 'scratch', 'usr', 'tmp', 'var', 'opt', 'etc', 'proc', 'root', 'mnt', 'srv',
  'data', 'nfs', 'lib', 'lib64', 'bin', 'sbin', 'dev', 'run', 'sys', 'private', 'Users', 'Volumes'];
const PATH_RE = new RegExp('(^|[\\s\'"`(=,\\[{<])(/(?:' + FS_ROOTS.join('|') + ')(?:/[^\\s\'"`()<>,;:\\[\\]{}]+)+)', 'g');

function redactPaths(str) {
  if (typeof str !== 'string') return str;
  return str.replace(PATH_RE, function (m, lead, p) { return lead + path.basename(p); });
}

function redactDeep(value, depth) {
  depth = depth || 0;
  if (typeof value === 'string') return redactPaths(value);
  if (depth > 6 || value === null || typeof value !== 'object') return value;
  if (Array.isArray(value)) return value.map(function (v) { return redactDeep(v, depth + 1); });
  const out = {};
  Object.keys(value).forEach(function (k) { out[k] = redactDeep(value[k], depth + 1); });
  return out;
}

function normalizeStatus(status) {
  const s = Number(status);
  return Number.isInteger(s) && s >= 400 && s <= 599 ? s : 500;
}

class PrimerHttpError extends Error {
  constructor(status, code, message, details) {
    super(message || code || 'error');
    this.name = 'PrimerHttpError';
    this.status = normalizeStatus(status);
    this.statusCode = this.status;
    this.code = typeof code === 'string' && code ? code : 'INTERNAL';
    let d = details !== null && typeof details === 'object' && !Array.isArray(details) ? Object.assign({}, details) : {};
    if (this.status === 503 && !Number.isFinite(Number(d.retry_after_s))) d.retry_after_s = DEFAULT_RETRY_AFTER_S;
    this.details = d;
  }

  toJSON() {
    return { message: this.message, code: this.code, details: this.details };
  }
}

function isPrimerHttpError(err) {
  return err instanceof PrimerHttpError ||
    (!!err && err.name === 'PrimerHttpError' && Number.isInteger(err.status) && typeof err.code === 'string');
}

// Map any thrown value to a client-safe PrimerHttpError.
function toPublicError(err, log) {
  log = log || console;
  if (isPrimerHttpError(err)) {
    if (err.status >= 500) logSafe(log, 'primers ' + err.status + ' ' + err.code + ': ' + err.message);
    return new PrimerHttpError(err.status, err.code, redactPaths(err.message), redactDeep(err.details));
  }
  // http-errors style client errors (body-parser: 400 entity.parse.failed, 413 entity.too.large).
  const st = err && Number(err.status || err.statusCode);
  if (err && err.expose === true && Number.isInteger(st) && st >= 400 && st < 500) {
    const code = st === 413 ? 'PAYLOAD_TOO_LARGE' : 'INVALID_REQUEST';
    return new PrimerHttpError(st, code, redactPaths(String(err.message || code)), {});
  }
  logSafe(log, 'primers unhandled error: ' + (err && err.stack ? err.stack : String(err)));
  return new PrimerHttpError(500, 'INTERNAL', 'internal error', {});
}

function logSafe(log, msg) {
  try {
    (log.error || log.log || function () {}).call(log, msg);
  } catch (e) { /* logging must never throw */ }
}

function setHeader(res, name, value) {
  if (typeof res.set === 'function') res.set(name, value);
  else if (typeof res.setHeader === 'function') res.setHeader(name, value);
}

// Send err as a JSON error response. Returns false when the response was already started.
function sendError(res, err, opts) {
  const log = (opts && opts.log) || console;
  const e = toPublicError(err, log);
  if (!res) return false;
  if (res.headersSent) {
    logSafe(log, 'primers error after headers were sent: ' + e.code);
    try {
      if (!res.writableEnded && typeof res.end === 'function') res.end();
    } catch (x) { /* ignore */ }
    return false;
  }
  const body = { message: e.message, code: e.code, details: e.details };
  setHeader(res, 'Cache-Control', 'no-store');
  if (e.status === 503) {
    setHeader(res, 'Retry-After', String(Math.max(0, Math.ceil(Number(e.details.retry_after_s)))));
  }
  if (typeof res.status === 'function' && typeof res.json === 'function') {
    res.status(e.status).json(body);
  } else {
    res.statusCode = e.status;
    setHeader(res, 'Content-Type', 'application/json; charset=utf-8');
    res.end(JSON.stringify(body));
  }
  return true;
}

// Promise-wrap an (async) express handler so a throw or rejection becomes a JSON error
// and never an unhandled rejection. Resolves when the handler (or the error send) is done.
function wrap(handler, opts) {
  return function (req, res, next) {
    return Promise.resolve()
      .then(function () { return handler(req, res, next); })
      .catch(function (e) { sendError(res, e, opts); });
  };
}

module.exports = {
  PrimerHttpError,
  isPrimerHttpError,
  toPublicError,
  sendError,
  wrap,
  redactPaths,
  DEFAULT_RETRY_AFTER_S
};
