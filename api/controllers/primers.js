'use strict';

// /primers endpoints (docs/primer_design_api.md; spec §A.2, §A.9).
//
//   POST /primers/design         designPrimers      synchronous Primer3 design
//   GET  /primers/genomes        primerGenomes      same-species genomes for a system_name
//   POST /primers/check          submitPrimerCheck  queue (or find) a specificity / pan-genome check job
//   GET  /primers/check/{job_id} getPrimerCheck     job status, progress and (partial) results
//   GET  /primers/variants       listPrimerVariants known variants (Ensembl) in a window, normalized
//   GET  /primers/variants/{variant_id} getPrimerVariant  one Ensembl variation id, normalized
//
// Every handler is promise-wrapped: a synchronous throw or a rejected promise becomes a JSON error
// {message, code, details} via errors.sendError, never an unhandled rejection (Node 24 would crash the
// whole API on one). Every response carries Cache-Control: no-store.
//
// The helper modules are required lazily inside the handlers, so loading this controller (app.js does it
// at startup for the /primers middleware) cannot fail because of a helper problem, and the endpoints that
// do not use Redis (design, genomes) never load the job store.
//
// app.js also uses the middleware exported here: noStore and notFound before swaggerExpress.register,
// notFound and errorHandler after it (JSON 404/405 for /primers/*), and startup(basePath).

const errors = require('../helpers/primers/errors');

const { PrimerHttpError } = errors;

const FEATURE_DISABLED_RETRY_AFTER_S = 300;
const MAX_ECHOED_PATH = 200;

function setNoStore(res) {
  if (!res || res.headersSent) return;
  if (typeof res.setHeader === 'function') res.setHeader('Cache-Control', 'no-store');
}

// req.swagger.params[name].value, with a plain-express fallback (body / query / path) so a handler
// mounted outside the swagger pipe still reads its input.
function paramValue(req, name, location) {
  const p = req && req.swagger && req.swagger.params && req.swagger.params[name];
  if (p && typeof p === 'object' && 'value' in p) return p.value;
  if (!req) return undefined;
  if (location === 'body') return req.body;
  if (location === 'query') return req.query ? req.query[name] : undefined;
  if (location === 'path') return req.params ? req.params[name] : undefined;
  return undefined;
}

function echoPath(req) {
  const p = String((req && (req.originalUrl || req.url)) || '').split('?')[0];
  return p.length > MAX_ECHOED_PATH ? p.slice(0, MAX_ECHOED_PATH) + '...' : p;
}

// createController(deps) -> handlers and middleware.
//   deps (tests): {config, design, genomes, jobs, variation (module-like objects), log}
function createController(deps) {
  deps = deps || {};
  const log = deps.log || console;

  const modules = {
    config: function () { return deps.config || require('../helpers/primers/config'); },
    design: function () { return deps.design || require('../helpers/primers/design'); },
    genomes: function () { return deps.genomes || require('../helpers/primers/genomes'); },
    jobs: function () { return deps.jobs || require('../helpers/primers/jobs'); },
    variation: function () { return deps.variation || require('../helpers/primers/variation'); }
  };

  function ensureEnabled() {
    const cfg = modules.config().get();
    if (cfg && cfg.enabled === false) {
      throw new PrimerHttpError(503, 'FEATURE_DISABLED', 'primer design is disabled on this server',
        { retry_after_s: FEATURE_DISABLED_RETRY_AFTER_S });
    }
    return cfg;
  }

  function sendJson(res, status, body) {
    if (res.headersSent || res.destroyed) return;
    setNoStore(res);
    res.status(status).json(body);
  }

  // (req, res[, next]) => Promise that always resolves. The swagger router calls controllers as
  // middleware (req, res, cb); cb is never called because every handler ends the response itself.
  function wrap(name, handler) {
    const wrapped = function (req, res) {
      return Promise.resolve()
        .then(function () {
          setNoStore(res);
          return handler(req, res);
        })
        .catch(function (err) {
          if (res && res.destroyed && !res.headersSent) {
            // The client went away; there is nobody to answer. Log server-side failures only.
            if (!errors.isPrimerHttpError(err) || err.status >= 500) {
              log.error('primers ' + name + ': failed after the client disconnected: ' + (err && err.message));
            }
            return;
          }
          errors.sendError(res, err, { log: log });
        });
    };
    Object.defineProperty(wrapped, 'name', { value: name });
    return wrapped;
  }

  const designPrimers = wrap('designPrimers', async function (req, res) {
    ensureEnabled();
    const body = paramValue(req, 'body', 'body');
    // A client that disconnects aborts the design: Primer3 and any megablast mask are killed.
    const ac = new AbortController();
    const onClose = function () {
      if (!res.writableFinished) {
        ac.abort(new PrimerHttpError(400, 'CLIENT_CLOSED_REQUEST', 'the client closed the connection', {}));
      }
    };
    res.once('close', onClose);
    try {
      const result = await modules.design().design(body, { signal: ac.signal, log: log });
      if (ac.signal.aborted) return;
      sendJson(res, 200, result);
    } finally {
      res.removeListener('close', onClose);
    }
  });

  const primerGenomes = wrap('primerGenomes', async function (req, res) {
    ensureEnabled();
    const systemName = paramValue(req, 'system_name', 'query');
    const result = await modules.genomes().genomesResponse(systemName, { log: log });
    sendJson(res, 200, result);
  });

  const submitPrimerCheck = wrap('submitPrimerCheck', async function (req, res) {
    ensureEnabled();
    const body = paramValue(req, 'body', 'body');
    const submitted = await modules.jobs().submit(body, { log: log });
    const out = Object.assign({}, submitted);
    const created = out.created === true;
    delete out.created;
    sendJson(res, created ? 202 : 200, out);
  });

  const getPrimerCheck = wrap('getPrimerCheck', async function (req, res) {
    ensureEnabled();
    const jobId = paramValue(req, 'job_id', 'path');
    const result = await modules.jobs().status(jobId, { log: log });
    sendJson(res, 200, result);
  });

  // The abort-on-close pattern of designPrimers as a helper: a client that disconnects before the response is
  // written aborts `signal` (reason 400 CLIENT_CLOSED_REQUEST); dispose() removes the listener.
  function clientAbort(res) {
    const ac = new AbortController();
    const onClose = function () {
      if (!res.writableFinished) {
        ac.abort(new PrimerHttpError(400, 'CLIENT_CLOSED_REQUEST', 'the client closed the connection', {}));
      }
    };
    res.once('close', onClose);
    return { signal: ac.signal, dispose: function () { res.removeListener('close', onClose); } };
  }

  // Ensembl waits end when the client goes away; the shared outbound request and its cache entry do not (§3.3).
  const listPrimerVariants = wrap('listPrimerVariants', async function (req, res) {
    ensureEnabled();
    const query = {};
    ['system_name', 'region', 'start', 'end', 'types', 'include_ems', 'limit'].forEach(function (name) {
      const value = paramValue(req, name, 'query');
      if (value !== undefined) query[name] = value;
    });
    const abort = clientAbort(res);
    try {
      const result = await modules.variation().listVariants(query, { signal: abort.signal, log: log });
      if (abort.signal.aborted) return;
      sendJson(res, 200, result);
    } finally {
      abort.dispose();
    }
  });

  const getPrimerVariant = wrap('getPrimerVariant', async function (req, res) {
    ensureEnabled();
    const query = { variant_id: paramValue(req, 'variant_id', 'path'), system_name: paramValue(req, 'system_name', 'query') };
    const abort = clientAbort(res);
    try {
      const result = await modules.variation().lookupVariant(query, { signal: abort.signal, log: log });
      if (abort.signal.aborted) return;
      sendJson(res, 200, result);
    } finally {
      abort.dispose();
    }
  });

  // Before register: every /primers response (validator 400s and body-parser 413s included) is no-store.
  function noStore(req, res, next) {
    setNoStore(res);
    next();
  }

  // Unknown /primers paths (and the bare /primers, which would otherwise reach the /{collection}
  // catch-all). The swagger middleware calls next() after it has already written a validator error,
  // so an answered request is left alone.
  function notFound(req, res, next) {
    if (res.headersSent || res.writableEnded) return;
    errors.sendError(res, new PrimerHttpError(404, 'NOT_FOUND', 'unknown primers endpoint: ' + req.method + ' ' + echoPath(req),
      { method: req.method, path: echoPath(req) }), { log: log });
  }

  // After register: errors the swagger middleware passes to express for /primers paths, chiefly
  // 405 "Path defined in Swagger, but <METHOD> operation is not".
  function errorHandler(err, req, res, next) {
    if (res.headersSent) return next(err);
    const status = Number(err && (err.statusCode || err.status));
    if (status === 405) {
      const allowed = Array.isArray(err.allowedMethods) ? err.allowedMethods.slice().sort() : [];
      if (allowed.length && typeof res.getHeader === 'function' && !res.getHeader('Allow')) {
        res.setHeader('Allow', allowed.join(', '));
      }
      return errors.sendError(res, new PrimerHttpError(405, 'METHOD_NOT_ALLOWED',
        req.method + ' is not allowed on ' + echoPath(req), { method: req.method, allowed_methods: allowed }), { log: log });
    }
    if (err && Array.isArray(err.errors) && status >= 400 && status < 500) {
      // Swagger validator shape, kept as the validator writes it: {message, code?, errors[]}.
      setNoStore(res);
      return res.status(status).json({ message: String(err.message || 'Validation errors'), code: err.code || 'INVALID_REQUEST', errors: err.errors });
    }
    return errors.sendError(res, err, { log: log });
  }

  // app.js, after register: adopt the runtime basePath and log the site key in the same format as the
  // worker's first line, so `grep 'primers site_key='` over both logs must show one value. Never throws.
  function startup(basePath) {
    try {
      const config = modules.config();
      if (basePath) config.setBasePath(basePath);
      const cfg = config.get();
      const key = config.siteKey();
      log.log(modules.jobs().siteKeyLogLine(key));
      log.log('primers api: enabled ' + cfg.enabled + ', basePath ' + config.basePath() + ', store ' + cfg.check.store +
        ', global_prefix ' + cfg.check.global_prefix);
      return key;
    } catch (err) {
      log.error('primers api: startup check failed (endpoints will answer with JSON errors): ' + (err && err.message));
      return null;
    }
  }

  return {
    designPrimers: designPrimers,
    primerGenomes: primerGenomes,
    submitPrimerCheck: submitPrimerCheck,
    getPrimerCheck: getPrimerCheck,
    listPrimerVariants: listPrimerVariants,
    getPrimerVariant: getPrimerVariant,
    noStore: noStore,
    notFound: notFound,
    errorHandler: errorHandler,
    startup: startup
  };
}

module.exports = createController();
module.exports.createController = createController;
module.exports.paramValue = paramValue;
