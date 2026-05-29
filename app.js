'use strict';

// Node 23 removed the whole `util.isX` family of deprecated type-checks
// (isArray, isBoolean, isBuffer, isDate, isError, isFunction, isNull,
// isNullOrUndefined, isNumber, isObject, isPrimitive, isRegExp, isString,
// isSymbol, isUndefined). The unmaintained `config@1.x` — pulled in
// transitively by swagger-node-runner — calls several of them inside
// cloneDeep. Restore them before anything else loads; the require cache
// hands the same `util` module to every later require('util'). Each
// shim is a no-op on Node <= 22 where the original still exists.
{
  const util = require('util');
  const shims = {
    isArray: (x) => Array.isArray(x),
    isBoolean: (x) => typeof x === 'boolean',
    isBuffer: (x) => Buffer.isBuffer(x),
    isDate: (x) => x instanceof Date,
    isError: (x) => x instanceof Error,
    isFunction: (x) => typeof x === 'function',
    isNull: (x) => x === null,
    isNullOrUndefined: (x) => x == null,
    isNumber: (x) => typeof x === 'number',
    isObject: (x) => typeof x === 'object' && x !== null,
    isPrimitive: (x) =>
      x === null ||
      (typeof x !== 'object' && typeof x !== 'function'),
    isRegExp: (x) => x instanceof RegExp,
    isString: (x) => typeof x === 'string',
    isSymbol: (x) => typeof x === 'symbol',
    isUndefined: (x) => x === undefined,
  };
  for (const [name, fn] of Object.entries(shims)) {
    if (typeof util[name] !== 'function') util[name] = fn;
  }
}

// Loaded before anything that reads process.env so a local .env (gitignored)
// can populate SESSION_SECRET, FIREBASE_CREDENTIALS_PATH, GOOGLE_CLIENT_ID
// etc. without requiring an inline `KEY=val node app.js` invocation.
// Production deploys that already inject these via systemd / docker-compose
// are unaffected — dotenv only fills in keys that aren't already set.
require('dotenv').config();

var SwaggerExpress = require('swagger-express-mw');
var express = require('express');
const passport = require('passport');
var cors = require('cors');

var app = express();
app.use(cors());
app.use(require('express-session')({
  secret: process.env.SESSION_SECRET,
  // resave/saveUninitialized: false so anonymous API reads don't get a
  // Set-Cookie response. Without this, every cacheable GET would emit a
  // session id, and an upstream cache would hand the same id to every
  // client on a HIT. Sessions are still created when something writes to
  // req.session (e.g. passport.authenticate during the OAuth flow).
  resave: false,
  saveUninitialized: false,
}));
app.use(passport.initialize());
app.use(passport.session());

// ---------------------------------------------------------------------------
// Cache-Control middleware
//
// Stamps Cache-Control headers so an upstream reverse-proxy cache (Apache
// mod_cache, Varnish, etc.) can safely cache idempotent read endpoints and
// must-not-cache user/auth endpoints.
//
// The default TTL is 1 hour. The cache key upstream is the full URL
// (path + query string), so different facet/filter combinations are cached
// independently. A new data release should trigger an upstream cache flush
// (e.g. `htcacheclean -t`) since the cache is not release-aware.
// ---------------------------------------------------------------------------
const CACHE_TTL_SECONDS = 3600;
const NO_CACHE_PATTERNS = [
  /(^|\/)auth(\/|$)/,
  /(^|\/)logout(\/|$)/,
  /(^|\/)protected-api(\/|$)/,
  /(^|\/)gene_lists(\/|$)/,
  /(^|\/)gxa(\/|$)/,
];
function isCacheableRequest(req) {
  // swagger.yaml only declares `get:` operations; HEAD is rejected by
  // swagger-express-mw with a 500 before the response is built, so there's
  // nothing to cache.
  if (req.method !== 'GET') return false;
  if (NO_CACHE_PATTERNS.some(re => re.test(req.path))) return false;
  // BED exports stream whole genomes — too large for a useful cache,
  // and the wt=bed key would crowd out smaller, hotter responses.
  if (req.query && req.query.wt === 'bed') return false;
  return true;
}
app.use(function cacheControl(req, res, next) {
  if (isCacheableRequest(req)) {
    res.setHeader('Cache-Control', 'public, max-age=' + CACHE_TTL_SECONDS);
    res.setHeader('Vary', 'Accept-Encoding');
  } else {
    res.setHeader('Cache-Control', 'no-store');
  }
  next();
});

app.get('/auth/google', passport.authenticate('google', { scope: ['profile', 'email'] }));
app.get('/auth/google/callback', passport.authenticate('google', { failureRedirect: '/' }),
  (req, res) => {
    res.redirect('/');
  });
app.get('/logout', (req, res) => {
  req.logout();
  res.redirect('/');
});

function isAuthenticated(req, res, next) {
  if (req.isAuthenticated()) {
    return next();
  }
  res.status(401).send('Unauthorized');
}

app.get('/protected-api', isAuthenticated, (req, res) => {
  res.send('This is protected data');
});

module.exports = app;

var config = {
  appRoot: __dirname // required config
};

SwaggerExpress.create(config, function (err, swaggerExpress) {
  if (err) { throw err; }

  var basePath = swaggerExpress.runner.swagger.basePath;

  app.get('/', function(req, res) {
    res.redirect(basePath);
  });

  // define routes for documentation
  app.get(basePath, function (req, res, next) { // redirect to /docs with the correct schema
    res.redirect(basePath+'/docs?url='+basePath+'/swagger');
  });

  app.use(basePath+'/docs', express.static('node_modules/swagger-ui/dist'));

  // install swagger server middleware
  swaggerExpress.register(app);

  // Redirect unimplemented routes to EBI Atlas. The swagger middleware above
  // doesn't call next() after handling its routes in the happy path, but
  // res.headersSent === true is the only reliable signal — Express still
  // runs this catch-all on some swagger paths (e.g. /swagger, which goes
  // through the x-swagger-pipe), and re-redirecting after headers are out
  // throws ERR_HTTP_HEADERS_SENT.
  app.all(`${basePath}/*`, (req, res) => {
    if (res.headersSent) return;
    // Construct the external URL, preserving the original path and query
    const ebiBase = 'https://www.ebi.ac.uk';
    const gxaUrl = ebiBase + req.originalUrl.replace(basePath,'');

    console.log(`Redirecting unhandled route ${req.originalUrl} to ${gxaUrl}`);
    res.redirect(301, gxaUrl); // 301 for permanent, 302 for temporary redirect
  });
  // start it up
  var version = +basePath.match(/\d+/);
  var port = 50008; // process.env.PORT || 10000 + version;
  app.listen(port);

  console.log('Listening on', port);
});

