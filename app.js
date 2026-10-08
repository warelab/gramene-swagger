'use strict';

// Node 24 removed the legacy util.is* family that config@1.x (pulled in by
// swagger-express-mw) still calls. Polyfill the whole set before anything loads
// `config`. (These were deprecated since Node 4 and dropped in Node 23.)
var util = require('util');
var tag = function (t) { return function (v) { return Object.prototype.toString.call(v) === '[object ' + t + ']'; }; };
var polyfills = {
  isArray: Array.isArray,
  isBoolean: function (v) { return typeof v === 'boolean'; },
  isBuffer: function (v) { return Buffer.isBuffer(v); },
  isDate: tag('Date'),
  isError: function (v) { return tag('Error')(v) || v instanceof Error; },
  isFunction: function (v) { return typeof v === 'function'; },
  isNull: function (v) { return v === null; },
  isNullOrUndefined: function (v) { return v == null; },
  isNumber: function (v) { return typeof v === 'number'; },
  isObject: function (v) { return v !== null && typeof v === 'object'; },
  isPrimitive: function (v) { return v === null || (typeof v !== 'object' && typeof v !== 'function'); },
  isRegExp: tag('RegExp'),
  isString: function (v) { return typeof v === 'string'; },
  isSymbol: function (v) { return typeof v === 'symbol'; },
  isUndefined: function (v) { return v === undefined; }
};
Object.keys(polyfills).forEach(function (k) { if (typeof util[k] !== 'function') util[k] = polyfills[k]; });

var SwaggerExpress = require('swagger-express-mw');
var express = require('express');
var cors = require('cors');

var app = express();
app.use(cors());

//var fs = require('fs');
//var yaml = require('js-yaml');

module.exports = app;

var config = {
  appRoot: __dirname // required config
};

// Optional public-host override for a dev instance (e.g. SWAGGER_HOST=localhost:50111 SWAGGER_SCHEMES=http),
// so the docs page targets that instance. Unset: the runner loads api/swagger/swagger.yaml itself, as before.
if (process.env.SWAGGER_HOST) {
  var swaggerDoc = require('js-yaml').safeLoad(
    require('fs').readFileSync(require('path').join(__dirname, 'api', 'swagger', 'swagger.yaml'), 'utf8'));
  swaggerDoc.host = process.env.SWAGGER_HOST;
  swaggerDoc.schemes = (process.env.SWAGGER_SCHEMES || 'http').split(',')
    .map(function (s) { return s.trim(); })
    .filter(function (s) { return s; });
  config.swagger = swaggerDoc;
}

// Primer design endpoints (/primers/*, api/controllers/primers.js). Loaded here only for their
// middleware; a load failure disables that middleware and is logged, the rest of the API is unaffected.
var primersApi = null;
try {
  primersApi = require('./api/controllers/primers');
} catch (e) {
  console.error('primers: controller failed to load:', e && e.stack ? e.stack : e);
}

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

  if (primersApi) {
    // Before register, so validator 400s and body-parser 413s are covered too: no-store on every
    // /primers response, Retry-After readable by browsers, and the bare /primers path answered with a
    // JSON 404 instead of reaching the /{collection} catch-all.
    app.use(basePath + '/primers', primersApi.noStore);
    app.use(basePath + '/primers', cors({ exposedHeaders: ['Retry-After'] }));
    app.all(basePath + '/primers', primersApi.notFound);
  }

  // install swagger server middleware
  swaggerExpress.register(app);

  if (primersApi) {
    // JSON 404 for unknown /primers paths and JSON 405 for wrong methods (instead of express's HTML).
    app.use(basePath + '/primers', primersApi.notFound);
    app.use(basePath + '/primers', primersApi.errorHandler);
    primersApi.startup(basePath); // logs "primers site_key=<key>", which must match the worker's line
  }

  // start it up
  var version = +basePath.match(/\d+/);
  var port = process.env.PORT || 50011;
  app.listen(port, process.env.HOST); // HOST unset: all interfaces, as before

  console.log('Listening on', port, process.env.HOST ? 'host ' + process.env.HOST : '');
});
