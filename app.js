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

  // start it up
  var version = +basePath.match(/\d+/);
  var port = process.env.PORT || 50011;
  app.listen(port);

  console.log('Listening on', port);
});

