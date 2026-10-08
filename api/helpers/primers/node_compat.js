'use strict';

// Node 23+ removed the legacy util.is* family that config@1.x (and other old
// dependencies such as swagger-express-mw) still call. This is the polyfill from
// app.js, packaged so standalone entry points (the check worker, unit tests) can
// load it before `config` or any swagger library. Idempotent: existing functions
// are never replaced, so requiring it after app.js has patched util is harmless.

const util = require('util');

const tag = function (t) {
  return function (v) { return Object.prototype.toString.call(v) === '[object ' + t + ']'; };
};

const polyfills = {
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

const installed = [];
Object.keys(polyfills).forEach(function (k) {
  if (typeof util[k] !== 'function') {
    util[k] = polyfills[k];
    installed.push(k);
  }
});

module.exports = {
  names: Object.freeze(Object.keys(polyfills)),
  installed: Object.freeze(installed)
};
