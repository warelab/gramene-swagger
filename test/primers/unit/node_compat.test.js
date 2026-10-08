'use strict';

const test = require('node:test');
const should = require('should');
const util = require('util');

const MODULE = require.resolve('../../../api/helpers/primers/node_compat');

test('node_compat installs the legacy util.is* family', function () {
  const nc = require(MODULE);
  nc.names.length.should.equal(15);
  nc.names.forEach(function (n) {
    (typeof util[n]).should.equal('function', n);
  });
  util.isArray([]).should.equal(true);
  util.isBoolean(false).should.equal(true);
  util.isBuffer(Buffer.alloc(1)).should.equal(true);
  util.isDate(new Date()).should.equal(true);
  util.isError(new TypeError('x')).should.equal(true);
  util.isFunction(function () {}).should.equal(true);
  util.isNull(null).should.equal(true);
  util.isNullOrUndefined(undefined).should.equal(true);
  util.isNumber(1).should.equal(true);
  util.isObject({}).should.equal(true);
  util.isObject(null).should.equal(false);
  util.isPrimitive('s').should.equal(true);
  util.isPrimitive({}).should.equal(false);
  util.isRegExp(/x/).should.equal(true);
  util.isString('s').should.equal(true);
  util.isSymbol(Symbol('s')).should.equal(true);
  util.isUndefined(undefined).should.equal(true);
});

test('node_compat is idempotent and never replaces an existing function', function () {
  require(MODULE);
  const saved = util.isNumber;
  const custom = function () { return 'custom'; };
  util.isNumber = custom;
  try {
    delete require.cache[MODULE];
    const again = require(MODULE);
    util.isNumber.should.equal(custom);
    again.installed.should.not.containEql('isNumber');
  } finally {
    util.isNumber = saved;
    delete require.cache[MODULE];
    require(MODULE);
  }
});

test('config@1.x loads after node_compat', function () {
  require(MODULE);
  const config = require('config');
  should(config.util.toObject({ a: [1, { b: 2 }] })).eql({ a: [1, { b: 2 }] });
});
