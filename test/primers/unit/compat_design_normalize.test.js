'use strict';

// Compatibility C2 (genotyping spec §7.2): design.normalize, with the built-in config, of every design request fixture
// that predates the feature deep-equals the snapshot recorded before any genotyping code, on primer-design d4af3a1 code
// (fixtures/contract/compat_design_normalize.json). /primers/design requests normalize exactly as before: no new key,
// no changed default or preset value.
// Re-record only deliberately: node test/primers/tools/record_compat.js --write

const test = require('node:test');
const should = require('should');

const design = require('../../../api/helpers/primers/design');
const compat = require('../tools/record_compat');
const SNAPSHOT = require('../fixtures/contract/compat_design_normalize.json').normalized;

const NAMES = compat.designFixtureNames();

test('the snapshot covers exactly the 11 existing design-*.json fixtures', function () {
  NAMES.length.should.equal(11);
  Object.keys(SNAPSHOT).sort().should.eql(NAMES);
});

NAMES.forEach(function (name) {
  test('C2 ' + name + ': design.normalize deep-equals the snapshot', function () {
    should(design.normalize(compat.readRequest(name), compat.CFG)).eql(SNAPSHOT[name]);
  });
});
