'use strict';

// Compatibility C1 (genotyping spec §7.2, §5.3): the job ids of the check request fixtures that predate the feature are
// pinned. Each check-*.json body in fixtures/contract/requests is normalized offline (sorghum_v11 catalog fixture, stub
// assemblies with fingerprint sha1('fp:' + system_name), stub gene SORBI_3004G087700; the stubs live in
// test/primers/tools/record_compat.js) and hashed with jobs.jobId(norm.request, norm.dbs, check.ALGORITHM_VERSION).
// fixtures/contract/compat_job_ids.json was recorded before any genotyping code, with api/helpers/primers/check and
// jobs/index.js identical to primer-design d4af3a1. A request without `genotyping` keeps its id after the feature lands,
// from jobId and through jobs.submit (whose algorithm version M7 routes through check.algorithmVersionFor).
// The feature's own check-genotyping-*.json fixtures are not pinned here.
// Re-record only deliberately: node test/primers/tools/record_compat.js --write

const test = require('node:test');
const should = require('should');

const jobs = require('../../../api/helpers/primers/jobs');
const check = require('../../../api/helpers/primers/check');
const { normalize } = require('../../../api/helpers/primers/check/normalize');
const { createMemoryStore } = require('../../../api/helpers/primers/jobs/memory_store');
const compat = require('../tools/record_compat');
const PINNED = require('../fixtures/contract/compat_job_ids.json');

const NAMES = compat.checkFixtureNames();

test('compat_job_ids.json holds exactly one 32-hex id per existing check fixture (11), all distinct', function () {
  NAMES.length.should.equal(11);
  Object.keys(PINNED.ids).sort().should.eql(NAMES);
  Object.keys(PINNED.request_sha256).sort().should.eql(NAMES);
  NAMES.forEach(function (name) {
    PINNED.ids[name].should.be.a.String();
    PINNED.ids[name].should.match(jobs.JOB_ID_RE, name);
  });
  new Set(Object.values(PINNED.ids)).size.should.equal(NAMES.length);
});

test('ALGORITHM_VERSION is still the version the ids were recorded with (\'2\')', function () {
  PINNED.algorithm_version.should.equal('2');
  check.ALGORITHM_VERSION.should.equal(PINNED.algorithm_version);
});

NAMES.forEach(function (name) {
  test('C1 ' + name + ': the job id is unchanged, from jobs.jobId and through jobs.submit', async function () {
    const body = compat.readRequest(name);
    compat.requestDigest(body).should.equal(PINNED.request_sha256[name], name + ' changed since its id was recorded, so the pin no longer applies');

    const norm = await normalize(JSON.parse(JSON.stringify(body)), compat.checkDeps());
    should(norm.request).not.have.property('genotyping');
    jobs.jobId(norm.request, norm.dbs, check.ALGORITHM_VERSION).should.equal(PINNED.ids[name]);

    // The API path with the real check module, as the controller calls it.
    const store = createMemoryStore({ cfg: compat.CFG, siteKey: 'compat_ids_test', shared: { slots: new Map(), panSlots: new Map() } });
    const submitted = await jobs.submit(JSON.parse(JSON.stringify(body)), compat.checkDeps({ store: store, check: check }));
    submitted.job_id.should.equal(PINNED.ids[name]);
    submitted.created.should.be.true();
  });
});
