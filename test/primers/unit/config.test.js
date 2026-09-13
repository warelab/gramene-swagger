'use strict';

require('../../../api/helpers/primers/node_compat');

const test = require('node:test');
const should = require('should');
const fs = require('fs');
const os = require('os');
const path = require('path');
const { execFileSync } = require('child_process');
const yaml = require('js-yaml');

const config = require('../../../api/helpers/primers/config');

const REPO = path.resolve(__dirname, '..', '..', '..');

function yamlPrimersBlock() {
  return yaml.safeLoad(fs.readFileSync(path.join(REPO, 'config', 'default.yaml'), 'utf8')).primers;
}

test('get() merges config/default.yaml primers over DEFAULTS and deep-freezes the result', function () {
  config._setForTests(null);
  const cfg = config.get();
  cfg.enabled.should.equal(true);
  should(cfg.site_key).equal(null);
  cfg.primer3_core.should.equal('/home/olson/bin/primer3_core');
  cfg.blastdbcmd.should.equal('/home/olson/bin/blastdbcmd');
  cfg.check.local_max_jobs.should.equal(2);
  cfg.check.global_prefix.should.equal('primers:global:');
  cfg.check.global_max_jobs.should.equal(2);
  cfg.check.pangenome_max_jobs.should.equal(1);
  cfg.check.defaults.max_product_size.should.equal(4000);
  cfg.design.max_template_length.should.equal(50000);
  cfg.design.max_fetch_length.should.equal(2000000);
  cfg.repeat_mask.evalue.should.equal(1e-10);
  cfg.check.cpu_s_per_primer_gb.should.eql({ ws5: 5.2, ws6: 2.2, ws7: 1.2 });
  Object.prototype.hasOwnProperty.call(cfg.check, 'in_process').should.equal(false); // no supervisor-only keys
  Object.isFrozen(cfg).should.equal(true);
  Object.isFrozen(cfg.check.defaults).should.equal(true);
  (function () { cfg.check.local_max_jobs = 9; }).should.throw(TypeError);
  config.get().should.equal(cfg); // cached
});

test('default.yaml primers block and DEFAULTS agree (a process without the yaml behaves the same)', function () {
  const block = yamlPrimersBlock();
  should.exist(block);
  const fromYaml = config._build({ env: {}, fileConfig: block }).config;
  const fromDefaults = config._build({ env: {}, fileConfig: {} }).config;
  JSON.parse(JSON.stringify(fromYaml)).should.eql(JSON.parse(JSON.stringify(fromDefaults)));
  JSON.parse(JSON.stringify(block)).should.eql(JSON.parse(JSON.stringify(config.DEFAULTS)));
});

test('environment overrides (A.7) with typed parsing', function () {
  const env = {
    PRIMERS_ENABLED: '0',
    PRIMERS_SITE_KEY: 'sorghum_v11_dev',
    PRIMER3_CORE: '/opt/p3/primer3_core',
    BLASTN: '/opt/blast/blastn',
    BLASTDBCMD: '/opt/blast/blastdbcmd',
    PRIMERS_FASTA_ROOT: '/data/fasta',
    PRIMERS_JOB_STORE: 'memory',
    PRIMERS_REDIS_URL: 'redis://localhost:6399/1',
    PRIMERS_GLOBAL_PREFIX: 'primers:test:global:',
    PRIMERS_GLOBAL_MAX_JOBS: '1',
    PRIMERS_MAX_QUEUED: '5'
  };
  const r = config._build({ env: env, fileConfig: {} });
  r.warnings.should.eql([]);
  const c = r.config;
  c.enabled.should.equal(false);
  c.site_key.should.equal('sorghum_v11_dev');
  c.primer3_core.should.equal('/opt/p3/primer3_core');
  c.blastn.should.equal('/opt/blast/blastn');
  c.blastdbcmd.should.equal('/opt/blast/blastdbcmd');
  c.fasta_root.should.equal('/data/fasta');
  c.check.store.should.equal('memory');
  c.check.redis_url.should.equal('redis://localhost:6399/1');
  c.check.global_prefix.should.equal('primers:test:global:');
  c.check.global_max_jobs.should.equal(1);
  c.check.max_queued.should.equal(5);
  c.check.local_max_jobs.should.equal(2); // untouched siblings survive

  ['1', 'true', 'TRUE', 'yes', 'Yes'].forEach(function (v) {
    config._build({ env: { PRIMERS_ENABLED: v }, fileConfig: {} }).config.enabled.should.equal(true, v);
  });
  ['0', 'false', 'no', 'off', 'enabled'].forEach(function (v) {
    config._build({ env: { PRIMERS_ENABLED: v }, fileConfig: { enabled: true } }).config.enabled.should.equal(false, v);
  });
});

test('malformed or empty env values are ignored (with warnings), never fatal', function () {
  const r = config._build({
    env: { PRIMERS_GLOBAL_MAX_JOBS: '-1', PRIMERS_MAX_QUEUED: '1.5', PRIMERS_JOB_STORE: 'mongo', PRIMERS_SITE_KEY: '', PRIMER3_CORE: '' },
    fileConfig: {}
  });
  r.config.check.global_max_jobs.should.equal(2);
  r.config.check.max_queued.should.equal(50);
  r.config.check.store.should.equal('redis');
  should(r.config.site_key).equal(null);
  r.config.primer3_core.should.equal('/home/olson/bin/primer3_core');
  r.warnings.length.should.equal(3);
  r.warnings.join('\n').should.match(/PRIMERS_GLOBAL_MAX_JOBS/).and.match(/PRIMERS_MAX_QUEUED/).and.match(/PRIMERS_JOB_STORE/);
  config._build({ env: { PRIMERS_GLOBAL_MAX_JOBS: '0' }, fileConfig: {} }).config.check.global_max_jobs.should.equal(0);
});

test('file config merges deeply, arrays and null replace, prototype keys are ignored', function () {
  const file = JSON.parse('{"design":{"max_concurrent":2},"assembly_overrides":{"vitis_vinifera":"Vitis_vinifera.12X"},' +
    '"__proto__":{"polluted":1},"check":{"defaults":{"include_unlikely":true}},"site_key":"from_yaml"}');
  const c = config._build({ env: {}, fileConfig: file }).config;
  c.design.max_concurrent.should.equal(2);
  c.design.max_waiting.should.equal(16);
  c.assembly_overrides.should.eql({ vitis_vinifera: 'Vitis_vinifera.12X' });
  c.check.defaults.include_unlikely.should.equal(true);
  c.check.defaults.max_product_size.should.equal(4000);
  c.site_key.should.equal('from_yaml');
  should(({}).polluted).be.undefined();
  should(c.polluted).be.undefined();
  config._build({ env: { PRIMERS_SITE_KEY: 'env_wins' }, fileConfig: file }).config.site_key.should.equal('env_wins');
  config._build({ env: {}, fileConfig: 'oops' }).warnings.length.should.equal(1);
});

test('basePath() is read from api/swagger/swagger.yaml', function () {
  config._setForTests(null);
  const doc = yaml.safeLoad(fs.readFileSync(path.join(REPO, 'api', 'swagger', 'swagger.yaml'), 'utf8'));
  config.basePath().should.equal(doc.basePath);
  config.basePath().should.equal('/sorghum_v11');

  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'primers-config-test-'));
  try {
    fs.writeFileSync(path.join(dir, 'ok.yaml'), 'swagger: "2.0"\nbasePath: /foo_v1\npaths: {}\n');
    fs.writeFileSync(path.join(dir, 'none.yaml'), 'swagger: "2.0"\npaths: {}\n');
    fs.writeFileSync(path.join(dir, 'bad.yaml'), 'basePath: "no-slash"\n');
    config.readBasePath(path.join(dir, 'ok.yaml')).should.equal('/foo_v1');
    (function () { config.readBasePath(path.join(dir, 'none.yaml')); }).should.throw(/basePath/);
    (function () { config.readBasePath(path.join(dir, 'bad.yaml')); }).should.throw(/basePath/);
  } finally {
    fs.rmSync(dir, { recursive: true, force: true });
  }
});

test('siteKey(): PRIMERS_SITE_KEY || primers.site_key || "<basePath w/o slash>:<mongo db>"', function () {
  const mongo = require('gramene-mongodb-config');
  (typeof mongo.getMongoConfig).should.equal('function');
  const db = mongo.getMongoConfig().db;
  db.should.equal('sorghum11');
  const savedEnv = process.env.PRIMERS_SITE_KEY;
  delete process.env.PRIMERS_SITE_KEY;
  try {
    config._setForTests(null);
    config.siteKey().should.equal('sorghum_v11:sorghum11');
    config.setBasePath('/sorghum_v12');
    config.siteKey().should.equal('sorghum_v12:' + db);
    config.basePath().should.equal('/sorghum_v12');
    (function () { config.setBasePath('sorghum'); }).should.throw(TypeError);
    config._setForTests({ site_key: 'dev_key' });
    config.siteKey().should.equal('dev_key');
    config._setForTests(null);
    config.basePath().should.equal('/sorghum_v11'); // override dropped
    process.env.PRIMERS_SITE_KEY = 'sorghum_v11_dev';
    config._setForTests(null);
    config.siteKey().should.equal('sorghum_v11_dev');
  } finally {
    if (savedEnv === undefined) delete process.env.PRIMERS_SITE_KEY;
    else process.env.PRIMERS_SITE_KEY = savedEnv;
    config._setForTests(null);
  }
});

test('a separate process started elsewhere (like the pm2 worker) derives the same config and site_key', function () {
  const env = {};
  Object.keys(process.env).forEach(function (k) {
    if (!/^(PRIMERS_|PRIMER3_CORE$|BLASTN$|BLASTDBCMD$|NODE_CONFIG|NODE_ENV$|NODE_APP_INSTANCE$|NODE_TEST)/.test(k)) env[k] = process.env[k];
  });
  const code = "const c = require(" + JSON.stringify(path.join(REPO, 'api', 'helpers', 'primers', 'config.js')) + ");" +
    "process.stdout.write(JSON.stringify({cwd: process.cwd(), dir: process.env.NODE_CONFIG_DIR, basePath: c.basePath(), " +
    "siteKey: c.siteKey(), local: c.get().check.local_max_jobs, blastdbcmd: c.get().blastdbcmd}));";
  const out = execFileSync(process.execPath, ['--no-deprecation', '-e', code], { cwd: '/', env: env, encoding: 'utf8', timeout: 30000 });
  const got = JSON.parse(out);
  got.cwd.should.equal('/');
  got.dir.should.equal(path.join(REPO, 'config'));
  got.basePath.should.equal('/sorghum_v11');
  got.siteKey.should.equal('sorghum_v11:sorghum11');
  got.local.should.equal(2);
  got.blastdbcmd.should.equal('/home/olson/bin/blastdbcmd');
  // and it equals what this process derives from the same sources
  const savedEnv = process.env.PRIMERS_SITE_KEY;
  delete process.env.PRIMERS_SITE_KEY;
  try {
    config._setForTests(null);
    got.siteKey.should.equal(config.siteKey());
    got.basePath.should.equal(config.basePath());
  } finally {
    if (savedEnv !== undefined) process.env.PRIMERS_SITE_KEY = savedEnv;
    config._setForTests(null);
  }
});
