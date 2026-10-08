'use strict';

require('../../../api/helpers/primers/node_compat');

const test = require('node:test');
const should = require('should');
const fs = require('fs');
const os = require('os');
const path = require('path');

const primer3 = require('../../../api/helpers/primers/primer3');
const boulder = require('../../../api/helpers/primers/boulder');
const config = require('../../../api/helpers/primers/config');
const { PrimerHttpError } = require('../../../api/helpers/primers/errors');

const REAL_BIN = config.get().primer3_core;
const HAVE_REAL = fs.existsSync(REAL_BIN);
const FX = path.join(__dirname, '..', 'fixtures', 'primer3');
const read = function (name) { return fs.readFileSync(path.join(FX, name), 'utf8'); };

const TMP = fs.mkdtempSync(path.join(os.tmpdir(), 'primers-primer3-test-'));
const pidsToReap = [];

function script(name, body, mode) {
  const p = path.join(TMP, name);
  fs.writeFileSync(p, '#!/bin/sh\n' + body + '\n');
  fs.chmodSync(p, mode === undefined ? 0o755 : mode);
  return p;
}

function quietLog() {
  const lines = [];
  return { lines: lines, error: function (m) { lines.push(String(m)); } };
}

async function rejection(p) {
  try {
    await p;
  } catch (e) {
    return e;
  }
  throw new Error('expected a rejection');
}

function alive(pid) {
  try {
    process.kill(pid, 0);
    return true;
  } catch (e) {
    return e.code !== 'ESRCH';
  }
}

async function waitDead(pid, ms) {
  const until = Date.now() + ms;
  while (alive(pid) && Date.now() < until) await new Promise(function (r) { setTimeout(r, 20); });
  return !alive(pid);
}

test.after(function () {
  pidsToReap.forEach(function (pid) {
    try { process.kill(pid, 'SIGKILL'); } catch (e) { /* gone */ }
  });
  config._setForTests(null);
  primer3._clearVersionCache();
  fs.rmSync(TMP, { recursive: true, force: true });
});

const RECORD = 'SEQUENCE_ID=t1\nPRIMER_TASK=generic\n=\n';

test('spawns with -strict_tags, no shell, env {PATH:/usr/bin:/bin}, a private 0700 cwd under tmpDir (removed afterwards), record on stdin', async function () {
  const bin = script('echo.sh', [
    '[ "$1" = "-strict_tags" ] || exit 9',
    '[ "$#" = "1" ] || exit 8',
    'echo "ARGV_OK=1"',
    'echo "PATH_SEEN=$PATH"',
    'echo "HOME_SEEN=${HOME:-unset}"',
    'echo "CWD_SEEN=$(pwd -P)"',
    'echo "CWD_MODE=$(stat -c %a .)"',
    'cat'
  ].join('\n'));
  const r = await primer3.run(RECORD, { bin: bin, timeoutMs: 5000, log: quietLog() });
  r.tags.ARGV_OK.should.equal('1');
  r.tags.PATH_SEEN.should.equal('/usr/bin:/bin');
  r.tags.HOME_SEEN.should.equal('unset');
  // never the shared os.tmpdir() itself: a private mkdtemp directory under it, gone after the run
  r.tags.CWD_SEEN.should.not.equal(fs.realpathSync(os.tmpdir()));
  path.dirname(r.tags.CWD_SEEN).should.equal(fs.realpathSync(os.tmpdir()));
  path.basename(r.tags.CWD_SEEN).indexOf('primers-primer3-').should.equal(0);
  r.tags.CWD_MODE.should.equal('700');
  fs.existsSync(r.tags.CWD_SEEN).should.equal(false);
  const base = fs.mkdtempSync(path.join(TMP, 'base-'));
  const rb = await primer3.run(RECORD, { bin: bin, timeoutMs: 5000, tmpDir: base });
  path.dirname(rb.tags.CWD_SEEN).should.equal(fs.realpathSync(base));
  fs.readdirSync(base).should.eql([]);
  r.tags.SEQUENCE_ID.should.equal('t1');
  should(r.error).equal(null);
  r.exitCode.should.equal(0);
  // a tags object is serialized (and guarded) first
  const r2 = await primer3.run({ SEQUENCE_ID: 't2', PRIMER_TASK: 'generic' }, { bin: bin, timeoutMs: 5000 });
  r2.tags.SEQUENCE_ID.should.equal('t2');
  const bad = await rejection(primer3.run({ SEQUENCE_ID: 'x\nPRIMER_TASK=evil' }, { bin: bin }));
  bad.code.should.equal('BOULDER_INVALID_TAG');
});

test('a working directory that cannot be created -> 500 PRIMER3_FAILED (logged), nothing spawned', async function () {
  const counter = path.join(TMP, 'nowd.count');
  const counting = script('nowd.sh', 'echo x >> ' + counter + '\ncat');
  const log = quietLog();
  const e = await rejection(primer3.run(RECORD, { bin: counting, timeoutMs: 5000, tmpDir: path.join(TMP, 'no-such-base'), log: log }));
  e.status.should.equal(500);
  e.code.should.equal('PRIMER3_FAILED');
  log.lines.join('\n').should.match(/working directory could not be created: ENOENT/);
  fs.existsSync(counter).should.equal(false);
});

test('missing or non-executable binary -> 503 PRIMER3_UNAVAILABLE', async function () {
  let e = await rejection(primer3.run(RECORD, { bin: path.join(TMP, 'does-not-exist'), timeoutMs: 5000, log: quietLog() }));
  e.should.be.instanceOf(PrimerHttpError);
  e.status.should.equal(503);
  e.code.should.equal('PRIMER3_UNAVAILABLE');
  e.details.retry_after_s.should.be.a.Number();
  JSON.stringify(e).should.not.match(/primers-primer3-test/);

  const noexec = script('noexec.sh', 'exit 0', 0o644);
  e = await rejection(primer3.run(RECORD, { bin: noexec, timeoutMs: 5000, log: quietLog() }));
  e.status.should.equal(503);
  e.code.should.equal('PRIMER3_UNAVAILABLE');

  // the configured binary is used by default
  config._setForTests({ primer3_core: path.join(TMP, 'also-missing') });
  try {
    e = await rejection(primer3.run(RECORD, { timeoutMs: 5000, log: quietLog() }));
    e.code.should.equal('PRIMER3_UNAVAILABLE');
  } finally {
    config._setForTests(null);
  }
});

test('non-zero exit -> 500 PRIMER3_FAILED; exit code and redacted stderr are logged only', async function () {
  const bin = script('fail.sh', 'cat >/dev/null\necho "fatal: cannot open /home/olson/secret/tables.ds" >&2\nexit 3');
  const log = quietLog();
  const e = await rejection(primer3.run(RECORD, { bin: bin, timeoutMs: 5000, log: log }));
  e.status.should.equal(500);
  e.code.should.equal('PRIMER3_FAILED');
  JSON.stringify(e).should.not.match(/secret|tables|fatal/);
  log.lines.join('\n').should.match(/code 3/);
  log.lines.join('\n').should.match(/tables\.ds/);
  log.lines.join('\n').should.not.match(/\/home\/olson/);
});

test('exit 0 without a record terminator -> 500 PRIMER3_FAILED', async function () {
  const bin = script('noterm.sh', 'cat >/dev/null\necho PRIMER_PAIR_NUM_RETURNED=0');
  const e = await rejection(primer3.run(RECORD, { bin: bin, timeoutMs: 5000, log: quietLog() }));
  e.code.should.equal('PRIMER3_FAILED');
  e.message.should.match(/incomplete/);
});

test('timeout -> SIGKILL + 504 DEADLINE_EXCEEDED', async function () {
  const pidFile = path.join(TMP, 'sleep.pid');
  const bin = script('sleep.sh', 'echo $$ > ' + pidFile + '\nexec sleep 30');
  const t0 = Date.now();
  const e = await rejection(primer3.run(RECORD, { bin: bin, timeoutMs: 300, log: quietLog() }));
  const pid = Number(fs.readFileSync(pidFile, 'utf8'));
  pidsToReap.push(pid);
  e.status.should.equal(504);
  e.code.should.equal('DEADLINE_EXCEEDED');
  (Date.now() - t0).should.be.below(5000);
  (await waitDead(pid, 2000)).should.equal(true);
});

test('a grandchild holding the pipes open does not delay the timeout', async function () {
  const pidFile = path.join(TMP, 'grandchild.pid');
  // background grandchild keeps stdout open; the foreground is exec'd so SIGKILL leaves no orphan
  const bin = script('grandchild.sh', 'sleep 30 &\necho $! > ' + pidFile + '\nexec sleep 30');
  const t0 = Date.now();
  const e = await rejection(primer3.run(RECORD, { bin: bin, timeoutMs: 300, log: quietLog() }));
  const gpid = Number(fs.readFileSync(pidFile, 'utf8'));
  pidsToReap.push(gpid);
  e.code.should.equal('DEADLINE_EXCEEDED');
  (Date.now() - t0).should.be.below(5000);
  try { process.kill(gpid, 'SIGKILL'); } catch (x) { /* gone */ }
  (await waitDead(gpid, 2000)).should.equal(true);
});

test('stdout over the cap -> SIGKILL + 500 PRIMER3_FAILED', async function () {
  const bin = script('yes.sh', 'exec yes PRIMER_LEFT_0_SEQUENCE=ACGTACGTACGTACGTACGT');
  const log = quietLog();
  const e = await rejection(primer3.run(RECORD, { bin: bin, timeoutMs: 5000, maxStdoutBytes: 10000, log: log }));
  e.status.should.equal(500);
  e.code.should.equal('PRIMER3_FAILED');
  e.details.limit_bytes.should.equal(10000);
  log.lines.join('\n').should.match(/exceeded 10000 bytes/);
});

test('abort signal kills the child; a pre-aborted signal or no time left never spawns', async function () {
  const pidFile = path.join(TMP, 'abort.pid');
  const bin = script('abort.sh', 'echo $$ > ' + pidFile + '\nexec sleep 30');
  const ac = new AbortController();
  const reason = new PrimerHttpError(504, 'DEADLINE_EXCEEDED', 'design deadline exceeded', { deadline_ms: 45000 });
  setTimeout(function () { ac.abort(reason); }, 200);
  const e = await rejection(primer3.run(RECORD, { bin: bin, timeoutMs: 10000, signal: ac.signal, log: quietLog() }));
  e.should.equal(reason);
  const pid = Number(fs.readFileSync(pidFile, 'utf8'));
  pidsToReap.push(pid);
  (await waitDead(pid, 2000)).should.equal(true);

  const counter = path.join(TMP, 'spawned.count');
  const counting = script('counting.sh', 'echo x >> ' + counter + '\ncat');
  const ac2 = new AbortController();
  ac2.abort();
  (await rejection(primer3.run(RECORD, { bin: counting, signal: ac2.signal }))).code.should.equal('DEADLINE_EXCEEDED');
  (await rejection(primer3.run(RECORD, { bin: counting, timeoutMs: 0 }))).code.should.equal('DEADLINE_EXCEEDED');
  fs.existsSync(counter).should.equal(false);
});

test('real primer3_core reproduces the captured spliced-cDNA fixture', { skip: !HAVE_REAL && 'primer3_core not installed' }, async function () {
  const want = boulder.parse(read('transcript_junction.output.txt'));
  const r = await primer3.run(read('transcript_junction.input.txt'), { timeoutMs: 30000 });
  r.exitCode.should.equal(0);
  should(r.error).equal(null);
  should(r.warning).equal(null);
  r.tags.should.eql(want);
  boulder.extractPairs(r.tags).length.should.equal(5);
  // same result from a tags object
  const r2 = await primer3.run(boulder.parse(read('transcript_junction.input.txt')), { timeoutMs: 30000 });
  r2.tags.should.eql(want);
});

test('real primer3_core: global PRIMER_ERROR (exit 252) is surfaced, not thrown', { skip: !HAVE_REAL && 'primer3_core not installed' }, async function () {
  const r = await primer3.run(read('error_junction_overlap.input.txt'), { timeoutMs: 30000 });
  r.exitCode.should.equal(252);
  r.error.should.equal('PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION > PRIMER_MAX_SIZE / 2');
  r.stderrTail.should.match(/^primer3_core: PRIMER_MIN_5_PRIME/); // path reduced to basename
  boulder.extractPairs(r.tags).should.eql([]);
});

test('real primer3_core: 50 kb template with 20 pairs well under the 10 s gate', { skip: !HAVE_REAL && 'primer3_core not installed' }, async function () {
  const r = await primer3.run(read('genomic_50kb.input.txt'), { timeoutMs: 30000 });
  r.elapsedMs.should.be.below(10000);
  const pairs = boulder.extractPairs(r.tags);
  pairs.length.should.equal(20);
  Object.keys(r.tags).filter(function (k) { return /^PRIMER_PAIR_\d+_PRODUCT_TM_OLIGO_TM_DIFF$/.test(k); }).length.should.equal(20);
  r.tags.should.eql(boulder.parse(read('genomic_50kb.output.txt')));
});

test('version() parses -about, caches success, and does not cache failure', async function () {
  primer3._clearVersionCache();
  const counter = path.join(TMP, 'about.count');
  const about = script('about.sh', 'echo x >> ' + counter + '\n[ "$1" = "-about" ] || exit 4\necho "libprimer3 release 9.9.9"');
  const p1 = primer3.version({ bin: about });
  (await p1).should.equal('9.9.9');
  primer3.version({ bin: about }).should.equal(p1);
  (await primer3.version({ bin: about })).should.equal('9.9.9');
  fs.readFileSync(counter, 'utf8').split('\n').filter(Boolean).length.should.equal(1);

  const failCounter = path.join(TMP, 'aboutfail.count');
  const broken = script('aboutfail.sh', 'echo x >> ' + failCounter + '\nexit 1');
  (await rejection(primer3.version({ bin: broken }))).code.should.equal('PRIMER3_UNAVAILABLE');
  (await rejection(primer3.version({ bin: broken }))).code.should.equal('PRIMER3_UNAVAILABLE');
  fs.readFileSync(failCounter, 'utf8').split('\n').filter(Boolean).length.should.equal(2);

  // -about also runs in a private working directory, removed afterwards
  const aboutCwd = path.join(TMP, 'about.cwd');
  const aboutPwd = script('aboutpwd.sh', 'pwd -P > ' + aboutCwd + '\necho "libprimer3 release 1.2.3"');
  (await primer3.version({ bin: aboutPwd })).should.equal('1.2.3');
  const seen = fs.readFileSync(aboutCwd, 'utf8').trim();
  seen.should.not.equal(fs.realpathSync(os.tmpdir()));
  path.dirname(seen).should.equal(fs.realpathSync(os.tmpdir()));
  fs.existsSync(seen).should.equal(false);

  const hang = script('abouthang.sh', 'exec sleep 30');
  const t0 = Date.now();
  (await rejection(primer3.version({ bin: hang, timeoutMs: 200 }))).code.should.equal('PRIMER3_UNAVAILABLE');
  (Date.now() - t0).should.be.below(5000);
  (await rejection(primer3.version({ bin: path.join(TMP, 'nope') }))).status.should.equal(503);
});

test('version() of the installed primer3_core is 2.6.1', { skip: !HAVE_REAL && 'primer3_core not installed' }, async function () {
  primer3._clearVersionCache();
  (await primer3.version()).should.equal('2.6.1');
});
