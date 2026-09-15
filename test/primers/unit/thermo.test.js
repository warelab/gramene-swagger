'use strict';

require('../../../api/helpers/primers/node_compat');

const test = require('node:test');
const should = require('should');
const fs = require('fs');

const thermo = require('../../../api/helpers/primers/thermo');
const decimal = require('../../../api/helpers/primers/genotyping/decimal');
const { PrimerHttpError } = require('../../../api/helpers/primers/errors');
const { revcomp } = require('../../../api/helpers/primers/sequence');
const cases = require('../fixtures/primer3/genotyping/cases');

// thermo.js (spec §4.10, §4.14, §6.1) and the thermo rows of §7.3-§7.4. ntthal itself is spawned only by the
// PRIMERS_REALDATA=1 tests; the rest inject a spawn (or run small system utilities through runBinary).

const REALDATA = process.env.PRIMERS_REALDATA === '1';
const NTTHAL = '/home/olson/primer3-2.6.1/bin/ntthal';
const FAM = 'GAAGGTGACCAAGTTCATGCT';
const HEX = 'GAAGGTCGGAGTCAACGGATT';
const S1_REF = 'ATCTTTGACTAGCGAGAAATTCAG';
const S1_ALT = 'ATCTTTGACTAGCGAGAAATTCAT';

function cfg(extra) {
  return cases.cfg(null, extra);
}

async function rejection(promise) {
  try {
    await promise;
  } catch (e) {
    return e;
  }
  throw new Error('expected a rejection');
}

// A spawn that answers `value(argv)` and records every call; opts.hold keeps calls pending until release().
function fakeSpawn(value, opts) {
  opts = opts || {};
  const calls = [];
  const pending = [];
  let active = 0;
  let maxActive = 0;
  const fn = function (bin, args, o) {
    calls.push({ bin: bin, args: args, opts: o });
    active++;
    maxActive = Math.max(maxActive, active);
    const answer = function () {
      active--;
      return { stdout: value(args) + '\n' };
    };
    if (!opts.hold) return Promise.resolve().then(answer);
    return new Promise(function (resolve) { pending.push(function () { resolve(answer()); }); });
  };
  fn.calls = calls;
  fn.maxActive = function () { return maxActive; };
  fn.releaseAll = function () { while (pending.length) pending.shift()(); };
  return fn;
}

test('argv: Primer3 salts (-mv 50 -dv 1.5 -n 0.6 -d 50), 37 °C, Tm only; effective params override the salts', function () {
  const a = thermo.argv(null, 'HAIRPIN', 'ACGT');
  a.join(' ').should.equal('-mv 50 -dv 1.5 -n 0.6 -d 50 -t 37 -r -a HAIRPIN -s1 ACGT');
  a.join(' ').should.containEql('-mv 50 -dv 1.5 -n 0.6 -d 50 -t 37 -r');
  thermo.argv({ salt_monovalent: 60, dna_conc: 250, opt_tm: 60 }, 'END1', 'AC', 'GT').join(' ')
    .should.equal('-mv 60 -dv 1.5 -n 0.6 -d 250 -t 37 -r -a END1 -s1 AC -s2 GT');
  thermo.PRIMER3_SALTS.should.eql({ salt_monovalent: 50, salt_divalent: 1.5, dntp_conc: 0.6, dna_conc: 50 });
});

test('parseOutput: "No secondary structure…" and negative values are 0; a value keeps its printed decimals (42.616970 -> 42.62)', function () {
  thermo.parseOutput('No secondary structure could be calculated\n').should.equal('0');
  thermo.parseOutput('-11.46\n').should.equal('0');
  thermo.parseOutput('-26.885942').should.equal('0');
  thermo.parseOutput('0.000000\n').should.equal('0.000000');
  thermo.parseOutput(' 42.616970\n').should.equal('42.616970');
  decimal.roundValue(thermo.parseOutput('42.616970'), 2).should.equal(42.62);
  decimal.roundValue(thermo.parseOutput('-11.46'), 2).should.equal(0);
  should(thermo.parseOutput('Tm: 42.6 degrees')).be.null();
  should(thermo.parseOutput('')).be.null();
});

test('memoization: one spawn per distinct argv; calls counts distinct requests; cross is ANY once and END1 both ways (the larger)', async function () {
  const spawn = fakeSpawn(function (args) {
    const key = args.slice(-4).join(' ');
    if (key === '-s1 ' + FAM + S1_REF + ' -s2 ' + HEX + S1_ALT) return /END1/.test(args.join(' ')) ? '-12.203375' : '6.889811';
    if (key === '-s1 ' + HEX + S1_ALT + ' -s2 ' + FAM + S1_REF) return '8.207452';
    return '42.616970';
  });
  const t = thermo.createThermo({ cfg: cfg(), spawn: spawn, log: { error: function () {} } });
  (await t.hairpin(FAM + S1_REF)).should.equal('42.616970');
  (await t.hairpin(FAM + S1_REF)).should.equal('42.616970');
  spawn.calls.should.have.length(1);
  t.calls.should.equal(1);
  const both = await Promise.all([t.selfAny(S1_REF), t.selfAny(S1_REF), t.selfEnd(S1_REF)]);
  both.should.have.length(3);
  spawn.calls.should.have.length(3);
  t.calls.should.equal(3);
  const x = await t.cross(FAM + S1_REF, HEX + S1_ALT);
  x.should.eql({ any: '6.889811', end: '8.207452' });
  spawn.calls.map(function (c) { return c.args[c.args.indexOf('-a') + 1]; }).slice(3).should.eql(['ANY', 'END1', 'END1']);
  t.calls.should.equal(6);
  spawn.calls[0].bin.should.equal(cfg().ntthal);
  spawn.calls[0].opts.timeoutMs.should.equal(5000);
});

test('pool: at most thermo_concurrency (4) ntthal processes at a time', async function () {
  const spawn = fakeSpawn(function () { return '1.000000'; }, { hold: true });
  const t = thermo.createThermo({ cfg: cfg(), spawn: spawn });
  const seqs = ['A', 'C', 'G', 'T', 'AC', 'AG', 'AT', 'CA', 'CG', 'CT'].map(function (s) { return s + 'ACGTACGTACGTACGT'; });
  const all = Promise.all(seqs.map(function (s) { return t.hairpin(s); }));
  for (let i = 0; i < 5; i++) await new Promise(function (resolve) { setImmediate(resolve); });
  spawn.calls.should.have.length(4);
  t.calls.should.equal(10);
  for (let i = 0; i < 20 && spawn.calls.length < 10; i++) {
    spawn.releaseAll();
    await new Promise(function (resolve) { setImmediate(resolve); });
  }
  spawn.releaseAll();
  (await all).should.have.length(10);
  spawn.maxActive().should.equal(4);

  const held = fakeSpawn(function () { return '0'; }, { hold: true });
  const two = thermo.createThermo({ cfg: cases.cfg({ thermo_concurrency: 2 }), spawn: held });
  const five = Promise.all(seqs.slice(0, 5).map(function (s) { return two.hairpin(s); }));
  for (let i = 0; i < 5; i++) await new Promise(function (resolve) { setImmediate(resolve); });
  held.calls.should.have.length(2);
  for (let i = 0; i < 20 && held.calls.length < 5; i++) {
    held.releaseAll();
    await new Promise(function (resolve) { setImmediate(resolve); });
  }
  held.releaseAll();
  (await five).should.eql(['0', '0', '0', '0', '0']);
  held.maxActive().should.equal(2);
});

test('invalid input throws before anything is spawned: N, lowercase, over 60 nt, missing partner', function () {
  const spawn = fakeSpawn(function () { return '0'; });
  const t = thermo.createThermo({ cfg: cfg(), spawn: spawn });
  ['ACGN' + 'ACGTACGTACGT', 'acgtacgtacgtacgt', 'A'.repeat(61), '', null].forEach(function (s) {
    (function () { t.hairpin(s); }).should.throw({ code: 'THERMO_INVALID_SEQUENCE' });
  });
  (function () { t.duplex(S1_REF, 'ACGN'); }).should.throw(TypeError);
  (function () { t.run('ANY', S1_REF); }).should.throw(TypeError);
  (function () { t.run('END2', S1_REF, S1_REF); }).should.throw(TypeError);
  (function () { t.run('HAIRPIN', S1_REF, S1_REF); }).should.throw(TypeError);
  (typeof t.hairpin('A'.repeat(60)).then).should.equal('function');
  spawn.calls.should.have.length(0); // the pool slot is taken on the next tick
});

test('errors: unparseable output is 500 THERMO_FAILED; spawn errors pass through; an aborted signal spawns nothing', async function () {
  const log = { lines: [], error: function (m) { this.lines.push(m); } };
  const garbled = thermo.createThermo({ cfg: cfg(), spawn: fakeSpawn(function () { return 'oops'; }), log: log });
  (await rejection(garbled.hairpin(S1_REF))).should.match({ status: 500, code: 'THERMO_FAILED', details: { binary: 'ntthal' } });
  log.lines.join('\n').should.match(/unparseable/);

  const missing = thermo.createThermo({ cfg: cfg(), spawn: async function () { throw new PrimerHttpError(503, 'THERMO_UNAVAILABLE', 'x', { binary: 'ntthal', retry_after_s: 60 }); } });
  (await rejection(missing.hairpin(S1_REF))).should.match({ status: 503, code: 'THERMO_UNAVAILABLE' });

  const ac = new AbortController();
  ac.abort(new PrimerHttpError(504, 'DEADLINE_EXCEEDED', 'late', {}));
  const spawn = fakeSpawn(function () { return '0'; });
  const aborted = thermo.createThermo({ cfg: cfg(), spawn: spawn, signal: ac.signal });
  (await rejection(aborted.hairpin(S1_REF))).should.match({ code: 'DEADLINE_EXCEEDED' });
  spawn.calls.should.have.length(0);

  const late = thermo.createThermo({ cfg: cfg(), spawn: spawn, deadline: { remaining: function () { return 0; } } });
  (await rejection(late.hairpin(S1_REF))).should.match({ code: 'DEADLINE_EXCEEDED' });
  const soon = thermo.createThermo({ cfg: cfg(), spawn: spawn, deadline: { remaining: function () { return 120.7; } } });
  await soon.hairpin(S1_REF);
  spawn.calls[0].opts.timeoutMs.should.equal(120);
});

test('runBinary (no ntthal): a missing binary is 503 THERMO_UNAVAILABLE, a non-zero exit 500 THERMO_FAILED, a slow one 504', async function () {
  const quiet = { error: function () {} };
  (await rejection(thermo.runBinary('/nonexistent/ntthal', ['-r'], { timeoutMs: 1000, log: quiet })))
    .should.match({ status: 503, code: 'THERMO_UNAVAILABLE', details: { binary: 'ntthal', retry_after_s: 60 } });
  if (fs.existsSync('/usr/bin/false')) {
    (await rejection(thermo.runBinary('/usr/bin/false', [], { timeoutMs: 1000, log: quiet }))).should.match({ status: 500, code: 'THERMO_FAILED' });
  }
  if (fs.existsSync('/bin/echo')) {
    (await thermo.runBinary('/bin/echo', ['42.616970'], { timeoutMs: 1000, log: quiet })).stdout.should.equal('42.616970\n');
  }
  if (fs.existsSync('/bin/sleep')) {
    (await rejection(thermo.runBinary('/bin/sleep', ['5'], { timeoutMs: 50, log: quiet }))).should.match({ status: 504, code: 'DEADLINE_EXCEEDED' });
  }
});

test('recorded fixture: every entry is an argv thermo.js builds and parses to a value', function () {
  const fixture = cases.loadThermo();
  const keys = Object.keys(fixture.calls);
  keys.length.should.be.above(400);
  keys.forEach(function (k) {
    k.should.startWith('-mv 50 -dv 1.5 -n 0.6 -d 50 -t 37 -r -a ');
    const args = k.split(' ');
    const s2 = args.indexOf('-s2') >= 0 ? args[args.indexOf('-s2') + 1] : undefined;
    thermo.argv(null, args[args.indexOf('-a') + 1], args[args.indexOf('-s1') + 1], s2).join(' ').should.equal(k);
    should(thermo.parseOutput(fixture.calls[k])).not.be.null(k);
  });
});

const realOpts = { skip: REALDATA && fs.existsSync(NTTHAL) ? false : 'set PRIMERS_REALDATA=1 (needs ' + NTTHAL + ')' };

test('real ntthal: FAM+REF hairpin 42.62, HEX+ALT 51.50, tails alone 0, FAM-REF x HEX-ALT ANY 6.89, -2 mismatch duplex 54.69', realOpts, async function () {
  const t = thermo.createThermo({ cfg: cfg() });
  const r2 = async function (p) { return decimal.roundValue(await p, 2); };
  (await r2(t.hairpin(FAM + S1_REF))).should.equal(42.62);
  (await r2(t.hairpin(HEX + S1_ALT))).should.equal(51.5);
  (await r2(t.hairpin(FAM))).should.equal(0);
  (await r2(t.hairpin(HEX))).should.equal(0);
  (await r2(t.run('ANY', FAM + S1_REF, HEX + S1_ALT))).should.equal(6.89);
  (await r2(t.duplex('ATCTTTGACTAGCGAGAAATTCGG', revcomp(S1_REF)))).should.equal(54.69);
});

test('real ntthal reproduces 50 recorded entries (§7.4 equivalence)', realOpts, async function () {
  const fixture = cases.loadThermo();
  const keys = Object.keys(fixture.calls);
  const step = Math.max(1, Math.floor(keys.length / 50));
  let checked = 0;
  for (let i = 0; i < keys.length && checked < 50; i += step) {
    const r = await thermo.runBinary(NTTHAL, keys[i].split(' '), { timeoutMs: 5000 });
    thermo.parseOutput(r.stdout).should.equal(thermo.parseOutput(fixture.calls[keys[i]]), keys[i]);
    checked++;
  }
  checked.should.equal(50);
});
