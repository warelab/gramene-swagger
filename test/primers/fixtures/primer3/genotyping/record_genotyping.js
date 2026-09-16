'use strict';

// Records the Primer3 runs and ntthal calls of the genotyping examples (spec §7.4): runs genotyping/design.js
// designGenotyping on every case of cases.js with the real primer3_core and ntthal, and writes
//   - one <case>_<orientation>_L<level>[_p<pair>_<alt|mmref|mmalt>].input.txt / .output.txt pair per Primer3 run (design
//     runs, and the check_primers scoring runs of §4.10) plus index.json, keyed by the sha256 of the exact serialized
//     input;
//   - ../../thermo/genotyping.json: every ntthal result, keyed by its argv (' '-joined).
// The inputs are therefore exactly what the code generates; the unit tests replay them and fail on any other input.
//
//   node test/primers/fixtures/primer3/genotyping/record_genotyping.js [--bin /home/olson/bin/primer3_core]
//        [--ntthal /home/olson/primer3-2.6.1/bin/ntthal]
//
// Run manually on squam, and only when a change to the records is intended. Earlier recordings are replaced.
// primer3_core runs as primer3.js runs it (-strict_tags, a private working directory, PATH=/usr/bin:/bin); ntthal runs
// through thermo.js runBinary itself.

require('../../../../../api/helpers/primers/node_compat');

const fs = require('fs');
const os = require('os');
const path = require('path');
const { spawnSync } = require('child_process');

const cases = require('./cases');
const ROOT = path.resolve(__dirname, '..', '..', '..', '..', '..');
const boulder = require(path.join(ROOT, 'api/helpers/primers/boulder'));
const genotyping = require(path.join(ROOT, 'api/helpers/primers/genotyping/design'));
const thermoModule = require(path.join(ROOT, 'api/helpers/primers/thermo'));

const SPAWN_ENV = Object.freeze({ PATH: '/usr/bin:/bin' });
const ID_RE = /_(forward|reverse)_L(\d)(?:_p(\d+)_(alt|mmref|mmalt))?$/;

function argValue(name, dflt) {
  const i = process.argv.indexOf(name);
  return i >= 0 && process.argv[i + 1] ? process.argv[i + 1] : dflt;
}

function spawnIn(bin, args, input) {
  const cwd = fs.mkdtempSync(path.join(os.tmpdir(), 'primers-record-genotyping-'));
  try {
    return spawnSync(bin, args, { input: input, cwd: cwd, env: SPAWN_ENV, encoding: 'utf8', maxBuffer: 16 * 1024 * 1024 });
  } finally {
    fs.rmSync(cwd, { recursive: true, force: true });
  }
}

function primer3Version(bin) {
  const r = spawnIn(bin, ['-about'], '');
  const m = /release\s+([0-9][^\s]*)/.exec(r.stdout || '');
  if (r.status !== 0 || !m) throw new Error(bin + ' -about failed');
  return m[1];
}

async function main() {
  const bin = argValue('--bin', '/home/olson/bin/primer3_core');
  const ntthal = argValue('--ntthal', '/home/olson/primer3-2.6.1/bin/ntthal');
  const version = primer3Version(bin);
  fs.readdirSync(cases.DIR).filter(function (f) { return /\.(input|output)\.txt$/.test(f); }).forEach(function (f) {
    fs.unlinkSync(path.join(cases.DIR, f));
  });
  const index = {
    generated_by: 'test/primers/fixtures/primer3/genotyping/record_genotyping.js',
    primer3: bin + ' -strict_tags',
    primer3_version: version,
    sequence: 'test/primers/fixtures/design/' + cases.WINDOW.file,
    cases: {},
    runs: {}
  };
  const thermoCalls = {};
  for (const name of Object.keys(cases.CASES)) {
    const names = [];
    const primer3 = {
      run: async function (tags) {
        const input = boulder.serialize(tags);
        const m = ID_RE.exec(tags.SEQUENCE_ID);
        if (!m) throw new Error('unexpected SEQUENCE_ID ' + tags.SEQUENCE_ID);
        const runName = name + '_' + m[1] + '_L' + m[2] + (m[3] !== undefined ? '_p' + m[3] + '_' + m[4] : '');
        const hash = cases.sha256(input);
        if (!Object.prototype.hasOwnProperty.call(index.runs, hash)) {
          const r = spawnIn(bin, ['-strict_tags'], input);
          if (r.status !== 0 && !/PRIMER_ERROR=/.test(r.stdout || '')) {
            throw new Error(runName + ': primer3_core exited with ' + r.status + ': ' + String(r.stderr).slice(0, 500));
          }
          fs.writeFileSync(path.join(cases.DIR, runName + '.input.txt'), input);
          fs.writeFileSync(path.join(cases.DIR, runName + '.output.txt'), r.stdout);
          index.runs[hash] = { name: runName, input: runName + '.input.txt', output: runName + '.output.txt' };
        }
        names.push(index.runs[hash].name);
        return cases.resultOf(cases.readRecording(index.runs[hash].output));
      },
      version: async function () { return version; }
    };
    const config = cases.cfg(null, { ntthal: ntthal });
    const thermo = thermoModule.createThermo({
      cfg: config,
      log: console,
      spawn: async function (b, args, opts) {
        const r = await thermoModule.runBinary(b, args, opts);
        thermoCalls[args.join(' ')] = String(r.stdout).trim();
        return r;
      }
    });
    const res = await genotyping.designGenotyping(cases.body(name), cases.deps(name, { cfg: config, primer3: primer3, thermo: thermo }));
    index.cases[name] = names;
    console.log(name + ': ' + names.length + ' Primer3 runs, ' + res.budget.thermo_calls + ' ntthal calls | ' + ['forward', 'reverse'].map(function (o) {
      const x = res.orientations[o];
      return o + ' ' + x.status + ' L' + x.relaxation_level + ' ' + x.attempts.map(function (a) {
        return 'returned ' + a.pairs_returned + ' not_scored ' + a.not_scored + ' sets ' + a.sets;
      }).join(' / ');
    }).join(' | ') + ' | sets ' + res.sets.map(function (s) { return s.id + ' ' + s.key + ' ' + s.score; }).join(', '));
  }
  fs.writeFileSync(cases.INDEX, JSON.stringify(index, null, 1) + '\n');
  const sorted = {};
  Object.keys(thermoCalls).sort().forEach(function (k) { sorted[k] = thermoCalls[k]; });
  fs.mkdirSync(path.dirname(cases.THERMO), { recursive: true });
  fs.writeFileSync(cases.THERMO, JSON.stringify({
    generated_by: 'test/primers/fixtures/primer3/genotyping/record_genotyping.js',
    ntthal: ntthal,
    note: 'ntthal stdout (trimmed) keyed by the argv thermo.js builds; Primer3 salts, -t 37 -r',
    calls: sorted
  }, null, 1) + '\n');
  console.log('wrote ' + Object.keys(index.runs).length + ' Primer3 runs to ' + cases.INDEX + ' and ' + Object.keys(sorted).length +
    ' ntthal calls to ' + cases.THERMO);
}

main().catch(function (err) {
  console.error(err && err.stack ? err.stack : err);
  process.exit(1);
});
