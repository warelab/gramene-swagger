'use strict';

// Records the Primer3 design runs of the genotyping examples (spec §7.4): runs genotyping/design.js designGenotyping
// on every case of cases.js with the real primer3_core, and writes one <case>_<orientation>_L<level>.input.txt /
// .output.txt pair per run plus index.json, keyed by the sha256 of the exact serialized input. The inputs are
// therefore exactly what the code generates; the unit tests replay them and fail on any other input.
//
//   node test/primers/fixtures/primer3/genotyping/record_genotyping.js [--bin /home/olson/bin/primer3_core]
//
// Run manually on squam, and only when a change to the records is intended. Earlier recordings are replaced.
// primer3_core runs as primer3.js runs it: -strict_tags, a private working directory and PATH=/usr/bin:/bin.

require('../../../../../api/helpers/primers/node_compat');

const fs = require('fs');
const os = require('os');
const path = require('path');
const { spawnSync } = require('child_process');

const cases = require('./cases');
const ROOT = path.resolve(__dirname, '..', '..', '..', '..', '..');
const boulder = require(path.join(ROOT, 'api/helpers/primers/boulder'));
const genotyping = require(path.join(ROOT, 'api/helpers/primers/genotyping/design'));

const SPAWN_ENV = Object.freeze({ PATH: '/usr/bin:/bin' });
const ID_RE = /_(forward|reverse)_L(\d)$/;

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
  for (const name of Object.keys(cases.CASES)) {
    const names = [];
    const primer3 = {
      run: async function (tags) {
        const input = boulder.serialize(tags);
        const m = ID_RE.exec(tags.SEQUENCE_ID);
        if (!m) throw new Error('unexpected SEQUENCE_ID ' + tags.SEQUENCE_ID);
        const runName = name + '_' + m[1] + '_L' + m[2];
        const r = spawnIn(bin, ['-strict_tags'], input);
        if (r.status !== 0) throw new Error(runName + ': primer3_core exited with ' + r.status + ': ' + String(r.stderr).slice(0, 500));
        fs.writeFileSync(path.join(cases.DIR, runName + '.input.txt'), input);
        fs.writeFileSync(path.join(cases.DIR, runName + '.output.txt'), r.stdout);
        const hash = cases.sha256(input);
        if (!Object.prototype.hasOwnProperty.call(index.runs, hash)) {
          index.runs[hash] = { name: runName, input: runName + '.input.txt', output: runName + '.output.txt' };
        }
        names.push(runName);
        return cases.resultOf(r.stdout);
      },
      version: async function () { return version; }
    };
    const res = await genotyping.designGenotyping(cases.body(name), cases.deps(name, { primer3: primer3 }));
    index.cases[name] = names;
    console.log(name + ': ' + names.join(', ') + ' | ' + ['forward', 'reverse'].map(function (o) {
      const x = res.orientations[o];
      return o + ' ' + x.status + ' L' + x.relaxation_level + ' ' + x.attempts.map(function (a) {
        return 'returned ' + a.pairs_returned + ' not_scored ' + a.not_scored + ' sets ' + a.sets;
      }).join(' / ');
    }).join(' | '));
  }
  fs.writeFileSync(cases.INDEX, JSON.stringify(index, null, 1) + '\n');
  console.log('wrote ' + Object.keys(index.runs).length + ' runs to ' + cases.INDEX);
}

main().catch(function (err) {
  console.error(err && err.stack ? err.stack : err);
  process.exit(1);
});
