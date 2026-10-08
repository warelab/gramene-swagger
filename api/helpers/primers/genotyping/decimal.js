'use strict';

// Exact decimal arithmetic for the genotyping design (spec §2.1, the rounding rule).
//
// Every input is an exact decimal as a tool prints it (Primer3 and ntthal print at most 6 decimal places), a short
// config decimal, or a ratio of integers. Values are carried as BigInt counts of 10^-6 ("micro units"); a derived
// number (a difference, a maximum, a sum such as a score) is computed on those integers, never in binary floating
// point, and rounded once, half away from zero. A ratio is rounded directly from its integer numerator and
// denominator. JS numbers are accepted as inputs: String(n) of a number parsed from a printed decimal is that decimal.

const SCALE = 6;
const UNIT = 1000000n;
const DECIMAL_RE = /^([-+]?)(\d+)(?:\.(\d+))?$/;

// micro(x) -> BigInt micro units of a finite decimal number or string ('57.099', 4.55, '-1.0'). Digits beyond the
// sixth decimal place must be zeros. Throws TypeError otherwise (an exponent, NaN, a non-number).
function micro(x) {
  const s = typeof x === 'number' ? (Number.isFinite(x) ? String(x) : '') : typeof x === 'string' ? x.trim() : '';
  const m = DECIMAL_RE.exec(s);
  if (!m || (m[3] && /[1-9]/.test(m[3].slice(SCALE)))) {
    throw new TypeError('not an exact decimal with at most ' + SCALE + ' decimal places: ' + String(x));
  }
  const v = BigInt(m[2]) * UNIT + BigInt((m[3] || '').slice(0, SCALE).padEnd(SCALE, '0'));
  return m[1] === '-' ? -v : v;
}

function fromInt(n) {
  if (!Number.isSafeInteger(n)) throw new TypeError('not a safe integer: ' + String(n));
  return BigInt(n) * UNIT;
}

function abs(v) {
  return v < 0n ? -v : v;
}

function max() {
  let best = null;
  for (let i = 0; i < arguments.length; i++) if (best === null || arguments[i] > best) best = arguments[i];
  return best;
}

function compare(a, b) {
  return a < b ? -1 : a > b ? 1 : 0;
}

function pow10(n) {
  return 10n ** BigInt(n);
}

// round(v, dp) -> Number: micro units rounded once, half away from zero, to dp decimal places (0-6).
function round(v, dp) {
  if (typeof v !== 'bigint') throw new TypeError('round expects micro units (BigInt)');
  if (!Number.isInteger(dp) || dp < 0 || dp > SCALE) throw new RangeError('dp must be an integer from 0 to ' + SCALE);
  const unit = pow10(SCALE - dp);
  const q = (abs(v) + unit / 2n) / unit;
  const n = Number(q) / Math.pow(10, dp);
  return v < 0n && q !== 0n ? -n : n;
}

// A tool value rounded to dp decimals: round(micro(x), dp).
function roundValue(x, dp) {
  return round(micro(x), dp);
}

// ratio(num, den, dp) -> Number: the exact ratio of two integers rounded half away from zero to dp decimal places.
function ratio(num, den, dp) {
  if (!Number.isSafeInteger(num) || !Number.isSafeInteger(den) || den <= 0) throw new TypeError('ratio expects integers with den > 0');
  const unit = pow10(dp);
  const a = BigInt(Math.abs(num)) * unit;
  const d = BigInt(den);
  const q = (2n * a + d) / (2n * d);
  const n = Number(q) / Math.pow(10, dp);
  return num < 0 && q !== 0n ? -n : n;
}

// 100 x GC / length of a sequence, rounded to 2 decimals from the exact ratio (§2.1, §2.8 gc).
function gcPercent(seq) {
  const s = String(seq);
  if (s.length === 0) throw new TypeError('empty sequence');
  let gc = 0;
  for (let i = 0; i < s.length; i++) if (s[i] === 'G' || s[i] === 'C') gc++;
  return ratio(100 * gc, s.length, 2);
}

// Exact test 100 x GC / length >= pct (a decimal), for the §4.6 GC floor.
function gcAtLeast(seq, pct) {
  const s = String(seq);
  let gc = 0;
  for (let i = 0; i < s.length; i++) if (s[i] === 'G' || s[i] === 'C') gc++;
  return fromInt(100 * gc) >= micro(pct) * BigInt(s.length);
}

// '9.179713' for 9179713n: the exact value, for logs and tests.
function format(v) {
  const a = abs(v);
  return (v < 0n ? '-' : '') + String(a / UNIT) + '.' + String(a % UNIT).padStart(SCALE, '0');
}

module.exports = {
  SCALE,
  UNIT,
  micro,
  fromInt,
  abs,
  max,
  compare,
  round,
  roundValue,
  ratio,
  gcPercent,
  gcAtLeast,
  format
};
