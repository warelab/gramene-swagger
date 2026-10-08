'use strict';

const { describe, it } = require('node:test');
const should = require('should');

const sensitivity = require('../../../api/helpers/primers/check/sensitivity');

describe('sensitivity.guaranteedMaxMismatches', () => {
  it('matches the spec examples (§B.3)', () => {
    should(sensitivity.guaranteedMaxMismatches(20, 5)).equal(3);
    should(sensitivity.guaranteedMaxMismatches(20, 6)).equal(2);
    should(sensitivity.guaranteedMaxMismatches(22, 6)).equal(2);
    should(sensitivity.guaranteedMaxMismatches(25, 5)).equal(4);
  });

  it('agrees with the measured word-size behaviour', () => {
    // m20_7_14 (2 mm, runs 6/6/6) is found at ws6 but missed at ws7.
    should(sensitivity.guaranteedMaxMismatches(20, 6)).be.aboveOrEqual(2);
    should(sensitivity.guaranteedMaxMismatches(20, 7)).be.below(2);
    // m22_6_12_18 (3 mm) is found at ws5.
    should(sensitivity.guaranteedMaxMismatches(22, 5)).be.aboveOrEqual(3);
  });

  it('is the largest k satisfying both bounds', () => {
    for (let L = 15; L <= 36; L++) {
      for (let W = 4; W <= 8; W++) {
        const k = sensitivity.guaranteedMaxMismatches(L, W);
        const ok = (j) => Math.ceil((L - j) / (j + 1)) >= W && L - 2 * j >= 12;
        should(ok(k)).be.true();
        should(ok(k + 1)).be.false();
      }
    }
  });

  it('never increases with word size and never decreases with length', () => {
    for (let L = 15; L <= 36; L++) {
      for (let W = 4; W < 8; W++) {
        should(sensitivity.guaranteedMaxMismatches(L, W + 1)).be.belowOrEqual(sensitivity.guaranteedMaxMismatches(L, W));
      }
    }
    for (let W = 4; W <= 8; W++) {
      for (let L = 15; L < 36; L++) {
        should(sensitivity.guaranteedMaxMismatches(L + 1, W)).be.aboveOrEqual(sensitivity.guaranteedMaxMismatches(L, W));
      }
    }
  });

  it('returns 0 when no guarantee holds and rejects bad input', () => {
    should(sensitivity.guaranteedMaxMismatches(10, 5)).equal(0);
    should(sensitivity.guaranteedMaxMismatches(4, 5)).equal(0);
    should(() => sensitivity.guaranteedMaxMismatches(0, 5)).throw(TypeError);
    should(() => sensitivity.guaranteedMaxMismatches(20, 0)).throw(TypeError);
    should(() => sensitivity.guaranteedMaxMismatches(20.5, 5)).throw(TypeError);
  });
});

describe('sensitivity.primerSensitivity and note', () => {
  it('builds the results.primers[].sensitivity block', () => {
    should(sensitivity.primerSensitivity(22, { referenceWordSize: 5, pangenomeWordSize: 6 })).eql({
      reference: { word_size: 5, guaranteed_max_mismatches: 3 },
      pangenome: { word_size: 6, guaranteed_max_mismatches: 2 }
    });
    should(sensitivity.primerSensitivity(20, { referenceWordSize: 5 })).eql({
      reference: { word_size: 5, guaranteed_max_mismatches: 3 }
    });
  });

  it('carries the §B.12 sensitivity note, including the indel caveat', () => {
    const note = sensitivity.sensitivityNote();
    should(note).startWith(
      'Sites are detected only if they contain an exact match of at least the word size and an ungapped score of about 11 or more. ');
    should(note).match(/indels near the middle of primers shorter than about 22 nt may be missed/);
    should(note).endWith('See primers[].sensitivity.');
    should(sensitivity.SENSITIVITY_NOTE).equal(note);
    // Defaults: cap 3, reference word size 5, 20-22 nt primers, no pan-genome search.
    should(note).match(/A product amplifies only when each primer has at most 3 mismatches \(max_amplifying_mismatches\)\./);
    should(note).match(/The reference search \(word size 5\) finds every site with up to 3 mismatches for these 20-22 nt primers, which covers the cap\./);
    should(note).not.match(/pan-genome/);
  });

  it('is request-specific: the pan-genome sentence says when its guarantee is below the cap', () => {
    for (const L of [20, 21, 22]) {
      should(sensitivity.guaranteedMaxMismatches(L, 6)).equal(2);
      should(sensitivity.guaranteedMaxMismatches(L, 5)).equal(3);
    }
    should(require('../../../api/helpers/primers/check/classify').DEFAULT_PARAMS.max_amplifying_mismatches).equal(3);
    const pan = sensitivity.sensitivityNote({ maxAmplifyingMismatches: 3, referenceWordSize: 5, pangenomeWordSize: 6, primerLengths: [22, 20, 21] });
    should(pan).match(/The reference search \(word size 5\) finds every site with up to 3 mismatches for these 20-22 nt primers, which covers the cap\./);
    should(pan).match(/The pan-genome search \(word size 6\) finds every site with up to 2 mismatches for these 20-22 nt primers, so products with 3 mismatches per primer can be missed in other genomes\./);
    should(pan).endWith('See primers[].sensitivity.');

    // P3 (22 and 24 nt): the pan-genome guarantee reaches the cap only for the longer primer.
    should(sensitivity.guaranteedMaxMismatches(24, 6)).equal(3);
    const p3 = sensitivity.sensitivityNote({ maxAmplifyingMismatches: 3, referenceWordSize: 5, pangenomeWordSize: 6, primerLengths: [24, 22] });
    should(p3).match(/pan-genome search \(word size 6\) finds every site with up to 2-3 mismatches for these 22-24 nt primers, so products with 3 mismatches per primer can be missed in other genomes for the shorter primers\./);

    // Other caps: 1 is covered everywhere; 5 is beyond the reference guarantee for 20 nt primers.
    const one = sensitivity.sensitivityNote({ maxAmplifyingMismatches: 1, referenceWordSize: 5, pangenomeWordSize: 6, primerLengths: [20] });
    should(one).match(/at most 1 mismatch \(max_amplifying_mismatches\)/);
    should(one).match(/word size 6\) finds every site with up to 2 mismatches for these 20 nt primers, which covers the cap\./);
    const five = sensitivity.sensitivityNote({ maxAmplifyingMismatches: 5, referenceWordSize: 5, primerLengths: [20] });
    should(five).match(/reference search \(word size 5\) finds every site with up to 3 mismatches for these 20 nt primers, so products with 4-5 mismatches per primer can be missed\./);
    should(five).not.match(/pan-genome/);
  });
});
