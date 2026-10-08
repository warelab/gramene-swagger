'use strict';

require('../../../api/helpers/primers/node_compat');

const test = require('node:test');
const should = require('should');
const fs = require('fs');
const path = require('path');

const guard = require('../../../api/helpers/primers/genotyping/guard');
const normalize = require('../../../api/helpers/primers/variation/normalize');
const cases = require('../fixtures/primer3/genotyping/cases');

// Spec §4.4 and the genotyping_guard rows of §7.3. Canonical entries come from the committed genome windows and Ensembl
// overlap fixtures (provenance in variation_normalize.test.js); template windows follow §4.3 (flank 400).

const FIX = path.join(__dirname, '..', 'fixtures');
const G1 = cases.genome();
const G2 = normalize.sequenceWindow(fs.readFileSync(path.join(FIX, 'design', 'sorghum_bicolor_1_13400-14000.plus.txt'), 'utf8').trim(),
  13400, 80884392);

function entry(key, overlap, genome) {
  const records = JSON.parse(fs.readFileSync(path.join(FIX, 'variation', 'overlap_1_' + overlap + '.json'), 'utf8'));
  const e = normalize.recordsToEntries(records, genome, { region: '1' }).entries.find(function (x) { return x.key === key; });
  should.exist(e, key);
  return e;
}

// The §4.3 template of an entry, arithmetic only: {start, length, zone (template), t(g)}.
function templateOf(e) {
  const fwd = e.discriminating.forward.position;
  const rev = e.discriminating.reverse.position;
  const start = Math.min(fwd, rev) - 400;
  const end = Math.max(fwd, rev) + 400;
  const t = function (g) { return g - start + 1; };
  return { start: start, length: end - start + 1, t: t, zone: { start: t(e.zone.start), end: t(e.zone.end) }, delta: e.vcf.alt.length - e.vcf.ref.length };
}

const RS871 = templateOf(entry('1:11109:C:A', '10709-11509', G1));
const EMS11193 = templateOf(entry('1:11193:C:T', '10793-11593', G1));
const DEL = templateOf(entry('1:11282:CA:C', '10883-11685', G1));
const INS = templateOf(entry('1:11502:C:CGT', '11102-11903', G1));
const TGG = entry('1:13735:TTGG:T', '13600-13820', G2);
const TGGT = templateOf(TGG);

test('targetTag: SEQUENCE_TARGET for the five example variants (402,10 391,10 405,10 390,10; rs5413863413 {t(13753)},10)', function () {
  guard.targetTag(RS871.zone, 'forward', RS871, 10).should.equal('402,10');
  guard.targetTag(RS871.zone, 'reverse', RS871, 10).should.equal('391,10');
  guard.targetTag(EMS11193.zone, 'forward', EMS11193, 10).should.equal('402,10');
  guard.targetTag(EMS11193.zone, 'reverse', EMS11193, 10).should.equal('391,10');
  guard.targetTag(DEL.zone, 'forward', DEL, 10).should.equal('405,10');
  guard.targetTag(DEL.zone, 'reverse', DEL, 10).should.equal('390,10');
  guard.targetTag(INS.zone, 'forward', INS, 10).should.equal('403,10');
  guard.targetTag(INS.zone, 'reverse', INS, 10).should.equal('391,10');
  TGG.zone.should.eql({ start: 13735, end: 13752 });
  TGGT.start.should.equal(13338);
  guard.targetTag(TGGT.zone, 'forward', TGGT, 10).should.equal(TGGT.t(13753) + ',10');
  guard.targetTag(TGGT.zone, 'reverse', TGGT, 10).should.equal(TGGT.t(13725) + ',10');
});

test('targetTag: the gap shortens at a template edge, never below 1', function () {
  guard.targetTag({ start: 50, end: 97 }, 'forward', { length: 100 }, 10).should.equal('98,3');
  guard.targetTag({ start: 4, end: 20 }, 'reverse', { length: 100 }, 10).should.equal('1,3');
  should(guard.targetTag({ start: 50, end: 100 }, 'forward', { length: 100 }, 10)).be.null();
  should(guard.targetTag({ start: 1, end: 20 }, 'reverse', { length: 100 }, 10)).be.null();
  (function () { guard.targetTag({ start: 5, end: 4 }, 'forward', { length: 100 }, 10); }).should.throw(RangeError);
  (function () { guard.targetTag({ start: 5, end: 6 }, 'both', { length: 100 }, 10); }).should.throw(TypeError);
  (function () { guard.targetTag({ start: 5, end: 6 }, 'forward', { length: 100 }, 0); }).should.throw(RangeError);
  (function () { guard.targetTag({ start: 5, end: 106 }, 'forward', { length: 100 }, 10); }).should.throw(RangeError);
});

test('commonPrimerOk: the unguarded overlap of §4.4 (common primer at 1:11108-11131) is rejected for rs871475760', function () {
  const t = RS871.t;
  const as = { start: t(11081), end: t(11109) };
  guard.commonPrimerOk({ orientation: 'forward', as: as, common: { start: t(11108), end: t(11131) } }, RS871.zone, RS871).should.be.false();
  // The §2.9 sets: S2 forward common 1:11149-11172, S1 reverse common 1:11068-11090.
  guard.commonPrimerOk({ orientation: 'forward', as: as, common: { start: t(11149), end: t(11172) } }, RS871.zone, RS871).should.be.true();
  const rev = { start: t(11109), end: t(11132) };
  guard.commonPrimerOk({ orientation: 'reverse', as: rev, common: { start: t(11068), end: t(11090) } }, RS871.zone, { delta: 0 }).should.be.true();
  // [zone.start - 1, zone.end + 1]: ending on 11108 touches it, ending on 11107 does not.
  guard.commonPrimerOk({ orientation: 'reverse', as: rev, common: { start: t(11085), end: t(11108) } }, RS871.zone, { delta: 0 }).should.be.false();
  guard.commonPrimerOk({ orientation: 'reverse', as: rev, common: { start: t(11084), end: t(11107) } }, RS871.zone, { delta: 0 }).should.be.true();
  guard.commonPrimerOk({ orientation: 'forward', as: as, common: { start: t(11111), end: t(11134) } }, RS871.zone, { delta: 0 }).should.be.true();
  guard.commonPrimerOk({ orientation: 'forward', as: as, common: { start: t(11110), end: t(11134) } }, RS871.zone, { delta: 0 }).should.be.false();
});

test('commonPrimerOk: a common primer at 13751-13780 is rejected for rs5413863413, on the ALT haplotype even with a REF-only zone', function () {
  const t = TGGT.t;
  TGG.discriminating.forward.should.eql({ position: 13749, ref_base: 'G', alt_base: 'A', alt_maps_to: 13752 });
  const set = { orientation: 'forward', as: { start: t(13725), end: t(13749) }, common: { start: t(13751), end: t(13780) } };
  guard.commonPrimerOk(set, TGGT.zone, { delta: TGGT.delta }).should.be.false();
  // A guard built from the REF span alone ([13735, 13738]) lets the primer through rule 1; on ALT carriers the bases after the
  // deletion move by -3, so it covers the ALT primer's 3' base, and rule 2 still rejects it.
  const refOnly = { start: t(13735), end: t(13738) };
  guard.commonPrimerOk(set, refOnly, { delta: TGGT.delta }).should.be.false();
  guard.commonPrimerOk({ orientation: 'forward', as: set.as, common: { start: t(13755), end: t(13780) } }, refOnly, { delta: -3 }).should.be.true();
  guard.commonPrimerOk({ orientation: 'forward', as: set.as, common: { start: t(13755), end: t(13780) } }, TGGT.zone, { delta: -3 }).should.be.true();
});

test('commonPrimerOk: the §2.10 indel sets pass on both haplotypes; the variant template works as haplotypes', function () {
  let t = DEL.t;
  guard.commonPrimerOk({ orientation: 'forward', as: { start: t(11257), end: t(11285) }, common: { start: t(11323), end: t(11344) } },
    DEL.zone, { length: 803, alt_length: 802 }).should.be.true();
  guard.commonPrimerOk({ orientation: 'reverse', as: { start: t(11283), end: t(11307) }, common: { start: t(11179), end: t(11206) } },
    DEL.zone, { delta: -1 }).should.be.true();
  t = INS.t;
  guard.commonPrimerOk({ orientation: 'reverse', as: { start: t(11502), end: t(11529) }, common: { start: t(11458), end: t(11479) } },
    INS.zone, { length: 802, alt_length: 804 }).should.be.true();
  // A reverse allele-specific primer moves by delta on ALT; a common primer overlapping it there is rejected.
  guard.commonPrimerOk({ orientation: 'reverse', as: { start: 50, end: 70 }, common: { start: 20, end: 49 } }, { start: 48, end: 50 }, { delta: 0 })
    .should.be.false();
  (function () { guard.commonPrimerOk({ orientation: 'forward', as: { start: 1, end: 2 } }, { start: 5, end: 6 }, { delta: 0 }); })
    .should.throw(TypeError);
  (function () { guard.commonPrimerOk({ orientation: 'forward', as: { start: 1, end: 2 }, common: { start: 9, end: 20 } }, { start: 5, end: 6 }, {}); })
    .should.throw(TypeError);
});
