'use strict';

require('../../../api/helpers/primers/node_compat');

const test = require('node:test');
const should = require('should');

const mismatch = require('../../../api/helpers/primers/genotyping/mismatch');
const sets = require('../../../api/helpers/primers/genotyping/sets');
const { revcomp } = require('../../../api/helpers/primers/sequence');

// Spec §4.11 and the genotyping_mismatch rows of §7.3: the literal Little 1995 Table 9.8.1 (32 cells), the Table 9.8.3
// classes, and the pinned deliberate-mismatch primers of rs871475760 and tmp_1_11193_C_T at -2 and -3. Pure.

test('Table 9.8.1: all 32 cells, literally', function () {
  const want = {
    AA: { A: 'A', G: 'G', C: 'A', T: 'G' },
    GG: { A: 'A', G: 'G', C: 'A', T: 'G' },
    AG: { A: 'C', G: 'T', C: 'A', T: 'G' },
    TT: { A: 'C', G: 'T', C: 'A', T: 'G' },
    CT: { A: 'C', G: 'T', C: 'A', T: 'G' },
    CC: { A: 'C', G: 'T', C: 'A', T: 'G' },
    AC: { A: 'G', G: 'A', C: 'C', T: 'T' },
    GT: { A: 'G', G: 'A', C: 'T', T: 'C' }
  };
  mismatch.TABLE_9_8_1.should.eql(want);
  let cells = 0;
  Object.keys(want).forEach(function (row) {
    Object.keys(want[row]).forEach(function (col) {
      cells++;
      // every cell is a real mismatch against the template base
      mismatch.TABLE_9_8_1[row][col].should.not.equal(mismatch.COMPLEMENT[col], row + col);
    });
  });
  cells.should.equal(32);
  Object.isFrozen(mismatch.TABLE_9_8_1.GT).should.be.true();
  mismatch.CLASSES.should.eql({ AG: 'max', CT: 'max', TT: 'max', CC: 'strong', AA: 'medium', GG: 'medium', AC: 'weak', GT: 'weak' });
});

test('terminal pairs: unordered, alphabetical; the template base is the complement for a forward primer and the plus base for a reverse one', function () {
  mismatch.terminalPair('G', 'A').should.equal('AG');
  mismatch.terminalPair('T', 'C').should.equal('CT');
  mismatch.templateBase('forward', 'A').should.equal('T');
  mismatch.templateBase('reverse', 'A').should.equal('A');
  // rs871475760 C/A reverse: as_ref ...CAG against A -> AG max; as_alt ...CAT against C -> CT max
  mismatch.terminalClass('ATCTTTGACTAGCGAGAAATTCAG', 'reverse', 'A').should.equal('max');
  mismatch.terminalClass('ATCTTTGACTAGCGAGAAATTCAT', 'reverse', 'C').should.equal('max');
  // tmp_1_11193_C_T C/T sites are weak + weak in both orientations (§4.12)
  mismatch.terminalClass('ACAAAGATAGATAACAAAAATAGCTCTC', 'forward', 'T').should.equal('weak');
  mismatch.terminalClass('ACAAAGATAGATAACAAAAATAGCTCTT', 'forward', 'C').should.equal('weak');
  mismatch.terminalClass('ACAACTTTTTAATATATTGTGTATACTCTAG', 'reverse', 'T').should.equal('weak');
  (function () { mismatch.terminalPair('N', 'A'); }).should.throw(TypeError);
  (function () { mismatch.templateBase('up', 'A'); }).should.throw(TypeError);
});

function mm(primer, orientation, otherPlus, k) {
  const block = mismatch.chooseMismatch(primer, orientation, otherPlus, k);
  return { block: block, seq: mismatch.applyMismatch(primer, block) };
}

test('rs871475760 at -2, both orientations and both alleles (§4.11 table, §2.10(d) blocks)', function () {
  const rr = mm('ATCTTTGACTAGCGAGAAATTCAG', 'reverse', 'A', 2);
  rr.seq.should.equal('ATCTTTGACTAGCGAGAAATTCGG');
  rr.block.should.eql({ position: 2, original_base: 'A', new_base: 'G', template_base: 'T', terminal_pair: 'AG', terminal_mismatch_class: 'max',
    added_mismatch_class: 'weak', source: 'Little 1995 Table 9.8.1' });
  const ra = mm('ATCTTTGACTAGCGAGAAATTCAT', 'reverse', 'C', 2);
  ra.seq.should.equal('ATCTTTGACTAGCGAGAAATTCGT');
  ra.block.should.eql({ position: 2, original_base: 'A', new_base: 'G', template_base: 'T', terminal_pair: 'CT', terminal_mismatch_class: 'max',
    added_mismatch_class: 'weak', source: 'Little 1995 Table 9.8.1' });
  const fr = mm('TGGTTATCCGAATATAGTCATACTCTATTC', 'forward', 'A', 2);
  fr.seq.should.equal('TGGTTATCCGAATATAGTCATACTCTATCC');
  [fr.block.original_base, fr.block.new_base].should.eql(['T', 'C']);
  const fa = mm('TGGTTATCCGAATATAGTCATACTCTATTA', 'forward', 'C', 2);
  fa.seq.should.equal('TGGTTATCCGAATATAGTCATACTCTATCA');
  [fa.block.original_base, fa.block.new_base].should.eql(['T', 'C']);
  mismatch.note(rr.block).should.equal('deliberate mismatch −2 A→G (Little 1995 Table 9.8.1)');
});

test('tmp_1_11193_C_T at -2: reverse REF A->C and ALT A->T (where a class shortcut gives C for both); forward both T->G', function () {
  const rr = mm('ACAACTTTTTAATATATTGTGTATACTCTAG', 'reverse', 'T', 2);
  const ra = mm('ACAACTTTTTAATATATTGTGTATACTCTAA', 'reverse', 'C', 2);
  rr.seq.should.equal('ACAACTTTTTAATATATTGTGTATACTCTCG');
  ra.seq.should.equal('ACAACTTTTTAATATATTGTGTATACTCTTA');
  [rr.block.terminal_pair, rr.block.new_base, ra.block.terminal_pair, ra.block.new_base].should.eql(['GT', 'C', 'AC', 'T']);
  rr.block.new_base.should.not.equal(ra.block.new_base);
  const fr = mm('ACAAAGATAGATAACAAAAATAGCTCTC', 'forward', 'T', 2);
  const fa = mm('ACAAAGATAGATAACAAAAATAGCTCTT', 'forward', 'C', 2);
  fr.seq.should.equal('ACAAAGATAGATAACAAAAATAGCTCGC');
  fa.seq.should.equal('ACAAAGATAGATAACAAAAATAGCTCGT');
  [fr.block.original_base, fr.block.new_base, fa.block.original_base, fa.block.new_base].should.eql(['T', 'G', 'T', 'G']);
});

test('tmp_1_11193_C_T at -3: reverse both T->G, forward both C->A; the source names the -3 extrapolation', function () {
  const rr = mm('ACAACTTTTTAATATATTGTGTATACTCTAG', 'reverse', 'T', 3);
  const ra = mm('ACAACTTTTTAATATATTGTGTATACTCTAA', 'reverse', 'C', 3);
  rr.seq.should.equal('ACAACTTTTTAATATATTGTGTATACTCGAG');
  ra.seq.should.equal('ACAACTTTTTAATATATTGTGTATACTCGAA');
  const fr = mm('ACAAAGATAGATAACAAAAATAGCTCTC', 'forward', 'T', 3);
  const fa = mm('ACAAAGATAGATAACAAAAATAGCTCTT', 'forward', 'C', 3);
  fr.seq.should.equal('ACAAAGATAGATAACAAAAATAGCTATC');
  fa.seq.should.equal('ACAAAGATAGATAACAAAAATAGCTATT');
  [rr, ra, fr, fa].forEach(function (x) {
    x.block.position.should.equal(3);
    x.block.source.should.equal('Little 1995 Table 9.8.1 applied at -3 (extrapolated)');
    x.block.source.should.match(/-3 \(extrapolated\)/);
  });
  (function () { mismatch.chooseMismatch('ACGTACGTACGTACGTACGT', 'forward', 'A', 4); }).should.throw(RangeError);
  (function () { mismatch.applyMismatch('ACGTACGTACGTACGTACGT', { position: 2, original_base: 'A', new_base: 'C' }); }).should.throw(RangeError);
});

test('substituteTemplate: the mismatch primer is present verbatim on the substituted haplotype copy, forward and reverse', function () {
  // a plus-strand stretch with a forward primer at 5-24 and a reverse primer at 5-24 (revcomp)
  const seq = 'GGGG' + 'TGGTTATCCGAATATAGTCA' + 'CCCCCCCC';
  const fp = { start: 5, end: 24 };
  const fwd = 'TGGTTATCCGAATATAGTCA';
  const fblock = mismatch.chooseMismatch(fwd, 'forward', 'T', 2);
  const fsub = mismatch.substituteTemplate(seq, fp, 'forward', fblock);
  fsub.slice(4, 24).should.equal(mismatch.applyMismatch(fwd, fblock));
  const rev = revcomp('TGGTTATCCGAATATAGTCA');
  const rblock = mismatch.chooseMismatch(rev, 'reverse', 'G', 3);
  const rsub = mismatch.substituteTemplate(seq, fp, 'reverse', rblock);
  revcomp(rsub.slice(4, 24)).should.equal(mismatch.applyMismatch(rev, rblock));
  fsub.length.should.equal(seq.length);
  (function () { mismatch.substituteTemplate('ACGT', { start: 1, end: 30 }, 'forward', fblock); }).should.throw(RangeError);
});

test('an indel whose -k base differs between the haplotypes gets no mismatch and warning MISMATCH_NOT_APPLICABLE', function () {
  // SNV and the example indels share the -2 base; a primer pair whose -2 bases differ does not.
  mismatch.applicable('ATCTTTGACTAGCGAGAAATTCAG', 'ATCTTTGACTAGCGAGAAATTCAT', 2).should.be.true();
  mismatch.applicable('ACAGATGATTTTCCAAATGATGATTCAAA', 'ACAGATGATTTTCCAAATGATGATTCAAG', 2).should.be.true();
  mismatch.applicable('GCAGGAAAAGAAATCCTAACATCATATG', 'GCAGGAAAAGAAATCCTAACATCATATA', 3).should.be.true();
  mismatch.applicable('ACGTACGTACGTACGTTAG', 'ACGTACGTACGTACGTAGT', 2).should.be.false();

  const none = { as_ref: [], as_alt: [], common: [] };
  const disc = { own_allele: { mm_pos: [], likelihood: 'likely' }, other_allele: { mm_pos: [1], likelihood: 'likely_weak' }, terminal_mismatch_class: 'max', in_shift_tract: false };
  const issues = sets.issuesFor({
    g: { structure_warn_th: 47, structure_high_th: 55, as_tm_diff_warn: 1.0, common_tm_low: -1.0, common_tm_high: 3.0, neighbour_3p_window: 5 },
    variant: { shift: 0 }, orientation: 'forward', problems: { as_ref: null, as_alt: null },
    balance: { as_tm_diff: 0n, common_minus_as: 1000000n, common_out: 0n }, dyes: { as_ref: null, as_alt: null }, tailed: null, tailed_pairs: null,
    mismatch_not_applicable: 2, mismatch_structures: null, neighbours: none, discrimination: { as_ref: disc, as_alt: disc }
  });
  issues.map(function (i) { return [i.code, i.severity, i.details]; }).should.eql([['MISMATCH_NOT_APPLICABLE', 'warn', { position: 2 }]]);
  should(issues[0].warning).match(/no deliberate mismatch was applied/);
});
