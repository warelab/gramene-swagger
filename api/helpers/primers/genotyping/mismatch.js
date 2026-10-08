'use strict';

// The deliberate mismatch of gel allele-specific PCR (spec §4.11): a literal lookup of Little, S. (1995), Current
// Protocols in Human Genetics 9.8, Table 9.8.1, with the mismatch classes of Table 9.8.3. Pure.
//
// Row: the primer's own terminal mispair, the unordered pair {the primer's 3' base, the template base opposite it on
// the other allele}, written in alphabetical order. For a forward (plus-strand) primer that template base is the
// complement of the other allele's plus-strand base; for a reverse primer it is the plus-strand base itself.
// Column: the template base opposite the mismatch position, i.e. the complement of the primer's original base there.
// Cell: the base the primer carries at that position instead. Each primer is looked up on its own row, so the REF and
// ALT primers of one set can receive different bases (tmp_1_11193_C_T reverse at -2: A->C and A->T).

const COMPLEMENT = Object.freeze({ A: 'T', C: 'G', G: 'C', T: 'A' });
const BASE_RE = /^[ACGT]$/;
const PRIMER_RE = /^[ACGT]+$/;

// Table 9.8.1, 8 rows x 4 template columns = 32 cells. Little prints "C or T" in the GT row's T column; C is chosen:
// it is a maximum mismatch in Little and strong in Ye 2001, so it survives their disagreement.
const TABLE_9_8_1 = Object.freeze({
  AA: Object.freeze({ A: 'A', G: 'G', C: 'A', T: 'G' }),
  GG: Object.freeze({ A: 'A', G: 'G', C: 'A', T: 'G' }),
  AG: Object.freeze({ A: 'C', G: 'T', C: 'A', T: 'G' }),
  TT: Object.freeze({ A: 'C', G: 'T', C: 'A', T: 'G' }),
  CT: Object.freeze({ A: 'C', G: 'T', C: 'A', T: 'G' }),
  CC: Object.freeze({ A: 'C', G: 'T', C: 'A', T: 'G' }),
  AC: Object.freeze({ A: 'G', G: 'A', C: 'C', T: 'T' }),
  GT: Object.freeze({ A: 'G', G: 'A', C: 'T', T: 'C' })
});

// Table 9.8.3 mismatch classes, keyed like the rows.
const CLASSES = Object.freeze({ AG: 'max', CT: 'max', TT: 'max', CC: 'strong', AA: 'medium', GG: 'medium', AC: 'weak', GT: 'weak' });

const POSITIONS = Object.freeze([2, 3]);
const SOURCES = Object.freeze({ 2: 'Little 1995 Table 9.8.1', 3: 'Little 1995 Table 9.8.1 applied at -3 (extrapolated)' });

function assertBase(b, name) {
  if (typeof b !== 'string' || !BASE_RE.test(b)) throw new TypeError(name + ' must be one of A, C, G, T');
}

function assertPrimer(p) {
  if (typeof p !== 'string' || !PRIMER_RE.test(p)) throw new TypeError('primer must be an uppercase A/C/G/T sequence');
}

function assertOrientation(o) {
  if (o !== 'forward' && o !== 'reverse') throw new TypeError('orientation must be forward or reverse');
}

// The unordered pair of two bases, alphabetical ('GA' -> 'AG').
function pairKey(a, b) {
  assertBase(a, 'base');
  assertBase(b, 'base');
  return a <= b ? a + b : b + a;
}

// The template base a primer of this orientation pairs with, where the plus strand has plusBase.
function templateBase(orientation, plusBase) {
  assertOrientation(orientation);
  assertBase(plusBase, 'plusBase');
  return orientation === 'forward' ? COMPLEMENT[plusBase] : plusBase;
}

// terminalPair(primer3Base, templateBase) -> the table row, e.g. terminalPair('G', 'A') = 'AG'.
function terminalPair(primer3Base, template) {
  return pairKey(primer3Base, template);
}

// The Little class of a primer's 3' base against the other allele: `otherPlusBase` is the other allele's plus-strand
// base at the primer's 3' position (discriminating.<orientation>.alt_base for as_ref, .ref_base for as_alt).
function terminalClass(primer, orientation, otherPlusBase) {
  assertPrimer(primer);
  return CLASSES[terminalPair(primer[primer.length - 1], templateBase(orientation, otherPlusBase))];
}

// chooseMismatch(primer, orientation, otherPlusBase, position) -> PrimerDeliberateMismatch
//   {position, original_base, new_base, template_base, terminal_pair, terminal_mismatch_class, added_mismatch_class, source}
//   primer: the perfect-match allele-specific primer (matched_seq), 5'->3'; position: 2 or 3 from the 3' end.
function chooseMismatch(primer, orientation, otherPlusBase, position) {
  assertPrimer(primer);
  if (POSITIONS.indexOf(position) < 0) throw new RangeError('position must be 2 or 3');
  if (primer.length <= position) throw new RangeError('the primer is too short for a mismatch at -' + position);
  const row = terminalPair(primer[primer.length - 1], templateBase(orientation, otherPlusBase));
  const original = primer[primer.length - position];
  const column = COMPLEMENT[original];
  const replacement = TABLE_9_8_1[row][column];
  return {
    position: position,
    original_base: original,
    new_base: replacement,
    template_base: column,
    terminal_pair: row,
    terminal_mismatch_class: CLASSES[row],
    added_mismatch_class: CLASSES[pairKey(replacement, column)],
    source: SOURCES[position]
  };
}

// applyMismatch(primer, block) -> the primer with block.new_base at block.position from the 3' end.
function applyMismatch(primer, block) {
  assertPrimer(primer);
  const i = primer.length - block.position;
  if (primer[i] !== block.original_base) throw new RangeError('the primer does not carry ' + block.original_base + ' at -' + block.position);
  assertBase(block.new_base, 'block.new_base');
  return primer.slice(0, i) + block.new_base + primer.slice(i + 1);
}

// §4.11 guard: the mismatch applies only when the base at -position is the same in both allele-specific primers;
// otherwise (an indel) the "extra" mismatch would itself be allele-specific (MISMATCH_NOT_APPLICABLE).
function applicable(refPrimer, altPrimer, position) {
  assertPrimer(refPrimer);
  assertPrimer(altPrimer);
  if (refPrimer.length <= position || altPrimer.length <= position) return false;
  return refPrimer[refPrimer.length - position] === altPrimer[altPrimer.length - position];
}

// The copy of a haplotype template on which the mismatch primer is present verbatim (§4.10): the template base
// opposite the mismatch is replaced. footprint: the primer's template span {start, end} (1-based) on `seq`.
function substituteTemplate(seq, footprint, orientation, block) {
  assertOrientation(orientation);
  const idx = orientation === 'forward' ? footprint.end - (block.position - 1) : footprint.start + (block.position - 1);
  if (!Number.isSafeInteger(idx) || idx < 1 || idx > seq.length) throw new RangeError('the mismatch lies outside the template');
  const plus = orientation === 'forward' ? block.new_base : COMPLEMENT[block.new_base];
  return seq.slice(0, idx - 1) + plus + seq.slice(idx);
}

// The order-row note of a mismatch primer (§4.17), e.g. "deliberate mismatch −2 A→G (Little 1995 Table 9.8.1)".
function note(block) {
  return 'deliberate mismatch −' + block.position + ' ' + block.original_base + '→' + block.new_base + ' (' + block.source + ')';
}

module.exports = {
  TABLE_9_8_1,
  CLASSES,
  COMPLEMENT,
  POSITIONS,
  SOURCES,
  pairKey,
  templateBase,
  terminalPair,
  terminalClass,
  chooseMismatch,
  applyMismatch,
  applicable,
  substituteTemplate,
  note
};
