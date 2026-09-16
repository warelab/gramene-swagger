'use strict';

// Regenerates indel_split.json from real blastn-short output (read-only; a tiny -subject run at nice 10):
//   node test/primers/fixtures/check_core/make_indel_split.js
// Subjects are sorghum_bicolor 4:7423301-7424000 (the sorghum_windows.json window) with
//   ins22   one extra base (T) after base 11 of P2_L (22 nt)          → HSPs q1-11 / q12-22
//   del26   base 13 of a 26-mer (P2_L + 4 nt) deleted                  → HSPs q1-13 / q14-26
//   ins22rc the reverse complement of ins22 (P2_L site on the minus strand)
// Subject coordinates are 1-based within each synthetic subject.

const fs = require('fs');
const os = require('os');
const path = require('path');
const { execFileSync } = require('child_process');

const ROOT = path.resolve(__dirname, '..', '..', '..', '..');
const DOC = require('./sorghum_windows.json');
const blast = require(path.join(ROOT, 'api/helpers/primers/check/blast'));
const { revcomp } = require(path.join(ROOT, 'api/helpers/primers/check/realign'));
const BLASTN = process.env.BLASTN || '/home/olson/bin/blastn';

const w = DOC.windows['4:7423301-7424000'];
const at = (pos) => pos - w.start;
const seq = w.seq.toUpperCase();
const P2_L = DOC.primers.P2_L;
const P2_R = DOC.primers.P2_R;
if (seq.slice(at(7423537), at(7423537) + 22) !== P2_L) throw new Error('P2_L not at 7423537');
const L26 = seq.slice(at(7423537), at(7423537) + 26);

const ins22 = seq.slice(0, at(7423547) + 1) + 'T' + seq.slice(at(7423547) + 1);
const del26 = seq.slice(0, at(7423549)) + seq.slice(at(7423549) + 1);
const ins22rc = revcomp(ins22);

const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'primers-indel-'));
try {
  const subj = path.join(dir, 'subjects.fa');
  fs.writeFileSync(subj, '>ins22\n' + ins22 + '\n>del26\n' + del26 + '\n>ins22rc\n' + ins22rc + '\n');
  const primers = [P2_L, P2_R, L26];
  const common = ['-task', 'blastn-short', '-reward', '1', '-penalty', '-1', '-word_size', '5', '-ungapped',
    '-evalue', '30000', '-searchsp', '15000000000', '-dust', 'no', '-soft_masking', 'false',
    '-max_target_seqs', '5000', '-max_hsps', '100000', '-subject', subj, '-query', '-'];
  const run = (fmt) => execFileSync('/usr/bin/nice', ['-n', '10', BLASTN].concat(common, ['-outfmt', fmt]),
    { input: blast.queryFasta(primers), cwd: dir, env: blast.CHILD_ENV }).toString();
  const out = {
    source: 'blastn 2.13.0 blastn-short r1 p-1 ws5 ungapped e30000 searchsp1.5e10, -subject; sorghum_bicolor 4:7423301-7424000 with planted indels (make_indel_split.js)',
    primers: { P2_L, P2_R, L26 },
    query_order: ['P2_L', 'P2_R', 'L26'],
    subjects: { ins22, del26, ins22rc },
    notes: {
      ins22: 'T inserted after genome 4:7423547 (P2_L base 11); P2_L 5\' end at subject 237, P2_R 5\' end at 447',
      del26: 'genome 4:7423549 (L26 base 13) deleted; L26 5\' end at 237, P2_R 5\' end at 445',
      ins22rc: 'reverse complement of ins22 (length ' + ins22.length + '); P2_L is on the minus strand with 5\' end at ' + (ins22.length - 237 + 1)
    },
    genome_hits: run(blast.GENOME_OUTFMT).split('\n').filter(Boolean),
    cdna_hits: run(blast.CDNA_OUTFMT).split('\n').filter(Boolean)
  };
  fs.writeFileSync(path.join(__dirname, 'indel_split.json'), JSON.stringify(out, null, 1) + '\n');
} finally {
  fs.rmSync(dir, { recursive: true, force: true });
}
