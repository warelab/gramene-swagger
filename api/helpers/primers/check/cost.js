'use strict';

// Check cost estimate (spec §B.13). The server estimate is authoritative; the gramene-primers
// client mirrors estimate() for display (src/cost.ts).
//
//   est_cpu_s = uniq_primers × [ ref_Gb × c(ws_ref) + (transcript ? ref_cdna_Gb × c(ws_ref) : 0)
//                                + f_pan × Σ_pan (transcript ? cdna_Gb : genome_Gb) × c(ws_pan)
//                                + genome_tasks × c_realign ]
//   c(ws) = cfg.check.cpu_s_per_primer_gb['ws' + ws]   ({ws5: 5.2, ws6: 2.2, ws7: 1.2}, single-thread, measured on sorghum)
//   f_pan = cfg.check.pangenome_cpu_factor, default 2.0 (values below 1 fall back to the default). Pan-genome
//           genomes run as up to 8 concurrent single-thread blastn processes, which cost more CPU per primer·Gb
//           than the single-thread measurement. Full-panel run on squam (algorithm v2): pair P3 (2 primers) against
//           all 119 other sorghum genomes (83.8 Gb) used 758.2 CPU-s (733.2 in blastn children + 25.0 in the
//           worker) against an estimate of 520 without the factor; word size 6 measured 4.33 CPU-s per primer·Gb
//           (factor 1.97). With 2.0 the same job estimates about 889 CPU-s.
//   genome_Gb = total_bases / 1e9, or FALLBACK_GENOME_GB when unknown; cDNA size is not known in the API, so
//               CDNA_GB_ESTIMATE (0.15 Gb) is used.
//   genome_tasks = 1 (reference genome) + (transcript ? 0 : number of pan-genome genomes); cDNA tasks re-align
//                  from BLAST's qseq/sseq and are not charged.
//   c_realign = cfg.check.realign_cpu_s_per_primer_task, default 0.6 CPU-s per primer per genome task: FASTA
//               re-alignment in the worker once sites whose gap-aware bound exceeds max_amplifying_mismatches are
//               no longer re-aligned. Measured on squam (process.cpuUsage, warm cache, two runs each):
//               sorghum_bicolor ws5, 2 primers (PA) 0.51-0.62 s, 9 primers incl. the repetitive P5_L 4.75-5.18 s;
//               sorghum_353 ws6, 2 primers at max size 5932 0.91-1.13 s, 9 primers 4.93-5.05 s
//               → 0.26-0.58 CPU-s per primer (re-aligning every site: 0.76-2.07). The full-panel run above measured
//               0.105 CPU-s per primer·task for re-alignment plus annotation, so 0.6 is conservative.
//
// A job over cfg.check.max_job_cpu_s (6000) is refused with 422 JOB_TOO_LARGE {estimate_cpu_s, limit}. With the
// factor, 10 primers against all 119 other sorghum genomes estimate about 4,400 CPU-s (accepted) and 20 primers
// about 8,900 CPU-s (refused).

const { PrimerHttpError } = require('../errors');

const DEFAULT_COEFFICIENTS = Object.freeze({ ws5: 5.2, ws6: 2.2, ws7: 1.2 });
const DEFAULT_REFERENCE_WORD_SIZE = 5;
const DEFAULT_PANGENOME_WORD_SIZE = 6;
const DEFAULT_MAX_JOB_CPU_S = 6000;
// Spec §B.13: "use 0.15 Gb in the API estimate". Real sorghum cdna.all DBs are ~0.10 Gb
// (sorghum_bicolor: 103,922,012 bases), so this is deliberately conservative.
const CDNA_GB_ESTIMATE = 0.15;
// An assembly without a .fai (total_bases null) is charged as 1 Gb. None exist in sorghum_v11.
const FALLBACK_GENOME_GB = 1.0;
const DEFAULT_REALIGN_CPU_S_PER_PRIMER_TASK = 0.6;
const DEFAULT_PANGENOME_CPU_FACTOR = 2.0;

function realignCoefficient(cfg) {
  const v = Number(checkCfg(cfg).realign_cpu_s_per_primer_task);
  return Number.isFinite(v) && v >= 0 ? v : DEFAULT_REALIGN_CPU_S_PER_PRIMER_TASK;
}

// Multiplier on the pan-genome BLAST term; a missing, non-numeric or < 1 value uses the default.
function pangenomeFactor(cfg) {
  const raw = checkCfg(cfg).pangenome_cpu_factor;
  const v = raw === null || raw === undefined || raw === '' ? NaN : Number(raw);
  return Number.isFinite(v) && v >= 1 ? v : DEFAULT_PANGENOME_CPU_FACTOR;
}

function checkCfg(cfg) {
  return (cfg && cfg.check) || {};
}

// CPU-s per primer·Gb for a BLAST word size. Unknown word sizes use the coefficient of the
// largest configured word size that is <= ws (smaller word sizes cost more), or the most
// expensive configured coefficient when ws is below every configured one.
function coefficient(wordSize, cfg) {
  const table = checkCfg(cfg).cpu_s_per_primer_gb || DEFAULT_COEFFICIENTS;
  const direct = Number(table['ws' + wordSize]);
  if (Number.isFinite(direct) && direct > 0) return direct;
  const known = Object.keys(table)
    .map(function (k) { const m = /^ws(\d+)$/.exec(k); return m ? { ws: Number(m[1]), c: Number(table[k]) } : null; })
    .filter(function (e) { return e && Number.isFinite(e.c) && e.c > 0; })
    .sort(function (a, b) { return a.ws - b.ws; });
  if (known.length === 0) {
    const d = DEFAULT_COEFFICIENTS['ws' + wordSize];
    return d !== undefined ? d : DEFAULT_COEFFICIENTS.ws5;
  }
  const lower = known.filter(function (e) { return e.ws <= wordSize; });
  if (lower.length) return lower[lower.length - 1].c;
  return Math.max.apply(null, known.map(function (e) { return e.c; }));
}

function genomeGb(asm) {
  const tb = asm && Number(asm.total_bases);
  return Number.isFinite(tb) && tb > 0 ? tb / 1e9 : FALLBACK_GENOME_GB;
}

function round1(x) {
  return Math.round(x * 10) / 10;
}

// Integer CPU-s, rounded up (float noise below 1e-6 is ignored so 37.0000000001 stays 37).
function ceilCpu(x) {
  return Math.ceil(Math.round(x * 1e6) / 1e6);
}

// estimate({unique_primers, mode, reference, pangenome, cfg, cdna_gb?})
//   unique_primers: number of distinct uppercase primers (<= 20)
//   mode: 'gene' | 'transcript' | 'region' | 'sequence'
//   reference: resolved assembly ({total_bases}); pangenome: array of resolved assemblies ([] when not checked)
//   cdna_gb: per-assembly cDNA size in Gb (default CDNA_GB_ESTIMATE)
// -> {cpu_s (integer, rounded up), total (tasks), breakdown {reference, transcriptome, pangenome, realign} (CPU-s, 1 dp;
//     pangenome includes the pan-genome CPU factor), word_sizes {reference, pangenome}}
// total = 1 reference genome task + 1 reference cDNA task in transcript mode + 1 per pan-genome genome (§B.2).
function estimate(opts) {
  const o = opts || {};
  const ccfg = checkCfg(o.cfg);
  const primers = Math.max(0, Number(o.unique_primers) || 0);
  const transcript = o.mode === 'transcript';
  const pan = Array.isArray(o.pangenome) ? o.pangenome : [];
  const wsRef = Number.isInteger(ccfg.reference_word_size) ? ccfg.reference_word_size : DEFAULT_REFERENCE_WORD_SIZE;
  const wsPan = Number.isInteger(ccfg.pangenome_word_size) ? ccfg.pangenome_word_size : DEFAULT_PANGENOME_WORD_SIZE;
  const cRef = coefficient(wsRef, o.cfg);
  const cPan = coefficient(wsPan, o.cfg);
  const cdnaGb = Number.isFinite(Number(o.cdna_gb)) && Number(o.cdna_gb) > 0 ? Number(o.cdna_gb) : CDNA_GB_ESTIMATE;

  const reference = primers * genomeGb(o.reference) * cRef;
  const transcriptome = transcript ? primers * cdnaGb * cRef : 0;
  let panGb = 0;
  pan.forEach(function (asm) { panGb += transcript ? cdnaGb : genomeGb(asm); });
  const pangenome = primers * panGb * cPan * pangenomeFactor(o.cfg);
  const genomeTasks = 1 + (transcript ? 0 : pan.length);
  const realign = primers * genomeTasks * realignCoefficient(o.cfg);

  return {
    cpu_s: ceilCpu(reference + transcriptome + pangenome + realign),
    total: 1 + (transcript ? 1 : 0) + pan.length,
    breakdown: { reference: round1(reference), transcriptome: round1(transcriptome), pangenome: round1(pangenome), realign: round1(realign) },
    word_sizes: { reference: wsRef, pangenome: wsPan }
  };
}

function limitOf(cfg) {
  const v = Number(checkCfg(cfg).max_job_cpu_s);
  return Number.isFinite(v) && v >= 0 ? v : DEFAULT_MAX_JOB_CPU_S;
}

// Throws 422 JOB_TOO_LARGE {estimate_cpu_s, limit} when est.cpu_s exceeds cfg.check.max_job_cpu_s.
function assertWithinLimit(est, cfg) {
  const limit = limitOf(cfg);
  if (est.cpu_s > limit) {
    throw new PrimerHttpError(422, 'JOB_TOO_LARGE',
      'the check is estimated at ' + est.cpu_s + ' CPU-seconds, over the limit of ' + limit +
      '; check fewer primers or genomes', { estimate_cpu_s: est.cpu_s, limit: limit });
  }
  return est;
}

module.exports = {
  DEFAULT_COEFFICIENTS,
  CDNA_GB_ESTIMATE,
  FALLBACK_GENOME_GB,
  DEFAULT_REALIGN_CPU_S_PER_PRIMER_TASK,
  DEFAULT_PANGENOME_CPU_FACTOR,
  coefficient,
  realignCoefficient,
  pangenomeFactor,
  genomeGb,
  estimate,
  assertWithinLimit,
  limitOf
};
