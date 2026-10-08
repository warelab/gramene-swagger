'use strict';

// Check plug-in interface (spec §A.8.5).
//   ALGORITHM_VERSION: bump => every job id changes (old results are never reused).
//     '2': max_amplifying_mismatches stringency, overlapping primer footprints discarded, gapped sites from
//          split BLAST HSPs, identical-primer pairs labelled LR, pan-genome summary.truncated.
//   algorithmVersionFor(request) -> the algorithm version jobs.submit hashes into the job id (genotyping spec §5.3):
//     ALGORITHM_VERSION for a request without genotyping, so every existing job id is unchanged, and
//     ALGORITHM_VERSION + '+' + genotype.GENOTYPING_VERSION ('2+g1') for one with a genotyping block, so a change to the
//     allele caller or the prediction invalidates genotyping jobs alone.
//   normalize(body, deps) -> {request, resolved, kind, warnings, estimate {cpu_s, total}, dbs}   (./normalize.js)
//   run(request, ctx) -> results (§B.12)                                                         (./run.js)
// Modules are required lazily so the API process never loads the BLAST pipeline and tests can stub run().

const ALGORITHM_VERSION = '2';

function algorithmVersionFor(request) {
  if (!request || request.genotyping === undefined || request.genotyping === null) return ALGORITHM_VERSION;
  return ALGORITHM_VERSION + '+' + require('./genotype').GENOTYPING_VERSION;
}

module.exports = {
  ALGORITHM_VERSION: ALGORITHM_VERSION,
  algorithmVersionFor: algorithmVersionFor,
  normalize: (...a) => require('./normalize').normalize(...a),
  run: (...a) => require('./run').run(...a)
};
