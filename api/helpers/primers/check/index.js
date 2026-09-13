'use strict';

// Check plug-in interface (spec §A.8.5).
//   ALGORITHM_VERSION: bump => every job id changes (old results are never reused).
//     '2': max_amplifying_mismatches stringency, overlapping primer footprints discarded, gapped sites from
//          split BLAST HSPs, identical-primer pairs labelled LR, pan-genome summary.truncated.
//   normalize(body, deps) -> {request, resolved, kind, warnings, estimate {cpu_s, total}, dbs}   (./normalize.js)
//   run(request, ctx) -> results (§B.12)                                                         (./run.js)
// Both are required lazily so the API process never loads the BLAST pipeline and tests can stub run().

module.exports = {
  ALGORITHM_VERSION: '2',
  normalize: (...a) => require('./normalize').normalize(...a),
  run: (...a) => require('./run').run(...a)
};
