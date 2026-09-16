'use strict';

// Entry point of the check worker, a separate pm2 app (plan Overrides: no supervisor; app.js does not fork it).
//   pm2 start api/helpers/primers/jobs/worker_main.js --name sorghum_primers11 --cwd <checkout> --exp-backoff-restart-delay 5000
//   dev: PRIMERS_SITE_KEY=sorghum_v11_dev node api/helpers/primers/jobs/worker_main.js
//
// - site_key comes from config.siteKey(): PRIMERS_SITE_KEY, primers.site_key, else
//   "<swagger.yaml basePath without '/'>:<mongo db>", the same derivation the API uses. It is logged at startup
//   as `primers site_key=<key>` (jobs/index.js siteKeyLogLine) so the two processes can be compared.
// - SIGTERM / SIGINT: running jobs go back to the FRONT of their queue, the lock is released, exit 0 within ~1 s.
//   A second signal exits immediately.
// - Mongo unavailable (mongoCollection() resolves undefined or times out) at startup, or a job failing with
//   MONGO_UNAVAILABLE: exit 75 so pm2 restarts the process (gramene-mongodb-config never reconnects).

require('../node_compat'); // before `config` and any library that still calls util.is*

const EXIT_OK = 0;
const EXIT_CRASH = 1;
const EXIT_MONGO_UNAVAILABLE = 75;
const MONGO_PROBE_TIMEOUT_MS = 15000;
const SHUTDOWN_EXIT_MS = 950;

function withTimeout(promise, ms, message) {
  let timer;
  const timeout = new Promise(function (resolve, reject) {
    timer = setTimeout(function () { reject(new Error(message)); }, ms);
  });
  return Promise.race([promise, timeout]).finally(function () { clearTimeout(timer); });
}

// true when the genes collection is reachable.
async function probeMongo(mongo, timeoutMs) {
  const coll = await withTimeout(Promise.resolve().then(function () { return mongo.genes.mongoCollection(); }),
    timeoutMs || MONGO_PROBE_TIMEOUT_MS, 'mongo probe timed out');
  return !!coll;
}

function info(log, msg) {
  (log.info || log.log).call(log, msg);
}

function errorLog(log, msg) {
  (log.error || log.log).call(log, msg);
}

// main(opts) -> {worker, store, siteKey, shutdown} | {idle: true, stop} | undefined (after exit 75)
// opts (tests): {config, log, exit, siteKey, store, check, mongo, probeMongo (false to skip), signals (false to skip),
//                timing, mongoProbeTimeoutMs}
async function main(opts) {
  opts = opts || {};
  const config = opts.config || require('../config');
  const jobs = require('./index');
  const cfg = config.get();
  const log = opts.log || console;
  const exit = opts.exit || function (code) { process.exit(code); };
  const siteKey = opts.siteKey || config.siteKey();
  const signals = opts.signals !== false;

  info(log, jobs.siteKeyLogLine(siteKey));
  info(log, 'primers worker: pid ' + process.pid + ', store ' + cfg.check.store + ', basePath ' + config.basePath() +
    ', global_prefix ' + cfg.check.global_prefix + ', local_max_jobs ' + cfg.check.local_max_jobs);

  function idle(reason) {
    const timer = setInterval(function () {}, 1 << 30);
    const stopIdle = function () {
      clearInterval(timer);
      exit(EXIT_OK);
    };
    // Handlers go in before the log line: "idling" is the readiness signal, and a SIGTERM that
    // arrived between the two would kill the process by default action instead of exiting 0.
    if (signals) {
      process.once('SIGTERM', stopIdle);
      process.once('SIGINT', stopIdle);
    }
    errorLog(log, 'primers worker: ' + reason + '; idling (no jobs will run)');
    return { idle: true, stop: stopIdle };
  }

  if (!cfg.enabled) return idle('primers are disabled (primers.enabled / PRIMERS_ENABLED)');
  if (!opts.store && cfg.check.store !== 'redis') {
    return idle('check.store is "' + cfg.check.store + '"; a separate worker process needs the redis store');
  }

  const mongo = opts.mongo || require('gramene-mongodb-config');
  if (opts.probeMongo !== false) {
    let ok = false;
    try {
      ok = await probeMongo(mongo, opts.mongoProbeTimeoutMs);
    } catch (err) {
      errorLog(log, 'primers worker: mongo probe failed: ' + (err && err.message));
    }
    if (!ok) {
      errorLog(log, 'primers worker: mongo is unavailable; exiting with code ' + EXIT_MONGO_UNAVAILABLE);
      exit(EXIT_MONGO_UNAVAILABLE);
      return undefined;
    }
  }

  const store = opts.store || require('./redis_store').createRedisStore({ siteKey: siteKey, cfg: cfg, role: 'worker', log: log });
  let worker = null;
  let shuttingDown = false;

  async function shutdown(reason, code) {
    if (shuttingDown) {
      if (reason === 'SIGTERM' || reason === 'SIGINT') exit(code);
      return;
    }
    shuttingDown = true;
    info(log, 'primers worker: ' + reason + '; requeueing running jobs and exiting');
    const hard = setTimeout(function () { exit(code); }, SHUTDOWN_EXIT_MS);
    try {
      if (worker) await worker.stop();
    } catch (err) {
      errorLog(log, 'primers worker: stop failed: ' + (err && err.message));
    }
    try {
      await withTimeout(Promise.resolve(store.close()), 150, 'store close timed out');
    } catch (err) { /* exiting anyway */ }
    clearTimeout(hard);
    exit(code);
  }

  const { createWorker } = require('./worker');
  worker = createWorker({
    store: store,
    cfg: cfg,
    siteKey: siteKey,
    check: opts.check,
    mongo: mongo,
    log: log,
    timing: opts.timing,
    onFatal: function (err) { shutdown('fatal ' + (err && err.code), EXIT_MONGO_UNAVAILABLE); }
  });

  if (signals) {
    process.on('SIGTERM', function () { shutdown('SIGTERM', EXIT_OK); });
    process.on('SIGINT', function () { shutdown('SIGINT', EXIT_OK); });
  }
  worker.start();
  return { worker: worker, store: store, siteKey: siteKey, shutdown: shutdown };
}

if (require.main === module) {
  process.on('unhandledRejection', function (reason) {
    console.error('primers worker: unhandled rejection:', reason && reason.stack ? reason.stack : reason);
  });
  main().catch(function (err) {
    console.error('primers worker: startup failed:', err && err.stack ? err.stack : err);
    process.exit(EXIT_CRASH);
  });
}

module.exports = {
  main,
  probeMongo,
  EXIT_OK,
  EXIT_CRASH,
  EXIT_MONGO_UNAVAILABLE
};
