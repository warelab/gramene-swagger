'use strict';

// Starts this worktree's API (app.js) as a loopback-only child process, for integration tests that need a server of their
// own, such as one pointed at the fake Ensembl (genotyping spec §6.4, §7.6). The command is the §6.4 dev instance:
//   HOST=127.0.0.1 PORT=<port> SWAGGER_HOST=localhost:<port> SWAGGER_SCHEMES=http PRIMERS_SITE_KEY=sorghum_v11_geno
//   PRIMERS_GLOBAL_MAX_JOBS=1 node app.js
// with opts.env added. No worker is started.
//
//   const api = await startApi({port: 50112, env: {PRIMERS_VARIATION_URL: fake.url}});
//   api.base       http://127.0.0.1:50112/sorghum_v11
//   api.siteKey    from the "primers site_key=" startup line
//   api.log()      everything the child has printed
//   await api.stop()   SIGTERM (SIGKILL after 10 s), then waits until the port is closed
//
// It refuses a port that already accepts connections (a dev API is running there) and the shared ports 50011 and 50111. The
// child inherits the environment except NODE_TEST_CONTEXT, which would make it a test-runner child, and NODE_CONFIG and
// PRIMERS_VARIATION_ENABLED, because the defaults are under test. opts.env can set any of them again.

const { spawn } = require('child_process');
const net = require('net');
const path = require('path');

const ROOT = path.resolve(__dirname, '..', '..', '..', '..');
const HOST = '127.0.0.1';
const FORBIDDEN_PORTS = new Set([50011, 50111]);
const MAX_LOG = 500000;
const DROPPED_ENV = ['NODE_TEST_CONTEXT', 'NODE_CONFIG', 'PRIMERS_VARIATION_ENABLED'];

const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));

function portOpen(port) {
  return new Promise(function (resolve) {
    const socket = net.connect({ port: port, host: HOST });
    socket.once('connect', function () {
      socket.destroy();
      resolve(true);
    });
    socket.once('error', function () { resolve(false); });
  });
}

async function startApi(opts) {
  opts = opts || {};
  const port = Number(opts.port);
  const basePath = opts.basePath || '/sorghum_v11';
  const readyTimeoutMs = opts.readyTimeoutMs || 90000;
  if (!Number.isSafeInteger(port) || port < 1024 || port > 65535) throw new TypeError('startApi: port must be an integer from 1024 to 65535');
  if (FORBIDDEN_PORTS.has(port)) throw new Error('startApi: port ' + port + ' belongs to a shared instance and is never used by tests');
  if (await portOpen(port)) {
    throw new Error('startApi: 127.0.0.1:' + port + ' is already in use. This test starts its own API there; stop the dev API first');
  }

  const env = Object.assign({}, process.env);
  DROPPED_ENV.forEach(function (k) { delete env[k]; });
  Object.assign(env, {
    HOST: HOST,
    PORT: String(port),
    SWAGGER_HOST: 'localhost:' + port,
    SWAGGER_SCHEMES: 'http',
    PRIMERS_SITE_KEY: 'sorghum_v11_geno',
    PRIMERS_GLOBAL_MAX_JOBS: '1'
  }, opts.env || {});

  const started = Date.now();
  const child = spawn(process.execPath, ['app.js'], { cwd: ROOT, env: env, stdio: ['ignore', 'pipe', 'pipe'] });
  let log = '';
  const append = function (d) {
    log += d;
    if (log.length > MAX_LOG) log = log.slice(log.length - MAX_LOG);
  };
  child.stdout.on('data', append);
  child.stderr.on('data', append);
  let exit = null;
  const exited = new Promise(function (resolve) {
    child.once('exit', function (code, signal) {
      exit = { code: code, signal: signal };
      resolve(exit);
    });
  });
  const killOnExit = function () {
    try { child.kill('SIGKILL'); } catch (e) { /* already gone */ }
  };
  process.once('exit', killOnExit);

  const base = 'http://' + HOST + ':' + port + basePath;
  let stopping = null;
  const api = {
    port: port,
    base: base,
    pid: child.pid,
    siteKey: null,
    startup_ms: null,
    log: function () { return log; },
    get exited() { return exit; },
    stop: function () {
      if (!stopping) {
        stopping = (async function () {
          if (!exit) {
            child.kill('SIGTERM');
            const t = setTimeout(function () { if (!exit) child.kill('SIGKILL'); }, 10000);
            await exited;
            clearTimeout(t);
          }
          process.removeListener('exit', killOnExit);
          for (let i = 0; i < 100 && await portOpen(port); i++) await sleep(100);
        })();
      }
      return stopping;
    }
  };

  try {
    for (;;) {
      if (exit) throw new Error('the API exited during startup (' + JSON.stringify(exit) + ')');
      if (Date.now() - started > readyTimeoutMs) throw new Error('the API did not answer within ' + readyTimeoutMs + ' ms');
      if (/Listening on \d+/.test(log)) {
        try {
          const res = await fetch(base + '/primers/genomes?system_name=sorghum_bicolor', { signal: AbortSignal.timeout(5000) });
          await res.arrayBuffer();
          break;
        } catch (e) { /* not yet */ }
      }
      await sleep(100);
    }
  } catch (err) {
    await api.stop();
    err.message = 'startApi on port ' + port + ': ' + err.message + '\n--- API log (tail) ---\n' + log.slice(-3000);
    throw err;
  }
  const m = /primers site_key=(\S+)/.exec(log);
  api.siteKey = m ? m[1] : null;
  api.startup_ms = Date.now() - started;
  return api;
}

module.exports = { startApi, portOpen, ROOT };
