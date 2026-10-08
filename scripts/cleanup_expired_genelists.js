#!/usr/bin/env node
// cleanup_expired_genelists.js — hard-delete gene lists that were soft-deleted (deletedAt set) more
// than 30 days ago. MONGO ONLY: it does NOT touch the genes-core saved_search field — cross-version
// saved_search propagation/cleanup is handled separately by a different script. Run daily via cron.
//
//   node cleanup_expired_genelists.js
//
// Targets whatever db `genelists` resolves to via gramene-mongodb-config (userData1). `require`
// resolves gramene-mongodb-config from the sibling gramene-swagger/node_modules symlink, so this
// works regardless of the invoking cwd.
'use strict';
var collections = require('gramene-mongodb-config');

var GRACE_DAYS = 30;
var cutoff = new Date(Date.now() - GRACE_DAYS * 24 * 60 * 60 * 1000);

collections.genelists.mongoCollection().then(function (col) {
  return col.deleteMany({ deletedAt: { $lt: cutoff } });
}).then(function (result) {
  var n = (result && (result.deletedCount != null ? result.deletedCount
                    : (result.result && result.result.n))) || 0;
  console.log('[' + new Date().toISOString() + '] cleanup_expired_genelists: purged ' + n +
              ' gene list(s) soft-deleted before ' + cutoff.toISOString());
  collections.closeMongoDatabase();
  setTimeout(function () { process.exit(0); }, 300);
}).catch(function (err) {
  console.error('[' + new Date().toISOString() + '] cleanup_expired_genelists FAILED: ' +
                (err && err.message || err));
  process.exit(1);
});
