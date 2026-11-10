'use strict';

var _ = require('lodash');
var mongoHelper = require('../helpers/mongo');
var solrHelper = require('../helpers/solr');
var mongoCollections = require('gramene-mongodb-config');
var treeShape = require('gramene-mongodb-config/tree_shape');
const murmur = require('murmurhash3js');
const axios = require('axios');

// POST /gene_lists guards. Mirrors the limits enforced by /gene_lists/validate in solr.js, so a
// list that validates can always be saved.
const MAX_LIST_IDS = 6000;
const MAX_ID_LENGTH = 255;
var JSONStream = require('JSONStream');
var csvStringify = require('csv-stringify');
var bedify = require('gramene-bedify');
var through2 = require('through2');
const admin = require("firebase-admin");
const { getAuth } = require("firebase-admin/auth")

var path = require('path');
var firebaseCredentialsPath = process.env.FIREBASE_CREDENTIALS_PATH
  ? path.resolve(process.cwd(), process.env.FIREBASE_CREDENTIALS_PATH)
  : "/usr/local/gramene/gramene-auth-firebase-adminsdk-c1sc0-263ff4cc4f.json";
var serviceAccount = require(firebaseCredentialsPath);

admin.initializeApp({
  credential: admin.credential.cert(serviceAccount)
});

(function init() {
  var toExport = {};

  // add a function to the controller for each mongo collection.
  _.forOwn(mongoCollections, function (coll, name) {
    toExport[name] = getFactory(coll.mongoCollection());
  });

  module.exports.get = function get(req, res) {
    var collection, toCall;
    collection = req.swagger.params.collection;
    toCall = toExport[collection];
    if (!toCall) {
      throw new Error("Required parameter `collection` not an allowed value");
    }
    return toCall(req, res);
  };
  // add a function to save a list
  toExport['saveList'] = saveList;
  toExport['deleteList'] = deleteList;
  toExport['updateList'] = updateList;
  toExport['restoreList'] = restoreList;
  // dedicated genelists GET: proper site/isPublic/uid filtering + soft-delete (active vs trash) views,
  // which the generic getFactory/buildQuery can't express ($exists / isPublic scoping).
  // /genetrees rebuilds the nested tree from the stored FLAT shape. Storage moved because the
  // nested form cannot hold a tree deeper than ~112 levels -- the server rejects it for
  // nested-object depth, which silently cost release 70 its 147 deepest trees. The rebuilt tree
  // is deep-equal to what this endpoint used to return: tree_shape's fromNested/toNested are
  // exact inverses, verified over real trees by a build gate. Only _id's position needs fixing
  // up below, and the response is byte-identical once it is.
  //
  // fetchWholeDoc is required, not an optimisation choice: `fl` used to be a mongo projection over
  // the nested document, but root fields such as taxon_name now live on nodes[0], so projecting in
  // mongo would return nothing for them. The whole document is fetched and `fl` applied to the
  // REBUILT tree instead, which is what preserves the old meaning of the parameter.
  toExport['genetrees'] = getFactory(mongoCollections.genetrees.mongoCollection(), {
    fetchWholeDoc: true,
    docTransform: function (fl) {
      return through2.obj(function (doc, enc, done) {
        var out = doc;
        if (doc && doc.nodes) {
          out = treeShape.toNested(doc);
          // Put _id back at the front. The server hoists _id to the first key of every stored
          // document, so the nested trees this endpoint used to return always led with it. A
          // rebuilt root instead carries _id wherever the builder happened to set it -- release
          // 70's builder sets it mid-object, just before `representative`. Deep-equal either
          // way, but this response IS the public document, so keep its old byte order. The copy
          // is one object deep; children are shared, not cloned.
          if (out._id !== undefined && Object.keys(out)[0] !== '_id') {
            var ordered = { _id: out._id };
            Object.keys(out).forEach(function (k) { if (k !== '_id') { ordered[k] = out[k]; } });
            out = ordered;
          }
        }
        if (fl && fl.length) {
          var picked = {};
          fl.forEach(function (f) { if (out[f] !== undefined) { picked[f] = out[f]; } });
          out = picked;   // as before, fl also drops _id unless it was asked for
        }
        this.push(out);
        done();
      });
    }
  });

  // `MSU_id` compatibility alias for the retained /oryza_v9 path.
  //
  // oryza9 served the rice RAP -> MSU/RGAP correspondence as a top-level scalar `MSU_id` on
  // 33,938 genes. This lineage deliberately keeps alternate ids in gene.alt_id instead (an
  // array, shared with Sobic./GRMZM/Vitvi ids for other species), so a straight rebuild drops
  // the old key from every response. The public path does not change across the v9 -> v10
  // rebuild, so the key is aliased back on here rather than being allowed to vanish -- the same
  // reasoning that kept the ORYZA9GT_ gene-tree ids.
  //
  // alt_id remains the real field; this is a read-time alias only, so nothing writes MSU_id back.
  var MSU_RE = /^LOC_Os/;
  toExport['genes'] = getFactory(mongoCollections.genes.mongoCollection(), {
    // a caller asking only for MSU_id still needs alt_id fetched to derive it
    projectionExtras: function (fl) {
      return (fl.indexOf('MSU_id') >= 0 && fl.indexOf('alt_id') < 0) ? ['alt_id'] : [];
    },
    docTransform: function (fl) {
      var all = !fl || !fl.length;
      var wantMsu = all || fl.indexOf('MSU_id') >= 0;
      var wantAlt = all || fl.indexOf('alt_id') >= 0;
      return through2.obj(function (doc, enc, done) {
        if (doc && wantMsu) {
          var alt = Array.isArray(doc.alt_id) ? doc.alt_id : (doc.alt_id ? [doc.alt_id] : []);
          var msu = alt.filter(function (x) { return MSU_RE.test(x); });
          // v9's MSU_id was a single string, and consumers parse it as one. RAP-DB maps 459 RAP
          // loci to more than one MSU gene -- which v9's scalar hand-made table could not
          // express at all -- so emit the first and leave the complete list in alt_id rather
          // than changing this field's type under existing callers.
          if (msu.length) { doc.MSU_id = msu[0]; }
        }
        if (doc && !wantAlt) { delete doc.alt_id; }
        this.push(doc);
        done();
      });
    }
  });

  toExport['genelists'] = genelistsHandler;

  // /saved_views GET has two modes: single-hash lookup vs listing. Route hash lookups to
  // getSavedViewByHash, and listings to the dedicated savedViewsListHandler (proper
  // site/isPublic/uid scoping, like genelistsHandler — the generic getFactory/buildQuery ANDs
  // uid=0 for anonymous requests, which hid all public saved views).
  toExport.savedviews = function (req, res) {
    if (req.swagger.params.hash && req.swagger.params.hash.value) {
      return getSavedViewByHash(req, res);
    }
    return savedViewsListHandler(req, res);
  };
  toExport['saveView'] = saveView;
  toExport['updateView'] = updateView;
  toExport['deleteView'] = deleteView;
  module.exports = toExport;
}());

// opts.docTransform(fl) -> a stream applied to each document before serialisation, and
// opts.fetchWholeDoc skips the mongo projection so that transform sees the complete document.
// opts.projectionExtras(fl) -> extra fields to FETCH (not to return) when fl is given, for a
// response field derived from another one.
// Used by /genetrees, which stores trees FLAT and rebuilds them on the way out.
function getFactory(collectionPromise, opts) {
  opts = opts || {};
  return function _get(req, res) {
    var params, nonSchemaParams, returnTsv,
      cursorPromise, transformer, mimetype;

    params = _.mapValues(req.swagger.params, 'value');
    nonSchemaParams = _.omit(req.query, _.keys(params));

    returnTsv = params.fl && params.wt == 'tab';

    if(returnTsv) {
      transformer = csvStringify({header:true, delimiter: '\t', columns: params.fl});
      mimetype = 'text/tab-separated-values';
    }
    else if (params.wt == 'bed') {
      if (params.bedFeature == 'gene') {
        params.fl = ['location','_id'];
      }
      else {
        params.fl = ['location','_id','gene_structure'];
      }
      params.wt='json';
      transformer = through2.obj(function(gene, enc, done) {
        this.push(bedify(gene,params));
        done();
      });
      mimetype = 'text/tab-separated-values';
    }
    else {
      transformer = JSONStream.stringify();
      mimetype = 'application/json';
    }
    if (req.swagger.operation.operationId === "genelists" ||
        req.swagger.operation.operationId === "savedviews") {
      nonSchemaParams.uid=0;
    }
    // opts.projectionExtras lets a collection fetch fields it DERIVES a response field from.
    // Applied here on purpose -- after the TSV/bed transformer has captured params.fl, and with
    // the caller's own fl kept in requestedFl for docTransform -- so an extra fetched field adds
    // neither a TSV column nor a key to the JSON the caller asked for.
    var requestedFl = params.fl;
    if (opts.projectionExtras && params.fl && params.fl.length) {
      var extras = opts.projectionExtras(params.fl);
      if (extras && extras.length) { params = _.assign({}, params, {fl: _.union(params.fl, extras)}); }
    }
    const authHeader = req.headers['authorization'];
    if (authHeader && authHeader.startsWith('Bearer ')) {
      const token = authHeader.split(' ')[1];
      getAuth()
      .verifyIdToken(token)
      .then((decodedToken) => {
        nonSchemaParams.uid = decodedToken.uid;
        cursorPromise = mongoHelper.cursorPromise(collectionPromise, opts.fetchWholeDoc ? _.omit(params, 'fl') : params, nonSchemaParams);
        cursorPromise.then(function(cursor) {
          res.contentType(mimetype);
          var src = cursor.stream();
          if (opts.docTransform) { src = src.pipe(opts.docTransform(requestedFl)); }
          src.pipe(transformer).pipe(res);
        });
      })
      .catch((error) => {
        res.status(401).send('Authorization failed');
      });
    } else {
      cursorPromise = mongoHelper.cursorPromise(collectionPromise, opts.fetchWholeDoc ? _.omit(params, 'fl') : params, nonSchemaParams);
      cursorPromise.then(function(cursor) {
        res.contentType(mimetype);
        var src2 = cursor.stream();
        if (opts.docTransform) { src2 = src2.pipe(opts.docTransform(requestedFl)); }
        src2.pipe(transformer).pipe(res);
      });
    }
  }
}

async function deleteList(req, res) {
  const authHeader = req.headers['authorization'];
  if (!authHeader || !authHeader.startsWith('Bearer ')) {
    return res.status(401).send('Authorization header missing or malformed');
  }
  const token = authHeader.split(' ')[1];
  let uid;
  try {
    const decoded = await getAuth().verifyIdToken(token);
    uid = decoded.uid;
  } catch (err) {
    return res.status(401).send('Authorization failed');
  }

  const listId = req.swagger.params.listId.value;
  if (!listId) {
    return res.status(400).send('listId query parameter is required');
  }
  const force = !!(req.swagger.params.force && req.swagger.params.force.value === true);
  try {
    const collection = await mongoCollections.genelists.mongoCollection();
    if (force) {
      // hard delete now (owner-forced). Does NOT touch the genes-core saved_search field —
      // cross-version saved_search propagation/cleanup is handled separately.
      const result = await collection.deleteOne({ _id: listId, uid: uid });
      if (result.deletedCount === 0) {
        return res.status(404).send('Gene list not found or not owned by user');
      }
      return res.json({ message: 'list permanently deleted' });
    }
    // soft delete: mark for deletion. Stays restorable by the owner for 30 days, then the daily
    // cleanup cron purges it (mongo only).
    const result = await collection.updateOne(
      { _id: listId, uid: uid, deletedAt: { $exists: false } },
      { $set: { deletedAt: new Date() } }
    );
    if (result.matchedCount === 0) {
      return res.status(404).send('Gene list not found, not owned by user, or already deleted');
    }
    res.json({ message: 'list marked for deletion (restorable for 30 days)' });
  } catch (err) {
    console.error('deleteList error:', err);
    res.status(500).send('Failed to delete gene list');
  }
}

// restore a soft-deleted list (un-set deletedAt) while it's still within the 30-day window.
async function restoreList(req, res) {
  const authHeader = req.headers['authorization'];
  if (!authHeader || !authHeader.startsWith('Bearer ')) {
    return res.status(401).send('Authorization header missing or malformed');
  }
  const token = authHeader.split(' ')[1];
  let uid;
  try {
    uid = (await getAuth().verifyIdToken(token)).uid;
  } catch (err) {
    return res.status(401).send('Authorization failed');
  }
  const listId = req.swagger.params.listId.value;
  if (!listId) {
    return res.status(400).send('listId query parameter is required');
  }
  try {
    const collection = await mongoCollections.genelists.mongoCollection();
    const result = await collection.updateOne(
      { _id: listId, uid: uid, deletedAt: { $exists: true } },
      { $unset: { deletedAt: "" } }
    );
    if (result.matchedCount === 0) {
      return res.status(404).send('No deleted gene list to restore (already purged or not owned by user)');
    }
    res.json({ message: 'list restored' });
  } catch (err) {
    console.error('restoreList error:', err);
    res.status(500).send('Failed to restore gene list');
  }
}

// dedicated GET handler for /gene_lists. Firebase Bearer auth (optional) -> uid; builds the mongo
// query explicitly (the generic buildQuery only does equality and ignores site/isPublic):
//   active view  -> {site, deletedAt:{$exists:false}} + (isPublic:true) OR (uid) scoping
//   trash view   -> {site, uid, deletedAt:{$exists:true}}   (owner only)
function genelistsHandler(req, res) {
  const params = _.mapValues(req.swagger.params, 'value');   // isPublic, site, rows, includeDeleted
  const trash = params.includeDeleted === 'trash';

  function run(uid) {
    const query = {};
    if (params.site) query.site = params.site;
    if (trash) {
      query.uid = uid;
      query.deletedAt = { $exists: true };
    } else {
      query.deletedAt = { $exists: false };
      if (params.isPublic === true) {
        query.isPublic = true;
      } else {
        query.uid = uid;
      }
    }
    const options = { limit: (params.rows && params.rows !== -1) ? params.rows : 20 };
    mongoCollections.genelists.mongoCollection().then(function (col) {
      res.contentType('application/json');
      col.find(query, options).stream().pipe(JSONStream.stringify()).pipe(res);
    });
  }

  const authHeader = req.headers['authorization'];
  if (authHeader && authHeader.startsWith('Bearer ')) {
    const token = authHeader.split(' ')[1];
    getAuth().verifyIdToken(token)
      .then(function (decoded) { run(decoded.uid); })
      .catch(function () { res.status(401).send('Authorization failed'); });
  } else {
    // anonymous: only public, non-trash lists are meaningful (uid 0 matches nothing private/trash)
    if (trash) return res.status(401).send('Authorization required for deleted lists');
    run(0);
  }
}

// dedicated GET-listing handler for /saved_views (non-hash mode), parallel to genelistsHandler.
// Firebase Bearer auth (optional) -> uid; scopes by site + (public OR owner). Saved views have no
// soft-delete, so there is no trash view. Reads userData1.savedviews (same shared db as genelists).
// Single-hash lookups are handled by getSavedViewByHash (routed by the wrapper in the IIFE above).
function savedViewsListHandler(req, res) {
  const params = _.mapValues(req.swagger.params, 'value');   // site, isPublic, rows
  function run(uid) {
    const query = {};
    if (params.site) query.site = params.site;
    if (params.isPublic === true) {
      query.isPublic = true;
    } else {
      query.uid = uid;
    }
    const options = { limit: (params.rows && params.rows !== -1) ? params.rows : 20 };
    mongoCollections.savedviews.mongoCollection().then(function (col) {
      res.contentType('application/json');
      col.find(query, options).stream().pipe(JSONStream.stringify()).pipe(res);
    });
  }

  const authHeader = req.headers['authorization'];
  if (authHeader && authHeader.startsWith('Bearer ')) {
    const token = authHeader.split(' ')[1];
    getAuth().verifyIdToken(token)
      .then(function (decoded) { run(decoded.uid); })
      .catch(function () { res.status(401).send('Authorization failed'); });
  } else {
    run(0);   // anonymous: only public views are meaningful (uid 0 matches nothing private)
  }
}

async function updateList(req, res) {
  const authHeader = req.headers['authorization'];
  if (!authHeader || !authHeader.startsWith('Bearer ')) {
    return res.status(401).send('Authorization header missing or malformed');
  }
  const token = authHeader.split(' ')[1];
  let uid;
  try {
    const decoded = await getAuth().verifyIdToken(token);
    uid = decoded.uid;
  } catch (err) {
    return res.status(401).send('Authorization failed');
  }

  const listId = req.swagger.params.listId.value;
  const body = req.swagger.params.updates.value || {};
  const updates = {};
  if (typeof body.label === 'string' && body.label.trim().length > 0) {
    updates.label = body.label.trim();
  }
  if (typeof body.isPublic === 'boolean') {
    updates.isPublic = body.isPublic;
  }
  if (Object.keys(updates).length === 0) {
    return res.status(400).send('No valid fields to update (allowed: label, isPublic)');
  }

  try {
    const collection = await mongoCollections.genelists.mongoCollection();
    const result = await collection.updateOne({ _id: listId, uid: uid }, { $set: updates });
    if (result.matchedCount === 0) {
      return res.status(404).send('Gene list not found or not owned by user');
    }
    res.json({ message: 'list updated', updated: updates });
  } catch (err) {
    console.error('updateList error:', err);
    res.status(500).send('Failed to update gene list');
  }
}

// Save a gene list. Takes the resolved identifiers in a JSON body and owns the two side effects
// that /gene_lists/validate used to perform: computing the hash, and tagging the genes core.
// Validate is now a pure lookup, so nothing is written until the user actually saves.
async function saveList(req, res) {
  var authHeader = req.headers['authorization'];
  if (!authHeader || !authHeader.startsWith('Bearer ')) {
    return res.status(401).send('Authorization header missing or malformed');
  }
  var token = authHeader.split(' ')[1];
  var decoded;
  try {
    decoded = await getAuth().verifyIdToken(token);
  } catch (err) {
    return res.status(401).send('Authorization failed');
  }
  var uid = decoded.uid;
  var owner = decoded.name || decoded.email || uid;

  var body = req.body || {};

  // The hash is derived here now. A client still sending one is on the old contract and would
  // otherwise have it silently ignored, so fail loudly instead.
  if (body.hash !== undefined) {
    return res.status(400).send('hash is computed by the server; remove it and post `ids` instead');
  }
  if (typeof body.label !== 'string' || !body.label.trim()) {
    return res.status(400).send('label is required');
  }
  if (typeof body.site !== 'string' || !body.site.trim()) {
    return res.status(400).send('site is required');
  }
  if (!Array.isArray(body.ids) || !body.ids.length) {
    return res.status(400).send('ids must be a non-empty array of gene identifiers');
  }
  if (body.ids.length > MAX_LIST_IDS) {
    return res.status(400).send('too many ids: ' + body.ids.length + ' (limit ' + MAX_LIST_IDS + ')');
  }
  for (var i = 0; i < body.ids.length; i++) {
    if (typeof body.ids[i] !== 'string' || !body.ids[i].length) {
      return res.status(400).send('every id must be a non-empty string');
    }
    if (body.ids[i].length > MAX_ID_LENGTH) {
      return res.status(400).send('id longer than ' + MAX_ID_LENGTH + ' characters');
    }
  }

  var ids = [...new Set(body.ids)].sort();
  // Joined with a separator, unlike the old ids.join('') which collided: ["AB","C"] and
  // ["A","BC"] hashed identically. Existing lists are unaffected — their hashes are stored in
  // genelists and in the cores, and sync_saved_search.js only ever compares stored values.
  var hash = murmur.x86.hash32(ids.join(','));
  var _id = hash + ' ' + uid;

  // Explicit allowlist. The previous version did `$set: params` over every swagger param, so a
  // missing hash produced `_id: "undefined <uid>"` and any future param silently joined the doc.
  var doc = {
    hash: hash,
    label: body.label.trim(),
    site: body.site.trim(),
    isPublic: body.isPublic === true,
    n_genes: ids.length,           // derived, never taken from the client
    uid: uid,
    owner: owner
  };

  try {
    var coll = await mongoCollections.genelists.mongoCollection();
    var result = await coll.updateOne(
      {_id: _id},
      {$set: doc, $setOnInsert: {createdAt: new Date()}},
      {upsert: true}
    );

    // Tag the member genes. Awaited and surfaced, unlike the old fire-and-forget POST that ran
    // after the response had already been sent. add-distinct keeps re-saves idempotent.
    var updates = ids.map(function (id) {
      return {id: id, saved_search: {'add-distinct': hash}};
    });
    await axios.post(solrHelper.genesURL + '/update?commit=true', updates,
                     {headers: {'Content-Type': 'application/json'}});

    res.json({
      message: 'list saved',
      hash: hash,
      _id: _id,
      n_genes: ids.length,
      upserted: !!(result && result.upsertedCount)
    });
  } catch (err) {
    console.error('saveList failed:', err && err.message || err);
    res.status(500).send('Failed to save gene list');
  }
}

// ── saved_views handlers ────────────────────────────────────────────────
//
// Mirror the gene_lists pattern (auth via Firebase Bearer ID token, owner
// scoping via `{_id: viewId, uid}`), but POST takes a JSON body because the
// snapshot blob doesn't fit in query params. The composite _id is
// `${hash} ${uid}` — matching genelists — so re-saving the same content by
// the same user upserts in place.

async function getSavedViewByHash(req, res) {
  var hash = req.swagger.params.hash.value;
  if (!hash) return res.status(400).send('hash query parameter is required');

  // Optional Bearer: anonymous is fine for public views. If a token is
  // present we honor it for private-view access; if invalid we still allow
  // public access (don't 401 a public lookup just because the caller's
  // token expired).
  var uid = null;
  var authHeader = req.headers['authorization'];
  if (authHeader && authHeader.startsWith('Bearer ')) {
    var token = authHeader.split(' ')[1];
    try {
      var decoded = await getAuth().verifyIdToken(token);
      uid = decoded.uid;
    } catch (_) { uid = null; }
  }

  try {
    var collection = await mongoCollections.savedviews.mongoCollection();
    var row = await collection.findOne({hash: hash});
    if (!row) return res.status(404).send('Saved view not found');
    if (!row.isPublic && row.uid !== uid) {
      return res.status(401).send('Private saved view — not authorized');
    }
    res.json(row);
  } catch (err) {
    console.error('getSavedViewByHash error:', err);
    res.status(500).send('Failed to fetch saved view');
  }
}

async function saveView(req, res) {
  var authHeader = req.headers['authorization'];
  if (!authHeader || !authHeader.startsWith('Bearer ')) {
    return res.status(401).send('Authorization header missing or malformed');
  }
  var token = authHeader.split(' ')[1];
  var decoded;
  try {
    decoded = await getAuth().verifyIdToken(token);
  } catch (err) {
    return res.status(401).send('Authorization failed');
  }

  var body = req.swagger.params.body.value || {};
  var hash = body.hash;
  if (!hash || !body.site || !body.label || !body.state) {
    return res.status(400).send('hash, site, label, and state are required');
  }

  var uid = decoded.uid;
  var owner = decoded.name || decoded.email || uid;
  var doc = {
    hash: hash,
    label: body.label,
    description: body.description || '',
    site: body.site,
    isPublic: !!body.isPublic,
    state: body.state,
    uid: uid,
    owner: owner
  };

  try {
    var collection = await mongoCollections.savedviews.mongoCollection();
    var id = hash + ' ' + uid;
    var result = await collection.updateOne(
      {_id: id},
      {$set: doc, $setOnInsert: {createdAt: new Date()}},
      {upsert: true}
    );
    res.json({message: 'view saved', hash: hash, _id: id, upserted: !!result.upsertedCount});
  } catch (err) {
    console.error('saveView error:', err);
    res.status(500).send('Failed to save view');
  }
}

async function updateView(req, res) {
  var authHeader = req.headers['authorization'];
  if (!authHeader || !authHeader.startsWith('Bearer ')) {
    return res.status(401).send('Authorization header missing or malformed');
  }
  var token = authHeader.split(' ')[1];
  var uid;
  try {
    uid = (await getAuth().verifyIdToken(token)).uid;
  } catch (err) {
    return res.status(401).send('Authorization failed');
  }

  var viewId = req.swagger.params.viewId.value;
  var body = req.swagger.params.updates.value || {};
  var updates = {};
  if (typeof body.label === 'string' && body.label.trim().length > 0) {
    updates.label = body.label.trim();
  }
  if (typeof body.isPublic === 'boolean') {
    updates.isPublic = body.isPublic;
  }
  if (Object.keys(updates).length === 0) {
    return res.status(400).send('No valid fields to update (allowed: label, isPublic)');
  }

  try {
    var collection = await mongoCollections.savedviews.mongoCollection();
    var result = await collection.updateOne({_id: viewId, uid: uid}, {$set: updates});
    if (result.matchedCount === 0) {
      return res.status(404).send('Saved view not found or not owned by user');
    }
    res.json({message: 'view updated', updated: updates});
  } catch (err) {
    console.error('updateView error:', err);
    res.status(500).send('Failed to update saved view');
  }
}

async function deleteView(req, res) {
  var authHeader = req.headers['authorization'];
  if (!authHeader || !authHeader.startsWith('Bearer ')) {
    return res.status(401).send('Authorization header missing or malformed');
  }
  var token = authHeader.split(' ')[1];
  var uid;
  try {
    uid = (await getAuth().verifyIdToken(token)).uid;
  } catch (err) {
    return res.status(401).send('Authorization failed');
  }

  var viewId = req.swagger.params.viewId.value;
  if (!viewId) return res.status(400).send('viewId query parameter is required');

  try {
    var collection = await mongoCollections.savedviews.mongoCollection();
    var result = await collection.deleteOne({_id: viewId, uid: uid});
    if (result.deletedCount === 0) {
      return res.status(404).send('Saved view not found or not owned by user');
    }
    res.json({message: 'view deleted'});
  } catch (err) {
    console.error('deleteView error:', err);
    res.status(500).send('Failed to delete saved view');
  }
}