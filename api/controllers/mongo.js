'use strict';

var _ = require('lodash');
var mongoHelper = require('../helpers/mongo');
var mongoCollections = require('gramene-mongodb-config');
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

function getFactory(collectionPromise) {
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
    const authHeader = req.headers['authorization'];
    if (authHeader && authHeader.startsWith('Bearer ')) {
      const token = authHeader.split(' ')[1];
      getAuth()
      .verifyIdToken(token)
      .then((decodedToken) => {
        nonSchemaParams.uid = decodedToken.uid;
        cursorPromise = mongoHelper.cursorPromise(collectionPromise, params, nonSchemaParams);
        cursorPromise.then(function(cursor) {
          res.contentType(mimetype);
          cursor.stream().pipe(transformer).pipe(res);
        });
      })
      .catch((error) => {
        res.status(401).send('Authorization failed');
      });
    } else {
      cursorPromise = mongoHelper.cursorPromise(collectionPromise, params, nonSchemaParams);
      cursorPromise.then(function(cursor) {
        res.contentType(mimetype);
        cursor.stream().pipe(transformer).pipe(res);
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

async function saveList(req, res) {
  let params = _.mapValues(req.swagger.params, 'value');

  const authHeader = req.headers['authorization'];
  if (authHeader && authHeader.startsWith('Bearer ')) {
    const token = authHeader.split(' ')[1];
    getAuth()
    .verifyIdToken(token)
    .then((decodedToken) => {
      params.uid = decodedToken.uid;
      params.owner = decodedToken.name || decodedToken.email || decodedToken.uid;
      // upsert to mongo collection if we have all the params
      mongoCollections.genelists.mongoCollection().then(function(mongo) {
        const id = `${params.hash} ${params.uid}`;
        mongo.updateOne(
          { _id: id },
          { $set: params, $setOnInsert: { createdAt: new Date() } },
          { upsert: true }
        ).then(function(result) {
          res.json({message:'list saved'});
        })
      })
    })
    .catch((error) => {
      res.status(401).send('Authorization failed');
    });
  } else {
    res.status(401).send('Authorization header missing or malformed');
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