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
    if (req.swagger.operation.operationId === "genelists") {
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
  try {
    const collection = await mongoCollections.genelists.mongoCollection();
    const result = await collection.deleteOne({ _id: listId, uid: uid });
    if (result.deletedCount === 0) {
      return res.status(404).send('Gene list not found or not owned by user');
    }
    res.json({ message: 'list deleted' });
  } catch (err) {
    console.error('deleteList error:', err);
    res.status(500).send('Failed to delete gene list');
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