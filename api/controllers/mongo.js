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

var serviceAccount = require("/usr/local/gramene/gramene-auth-firebase-adminsdk-c1sc0-263ff4cc4f.json");

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

    cursorPromise = mongoHelper.cursorPromise(collectionPromise, params, nonSchemaParams);
    cursorPromise.then(function(cursor) {
      res.contentType(mimetype);
      cursor.stream().pipe(transformer).pipe(res);
    });
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
      // upsert to mongo collection if we have all the params
      mongoCollections.genelists.mongoCollection().then(function(mongo) {
        const id = `${params.hash} ${params.uid}`;
        mongo.updateOne({ _id: id }, { $set:params }, { upsert: true }).then(function(result) {
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