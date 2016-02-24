'use strict';

var _ = require('lodash');
var mongoHelper = require('../helpers/mongo');
var mongoCollections = require('gramene-mongodb-config');
var JSONStream = require('JSONStream');
var csvStringify = require('csv-stringify');

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

  module.exports = toExport;
}());

function getFactory(collectionPromise) {
  return function _get(req, res) {
    var params, nonSchemaParams, returnCsv,
      cursorPromise, transformer, mimetype;

    params = _.mapValues(req.swagger.params, 'value');
    nonSchemaParams = _.omit(req.query, _.keys(params));
    returnCsv = params.fl && params.wt == 'tab';

    cursorPromise = mongoHelper.cursorPromise(collectionPromise, params, nonSchemaParams);
    if(returnCsv) {
      transformer = csvStringify({header:true, delimiter: '\t', columns: params.fl});
      mimetype = 'text/tab-separated-values';
    }
    else {
      transformer = JSONStream.stringify();
      mimetype = 'application/json';
    }

    cursorPromise.then(function(cursor) {
      res.contentType(mimetype);
      cursor.stream().pipe(transformer).pipe(res);
    });
  }
}