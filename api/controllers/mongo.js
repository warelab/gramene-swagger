'use strict';

var _ = require('lodash');
var mongoHelper = require('../helpers/mongo');
var mongoCollections = require('gramene-mongodb-config');
var JSONStream = require('JSONStream');
var csvStringify = require('csv-stringify');
var bedify = require('gramene-bedify');
var through2 = require('through2');

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