'use strict';

var _ = require('lodash');
var mongoHelper = require('../helpers/mongo');
var mongoCollections = require('gramene-mongodb-config');
var JSONStream = require('JSONStream');
var csvStringify = require('csv-stringify');

//module.exports = {
//  get: get
//};

// add a function to the controller for each mongo collection.
_.forOwn(mongoCollections, function(coll, name) {
  module.exports[name] = getFactory(coll.mongoCollection());
});

function getFactory(collectionPromise) {
  return function get(req, res) {
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