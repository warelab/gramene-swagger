'use strict';

var sanitize = require('mongo-sanitize')
  , isNumber = require("isnumber")
  , collections = require('gramene-mongodb-config')
  , _ = require('lodash');

module.exports = {
  cursorPromise: cursorPromise
};

function cursorPromise(params, nonSchemaParams) {
  var mongoCollectionPromise = collections[params.collection].mongoCollection();
  var query = buildQuery(params, nonSchemaParams);
  var options = {};

  if (params.rows) {
    if (params.rows !== -1) options.limit = params.rows;
  }
  else {
    options.limit = 20;
  }

  if (params.start) options.skip = params.start;
  if (params.fl) {
    options.fields = {_id:0};
    params.fl.forEach(function (f) {
      options.fields[f] = 1;
    });
  }

  return mongoCollectionPromise.then(function(collection) {
    return collection.find(query, options);
  });
}

function buildQuery(params, nonSchemaParams) {
  var qExprs = [], result;

  nonSchemaParams = _.mapValues(nonSchemaParams, sanitize);

  // free text query
  if (params.q) qExprs.push({'$text': {'$search': params['q']}});

  // location query
  if (params.l) {
    var a = params.l.split(':');
    var t = a[0], r = a[1], s = a[2], e = a[3];
    qExprs.push(
      {'location.taxon_id': t},
      {'location.region': r},
      {'location.start': {'$lte': +e}},
      {'location.end': {'$gte': +s}});
  }

  if (params.idList) {
    var list = params.idList.map(function (x) {
        return isNumber(x) ? +x : x;
      });
    qExprs.push({'_id': {'$in': list}});
  }

  function maybeCastToNumber(x) {
    return isNumber(x) ? +x : x;
  }

  _.forEach(nonSchemaParams, function(value, key) {
    var expr = {};
    if (Array.isArray(value)) {
      expr[key] = {'$in': value.map(maybeCastToNumber)};
    }
    else {
      expr[key] = maybeCastToNumber(value);
    }
    qExprs.push(expr);
  });

  if (qExprs.length > 1) {
    result = {'$and': qExprs};
  }
  else if (qExprs.length == 1) {
    result = qExprs[0];
  }
  else {
    result = {};
  }

  return result;
}


