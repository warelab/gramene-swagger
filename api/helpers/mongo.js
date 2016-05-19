'use strict';

var sanitize = require('mongo-sanitize')
  , isNumber = require("isnumber")
  , _ = require('lodash');

module.exports = {
  cursorPromise: cursorPromise
};

function cursorPromise(mongoCollectionPromise, params, nonSchemaParams) {
  var query, options;
  query = buildQuery(params, nonSchemaParams);
  options = {};

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
      {'taxon_id': t},
      {'location.region': r},
      {'location.start': {'$lte': +e}},
      {'location.end': {'$gte': +s}});
  }

  if (params.idList) {
    var idSet={};
    params.idList.forEach(function(x) {
      idSet[x] = 1;
    });
    var list = params.idList.map(function (x) {
      if (idSet[x]++ === 1) {
        return isNumber(x) ? +x : x;
      }
    });
    qExprs.push({'_id': {'$in': list}});
  }

  function maybeCastToNumber(x) {
    return isNumber(x) ? +x : x;
  }

  _.merge(nonSchemaParams,_.pick(params,['taxon_id','db_type', 'subset']));

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


