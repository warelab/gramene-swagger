'use strict';

var _ = require('lodash');
var mongoHelper = require('../helpers/mongo');
var JSONStream = require('JSONStream');
var csvStringify = require('csv-stringify');

module.exports = {
  get: get
};

function get(req, res) {
  var params = _.mapValues(req.swagger.params, 'value');
  var nonSchemaParams = _.omit(req.query, _.keys(params));
  var returnCsv = params.fl && params.wt == 'tab';

  var cursorPromise = mongoHelper.cursorPromise(params, nonSchemaParams);
  var transformer, mimetype;
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