'use strict';

var _ = require('lodash');
var solrHelper = require('../helpers/solr');

module.exports = {
  genes: genes,
  suggestions: suggestions
};

function genes(req, res) {
  solrRequest(req, res, solrHelper.streamGenes);
}

function suggestions(req, res) {
  solrRequest(req, res, solrHelper.streamSuggestions);
}

function solrRequest(req, res, streamMethod) {
  var params = _.mapValues(req.swagger.params, 'value');
  var solrStream = streamMethod(params);

  solrStream.pipe(res);
}

