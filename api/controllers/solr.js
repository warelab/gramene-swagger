'use strict';

var _ = require('lodash');
var solrHelper = require('../helpers/solr');

module.exports = {
  genes: genes,
  suggestions: suggestions
};

function genes(req, res) {
  console.error("genes",req.swagger.params);
  solrRequest(req, res, solrHelper.streamGenes, {facet: true, stats: true});
}

function suggestions(req, res) {
  solrRequest(req, res, solrHelper.streamSuggestions);
}

function solrRequest(req, res, streamMethod, additionalParams) {
  var params = _.mapValues(req.swagger.params, 'value');
  if(additionalParams) {
    _.assign(params, additionalParams);
  }
  var solrStream = streamMethod(params);

  solrStream.pipe(res);
}

