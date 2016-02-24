'use strict';

var _ = require('lodash');
var solrHelper = require('../helpers/solr');

module.exports = {
  genes: genes,
  suggestions: suggestions
};

function genes(req, res) {
  solrRequest(req, res, solrHelper.streamGenes, {facet: true});
}

function suggestions(req, res) {
  solrRequest(req, res, solrHelper.streamSuggestions);
}

function solrRequest(req, res, streamMethod, additionalParams) {
  var params = _.mapValues(req.swagger.params, changeVerticalBarToComma);

  if(additionalParams) {
    _.assign(params, additionalParams);
  }
  var solrStream = streamMethod(params);

  solrStream.pipe(res);
}

function changeVerticalBarToComma(param) {
  var result;

  function _changeVerticalBarToComma(_param) {
    if(_.isString(_param)) {
      return _param.replace(/\|/g, ',');
    }
  }

  if(_.isArray(param.value)) {
    result = param.value.map(_changeVerticalBarToComma);
  }
  else {
    result = _changeVerticalBarToComma(param.value);
  }

  return result;
}



