'use strict';

var _ = require('lodash');

function markers(req, res) {
  var params = _.mapValues(req.swagger.params, 'value');
  console.error("markers",params);
  res.json(params);
}

function genotypes(req, res) {
  var params = _.mapValues(req.swagger.params, 'value');
  console.error("genotypes",params);
  res.json(params);
}

function exons(req, res) {
  var params = _.mapValues(req.swagger.params, 'value');
  console.error("exons",params);
  res.json(params);
}

module.exports = {
  map: markers,
  genotypes: genotypes,
  features: exons
};
