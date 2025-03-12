'use strict';

const solrHelper = require("../helpers/solr");
const axios = require("axios");
const murmur = require("murmurhash3js");
module.exports = {
  baseline_refexperiment: baseline_refexperiment,
  baseline_experiments: baseline_experiments
};

async function baseline_experiments(req, res) {
  // get list of ids from request body
  // uniqify and sort ids
  const solrURL = solrHelper.genesURL + '/select';

  var ids = req.body.geneQuery;
  let uniqueIdentifiers = [...new Set(ids)].sort();

}
