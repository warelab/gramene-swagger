'use strict';

var _ = require('lodash');

const urlBase = 'https://data.sorghumbase.org/sorghum_v6';
var request = require('request');
var JSONStream = require('JSONStream');
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
  // https://data.sorghumbase.org/sorghum_v6/search?fl=id,region,start,end,system_name&q=system_name:sorghum_bicolor%20AND%20region:4%20AND%20start:[1%20TO%20100000]
  console.error("exons",params);
  let url = `${urlBase}/search?fl=id,region,start,end,system_name&q=`;
  url = `${url}system_name:${params.genome}`;
  url = `${url} AND region:${params.region}`;
  url = `${url} AND start:[${params.start} TO ${params.end}] AND end:[${params.start} TO ${params.end}]`;

  let searchParams = {
    system_name: params.genome,
    region: params.region,
    start: `[${params.start} TO ${params.end}]`,
    end: `[${params.start} TO ${params.end}]`
  };

  var stream = solrStream(url, {});
  stream.pipe(res);
//  res.json(params);
}
function solrStream(uri, params) {
  var stream = request({uri: uri, qs: params, qsStringifyOptions: { arrayFormat: 'repeat' }});

  // modify content-type of incoming stream header;
  // A call to res.contentType() in the controller
  // is overwritten by response stream from SOLR
  stream.on('response', function stripHeaders(r) {
    var type = 'text/plain';

    r.headers['content-type'] = type;
  });

  stream.on('error', function error(err) {
    console.log(err);
  });
  
  return stream.pipe(JSONStream.parse('response.docs.*',function(doc) {
    console.error(doc);
    var features = "# fjFile = FEATURES\n";
    var featureCols = [doc.id, doc.region, doc.start, doc.end];
    return featureCols.join('\t')+'\n';
  }));

  return stream;
}

module.exports = {
  map: markers,
  genotypes: genotypes,
  features: exons
};
