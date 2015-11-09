'use strict';

var request = require('request');
var version = require('gramene-mongodb-config').getMongoConfig().version;

var urlBase = 'http://brie:8983/solr/';

module.exports = {
  streamGenes: streamGenes,
  streamSuggestions: streamSuggestions
};

function streamGenes(params) {
  return solrStream(urlBase + 'genes' + version + '/query', params);
}

function streamSuggestions(params) {
  return solrStream(urlBase + 'suggestions' + version + '/try', params);
}

function solrStream(uri, params) {
  var stream = request({uri: uri, qs: params, qsStringifyOptions: { arrayFormat: 'repeat' }});

  // modify content-type of incoming stream header;
  // A call to res.contentType() in the controller
  // is overwritten by response stream from SOLR
  stream.on('response', function stripHeaders(r) {
    r.headers['content-type'] = 'application/json';
  });

  stream.on('error', function error(err) {
    console.log(err);
  });

  return stream;
}

