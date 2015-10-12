'use strict';

var request = require('request');

var urlBase = 'http://brie:8946/solr/';

module.exports = {
  streamGenes: streamGenes,
  streamSuggestions: streamSuggestions
};

function streamGenes(params) {
  return solrStream(urlBase + 'newgenes/query', params);
}

function streamSuggestions(params) {
  return solrStream(urlBase + 'suggestions/try', params);
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

