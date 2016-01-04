'use strict';

var request = require('request');
var version = require('gramene-mongodb-config').getMongoConfig().version;
var through2 = require('through2');
var csv2 = require('csv2');

var urlBase = 'http://brie:8983/solr/';

module.exports = {
  streamGenes: streamGenes,
  streamSuggestions: streamSuggestions
};

function streamGenes(params) {
  if (params.wt === 'bed') {
    params.fl = ['region,start,end,strand,id,taxon_id'];
    params.wt = 'csv';
    params.sort = 'gene_idx asc';
    params.isBed = true;
  }
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
    var type;
    switch(params.wt) {
      case 'xml':
        type = 'application/xml';
        break;
      case 'json':
        type = 'application/json';
        break;
      case 'bed':
        type = 'text/tab-separated-values';
        break;
      default:
        type = 'text/plain';
    }

    r.headers['content-type'] = type;
  });

  stream.on('error', function error(err) {
    console.log(err);
  });
  
  if (params.isBed) {
    return stream.pipe(csv2()).pipe(through2.obj(function (values, enc, done) {
      if (values[0] !== 'region') {
        var strand = (values[3] === '1') ? '+' : '-';
        var start = values[1] - 1;
        this.push(values[0] // chromosome (region)
          +'\t'+start      // start
          +'\t'+values[2]   // end
          +'\t'+values[4]   // name
          +'\t'+values[5]   // score (taxon_id)
          +'\t'+strand     // strand
          +'\n');
      }
      done();
    }));
  }
  return stream;
}

