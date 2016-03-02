'use strict';

var request = require('request');
var version = require('gramene-mongodb-config').getMongoConfig().version;
var through2 = require('through2');
var csv2 = require('csv2');
var JSONStream = require('JSONStream');

var urlBase = 'http://brie.cshl.edu:8983/solr/';

module.exports = {
  streamGenes: streamGenes,
  streamSuggestions: streamSuggestions
};

var bedFields = {
  region   : 1,
  start    : 1,
  end      : 1,
  id       : 1,
  gene_idx : 1,
  strand   : 1
};

function streamGenes(params) {
  if (params.wt === 'bed') {
    if (params.fl) {
      var addFields = [];
      params.fl.forEach(function(field) {
        if (!bedFields[field]) {
          addFields.push(field);
        }
      });
      params.fl = 'region,start,end,id,gene_idx,strand';
      if (addFields.length > 0) {
        params.fl += ',' + addFields.join(',');
      }
    }
    else {
      params.fl = 'region,start,end,id,gene_idx,strand';
    }
    params.sort = 'gene_idx asc';
    params.isBed = true;
    params.omitHeader=true;
    return solrStream(urlBase + 'genes' + version + '/export', params);
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
    return stream.pipe(JSONStream.parse('response.docs.*',function(doc) {
      doc.start--;
      doc.strand = (doc.strand === 1) ? '+' : '-';
      var bedCols = [doc.region, doc.start, doc.end, doc.gene_idx, doc.id, doc.strand];
      for (var field in doc) {
        if (!bedFields[field]) {
          bedCols.push(doc[field]);
        }
      }
      return bedCols.join('\t')+'\n';
    }));
  }
  return stream;
}

