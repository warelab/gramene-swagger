'use strict';

var _ = require('lodash');
var solrHelper = require('../helpers/solr');
const murmur = require('murmurhash3js');
const hash32 = murmur.x86.hash32('example');
// console.log(hash32); // Compact 32-bit hash
const axios = require('axios');

module.exports = {
  genes: genes,
  suggestions: suggestions,
  validate: validate
};

function genes(req, res) {
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

function validate(req, res) {
  var ids = req.body.ids;
  let uniqueIdentifiers = [...new Set(ids)].sort();
  console.log(uniqueIdentifiers);

  const chunkSize = 200;

  // Function to split the list into chunks
  function chunkArray(array, size) {
    const chunks = [];
    for (let i = 0; i < array.length; i += size) {
      chunks.push(array.slice(i, i + size));
    }
    return chunks;
  }

  // Split the unique identifiers into chunks
  const identifierChunks = chunkArray(uniqueIdentifiers, chunkSize);


  const solrURL = 'http://squam:8983/solr/sorghum_genes7/select';

  async function querySolr(chunk) {
    const params = {
      q: `id:(${chunk.join(' OR ')})`,
      fl: 'id',
      rows: chunk.length
    };
    
    try {
      const response = await axios.get(solrURL, {params});
      return response.data.response.docs.map(doc => doc.id);
    } catch (error) {
      console.error('error querying solr:', error);
      return [];
    }
  }
 
  async function checkIdentifiers() {
    let foundIds = [];
  
    for (const chunk of identifierChunks) {
      const result = await querySolr(chunk);
      foundIds = foundIds.concat(result);
    }

    // Compare the found IDs with the original list
    const missingIds = uniqueIdentifiers.filter(id => !foundIds.includes(id));

    console.log('Found IDs:', foundIds);
    console.log('Missing IDs:', missingIds);
    
    res.json({
      ids: foundIds,
      missing: missingIds,
      hash: murmur.x86.hash32(foundIds.join(''))
    });
  }

  checkIdentifiers();

  // console.log(solrURL, params);
  // axios.get(solrURL, { params })
  // .then(response => {
  //   const docs = response.data.response.docs;
  //   console.log('docs:', docs);
  //   const foundIds = response.data.response.docs.map(doc => doc.id);
  //   console.log('Found IDs:', foundIds);
  //
  //   // Compare the found IDs with the original list
  //   const missingIds = uniqueSortedIds.filter(id => !foundIds.includes(id));
  //   console.log('Missing IDs:', missingIds);
  //   res.json("done");
  // })
  // .catch(error => {
  //   console.error('Error querying Solr:', error);
  // });
  //

}