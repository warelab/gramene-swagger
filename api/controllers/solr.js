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

async function validate(req, res) {
  // get list of ids from request body
  // uniqify and sort ids
  // check if this set of ids has been validated before
  // split list into manageable chunks
  // validate the ids
  // create a hash key for the validated ids
  const solrURL = solrHelper.genesURL + '/select';

  var ids = req.body.ids;
  let uniqueIdentifiers = [...new Set(ids)].sort();

  // check if this set of ids has been validated already
  async function checkHash(hashed_input) {
    try {
      const params = {q:`saved_search:${hashed_input}`,rows:0};
      const response = await axios.get(solrURL, {params});
      // console.log(response.data.response.numFound, uniqueIdentifiers.length);
      return response.data.response.numFound === uniqueIdentifiers.length;
    } catch (error) {
      console.error('error querying solr for saved_search:', error);
    }
  }
  
  const hashed_input = murmur.x86.hash32(uniqueIdentifiers.join(''));
  const isSaved = await checkHash(hashed_input);

  if (isSaved) {
    // console.log("saved search found!");
    res.json({
      ids: uniqueIdentifiers,
      missing: [],
      hash: hashed_input
    });
    return;
  }
  // console.log("saved search not found");

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

    // build a hash for the found ids
    const murmurhash = murmur.x86.hash32(foundIds.sort().join('')); // sorting because async
    res.json({
      ids: foundIds,
      missing: missingIds,
      hash: murmurhash
    });
    // also do an atomic update to the gene docs
    const updateData = foundIds.map(id => {
      return {
        id: id,
        saved_search: { add: murmurhash }
      }
    });
    // console.log("posting atomic updates", updateData);
    axios.post(solrHelper.genesURL + '/update?commit=true', updateData, {
      headers: { 'Content-Type': 'application/json' }
    })
    .then(response => {
      // console.log('Update successful:', response.data)
    })
    .catch(error => {
      console.error('Error while updating:', error);
    });
  }
  checkIdentifiers();
}

function saveSearch(req, res) {
  // given
  // hash: 32bit integer
  // name: string name for saved search
  // public: boolean
  // site: site where the saved search comes from
  // authenticated user
  //
  // insert a record to mongodb collection
}