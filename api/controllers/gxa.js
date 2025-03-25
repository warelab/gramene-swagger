'use strict';
var mongo = require('gramene-mongodb-config');

module.exports = {
  baseline_refexperiment: baseline_refexperiment,
  baseline_experiments: baseline_experiments
};

function fetchBy(collectionPromise, ids, idField) {
  return collectionPromise.then(coll => {
    let query = {};
    query[idField] = { '$in': ids };
    return coll.find(query, {}, { rows: -1 }).toArray();
  });
}

async function baseline_experiments(req, res) {
  try {
    // Get list of ids from request body
    const ids = req.body.replace("geneQuery=", "").split(" ");
    
    // Uniquify ids
    let uniqueIdentifiers = [...new Set(ids)];
    // Fetch gene expression from MongoDB
    let genes = await fetchBy(mongo.expression.mongoCollection(), uniqueIdentifiers, '_id');
    // get baseline experiment assays with non-zero tpm
    var tpm = {};
    genes.forEach(g => {
      Object.keys(g).forEach(exp => {
        if (exp !== "_id") {
          g[exp].forEach(a => {
            if (a.hasOwnProperty('value') && a.value >= 0) {
              if (!tpm.hasOwnProperty(exp)) {
                tpm[exp] = {};
              }
              tpm[exp][a.group] = a.value < 0.5 ? 0 : a.value;
            }
          })
        }
      })
    })
    let experimentIdentifiers = Object.keys(tpm);

    // Fetch experiments from MongoDB
    let experiments = await fetchBy(mongo.experiments.mongoCollection(), experimentIdentifiers, '_id');
    let exp_lut = {};
    experiments.forEach(e => {
      exp_lut[e._id] = e
    })
    // Fetch assays from MongoDB
    let assays = await fetchBy(mongo.assays.mongoCollection(), experimentIdentifiers, 'experiment');

    // get the organism_part ontology terms for the anatomogram
    let anatomogram_parts = new Set();
    let columns = {};
    let groups = {};
    let species;
    assays.forEach(assay => {
      if (!species) {
        assay.characteristic.forEach(c => {
          if (c.type === "organism") {
            species = c.label.replace(" ","_").toLowerCase();
          }
        })
      }
      assay.factor.forEach(factor => {
        if (factor.type === "organism part") {
          if (tpm[assay.experiment] && tpm[assay.experiment].hasOwnProperty(assay.group)) {
            if (! groups[assay.experiment]) {
              groups[assay.experiment] = {};
            }
            if (! groups[assay.experiment][factor.label]) {
              groups[assay.experiment][factor.label] = [];
            }
            groups[assay.experiment][factor.label].push(assay);
            if (! columns.hasOwnProperty(factor.label)) {
              columns[factor.label] = {};
            }
            if (factor.ontology) {
              const termId = factor.id.replace(":","_");
              anatomogram_parts.add(termId);
              columns[factor.label].factorValueOntologyTermId = termId;
            }
          }
        }
      })
    });
    // sort the columns
    const colNames = Object.keys(columns).sort();
    let rows = {};
    // look for experiments that need multiple rows
    Object.keys(groups).forEach(eid => {
      colNames.forEach((opart,index) => {
        if (groups[eid][opart]) {
          if (groups[eid][opart].length > 1) {
            groups[eid][opart].forEach(assay => {
              let flabels = assay.factor.filter(factor => factor.type !== "organism part").map(factor => factor.label);
              let ef = `${eid} - ${flabels.join(' - ')}`;
              if (exp_lut[eid].name) {
                ef = `${exp_lut[eid].name} - ${flabels.join(' - ')}`;
              }
              if (! rows.hasOwnProperty(ef)) {
                rows[ef] = {
                  id: ef,
                  name: ef,
                  experimentType: "RNASEQ_MRNA_BASELINE",
                  expressionUnit: "TPM",
                  expressions: Array.from({ length: colNames.length }, () => ({})),
                  uri: `experiments/${eid}?geneQuery=%5B%7B%22value%22%3A%22${uniqueIdentifiers[0]}%22%7D%5D`
                }
                // prepend source of experiment (EBI, JGI, etc.)
                rows[ef].name = `${exp_lut[eid].source} - ${rows[ef].name}`;
              }
              rows[ef]['expressions'][index] = {value: tpm[assay.experiment][assay.group]}
            })
          }
          else {
            const assay = groups[eid][opart][0];
            if (! rows.hasOwnProperty(eid)) {
              rows[eid] = {
                id: eid,
                name: exp_lut[eid].name || eid,
                experimentType: "RNASEQ_MRNA_BASELINE",
                expressionUnit: "TPM",
                expressions: Array.from({ length: colNames.length }, () => ({})),
                uri: `experiments/${eid}?geneQuery=%5B%7B%22value%22%3A%22${uniqueIdentifiers[0]}%22%7D%5D`
              }
                // prepend source of experiment (EBI, JGI, etc.)
              rows[eid].name = `${exp_lut[eid].source} - ${rows[eid].name}`;
            }
            rows[eid]['expressions'][index] = {value: tpm[assay.experiment][assay.group]}
          }
        }
      })
    });
    let sortedRows = Object.values(rows).sort((a,b) => a.id.localeCompare(b.id));
    let result = {
      anatomogram: {
        species: species,
        allSvgPathIds: [...anatomogram_parts]
      },
      columnHeaders: colNames.map(k => {
        let v = columns[k];
        let ch = {
          assayGroupId: k,
          factorValue: k
        };
        if (v.factorValueOntologyTermId) {
          ch.factorValueOntologyTermId = v.factorValueOntologyTermId
        }
        return ch;
      }),
      columnGroupings:[],
      profiles: {
        rows: sortedRows,
        searchResultTotal: `${sortedRows.length}`
      },
      config: {
        columnType: "",
        conditionQuery: "%5B%5D",
        disclaimer: "",
        expressionUnit: "",
        geneQuery: "here",
        genomeBrowsers: [],
        species: "sorghum bicolor"
      }
    };
    // Send response
    res.json(result);
  } catch (error) {
    console.error("Error in baseline_experiments:", error);
    res.status(500).json({ error: "Internal Server Error" });
  }
}

async function baseline_refexperiment(req, res) {
  // get list of ids from request body
  var ids = req.body.geneQuery;
  // uniqify ids
  let uniqueIdentifiers = [...new Set(ids)];

  // fetch expression data from mongodb
  mongo.expression.mongoCollection().then(function(expr) {
    var options = {rows:-1};
    var query = {'_id': {'$in': uniqueIdentifiers}};
    expr.find(query,options).then(function(genes) {
    })
  });
}
