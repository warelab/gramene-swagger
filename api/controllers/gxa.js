'use strict';
var mongo = require('gramene-mongodb-config');

module.exports = {
  one_baseline_experiment: baseline_experiment,
  baseline_experiments: baseline_experiments,
  redirectToEbi: function(req, res) {
    const target = 'https://www.ebi.ac.uk' + req.originalUrl.replace('/auth_testing', '');
    res.redirect(301, target);
  }
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
    const searchArray = uniqueIdentifiers.map(id => {value: id});
    const encodedSearch = encodeURIComponent(JSON.stringify(searchArray));
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
            species = c.label.split(" ").slice(0,2).join("_").toLowerCase();
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
                  uri: `${atlasURL}/experiments/${eid}?geneQuery=%5B%7B%22value%22%3A%22${uniqueIdentifiers[0]}%22%7D%5D`
                }
                // prepend source of experiment (EBI, JGI, etc.)
                // rows[ef].name = `${exp_lut[eid].source} - ${rows[ef].name}`;
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
                uri: `${atlasURL}/experiments/${eid}?geneQuery=%5B%7B%22value%22%3A%22${uniqueIdentifiers[0]}%22%7D%5D`
              }
                // prepend source of experiment (EBI, JGI, etc.)
              // rows[eid].name = `${exp_lut[eid].source} - ${rows[eid].name}`;
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
        geneQuery: encodedSearch,
        genomeBrowsers: [],
        species: species.replace("_"," ")
      }
    };
    // Send response
    res.json(result);
  } catch (error) {
    console.error("Error in baseline_experiments:", error);
    res.status(500).json({ error: "Internal Server Error" });
  }
}
const atlasURL = 'https://www.ebi.ac.uk/gxa';
function differential_experiment(res, experiment, assays, genes, species) {
  const assayMap = new Map(assays.map(assay => [assay.group, assay]));
  let expressions = {};
  let sigGroupsSet = new Set();
  genes.forEach(g => {
    g[experiment._id].forEach(a => {
      if (a.hasOwnProperty('p_value') && a.p_value < 0.05) {
        sigGroupsSet.add(a.group);
        if(!expressions.hasOwnProperty(g._id)) {
          expressions[g._id]={};
        }
        expressions[g._id][a.group] = {
          pValue: a.p_value,
          foldChange: a.l2fc
        }
      }
    })
  })
  function pickName(ref,test) {
    let factor = {};
    ref.factor.forEach(f => {
      (factor[f.type] ||= {}).ref = f.label;
    });
    test.factor.forEach(f => {
      (factor[f.type] ||= {}).test = f.label;
    });
    let characteristic = {};
    ref.characteristic.forEach(f => {
      if (!factor[f.type]) {
        (characteristic[f.type] ||= {}).ref = f.label;
      }
    });
    test.characteristic.forEach(f => {
      if (!factor[f.type]) {
        (characteristic[f.type] ||= {}).test = f.label;
      }
    });
    // const fs = Object.keys(factor).filter(t => factor[t].ref && factor[t].test && factor[t].ref !== factor[t].test);
    // const cs = Object.keys(characteristic).filter(t => characteristic[t].ref && characteristic[t].test && characteristic[t].ref !== characteristic[t].test);
    const fs = Object.keys(factor).filter(t => factor[t].ref !== factor[t].test);
    const cs = Object.keys(characteristic).filter(t => characteristic[t].ref !== characteristic[t].test);
    const contrasts = fs.concat(cs).map(t => {
      const refLabel = factor[t] ? factor[t].ref : characteristic[t].ref;
      const testLabel = factor[t] ? factor[t].test : characteristic[t].test;
      return `${t}: ${refLabel} vs ${testLabel}`
    }).join(' and ');
    const sames = Object.keys(factor).filter(t => factor[t].ref && factor[t].test && factor[t].ref === factor[t].test).map(t => `${t}: ${factor[t].ref}`).join(' and ');
    return sames !== "" ? [contrasts,sames].join(" in ") : contrasts;
  }
  const sigGroups = [...sigGroupsSet];
  let de = {
    experiment: {
      accession: experiment._id,
      type: "rnaseq_mrna_differential",
      urls: {
        main_page: `${atlasURL}/experiments/${experiment._id}?geneQuery=%255B%255D`,
        genome_browsers: `${atlasURL}/experiments/${experiment._id}/redirect/genome-browsers`,
        download: `${atlasURL}/experiments-content/${experiment._id}/download/RNASEQ_MRNA_DIFFERENTIAL?geneQuery=%5B%5D&unit=FOLD_CHANGE&cutoff=0.05&heatmapMatrixSize=50&selectedColumnIds=&type=RNASEQ_MRNA_DIFFERENTIAL`
      },
      description: experiment.description,
      species: species
    },
    config: {
      geneQuery: "%5B%5D",
      species: species,
      genomeBrowsers: ["Ensembl Genomes"],
      disclaimer: '',
      columnType: '',
      conditionQuery: ''
    },
    columnGroupings: [],
    columnHeaders: sigGroups.map(sg => {
      const [refGroup, testGroup] = sg.split("_");
      const ref = assayMap.get(refGroup);
      const test = assayMap.get(testGroup);
      const isFactor = new Map(ref.factor.map(f => [f.type, true]));
      const contrastName = pickName(ref,test);
      let properties = ref.characteristic.map((prop,i) => {
        return {
          propertyName: prop.type,
          testValue: test.characteristic[i].label,
          contrastPropertyType: isFactor.get(prop.type) ? "FACTOR" : "SAMPLE",
          referenceValue: prop.label
        }
      });
      ref.factor.forEach((prop,i) => {
        properties.push({
          propertyName: prop.type,
          testValue: test.factor[i].label,
          contrastPropertyType: 'FACTOR',
          referenceValue: prop.label
        });
      })
      return {
        id: sg,
        displayName: contrastName,
        referenceAssayGroup: {
          id: refGroup,
          assayAccessions:['a','b','c'],
          replicates: 3
        },
        testAssayGroup: {
          id: testGroup,
          assayAccessions: ['a','b','c'],
          replicates: 3
        },
        contrastSummary: {
          properties: properties,
          experimentDescription: experiment.description,
          contrastDescription: contrastName,
          testReplicates: 3,
          referenceReplicates: 3
        },
        resources: [
          {
            type: "gsea_go",
            uri: `${atlasURL}/external-resources/${experiment._id}/${sg}/gsea_go.png`
          },
          {
            type: "gsea_interpro",
            uri: `${atlasURL}/external-resources/${experiment._id}/${sg}/gsea_interpro.png`
          },
          {
            type: "ma-plot",
            uri: `${atlasURL}/external-resources/${experiment._id}/${sg}/ma-plot.png`
          }
        ]
      }
    }),
    profiles: {
      searchResultTotal: Object.keys(expressions).length,
      rows: Object.keys(expressions).map(g => {
        return {
          id: g,
          name: g, // possibly lookup the name
          uri: `${atlasURL}/genes/${g}`,
          expressionUnit: "Log2 fold change",
          expressions: sigGroups.map(sg => expressions[g][sg] || {})
        }
      })
    }
  };
  res.json(de);
}

async function baseline_experiment(req, res) {
  try {
    // get experiment_id param
    var experiment_id = req.swagger.params.experiment_id.raw;
    // get list of ids from request body
    const ids = req.body.replace("geneQuery=", "").split(" ");
    // uniqify ids
    let uniqueIdentifiers = [...new Set(ids)];
    // Fetch experiments from MongoDB
    let experiments = await fetchBy(mongo.experiments.mongoCollection(), [experiment_id], '_id');
    const experiment = experiments[0];
    // Fetch taxonomy from MongoDB
    let taxa = await fetchBy(mongo.taxonomy.mongoCollection(), [experiment.taxon_id], '_id');
    const species = taxa[0].name;
    // Fetch assays from MongoDB
    let assays = await fetchBy(mongo.assays.mongoCollection(), [experiment_id], 'experiment');
    const searchArray = uniqueIdentifiers.map(id => {value: id});
    const encodedSearch = encodeURIComponent(JSON.stringify(searchArray));
    // Fetch gene expression from MongoDB
    let genes = await fetchBy(mongo.expression.mongoCollection(), uniqueIdentifiers, '_id');


    if (experiment.type === 'Differential') {
      return differential_experiment(res,experiment,assays,genes,species)
    }


    // get baseline experiment assays with non-zero tpm
    var tpm = {};
    genes.forEach(g => {
      tpm[g._id] = {};
      g[experiment_id].forEach(a => {
        if (a.hasOwnProperty('value') && a.value >= 0) {
          tpm[g._id][a.group] = a.value < 0.5 ? 0 : a.value;
        }
      })
    })
    // get the organism_part ontology terms for the anatomogram
    let anatomogram_parts = new Set();
    let columns = [];
    // let species;
    assays.forEach(assay => {
      let isFactor = {};
      assay.factor.forEach(factor => {
        isFactor[factor.label] = true;
      })
      const col = {
        assayGroupId: assay.group,
        factorValue: Object.keys(isFactor).join('; '),
        assayGroupSummary: {replicates:2, properties: []}
      }
      assay.characteristic.forEach(ch => {
        if (!species && ch.type === "organism") {
          species = ch.label.split(" ").slice(0,2).join("_").toLowerCase();
        }
        col.assayGroupSummary.properties.push({
          propertyName: ch.type,
          testValue: ch.label,
          contrastPropertyType: isFactor[ch.label] ? "FACTOR" : "SAMPLE"
        });
        if (isFactor[ch.label] && ch.ontology) {
          const termId = ch.id.replace(":","_");
          anatomogram_parts.add(termId);
          col.factorValueOntologyTermId = termId;
        }
      })
      columns.push(col);
      
    });
    // sort the columns
    const columnHeaders = columns.sort((a,b) => a.factorValue.localeCompare(b.factorValue));
    console.log("uids",uniqueIdentifiers);
    let rows = uniqueIdentifiers.map(gene => {
      return {
        id: gene,
        name: gene,
        uri: `genes/${gene}`,
        expressionUnit: "TPM",
        expressions: columnHeaders.map(ch => {
          if (tpm[gene][ch.assayGroupId]) {
            const v = tpm[gene][ch.assayGroupId];
            return {value: v, quartiles: {min: v, lower: v, median: v, upper: v, max: v}}
          }
          return {}
        })
      }
    });
    let result = {
      anatomogram: {
        species: species,
        allSvgPathIds: [...anatomogram_parts]
      },
      columnHeaders: columnHeaders,
      columnGroupings:[],
      profiles: {
        rows: rows,
        searchResultTotal: `${rows.length}`
      },
      config: {
        columnType: "",
        conditionQuery: "",
        disclaimer: "",
        geneQuery: encodedSearch,
        genomeBrowsers: [],
        species: species.replace("_"," ")
      },
      experiment: {
        accession: experiment_id,
        description: experiment.description,
        species: species.replace("_"," "),
        type: "rnaseq_mrna_baseline",
        urls: {
          download: `${atlasURL}/experiments-content/${experiment_id}/download/RNASEQ_MRNA_BASELINE?unit=TPM&cutoff=0.5&heatmapMatrixSize=50&selectedColumnIds=&type=RNASEQ_MRNA_BASELINE&geneQuery=${encodedSearch}`,
          genome_browsers: `${atlasURL}/experiments/${experiment_id}/redirect/genome-browsers`,
          main_page: `${atlasURL}/experiments/${experiment_id}?geneQuery=${encodedSearch}`
        }
      }
    };
    if (anatomogram_parts.size === 0) {
      delete result.anatomogram;
    }
    res.json(result);
  } catch (error) {
    console.error(`Error in baseline_experiment`, error);
    res.status(500).json({ error: "Internal Server Error" });
  }
}
