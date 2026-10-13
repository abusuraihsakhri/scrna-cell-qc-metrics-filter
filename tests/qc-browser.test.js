"use strict";
const test = require("node:test");
const assert = require("node:assert/strict");
const { parseCSV, calculateQC, filteredMatrix, toCSV } = require("../web/app.js");

const sample = "gene,good,high_mito,empty\nMT-CO1,10,80,0\nRPL13,10,10,0\nACTB,80,10,0\nGAPDH,100,0,0\n";
const settings = { minCounts: 200, minFeatures: 4, maxMito: 5, maxRibo: null, mitoPrefix: "MT-" };

test("correct counts, features and gene-type percentages", () => {
  const rows = calculateQC(parseCSV(sample), settings);
  assert.deepEqual(rows, [
    {cell_id:"good",n_counts:200,n_features:4,percent_mt:5,percent_ribo:5,pass_qc:true},
    {cell_id:"high_mito",n_counts:100,n_features:3,percent_mt:80,percent_ribo:10,pass_qc:false},
    {cell_id:"empty",n_counts:0,n_features:0,percent_mt:0,percent_ribo:0,pass_qc:false},
  ]);
});
test("filtered export retains only passing cell columns", () => {
  const table = parseCSV(sample);
  assert.deepEqual(filteredMatrix(table, calculateQC(table, settings)), [
    ["gene", "good"], ["MT-CO1", "10"], ["RPL13", "10"], ["ACTB", "80"], ["GAPDH", "100"]
  ]);
});
test("quoted CSV fields and roundtrip", () => {
  const rows = [["gene","cell,1"],["RPL\"3","4"]];
  assert.deepEqual(parseCSV(toCSV(rows)), rows);
});
test("optional ribosomal threshold", () => {
  const rows = calculateQC(parseCSV(sample), {...settings, maxRibo: 4});
  assert.equal(rows[0].pass_qc, false);
});
test("malformed inputs rejected", () => {
  for (const data of ["", "gene\n", "gene,a,a\nX,2,3\n", "gene,a\nX,-1\n",
    "gene,a\nX,NaN\n", "gene,a\nX,1.1\n", "gene,a\nX,1,2\n",
    "gene,=cmd\nX,1\n", "gene,a\n"]) {
    assert.throws(() => calculateQC(parseCSV(data), settings));
  }
  assert.throws(() => parseCSV('gene,a\n"unclosed'));
});
test("invalid thresholds rejected", () => {
  assert.throws(() => calculateQC(parseCSV(sample), {...settings,maxMito:101}));
  assert.throws(() => calculateQC(parseCSV(sample), {...settings,minCounts:-1}));
  assert.throws(() => calculateQC(parseCSV(sample), {...settings,maxRibo:NaN}));
});
