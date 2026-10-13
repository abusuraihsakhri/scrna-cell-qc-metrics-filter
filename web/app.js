"use strict";

/* Dense gene-by-cell matrix QC. Counts and calculations remain in the browser. */
const CSV_LIMIT_BYTES = 10 * 1024 * 1024;
const MAX_EXACT_COUNT = Number.MAX_SAFE_INTEGER;
const DEMO = "gene,good,high_mito,empty\nMT-CO1,10,80,0\nRPL13,10,10,0\nACTB,80,10,0\nGAPDH,100,0,0\n";

function parseCSV(source) {
  const text = source.replace(/^\uFEFF/, "");
  const rows = [];
  let row = [], field = "", quoted = false, closed = false;
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (quoted) {
      if (ch === '"' && text[i + 1] === '"') { field += '"'; i++; }
      else if (ch === '"') { quoted = false; closed = true; }
      else { field += ch; }
    } else if (ch === '"') {
      if (field !== "" || closed) throw new Error("Invalid CSV quoting");
      quoted = true;
    } else if (ch === "," || ch === "\n" || ch === "\r") {
      row.push(field);
      field = "";
      closed = false;
      if (ch !== ",") {
        if (row.some(x => x !== "")) rows.push(row);
        row = [];
        if (ch === "\r" && text[i + 1] === "\n") i++;
      }
    } else {
      if (closed) throw new Error("Invalid characters after closing quote");
      field += ch;
    }
  }
  if (quoted) throw new Error("Unterminated quoted CSV field");
  if (field !== "" || row.length > 0) {
    row.push(field);
    if (row.some(x => x !== "")) rows.push(row);
  }
  return rows;
}

function validLabel(label, kind) {
  const value = label.trim();
  if (!value || value.length > 200 || /[\x00-\x1f]/.test(value))
    throw new Error("Invalid " + kind + " label");
  if (/^[=+\-@]/.test(value))
    throw new Error("Unsafe " + kind + " label prefix");
  return value;
}

function calculateQC(table, settings) {
  if (!Array.isArray(table) || table.length < 2 || table[0].length < 2)
    throw new Error("Expected a CSV with a gene column, cell headers, and at least one gene");
  const { minCounts, minFeatures, maxMito, maxRibo, mitoPrefix } = settings;
  for (const [name, value, upper] of [
    ["Minimum counts", minCounts, Infinity],
    ["Minimum features", minFeatures, Infinity],
    ["Maximum mitochondrial %", maxMito, 100],
  ]) {
    if (!Number.isFinite(value) || value < 0 || value > upper)
      throw new Error(name + " must be a valid nonnegative number");
  }
  if (maxRibo !== null && (!Number.isFinite(maxRibo) || maxRibo < 0 || maxRibo > 100))
    throw new Error("Maximum ribosomal % must be between 0 and 100");
  if (typeof mitoPrefix !== "string" || !mitoPrefix.trim())
    throw new Error("Mitochondrial prefix must not be empty");

  const cells = table[0].slice(1).map(s => validLabel(s, "cell"));
  if (new Set(cells).size !== cells.length) throw new Error("Duplicate cell identifiers");
  const totals = cells.map(() => 0), features = cells.map(() => 0);
  const mito = cells.map(() => 0), ribo = cells.map(() => 0);
  const mp = mitoPrefix.trim().toUpperCase();
  for (let r = 1; r < table.length; r++) {
    const row = table[r];
    if (row.length !== cells.length + 1)
      throw new Error("Row " + (r + 1) + ": wrong number of columns");
    const gene = validLabel(row[0], "gene").toUpperCase();
    const isMito = gene.startsWith(mp);
    const isRibo = gene.startsWith("RPL") || gene.startsWith("RPS");
    for (let c = 0; c < cells.length; c++) {
      const raw = row[c + 1].trim();
      if (!/^\d+$/.test(raw)) throw new Error("Row " + (r + 1) + ", " + cells[c] + ": expected a nonnegative integer");
      const value = Number(raw);
      if (!Number.isSafeInteger(value) || !Number.isSafeInteger(totals[c] + value))
        throw new Error("Count exceeds safe integer range at row " + (r + 1));
      totals[c] += value;
      if (value > 0) features[c]++;
      if (isMito) mito[c] += value;
      if (isRibo) ribo[c] += value;
    }
  }
  const records = cells.map((cell, c) => {
    const total = totals[c];
    const pctMt = total ? mito[c] * 100 / total : 0;
    const pctRibo = total ? ribo[c] * 100 / total : 0;
    return {
      cell_id: cell, n_counts: total, n_features: features[c],
      percent_mt: Number(pctMt.toFixed(4)), percent_ribo: Number(pctRibo.toFixed(4)),
      pass_qc: total > 0 && total >= minCounts && features[c] >= minFeatures
        && pctMt <= maxMito && (maxRibo === null || pctRibo <= maxRibo),
    };
  });
  return records;
}

function escapeCSV(value) {
  const str = String(value);
  return /[",\r\n]/.test(str) ? '"' + str.replace(/"/g, '""') + '"' : str;
}
function toCSV(rows) {
  return rows.map(row => row.map(escapeCSV).join(",")).join("\r\n") + "\r\n";
}
function filteredMatrix(table, records) {
  const passed = new Set(records.filter(r => r.pass_qc).map(r => r.cell_id));
  const indices = table[0].map((c, i) => passed.has(c.trim()) ? i : -1).filter(i => i > 0);
  return table.map(row => [row[0], ...indices.map(i => row[i])]);
}

function init() {
  const $ = id => document.getElementById(id);
  let table = null, results = null;
  const status = $("status");
  function message(text, isError = false) {
    status.textContent = text;
    status.classList.toggle("error", isError);
  }
  function reset() {
    results = null;
    $("summary").hidden = true;
    $("results-section").hidden = true;
    $("export").disabled = true;
    $("export-matrix").disabled = true;
  }
  $("file").addEventListener("change", async event => {
    reset();
    const file = event.target.files[0];
    if (!file) { table = null; return; }
    if (file.size > CSV_LIMIT_BYTES) {
      table = null;
      message("CSV exceeds the 10 MiB browser limit; use the Python CLI for larger data.", true);
      return;
    }
    try {
      table = parseCSV(await file.text());
      message("Loaded " + file.name + ". Select thresholds and run QC.");
    } catch (err) {
      table = null;
      message(err.message, true);
    }
  });
  $("demo").addEventListener("click", () => {
    reset();
    $("file").value = "";
    table = parseCSV(DEMO);
    $("min-counts").value = 200;
    $("min-features").value = 4;
    $("max-mito").value = 5;
    $("max-ribo").value = "";
    $("mito-prefix").value = "MT-";
    message("Demo matrix loaded: 4 genes and 3 cells.");
  });
  $("form").addEventListener("submit", event => {
    event.preventDefault();
    reset();
    if (!table) { message("Choose a count matrix CSV or load the demo.", true); return; }
    try {
      const ribo = $("max-ribo").value.trim();
      results = calculateQC(table, {
        minCounts: Number($("min-counts").value),
        minFeatures: Number($("min-features").value),
        maxMito: Number($("max-mito").value),
        maxRibo: ribo === "" ? null : Number(ribo),
        mitoPrefix: $("mito-prefix").value,
      });
      const passed = results.filter(r => r.pass_qc).length;
      $("total").textContent = String(results.length);
      $("passed").textContent = String(passed);
      $("failed").textContent = String(results.length - passed);
      $("rate").textContent = (100 * passed / results.length).toFixed(1) + "%";
      $("summary").hidden = false;
      const tbody = $("result-body");
      tbody.replaceChildren();
      for (const r of results.slice(0, 150)) {
        const tr = document.createElement("tr");
        for (const value of [r.cell_id, r.n_counts, r.n_features,
          r.percent_mt.toFixed(2) + "%", r.percent_ribo.toFixed(2) + "%",
          r.pass_qc ? "PASS" : "FAIL"]) {
          const td = document.createElement("td");
          td.textContent = String(value);
          tr.appendChild(td);
        }
        tr.className = r.pass_qc ? "pass" : "fail";
        tbody.appendChild(tr);
      }
      $("preview-note").textContent = "Showing " + Math.min(results.length, 150) +
        " of " + results.length + " cells. Export includes all cells.";
      $("results-section").hidden = false;
      $("export").disabled = false;
      $("export-matrix").disabled = false;
      message("QC complete. Percentages are based on total counts per cell.");
    } catch (err) { message(err.message, true); }
  });
  function download(filename, csvText) {
    const blob = new Blob([csvText], { type: "text/csv;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = filename;
    document.body.appendChild(link);
    link.click();
    link.remove();
    URL.revokeObjectURL(url);
  }
  $("export").addEventListener("click", () => {
    if (!results) return;
    const fields = ["cell_id", "n_counts", "n_features", "percent_mt", "percent_ribo", "pass_qc"];
    download("cell_qc.csv", toCSV([fields, ...results.map(r => fields.map(f => r[f]))]));
  });
  $("export-matrix").addEventListener("click", () => {
    if (results && table) download("filtered_matrix.csv", toCSV(filteredMatrix(table, results)));
  });
}
if (typeof document !== "undefined") document.addEventListener("DOMContentLoaded", init);
if (typeof module !== "undefined" && module.exports) {
  module.exports = { parseCSV, calculateQC, toCSV, filteredMatrix };
}
