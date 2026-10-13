# scRNA-seq Cell QC Metrics & Filter

Calculate library size, detected features, mitochondrial and ribosomal count fractions for each cell in a dense gene-by-cell CSV. Apply adjustable QC cutoffs and export results or a count matrix containing only passing cells.

## Browser application

The standalone interface in [web/](web/) accepts a local CSV, reports per-cell QC metrics and filtering decisions, and exports two CSV files. Choose **Load example data** to test without uploading anything. Files are parsed locally by JavaScript; they are not sent to a server. The interface accepts CSV files up to 10 MiB; use the Python CLI for larger data. The browser does not require Python, Pyodide, external fonts or external packages.

After enabling **Settings → Pages → Build and deployment → Source: GitHub Actions** once, the deployment workflow publishes `web/` on pushes to `master`. The live URL must be verified after deployment.

## Input format

First column: gene identifiers. Remaining columns: cell barcodes. Values: **nonnegative integer** counts.

```csv
gene,good,high_mito,empty
MT-CO1,10,80,0
RPL13,10,10,0
ACTB,80,10,0
GAPDH,100,0,0
```

This is [example_matrix.csv](example_matrix.csv). Counts can be raw UMI or read counts, but must not be normalized or log-transformed. All data are interpreted as a **dense gene-by-cell matrix**, not 10x MTX, HDF5 or AnnData. Cell IDs must be unique. Unsafe spreadsheet-formula prefixes in gene/cell labels are rejected.

### Calculation and defaults

| Column | Definition |
| --- | --- |
| `n_counts` | Sum of counts across all gene rows for the cell |
| `n_features` | Number of gene rows with counts > 0 |
| `percent_mt` | 100 × counts for genes with prefix `MT-` / `n_counts` |
| `percent_ribo` | 100 × counts for genes with prefix `RPL` or `RPS` / `n_counts` |
| `pass_qc` | Nonzero total, counts ≥ minimum, features ≥ minimum, mitochondrial % ≤ maximum, and optional ribosomal % ≤ maximum |

Defaults: minimum 500 counts, minimum 200 detected features, maximum 20% mitochondrial counts, **no ribosomal exclusion**. The mitochondrial prefix is case-insensitive and configurable. Zero-count cells report 0% for both fractions but **always fail**. Cutoffs are inclusive. Prefix matching is a heuristic; inspect gene annotations, mitochondrial nomenclature and dataset distributions before using thresholds. Ribosomal pseudogene annotations may be included.

## Python CLI

Requires Python 3.10 or newer. The new matrix module uses only the standard library; the package also contains a separate legacy audit API requiring Pydantic.

```bash
python -m pip install -e ".[dev]"
python scrna_qc.py matrix -i example_matrix.csv -o cell_qc.csv \
  --filtered-matrix filtered_matrix.csv \
  --min-counts 200 --min-features 4 --max-mito 5
```

The example passes only cell `good`. Optional `--max-ribo 30` applies a ribosomal cutoff. For mouse data, use `--mito-prefix mt-`. All cell metrics are exported regardless of pass/fail status; the filtered matrix preserves gene rows and passing cell columns.

The earlier commands `python scrna_qc.py single` and `python scrna_qc.py batch` are retained for compatibility. **They calculate a generic weighted score and are not valid scRNA-seq QC metrics.** Likewise, `python cli.py audit|batch|chat|verify-audit|serve` exposes a separate legacy generic audit system, not single-cell QC or clinical decision support.

### API / Docker (legacy audit system)

```bash
python -m pip install -e ".[api]"
python cli.py serve --host 127.0.0.1 --port 8000
```

For Docker Compose, supply a strong `AUDIT_SECRET_KEY` from the deployment environment before starting. Never commit deployment secrets. Audit records are in-memory and are not persistent across restarts; the identifier screen is heuristic and **cannot guarantee PHI detection**. Do not submit identifiable patient data. GitHub Pages hosts only the static QC application, not FastAPI.

## Tests and development

```bash
python -m pip install -e ".[dev]"
pytest -q
node --check web/app.js
node --test tests/qc-browser.test.js
```

No frontend build is required. For local browsing run `python -m http.server 8000 --directory web` and open `http://localhost:8000/` in a current Chromium, Firefox or Safari browser; the UI uses standard File, Blob and DOM APIs. GitHub Actions tests the Python suite on 3.10–3.12 and the JavaScript QC logic with Node.

## Privacy and limitations

Browser input and downloaded results stay on the local device unless the user chooses to share them. JavaScript uses no external analytics, network requests or server calls. **Small dense CSV only** in the browser; huge single-cell matrices should be processed with established sparse-matrix workflows. This tool does not perform ambient RNA correction, doublet detection, automated outlier modeling, sample stratification or biological interpretation.

## License

MIT — see [LICENSE](LICENSE).
