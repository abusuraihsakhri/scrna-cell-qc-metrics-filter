# Scrna Cell Qc Metrics Filter

> **Domain:** Computational Biology & AI Drug Discovery
> **Reference Guidelines & Standards:** `wwPDB, IUPAC & CLSI Computational Guidelines`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

Single-Cell RNA-Seq (scRNA-Seq) Quality Filter
Calculates nCounts, nFeatures, and percent mitochondrial / ribosomal reads per cell with outlier filtering.

Zero-dependency Python implementation with single and batch evaluation.
Author: Dr. Abu Suraih Sakhri
License: MIT

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Analytical Functions

- **`calculate_metrics()`**: Core domain algorithm that computes a weighted score from numeric inputs and classifies results into Low/Moderate/High tiers.
- **`process_single()`** — Evaluates a single case with provided parameters.
- **`process_batch()`** — Processes a CSV file of records, appending QC metrics to each row.
- **`main()`** — CLI entry point supporting `single` and `batch` subcommands.

---

## 📐 Mathematical Formulation & Logic

```text
score = primary_val + sum(nv * (1/idx) for idx, nv in enumerate(remaining_vals, start=2))
rounded_score = round(score, 2)

Classification tiers:
  score < 10.0  → "Low / Standard"
  score < 25.0  → "Moderate / Intermediate"
  score >= 25.0 → "High / Severe"
```

---

## 💻 Installation

```bash
# Clone the repository
git clone https://github.com/abusuraihsakhri/scrna-cell-qc-metrics-filter.git
cd scrna-cell-qc-metrics-filter

# Install core dependencies (pydantic required)
pip install pydantic>=2.0

# Optional: Install API server dependencies
pip install fastapi>=0.110 "uvicorn[standard]>=0.27"

# Optional: Install development dependencies
pip install pytest>=8.0
```

---

## 💻 CLI Quickstart & Usage

### 1. Single Case Evaluation
```bash
python scrna_qc.py single --v1 10.0 --v2 5.0 --v3 2.0
```

### 2. Batch CSV Processing
```bash
python scrna_qc.py batch -i sample.csv -o results.csv
```

### 3. Enterprise CLI (requires pydantic)
```bash
# Audit a single task
python cli.py audit --task-id TASK-001 --primary 28.5 --secondary 14.2

# Batch process records
python cli.py batch -i sample.csv -o results.csv

# Verify audit trail integrity
python cli.py verify-audit

# Launch FastAPI REST server
python cli.py serve --host 127.0.0.1 --port 8000
```

### Parameter Reference
- `--v1`, `--v2`, `--v3`: Numeric parameters for single evaluation (default: 10.0, 5.0, 2.0)
- `-i`, `--input`: Input CSV file path (required for batch)
- `-o`, `--output`: Output CSV file path (default: results.csv)

### Input Data Schema

| Field | Description | Requirement |
|:------|:------------|:------------|
| `Patient_ID` | Patient identifier | Required |
| `v1` | Primary measurement | Required |
| `v2` | Secondary measurement | Required |
| `v3` | Tertiary measurement | Optional |

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Configurable Audit Secret:** Set `AUDIT_SECRET_KEY` environment variable for persistent audit integrity across restarts. Without it, a random key is generated (with a warning).
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Run specific test files:

```bash
pytest tests/ -v          # Enterprise module tests
pytest test_scrna_qc.py -v  # Core QC tests
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py 100
```

---

## 🐳 Container Deployment

```bash
docker build -t scrna-cell-qc-metrics-filter .
docker run -p 8000:8000 scrna-cell-qc-metrics-filter
```

---

## 📁 Project Structure

```
scrna-cell-qc-metrics-filter/
├── scrna_qc.py          # Core QC algorithm & CLI
├── cli.py               # Enterprise CLI with audit/PHI features
├── enrichment.py        # Enrichment analysis modules
├── simulator.py         # High-throughput simulation
├── agents/              # Enterprise agent framework
│   ├── base.py          # Security, PHI guard, audit trail
│   ├── models.py        # Pydantic data models
│   ├── supervisor.py    # Multi-agent orchestrator
│   ├── workers.py       # Specialized domain workers
│   ├── api.py           # FastAPI REST server
│   └── ...
├── tests/               # Test suite
├── sample.csv           # Sample input data
├── pyproject.toml       # Project configuration
└── Dockerfile           # Container build
```
