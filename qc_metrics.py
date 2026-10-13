"""Per-cell scRNA-seq QC metrics for a small dense gene-by-cell count matrix.

Input CSV: first column = gene identifier, subsequent columns = cell barcodes.
Values must be nonnegative integer counts (UMI counts or read counts).
No dependency on the unrelated legacy score or audit workflows.
"""

import csv
import math
from pathlib import Path

MAX_EXACT_COUNT = 2**53 - 1
RESULT_FIELDS = ["cell_id", "n_counts", "n_features", "percent_mt", "percent_ribo", "pass_qc"]


def _label(value, kind):
    value = value.strip()
    if not value or len(value) > 200 or any(ord(ch) < 32 for ch in value):
        raise ValueError(f"Invalid {kind}: empty, excessively long or control character")
    if value[0] in "=+-@":
        raise ValueError(f"Unsafe {kind} prefix: {value[0]!r}")
    return value


def _validate_threshold(value, name, percentage=False):
    if not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
        raise ValueError(f"{name} must be a finite nonnegative number")
    if percentage and value > 100:
        raise ValueError(f"{name} must not exceed 100")


def calculate_matrix_qc(
    input_csv,
    *,
    min_counts=500,
    min_features=200,
    max_mito=20.0,
    max_ribo=None,
    mito_prefix="MT-",
    ribo_prefixes=("RPL", "RPS"),
):
    """Calculate QC metrics and threshold-based flags for each cell.

    Thresholds are inclusive. Zero-count cells always fail. Ribosomal
    thresholds are disabled by default; ribosomal percentages remain reported.
    Percentages use total counts as denominator, not detected feature counts.
    """
    for value, name, pct in (
        (min_counts, "min_counts", False),
        (min_features, "min_features", False),
        (max_mito, "max_mito", True),
    ):
        _validate_threshold(value, name, pct)
    if max_ribo is not None:
        _validate_threshold(max_ribo, "max_ribo", True)
    if not isinstance(mito_prefix, str) or not mito_prefix:
        raise ValueError("mito_prefix must be a nonempty string")
    if not ribo_prefixes or any(not isinstance(p, str) or not p for p in ribo_prefixes):
        raise ValueError("ribo_prefixes must contain nonempty strings")

    with open(input_csv, encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, strict=True)
        header = next(reader, None)
        if not header or len(header) < 2:
            raise ValueError("Expected gene column followed by one or more cell columns")
        cells = [_label(c, "cell identifier") for c in header[1:]]
        if len(set(cells)) != len(cells):
            raise ValueError("Duplicate cell identifiers")
        totals = [0] * len(cells)
        features = [0] * len(cells)
        mito = [0] * len(cells)
        ribo = [0] * len(cells)
        genes = 0
        mp = mito_prefix.upper()
        rp = tuple(p.upper() for p in ribo_prefixes)

        for row_num, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(f"Row {row_num}: expected {len(header)} columns, got {len(row)}")
            gene = _label(row[0], "gene identifier")
            is_mito = gene.upper().startswith(mp)
            is_ribo = gene.upper().startswith(rp)
            genes += 1
            for i, raw in enumerate(row[1:]):
                try:
                    value = int(raw.strip())
                except ValueError as exc:
                    raise ValueError(f"Row {row_num}, cell {cells[i]}: expected integer count") from exc
                if value < 0 or value > MAX_EXACT_COUNT or totals[i] + value > MAX_EXACT_COUNT:
                    raise ValueError(f"Row {row_num}, cell {cells[i]}: count out of range")
                totals[i] += value
                features[i] += int(value > 0)
                if is_mito:
                    mito[i] += value
                if is_ribo:
                    ribo[i] += value

    if genes == 0:
        raise ValueError("Count matrix has no gene rows")

    results = []
    for i, cell in enumerate(cells):
        total = totals[i]
        pct_mt = (mito[i] * 100 / total) if total else 0.0
        pct_ribo = (ribo[i] * 100 / total) if total else 0.0
        passed = (
            total > 0
            and total >= min_counts
            and features[i] >= min_features
            and pct_mt <= max_mito
            and (max_ribo is None or pct_ribo <= max_ribo)
        )
        results.append({
            "cell_id": cell,
            "n_counts": total,
            "n_features": features[i],
            "percent_mt": round(pct_mt, 4),
            "percent_ribo": round(pct_ribo, 4),
            "pass_qc": passed,
        })
    return results


def write_qc_csv(rows, output_csv):
    """Write a per-cell QC table; reject unsafe filenames that overwrite input upstream."""
    path = Path(output_csv)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def write_filtered_matrix(input_csv, output_csv, passing_cells):
    """Export original genes and counts for passing cells, preserving source order."""
    source, dest = Path(input_csv), Path(output_csv)
    if source.resolve() == dest.resolve():
        raise ValueError("Output must not overwrite input")
    dest.parent.mkdir(parents=True, exist_ok=True)
    with source.open("r", encoding="utf-8-sig", newline="") as inp, dest.open(
        "w", encoding="utf-8", newline=""
    ) as out:
        reader = csv.reader(inp)
        writer = csv.writer(out)
        header = next(reader)
        indices = [i for i, cell in enumerate(header[1:], start=1) if cell.strip() in passing_cells]
        writer.writerow([header[0]] + [header[i] for i in indices])
        for row in reader:
            writer.writerow([row[0]] + [row[i] for i in indices])
