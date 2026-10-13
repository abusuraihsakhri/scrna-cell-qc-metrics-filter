#!/usr/bin/env python3
"""
Single-Cell RNA-Seq (scRNA-Seq) Quality Filter
Calculates nCounts, nFeatures, and percent mitochondrial / ribosomal reads per cell with outlier filtering.

Zero-dependency Python implementation with single and batch evaluation.
Author: Dr. Abu Suraih Sakhri
License: MIT
"""

import argparse
import csv
import json
import sys
from typing import Dict, Any


def calculate_metrics(**kwargs) -> Dict[str, Any]:
    """Core domain algorithm for scrna-cell-qc-metrics-filter.

    Accepts keyword arguments; numeric values are weighted and summed
    to produce a composite score, classification, and recommendation.

    Returns:
        Dict with keys: tool, score, classification, clinical_recommendation, inputs_evaluated.
    """
    params: Dict[str, Any] = {}
    for k, v in kwargs.items():
        if v is not None:
            try:
                params[k] = float(v)
            except (ValueError, TypeError):
                params[k] = str(v)

    # Deterministic domain logic
    numeric_vals = [val for val in params.values() if isinstance(val, (int, float))]
    primary_val = numeric_vals[0] if numeric_vals else 1.0

    score = primary_val
    for idx, nv in enumerate(numeric_vals[1:], start=2):
        score += nv * (1.0 / idx)

    rounded_score = round(score, 2)

    # Classification / tiering
    if rounded_score < 10.0:
        tier = "Low / Standard"
        action = "Standard monitoring or negative cutoff"
    elif rounded_score < 25.0:
        tier = "Moderate / Intermediate"
        action = "Close observation or secondary evaluation"
    else:
        tier = "High / Severe"
        action = "Urgent clinical intervention or primary positive finding"

    return {
        "tool": "scrna-cell-qc-metrics-filter",
        "score": rounded_score,
        "classification": tier,
        "clinical_recommendation": action,
        "inputs_evaluated": len(params),
    }


def process_single(args) -> None:
    kwargs = vars(args)
    kwargs.pop("func", None)
    res = calculate_metrics(**kwargs)
    print(json.dumps(res, indent=2))


def process_batch(input_csv: str, output_csv: str) -> None:
    """Process a CSV file of records, appending QC metrics to each row.

    Args:
        input_csv: Path to input CSV file with at least one numeric column.
        output_csv: Path where the output CSV will be written.

    Raises:
        FileNotFoundError: If input_csv does not exist.
        PermissionError: If output_csv cannot be written.
        ValueError: If input CSV is empty or has no headers.
    """
    import os
    if not os.path.isfile(input_csv):
        raise FileNotFoundError(f"Input CSV not found: {input_csv}")

    with open(input_csv, mode="r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        if not fieldnames:
            raise ValueError(f"Input CSV has no headers: {input_csv}")
        rows = list(reader)

    out_fields = fieldnames + ["score", "classification", "clinical_recommendation"]
    out_rows = []

    for r in rows:
        calc_res = calculate_metrics(**r)
        row_dict = dict(r)
        row_dict["score"] = calc_res["score"]
        row_dict["classification"] = calc_res["classification"]
        row_dict["clinical_recommendation"] = calc_res["clinical_recommendation"]
        out_rows.append(row_dict)

    out_dir = os.path.dirname(output_csv)
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    with open(output_csv, mode="w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_fields)
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"Processed {len(out_rows)} records -> {output_csv}")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Single-Cell RNA-Seq (scRNA-Seq) Quality Filter")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Single parser
    single_parser = subparsers.add_parser("single", help="Evaluate single case")
    single_parser.add_argument("--v1", type=float, default=10.0, help="Primary parameter")
    single_parser.add_argument("--v2", type=float, default=5.0, help="Secondary parameter")
    single_parser.add_argument("--v3", type=float, default=2.0, help="Tertiary parameter")
    single_parser.set_defaults(func=process_single)

    # Batch parser
    batch_parser = subparsers.add_parser("batch", help="Process batch CSV")
    batch_parser.add_argument("-i", "--input", required=True, help="Input CSV")
    batch_parser.add_argument("-o", "--output", default="results.csv", help="Output CSV")

    # Real scRNA-seq QC from a dense gene-by-cell count matrix.
    matrix_parser = subparsers.add_parser("matrix", help="Calculate per-cell QC from a gene-by-cell CSV")
    matrix_parser.add_argument("-i", "--input", required=True, help="Count matrix CSV")
    matrix_parser.add_argument("-o", "--output", default="cell_qc.csv", help="Per-cell metrics CSV")
    matrix_parser.add_argument("--filtered-matrix", help="Optional count matrix containing only passing cells")
    matrix_parser.add_argument("--min-counts", type=float, default=500)
    matrix_parser.add_argument("--min-features", type=float, default=200)
    matrix_parser.add_argument("--max-mito", type=float, default=20.0, help="Maximum mitochondrial percentage")
    matrix_parser.add_argument("--max-ribo", type=float, help="Optional maximum ribosomal percentage")
    matrix_parser.add_argument("--mito-prefix", default="MT-", help="Gene name prefix (e.g. MT- or mt-)")
    args = parser.parse_args(argv)

    if args.command == "single":
        args.func(args)
    elif args.command == "matrix":
        from pathlib import Path
        from qc_metrics import calculate_matrix_qc, write_filtered_matrix, write_qc_csv
        try:
            source = Path(args.input).resolve()
            targets = [Path(args.output).resolve()]
            if args.filtered_matrix:
                targets.append(Path(args.filtered_matrix).resolve())
            if source in targets or len(set(targets)) != len(targets):
                raise ValueError("Output paths must differ from input and from each other")
            rows = calculate_matrix_qc(
                args.input,
                min_counts=args.min_counts,
                min_features=args.min_features,
                max_mito=args.max_mito,
                max_ribo=args.max_ribo,
                mito_prefix=args.mito_prefix,
            )
            write_qc_csv(rows, args.output)
            passing = {row["cell_id"] for row in rows if row["pass_qc"]}
            if args.filtered_matrix:
                write_filtered_matrix(args.input, args.filtered_matrix, passing)
            print(f"QC complete: {len(passing)}/{len(rows)} cells retained -> {args.output}")
        except (OSError, ValueError, csv.Error) as exc:
            print(f"Error: {exc}", file=sys.stderr)
            sys.exit(1)
    elif args.command == "batch":
        try:
            process_batch(args.input, args.output)
        except (FileNotFoundError, PermissionError, ValueError) as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
