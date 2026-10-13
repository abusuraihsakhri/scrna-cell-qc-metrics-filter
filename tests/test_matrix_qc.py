"""Reference calculations and failure-mode tests for count matrix QC."""

import csv
import pytest

from qc_metrics import calculate_matrix_qc, write_qc_csv, write_filtered_matrix
from scrna_qc import main


@pytest.fixture
def matrix(tmp_path):
    path = tmp_path / "matrix.csv"
    path.write_text(
        "gene,good,high_mito,empty\n"
        "MT-CO1,10,80,0\n"
        "RPL13,10,10,0\n"
        "ACTB,80,10,0\n"
        "GAPDH,100,0,0\n",
        encoding="utf-8",
    )
    return path


def test_counts_features_percentages_and_inclusive_thresholds(matrix):
    rows = calculate_matrix_qc(matrix, min_counts=200, min_features=4, max_mito=5)
    assert rows[0] == {
        "cell_id": "good", "n_counts": 200, "n_features": 4,
        "percent_mt": 5.0, "percent_ribo": 5.0, "pass_qc": True,
    }
    assert rows[1]["n_counts"] == 100
    assert rows[1]["n_features"] == 3
    assert rows[1]["percent_mt"] == 80
    assert rows[1]["percent_ribo"] == 10
    assert rows[1]["pass_qc"] is False
    assert rows[2]["n_counts"] == 0
    assert rows[2]["pass_qc"] is False


def test_optional_ribosomal_filter(matrix):
    rows = calculate_matrix_qc(matrix, min_counts=1, min_features=1, max_mito=100, max_ribo=4)
    assert not rows[0]["pass_qc"]
    assert not rows[1]["pass_qc"]


def test_case_insensitive_mito_prefix(tmp_path):
    inp = tmp_path / "mouse.csv"
    inp.write_text("gene,c1\nmt-Nd1,9\nGAPDH,1\n", encoding="utf-8")
    result = calculate_matrix_qc(inp, min_counts=1, min_features=1, max_mito=90)
    assert result[0]["percent_mt"] == 90
    assert result[0]["pass_qc"] is True


@pytest.mark.parametrize("data", [
    "", "gene\n", "gene,a,a\nX,2,3\n", "gene,a\nX,-1\n",
    "gene,a\nX,1.2\n", "gene,a\nX,nan\n", "gene,a\nX,inf\n",
    "gene,a\nX,1,3\n", "gene,a\n", "gene,=HYPERLINK(1)\nX,2\n",
])
def test_invalid_matrices(tmp_path, data):
    inp = tmp_path / "bad.csv"
    inp.write_text(data, encoding="utf-8")
    with pytest.raises(ValueError):
        calculate_matrix_qc(inp)


def test_invalid_thresholds(matrix):
    for options in ({"min_counts": -1}, {"max_mito": 101}, {"max_ribo": float("nan")}):
        with pytest.raises(ValueError):
            calculate_matrix_qc(matrix, **options)


def test_csv_export_and_filtered_matrix(matrix, tmp_path):
    rows = calculate_matrix_qc(matrix, min_counts=100, min_features=4, max_mito=5)
    output = tmp_path / "qc.csv"
    filtered = tmp_path / "filtered.csv"
    write_qc_csv(rows, output)
    write_filtered_matrix(matrix, filtered, {r["cell_id"] for r in rows if r["pass_qc"]})
    with output.open(newline="") as handle:
        records = list(csv.DictReader(handle))
    assert len(records) == 3
    assert records[0]["pass_qc"] == "True"
    assert filtered.read_text().splitlines()[0] == "gene,good"
    assert "MT-CO1,10" in filtered.read_text()


def test_cli_matrix(matrix, tmp_path, capsys):
    output = tmp_path / "qc.csv"
    filtered = tmp_path / "filtered.csv"
    assert main([
        "matrix", "-i", str(matrix), "-o", str(output),
        "--filtered-matrix", str(filtered), "--min-counts", "200",
        "--min-features", "4", "--max-mito", "5",
    ]) is None
    assert output.exists() and filtered.exists()
    assert "1/3 cells retained" in capsys.readouterr().out


def test_cli_prevents_input_overwrite(matrix):
    with pytest.raises(SystemExit) as exc:
        main(["matrix", "-i", str(matrix), "-o", str(matrix)])
    assert exc.value.code == 1
