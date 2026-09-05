import pytest
from scrna_qc import calculate_metrics, process_batch, main


def test_scrna_qc_single():
    res = calculate_metrics(v1=12.0, v2=4.0)
    assert "score" in res
    assert "classification" in res
    assert res["score"] > 0


def test_scrna_qc_batch(tmp_path):
    csv_in = tmp_path / "in.csv"
    csv_out = tmp_path / "out.csv"
    csv_in.write_text("Patient,v1,v2\nPat_001,15.0,3.0\nPat_002,5.0,1.0\n", encoding="utf-8")

    process_batch(str(csv_in), str(csv_out))
    assert csv_out.exists()
    content = csv_out.read_text(encoding="utf-8")
    assert "Pat_001" in content
    assert "score" in content


def test_calculate_metrics_classification_low():
    res = calculate_metrics(v1=5.0, v2=2.0)
    assert res["classification"] == "Low / Standard"
    assert res["score"] < 10.0


def test_calculate_metrics_classification_moderate():
    res = calculate_metrics(v1=15.0)
    assert res["classification"] == "Moderate / Intermediate"
    assert 10.0 <= res["score"] < 25.0


def test_calculate_metrics_classification_high():
    res = calculate_metrics(v1=30.0)
    assert res["classification"] == "High / Severe"
    assert res["score"] >= 25.0


def test_calculate_metrics_empty_input():
    res = calculate_metrics()
    assert res["score"] == 1.0
    assert res["inputs_evaluated"] == 0


def test_calculate_metrics_string_values():
    res = calculate_metrics(v1="abc", v2="def")
    assert res["score"] == 1.0
    assert res["inputs_evaluated"] == 2


def test_calculate_metrics_mixed_types():
    res = calculate_metrics(v1=10.0, v2="abc")
    assert res["score"] == 10.0
    assert res["inputs_evaluated"] == 2


def test_process_batch_file_not_found(tmp_path):
    csv_out = tmp_path / "out.csv"
    with pytest.raises(FileNotFoundError):
        process_batch(str(tmp_path / "nonexistent.csv"), str(csv_out))


def test_process_batch_empty_headers(tmp_path):
    csv_in = tmp_path / "in.csv"
    csv_out = tmp_path / "out.csv"
    csv_in.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="no headers"):
        process_batch(str(csv_in), str(csv_out))


def test_main_single():
    result = main(["single", "--v1", "10.0", "--v2", "5.0"])
    assert result is None  # main() returns None for single command


def test_main_batch_file_not_found():
    with pytest.raises(SystemExit) as exc_info:
        main(["batch", "-i", "nonexistent.csv"])
    assert exc_info.value.code == 1


def test_scrna_qc_output_fields(tmp_path):
    csv_in = tmp_path / "in.csv"
    csv_out = tmp_path / "out.csv"
    csv_in.write_text("Patient,v1,v2\nPat_001,15.0,3.0\n", encoding="utf-8")

    process_batch(str(csv_in), str(csv_out))
    content = csv_out.read_text(encoding="utf-8")
    assert "clinical_recommendation" in content
    assert "classification" in content
