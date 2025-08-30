import argparse
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.append(str(Path(__file__).resolve().parents[1]))

import get_assay_data as gad
from library import chembl_library as cl


def _sample_assay_df() -> pd.DataFrame:
    data = {col: [""] for col in cl.ASSAY_COLUMNS}
    data["assay_chembl_id"] = ["CHEMBL123"]
    return pd.DataFrame(data)


def test_read_ids(tmp_path: Path) -> None:
    csv_file = tmp_path / "assays.csv"
    csv_file.write_text(
        "assay_chembl_id\nCHEMBL123\n#N/A\n\nCHEMBL456\n", encoding="utf8"
    )
    assert gad.read_ids(csv_file) == ["CHEMBL123", "CHEMBL456"]


def test_run_chembl(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(
        cl, "get_assays_all", lambda ids, chunk_size=5: _sample_assay_df()
    )
    input_csv = tmp_path / "assays.csv"
    input_csv.write_text("assay_chembl_id\nCHEMBL123\n", encoding="utf8")
    output_csv = tmp_path / "out.csv"
    args = argparse.Namespace(
        input_csv=input_csv,
        output_csv=output_csv,
        column="assay_chembl_id",
        sep=",",
        encoding="utf8",
        chunk_size=5,
    )
    assert gad.run_chembl(args) == 0
    df = pd.read_csv(output_csv, dtype=str)
    assert df.loc[0, "assay_chembl_id"] == "CHEMBL123"
