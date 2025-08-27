import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1]))

import get_testitem_data as gtd
from library import chembl_library as cl


def _sample_testitem_df() -> pd.DataFrame:
    data = {col: [""] for col in cl.TESTITEM_COLUMNS}
    data["molecule_chembl_id"] = ["CHEMBL25"]
    return pd.DataFrame(data)


def test_read_ids(tmp_path: Path) -> None:
    csv_file = tmp_path / "molecules.csv"
    csv_file.write_text(
        "molecule_chembl_id\nCHEMBL25\n#N/A\n\nCHEMBL12\n", encoding="utf8"
    )
    assert gtd.read_ids(csv_file) == ["CHEMBL25", "CHEMBL12"]


def test_run_chembl(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(
        cl, "get_testitem", lambda ids, chunk_size=5: _sample_testitem_df()
    )
    input_csv = tmp_path / "molecules.csv"
    input_csv.write_text("molecule_chembl_id\nCHEMBL25\n", encoding="utf8")
    output_csv = tmp_path / "out.csv"
    args = argparse.Namespace(
        input_csv=input_csv,
        output_csv=output_csv,
        column="molecule_chembl_id",
        sep=",",
        encoding="utf8",
        chunk_size=5,
    )
    assert gtd.run_chembl(args) == 0
    df = pd.read_csv(output_csv, dtype=str)
    assert df.loc[0, "molecule_chembl_id"] == "CHEMBL25"
