import argparse
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.append(str(Path(__file__).resolve().parents[1]))

import get_activity_data as gad
from library import chembl_library as cl


def _sample_activity_df() -> pd.DataFrame:
    data = {col: [""] for col in cl.ACTIVITY_COLUMNS}
    data["activity_id"] = ["1"]
    return pd.DataFrame(data)


def test_read_ids(tmp_path: Path) -> None:
    csv_file = tmp_path / "activities.csv"
    csv_file.write_text(
        "activity_id\n1\n#N/A\n\n2\n", encoding="utf8"
    )
    assert gad.read_ids(csv_file) == ["1", "2"]


def test_run_chembl(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(
        cl, "get_activities", lambda ids, chunk_size=5: _sample_activity_df()
    )
    input_csv = tmp_path / "activities.csv"
    input_csv.write_text("activity_id\n1\n", encoding="utf8")
    output_csv = tmp_path / "out.csv"
    args = argparse.Namespace(
        input_csv=input_csv,
        output_csv=output_csv,
        column="activity_id",
        sep=",",
        encoding="utf8",
        chunk_size=5,
    )
    assert gad.run_chembl(args) == 0
    df = pd.read_csv(output_csv, dtype=str)
    assert df.loc[0, "activity_id"] == "1"

