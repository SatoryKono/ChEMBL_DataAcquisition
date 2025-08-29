"""Unit tests for the target post-processing helpers."""

from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mylib.transforms import postprocess_targets


def load_raw() -> pd.DataFrame:
    path = Path(__file__).parent / "data" / "targets_raw.csv"
    return pd.read_csv(path, dtype=str)


def test_gene_name_and_synonyms() -> None:
    raw = load_raw()
    processed = postprocess_targets(raw)
    row1 = processed[processed["chembl_id"] == "CHEMBL1"].iloc[0]
    row2 = processed[processed["chembl_id"] == "CHEMBL2"].iloc[0]

    assert row1["gene_name"] == "GNA1"
    assert row2["gene_name"] == "N1L"

    assert "desc a" in row1["synonyms"].split("|")
    assert "desc b" in row2["synonyms"].split("|")

    assert set(row1["ec_number"].split("|")) == {"1.1.1.1", "2.2.2.2"}
    assert row2["ec_number"] == "3.3.3.3"


def test_identifier_columns() -> None:
    raw = load_raw()
    processed = postprocess_targets(raw)
    row = processed.iloc[0]
    assert row["uniprotkb_Id"] == "P12345"
    assert row["secondary_uniprot_id"] == "P12345"
