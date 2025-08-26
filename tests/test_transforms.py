"""Unit tests for IUPHAR transformations."""

from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1]))

from mylib.transforms import IUPHARData


def load_data() -> IUPHARData:
    base = Path(__file__).parent / "data"
    target = base / "target.csv"
    family = base / "family.csv"
    return IUPHARData.from_files(target, family)


def test_family_chain() -> None:
    data = load_data()
    assert data.family_chain("F2") == ["F2", "F1"]


def test_target_id_by_uniprot() -> None:
    data = load_data()
    assert data.target_id_by_uniprot("Q22222") == "T2"


def test_full_paths() -> None:
    """Ensure full ID and name paths are composed correctly."""
    data = load_data()
    assert data.all_id("T1") == "T1#F1"
    assert data.all_name("T1") == "Target One#Family One"


def test_uniprot_file_processing(tmp_path: Path) -> None:
    data = load_data()
    input_file = Path(__file__).parent / "data" / "uniprot_input.csv"
    output_file = tmp_path / "mapped.csv"
    data.map_uniprot_file(input_file, output_file)
    assert output_file.exists()
    out = pd.read_csv(output_file, dtype=str).fillna("")
    assert list(out["uniprot_id"]) == ["Q11111", "Q33333", "Q99999"]
    assert list(out["target_id"]) == ["T1", "T3", ""]
    assert list(out["full_id_path"]) == ["T1#F1", "T3#F2>F1", ""]
    assert list(out["full_name_path"]) == [
        "Target One#Family One",
        "Target Three#Family Two>Family One",
        "",
    ]
