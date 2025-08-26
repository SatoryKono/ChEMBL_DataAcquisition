"""Unit tests for IUPHAR transformations."""

from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1]))

from IUPHAR.script.transforms import IUPHARData, IUPHARClassifier


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
    assert list(out["uniprot_id"]) == [
        "Q11111",
        "Q33333",
        "Q99999",
        "Q88888",
        "Q77777",
        "Q66666",
        "UNKNOWN",
    ]
    assert list(out["target_id"]) == [
        "T1",
        "T3",
        "T2",
        "T1",
        "T3",
        "T3",
        "",
    ]
    assert list(out["full_id_path"]) == [
        "T1#F1",
        "T3#F2>F1",
        "T2#F1",
        "T1#F1",
        "T3#F2>F1",
        "T3#F2>F1",
        "",
    ]
    assert list(out["full_name_path"]) == [
        "Target One#Family One",
        "Target Three#Family Two>Family One",
        "Target Two#Family One",
        "Target One#Family One",
        "Target Three#Family Two>Family One",
        "Target Three#Family Two>Family One",
        "",
    ]
    assert "IUPHAR_class" in out.columns
    assert list(out["IUPHAR_chain"]) == [
        "F1",
        "F2>F1",
        "F1",
        "F1",
        "F2>F1",
        "F2>F1",
        "0690-1>0690",
    ]
    # ec_number fallback should classify unresolved targets
    ec_row = out[out["uniprot_id"] == "UNKNOWN"].iloc[0]
    assert ec_row["IUPHAR_class"] == "Enzyme"
    assert ec_row["IUPHAR_subclass"] == "Oxidoreductase"


def test_classify_by_target() -> None:
    data = load_data()
    classifier = IUPHARClassifier(data)
    rec = classifier.by_target_id("T1")
    assert rec.IUPHAR_target_id == "T1"
    assert rec.IUPHAR_family_id == "F1"
    assert rec.IUPHAR_tree == ["F1"]


def test_classify_by_uniprot() -> None:
    data = load_data()
    classifier = IUPHARClassifier(data)
    rec = classifier.by_uniprot_id("Q11111")
    assert rec.IUPHAR_target_id == "T1"
    assert rec.IUPHAR_family_id == "F1"


def test_classify_by_name_heuristic() -> None:
    data = load_data()
    classifier = IUPHARClassifier(data)
    rec = classifier.by_name("Example Kinase")
    assert rec.IUPHAR_type == "Enzyme.Transferase"
    assert rec.IUPHAR_class == "Enzyme"
    assert rec.IUPHAR_subclass == "Transferase"


def test_ec_number_mapping() -> None:
    data = load_data()
    classifier = IUPHARClassifier(data)
    rec = classifier.by_ec_number("1.2.3.4")
    assert rec.IUPHAR_type == "Enzyme.Oxidoreductase"
    assert rec.IUPHAR_tree == ["0690-1", "0690"]


def test_target_id_by_name_special_chars() -> None:
    """Special characters should be treated literally, not as regex."""
    data = load_data()
    # "[" would previously raise a regex error
    assert data.target_id_by_name("alp[ha") == ""
