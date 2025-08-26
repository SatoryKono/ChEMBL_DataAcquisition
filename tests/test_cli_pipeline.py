from __future__ import annotations

import argparse
from pathlib import Path
import sys
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1]))

import get_target_data as gtd
from library import chembl_library as cl

DATA_DIR = Path(__file__).resolve().parents[1] / "data"


def _sample_chembl_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "pref_name": ["MAP3K14"],
            "target_chembl_id": ["CHEMBL1075024"],
            "component_description": [
                "Mitogen-activated protein kinase kinase kinase 14"
            ],
            "component_id": ["123"],
            "relationship": ["single"],
            "gene": ["MAP3K14"],
            "uniprot_id": ["Q99558"],
            "chembl_alternative_name": ["NIK"],
            "ec_code": ["2.7.11.25"],
            "hgnc_name": ["MAP3K14"],
            "hgnc_id": ["6853"],
        }
    )


def test_run_chembl(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(cl, "get_targets", lambda ids: _sample_chembl_df())
    input_csv = tmp_path / "targets.csv"
    input_csv.write_text("chembl_id\nCHEMBL1075024\n", encoding="utf8")
    output_csv = tmp_path / "chembl.csv"
    args = argparse.Namespace(
        input_csv=input_csv,
        output_csv=output_csv,
        column="chembl_id",
        sep=",",
        encoding="utf8",
    )
    assert gtd.run_chembl(args) == 0
    df = pd.read_csv(output_csv, dtype=str)
    assert df.loc[0, "uniprot_id"] == "Q99558"


def test_run_uniprot(tmp_path: Path) -> None:
    input_csv = tmp_path / "uids.csv"
    input_csv.write_text("uniprot_id\nQ99558\n", encoding="utf8")
    output_csv = tmp_path / "uniprot.csv"
    args = argparse.Namespace(
        input_csv=input_csv,
        output_csv=output_csv,
        data_dir=DATA_DIR / "uniprot",
        sep=",",
        encoding="utf8",
    )
    assert gtd.run_uniprot(args) == 0
    df = pd.read_csv(output_csv, dtype=str)
    assert df.loc[0, "uniprot_id"] == "Q99558"


def test_run_iuphar(tmp_path: Path) -> None:
    input_csv = tmp_path / "uids.csv"
    input_csv.write_text("uniprot_id\nQ99558\n", encoding="utf8")
    output_csv = tmp_path / "iuphar.csv"
    args = argparse.Namespace(
        input_csv=input_csv,
        output_csv=output_csv,
        target_csv=DATA_DIR / "_IUPHAR_target.csv",
        family_csv=DATA_DIR / "_IUPHAR_family.csv",
        sep=",",
        encoding="utf8",
    )
    assert gtd.run_iuphar(args) == 0
    df = pd.read_csv(output_csv, dtype=str)
    assert df.loc[0, "target_id"] == "2074"


def test_run_all(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(cl, "get_targets", lambda ids: _sample_chembl_df())
    input_csv = tmp_path / "chembl_ids.csv"
    input_csv.write_text("chembl_id\nCHEMBL1075024\n", encoding="utf8")
    output_csv = tmp_path / "merged.csv"
    args = argparse.Namespace(
        input_csv=input_csv,
        output_csv=output_csv,
        chembl_out=None,
        uniprot_out=None,
        iuphar_out=None,
        data_dir=DATA_DIR / "uniprot",
        target_csv=DATA_DIR / "_IUPHAR_target.csv",
        family_csv=DATA_DIR / "_IUPHAR_family.csv",
        sep=",",
        encoding="utf8",
    )
    assert gtd.run_all(args) == 0
    df = pd.read_csv(output_csv, dtype=str)
    assert df.loc[0, "chembl_id"] == "CHEMBL1075024"
    assert df.loc[0, "uniprot_id"] == "Q99558"
    assert df.loc[0, "target_id"] == "2074"
    assert df.loc[0, "IUPHAR_family_id"] == "0624"
