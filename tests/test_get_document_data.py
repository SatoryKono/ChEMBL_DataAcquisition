from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1]))

import get_document_data as gdd
from library import chembl_library as cl


def _sample_doc_df() -> pd.DataFrame:
    data = {col: [""] for col in cl.DOCUMENT_COLUMNS}
    data["document_chembl_id"] = ["CHEMBL100"]
    data["pubmed_id"] = [12345678]  # int64 column
    return pd.DataFrame(data)


def _sample_pub_df() -> pd.DataFrame:
    return pd.DataFrame({"PubMed.PMID": ["12345678"], "PubMed.DOI": ["10.1000/test"]})


def test_run_all_merges_types(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(cl, "get_documents", lambda ids, chunk_size=5: _sample_doc_df())
    monkeypatch.setattr(gdd, "fetch_pubmed_records", lambda pmids, sleep: _sample_pub_df())

    input_csv = tmp_path / "docs.csv"
    input_csv.write_text("chembl_id\nCHEMBL100\n", encoding="utf8")
    output_csv = tmp_path / "out.csv"

    args = argparse.Namespace(
        input_csv=input_csv,
        output_csv=output_csv,
        column="chembl_id",
        sep=",",
        encoding="utf8",
        chunk_size=5,
        sleep=0.0,
    )
    assert gdd.run_all(args) == 0
    df = pd.read_csv(output_csv, dtype=str)
    assert df.loc[0, "PubMed.PMID"] == "12345678"
