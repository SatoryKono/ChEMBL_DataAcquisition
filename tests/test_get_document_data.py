from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1]))

import get_document_data as gdd
from library import chembl_library as cl
from library import pubmed_library as pl
from library import semantic_scholar_library as ssl
from library import openalex_crossref_library as ocl


def _sample_doc_df() -> pd.DataFrame:
    data = {col: [""] for col in cl.DOCUMENT_COLUMNS}
    data["document_chembl_id"] = ["CHEMBL100"]
    data["pubmed_id"] = [12345678]  # int64 column
    return pd.DataFrame(data)


def _sample_pub_df() -> pd.DataFrame:
    return pd.DataFrame({"PubMed.PMID": ["12345678"], "PubMed.DOI": ["10.1000/test"]})


def test_run_all_merges_types(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(cl, "get_documents", lambda ids, chunk_size=5: _sample_doc_df())
    monkeypatch.setattr(
        gdd,
        "fetch_pubmed_records",
        lambda pmids, sleep, workers=1, batch_size=200: _sample_pub_df(),
    )

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
        workers=1,
        batch_size=200,
    )
    assert gdd.run_all(args) == 0
    df = pd.read_csv(output_csv, dtype=str)
    assert df.loc[0, "PubMed.PMID"] == "12345678"


def test_fetch_pubmed_records_batches(monkeypatch) -> None:
    calls: list[list[str]] = []

    def fake_fetch_batch(session, pmids, sleep):
        calls.append(pmids)
        return [{"PubMed.PMID": pid, "PubMed.DOI": ""} for pid in pmids]

    monkeypatch.setattr(pl, "fetch_pubmed_batch", fake_fetch_batch)
    monkeypatch.setattr(
        ssl, "fetch_semantic_scholar", lambda s, pmid, sleep: {"scholar.DOI": ""}
    )
    monkeypatch.setattr(ocl, "fetch_openalex", lambda s, pmid, sleep: {})
    monkeypatch.setattr(ocl, "fetch_crossref", lambda s, doi, sleep: {})

    df = gdd.fetch_pubmed_records(["1", "2", "3"], sleep=0.0, max_workers=1, batch_size=2)
    assert calls == [["1", "2"], ["3"]]
    assert set(df["PubMed.PMID"]) == {"1", "2", "3"}
