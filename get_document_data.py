"""Command line interface for retrieving document metadata from external sources.

The tool integrates :mod:`library.pubmed_library` and
:mod:`library.chembl_library` to collect information about publications from
several public APIs.  The interface mirrors :mod:`get_target_data.py` and
provides three sub-commands:

``pubmed``
    Query PubMed, Semantic Scholar, OpenAlex and CrossRef for a list of PMIDs.
``chembl``
    Retrieve document information from the ChEMBL API.
``all``
    Run the ChEMBL and PubMed pipelines and merge the results.

Example
-------
Fetch PubMed metadata for identifiers listed in ``pmids.csv``::

    python get_document_data.py pubmed pmids.csv output.csv

The input file must contain a ``PMID`` column.
"""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path
from typing import Sequence

import pandas as pd
import requests

from library import chembl_library as cl
from library import pubmed_library as pl

logger = logging.getLogger(__name__)


def read_ids(
    path: str | Path,
    column: str = "document_chembl_id",
    sep: str = ",",
    encoding: str = "utf8",
) -> list[str]:
    """Read identifier values from a CSV file.

    Parameters
    ----------
    path:
        CSV file path.
    column:
        Name of the column containing identifiers.
    sep:
        Field delimiter used in the CSV file.
    encoding:
        File encoding of the CSV file.

    Returns
    -------
    list[str]
        Identifier values in the order they appear. Empty strings and ``"#N/A"``
        entries are discarded.
    """

    try:
        with Path(path).open("r", encoding=encoding, newline="") as fh:
            reader = csv.DictReader(fh, delimiter=sep)
            if reader.fieldnames is None or column not in reader.fieldnames:
                raise ValueError(f"column '{column}' not found in {path}")
            ids: list[str] = []
            for row in reader:
                value = (row.get(column) or "").strip()
                if value and value != "#N/A":
                    ids.append(value)
            return ids
    except FileNotFoundError as exc:  # pragma: no cover - trivial
        raise FileNotFoundError(f"input file not found: {path}") from exc
    except csv.Error as exc:  # pragma: no cover - malformed CSV
        raise ValueError(f"malformed CSV in file: {path}: {exc}") from exc


def fetch_pubmed_records(pmids: list[str], sleep: float) -> pd.DataFrame:
    """Retrieve metadata for a list of PubMed identifiers."""

    records: list[dict[str, str]] = []
    with requests.Session() as session:
        for pmid in pmids:
            pubmed = pl.fetch_pubmed(session, pmid, sleep)
            semsch = pl.fetch_semantic_scholar(session, pmid, sleep)
            openalex = pl.fetch_openalex(session, pmid, sleep)
            doi = pubmed.get("PubMed.DOI") or semsch.get("scholar.DOI") or ""
            crossref = pl.fetch_crossref(session, doi, sleep)
            combined: dict[str, str] = {}
            combined.update(pubmed)
            combined.update(semsch)
            combined.update(openalex)
            combined.update(crossref)
            records.append(combined)

    if not records:
        return pd.DataFrame()
    return pd.DataFrame(records)


def run_pubmed(args: argparse.Namespace) -> None:
    """Entry point for the ``pubmed`` sub-command.

    Parameters
    ----------
    args:
        Parsed command-line arguments.
    """

    pmids = pl.read_pmids(args.input_csv)
    df = fetch_pubmed_records(pmids, args.sleep)
    df.to_csv(args.output_csv, index=False)


def run_chembl(args: argparse.Namespace) -> None:
    """Entry point for the ``chembl`` sub-command.

    Parameters
    ----------
    args:
        Parsed command-line arguments.
    """

    ids = read_ids(
        args.input_csv, column=args.column, sep=args.sep, encoding=args.encoding
    )
    df = cl.get_documents(ids, chunk_size=args.chunk_size)
    df.to_csv(args.output_csv, index=False)


def run_all(args: argparse.Namespace) -> None:
    """Fetch data from ChEMBL and PubMed and merge the results.

    Parameters
    ----------
    args:
        Parsed command-line arguments.
    """

    ids = read_ids(
        args.input_csv, column=args.column, sep=args.sep, encoding=args.encoding
    )
    doc_df = cl.get_documents(ids, chunk_size=args.chunk_size)
    if doc_df.empty or "pubmed_id" not in doc_df:
        doc_df.to_csv(args.output_csv, index=False)
        return

    # Ensure PubMed IDs are numeric and collect non-null IDs for querying
    doc_df["pubmed_id"] = pd.to_numeric(doc_df["pubmed_id"], errors="coerce").astype(
        "Int64"
    )
    pmids = doc_df["pubmed_id"].dropna().astype(str).tolist()

    pub_df = fetch_pubmed_records(pmids, args.sleep)
    if pub_df.empty or "PubMed.PMID" not in pub_df.columns:
        doc_df.to_csv(args.output_csv, index=False)
        return

    pub_df = pub_df.add_prefix("pubmed.")
    pub_df["pubmed.PubMed.PMID"] = pd.to_numeric(
        pub_df["pubmed.PubMed.PMID"], errors="coerce"
    ).astype("Int64")

    merged = doc_df.merge(
        pub_df, how="left", left_on="pubmed_id", right_on="pubmed.PubMed.PMID"
    )
    merged.to_csv(args.output_csv, index=False)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Document data utilities")
    parser.add_argument("--log-level", default="INFO", help="Logging level")
    sub = parser.add_subparsers(dest="command", required=True)

    pubmed = sub.add_parser("pubmed", help="Fetch data from PubMed and related APIs")
    pubmed.add_argument("input_csv", type=Path, help="CSV with a PMID column")
    pubmed.add_argument("output_csv", type=Path, help="Destination CSV file")
    pubmed.add_argument(
        "--sleep", type=float, default=5.0, help="Seconds to sleep between requests"
    )
    pubmed.set_defaults(func=run_pubmed)

    chembl = sub.add_parser("chembl", help="Fetch document information from ChEMBL")
    chembl.add_argument("input_csv", type=Path, help="CSV with document_chembl_id column")
    chembl.add_argument("output_csv", type=Path, help="Destination CSV file")
    chembl.add_argument(
        "--column", default="document_chembl_id", help="Column name containing identifiers"
    )
    chembl.add_argument("--sep", default=",", help="CSV delimiter")
    chembl.add_argument("--encoding", default="utf8", help="File encoding")
    chembl.add_argument(
        "--chunk-size", type=int, default=5, help="Maximum number of IDs per request"
    )
    chembl.set_defaults(func=run_chembl)

    all_cmd = sub.add_parser("all", help="Run both ChEMBL and PubMed pipelines")
    all_cmd.add_argument("input_csv", type=Path, help="CSV with document_chembl_id column")
    all_cmd.add_argument("output_csv", type=Path, help="Destination CSV file")
    all_cmd.add_argument(
        "--column", default="document_chembl_id", help="Column in the input CSV"
    )
    all_cmd.add_argument("--sep", default=",", help="CSV delimiter")
    all_cmd.add_argument("--encoding", default="utf8", help="File encoding")
    all_cmd.add_argument(
        "--chunk-size", type=int, default=5, help="Maximum IDs per request"
    )
    all_cmd.add_argument(
        "--sleep", type=float, default=5.0, help="Seconds to sleep between PubMed requests"
    )
    all_cmd.set_defaults(func=run_all)

    return parser


def main(argv: Sequence[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO))
    args.func(args)


if __name__ == "__main__":
    main()

