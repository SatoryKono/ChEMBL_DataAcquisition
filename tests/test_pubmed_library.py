"""Unit tests for :mod:`library.pubmed_library`."""

from __future__ import annotations

import csv
from pathlib import Path

from library import pubmed_library as pl


def test_combine() -> None:
    assert pl.combine(["A", "", "B"]) == "A|B"


def test_read_pmids(tmp_path: Path) -> None:
    path = tmp_path / "pmids.csv"
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["PMID"])
        writer.writeheader()
        writer.writerows(
            [
                {"PMID": "123"},
                {"PMID": ""},
                {"PMID": "456"},
            ]
        )
    assert pl.read_pmids(path) == ["123", "456"]
