"""Tests for document-related helper libraries."""

from __future__ import annotations

import requests

from library import semantic_scholar_library as ss
from library import openalex_crossref_library as ocl


def test_semantic_scholar_error(monkeypatch) -> None:
    def fake_request(*args, **kwargs):
        return None, "boom"

    monkeypatch.setattr(ss, "_do_request", fake_request)
    with requests.Session() as session:
        res = ss.fetch_semantic_scholar(session, "123", 0)
    assert res["scholar.Error"] == "boom"


def test_openalex_error(monkeypatch) -> None:
    def fake_request(*args, **kwargs):
        return None, "oops"

    monkeypatch.setattr(ocl, "_do_request", fake_request)
    with requests.Session() as session:
        res = ocl.fetch_openalex(session, "123", 0)
    assert res["OpenAlex.Error"] == "oops"


def test_crossref_missing_doi() -> None:
    with requests.Session() as session:
        res = ocl.fetch_crossref(session, "", 0)
    assert res["crossref.Error"] == "Missing DOI"
