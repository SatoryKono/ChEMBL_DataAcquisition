from pathlib import Path
import sys

import requests
import pytest

sys.path.append(str(Path(__file__).resolve().parents[1]))

from library import openalex_crossref_library as ocl


def test_fetch_openalex(monkeypatch: pytest.MonkeyPatch) -> None:
    session = requests.Session()

    def fake_fetch(session_arg: requests.Session, pmid: str, sleep: float) -> dict[str, str]:
        assert session_arg is session
        assert pmid == "1"
        assert sleep == 0.1
        return {"result": "openalex"}

    monkeypatch.setattr(ocl._pl, "fetch_openalex", fake_fetch)
    result = ocl.fetch_openalex(session, "1", 0.1)
    assert result == {"result": "openalex"}


def test_fetch_crossref(monkeypatch: pytest.MonkeyPatch) -> None:
    session = requests.Session()

    def fake_fetch(session_arg: requests.Session, doi: str, sleep: float) -> dict[str, str]:
        assert session_arg is session
        assert doi == "10.1/abc"
        assert sleep == 0.2
        return {"result": "crossref"}

    monkeypatch.setattr(ocl._pl, "fetch_crossref", fake_fetch)
    result = ocl.fetch_crossref(session, "10.1/abc", 0.2)
    assert result == {"result": "crossref"}
