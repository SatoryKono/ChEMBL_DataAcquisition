from pathlib import Path
import sys

import requests
import pytest

sys.path.append(str(Path(__file__).resolve().parents[1]))

from library import semantic_scholar_library as ssl


def test_fetch_semantic_scholar(monkeypatch: pytest.MonkeyPatch) -> None:
    session = requests.Session()

    def fake_fetch(session_arg: requests.Session, pmid: str, sleep: float) -> dict[str, str]:
        assert session_arg is session
        assert pmid == "42"
        assert sleep == 0.5
        return {"ok": "1"}

    monkeypatch.setattr(ssl._pl, "fetch_semantic_scholar", fake_fetch)
    result = ssl.fetch_semantic_scholar(session, "42", 0.5)
    assert result == {"ok": "1"}


def test_fetch_semantic_scholar_batch(monkeypatch: pytest.MonkeyPatch) -> None:
    session = requests.Session()
    pmids = ["1", "2", "3"]

    def fake_fetch_batch(
        session_arg: requests.Session, pmid_list: list[str], sleep: float
    ) -> list[dict[str, str]]:
        assert session_arg is session
        assert pmid_list == pmids
        assert sleep == 0.5
        return [{"pmid": "1"}, {"pmid": "2"}, {"pmid": "3"}]

    monkeypatch.setattr(ssl._pl, "fetch_semantic_scholar_batch", fake_fetch_batch)
    result = ssl.fetch_semantic_scholar_batch(session, pmids, 0.5)
    assert result == [{"pmid": "1"}, {"pmid": "2"}, {"pmid": "3"}]
