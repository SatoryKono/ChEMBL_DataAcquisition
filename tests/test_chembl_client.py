import sys
import json
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
import requests

# Ensure project root is on the path for imports
sys.path.append(str(Path(__file__).resolve().parents[1]))

from mylib import chembl_client as cc

DATA_DIR = Path(__file__).parent / "data"


class MockResponse:
    def __init__(self, payload: Any, status: int = 200) -> None:
        self.payload = payload
        self.status = status

    def json(self) -> Any:  # pragma: no cover - simple delegate
        return self.payload

    def raise_for_status(self) -> None:
        if self.status >= 400:
            raise requests.HTTPError(f"status: {self.status}")


@pytest.fixture
def target_payload() -> dict[str, Any]:
    with (DATA_DIR / "target.json").open("r", encoding="utf8") as fh:
        return json.load(fh)


@pytest.fixture
def assay_payload() -> dict[str, Any]:
    with (DATA_DIR / "assays.json").open("r", encoding="utf8") as fh:
        return json.load(fh)


def test_parse_target_record(target_payload: dict[str, Any]) -> None:
    record = cc._parse_target_record(target_payload)
    assert record["pref_name"] == "Sample Target"
    assert record["gene"] == "GENE1|GENE_ONE"
    assert record["ec_code"] == "1.2.3.4"
    assert record["chembl_alternative_name"] == "ALT1"
    assert record["HGNC_id"] == "1"


def test_get_target(monkeypatch: pytest.MonkeyPatch, target_payload: dict[str, Any]) -> None:
    def fake_get(url: str, timeout: int = 30):
        return MockResponse(target_payload)

    monkeypatch.setattr(cc.requests, "get", fake_get)
    result = cc.get_target("CHEMBL123")
    assert result["pref_name"] == "Sample Target"


def test_get_targets(monkeypatch: pytest.MonkeyPatch, target_payload: dict[str, Any]) -> None:
    def fake_get(url: str, timeout: int = 30):
        return MockResponse({"targets": [target_payload]})

    monkeypatch.setattr(cc.requests, "get", fake_get)
    df = cc.get_targets(["CHEMBL123", "", "#N/A"])
    assert list(df.columns) == cc.TARGET_FIELDS
    assert df.iloc[0]["pref_name"] == "Sample Target"


def test_get_assay(monkeypatch: pytest.MonkeyPatch, assay_payload: dict[str, Any]) -> None:
    assay = assay_payload["assays"][0]

    def fake_get(url: str, timeout: int = 30):
        return MockResponse(assay)

    monkeypatch.setattr(cc.requests, "get", fake_get)
    df = cc.get_assay("CHEMBL1")
    assert df.iloc[0]["assay_chembl_id"] == "CHEMBL1"
    assert list(df.columns) == cc.ASSAY_COLUMNS


def test_get_assays(monkeypatch: pytest.MonkeyPatch, assay_payload: dict[str, Any]) -> None:
    def fake_get(url: str, timeout: int = 30):
        return MockResponse(assay_payload)

    monkeypatch.setattr(cc.requests, "get", fake_get)
    df = cc.get_assays(["CHEMBL1"])
    assert df.iloc[0]["assay_chembl_id"] == "CHEMBL1"


def test_extend_target(monkeypatch: pytest.MonkeyPatch, target_payload: dict[str, Any]) -> None:
    def fake_get_target(tid: str) -> dict[str, Any]:
        return cc._parse_target_record(target_payload)

    monkeypatch.setattr(cc, "get_target", fake_get_target)
    df = pd.DataFrame({"task_chembl_id": ["CHEMBL123"]})
    extended = cc.extend_target(df)
    assert "chembl_pref_name" in extended.columns
    assert extended.loc[0, "chembl_pref_name"] == "Sample Target"

