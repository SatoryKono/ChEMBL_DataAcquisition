import csv
import json
from pathlib import Path
from typing import Dict

import pytest

from uniprot import uniprot_utils as uu

SAMPLE_ENTRY: Dict[str, object] = {
    "proteinDescription": {
        "recommendedName": {
            "fullName": {"value": "Protein Kinase"},
            "shortName": {"value": "PK"},
        },
        "alternativeNames": [
            {
                "fullName": {"value": "Alt Protein"},
                "shortName": [{"value": "AP"}],
            }
        ],
    },
    "genes": [
        {"geneName": {"value": "GENE1"}, "synonyms": [{"value": "G1"}]}
    ],
    "organism": {
        "lineage": ["Eukaryota", "Metazoa", "Chordata", "Mammalia", "Homo"],
        "taxonId": 9606,
    },
}


def test_fetch_uniprot(monkeypatch: pytest.MonkeyPatch) -> None:
    expected = {"id": "P12345"}

    class MockResp:
        def raise_for_status(self) -> None:  # pragma: no cover - trivial
            pass

        def json(self) -> Dict[str, str]:
            return expected

    def mock_get(url: str) -> MockResp:
        assert url == uu.API_URL.format(id="P12345")
        return MockResp()

    monkeypatch.setattr(uu.requests, "get", mock_get)
    result = uu.fetch_uniprot("P12345")
    assert result == expected


def test_extract_names() -> None:
    names = uu.extract_names(SAMPLE_ENTRY)
    assert names == {
        "Protein Kinase",
        "PK",
        "Alt Protein",
        "AP",
        "GENE1",
        "G1",
    }


def test_extract_organism() -> None:
    org = uu.extract_organism(SAMPLE_ENTRY)
    assert org == {
        "genus": "Homo",
        "superkingdom": "Eukaryota",
        "phylum": "Chordata",
        "taxon_id": "9606",
    }


def test_iter_ids_and_validation(tmp_path: Path) -> None:
    good_csv = tmp_path / "ids.csv"
    good_csv.write_text("UniProt_id\nP12345\nQ8N158 \n", encoding="utf-8")
    ids = list(uu.iter_ids(str(good_csv)))
    assert ids == ["P12345", "Q8N158"]

    bad_csv = tmp_path / "bad.csv"
    bad_csv.write_text("id\nP12345\n", encoding="utf-8")
    with pytest.raises(ValueError):
        list(uu.iter_ids(str(bad_csv)))


def test_collect_info(tmp_path: Path) -> None:
    data_dir = tmp_path / "Data"
    data_dir.mkdir()
    (data_dir / "P12345.json").write_text(json.dumps(SAMPLE_ENTRY), encoding="utf-8")

    info = uu.collect_info("P12345", str(data_dir))
    assert info == {
        "UniProt_id": "P12345",
        "names": "AP|Alt Protein|G1|GENE1|PK|Protein Kinase",
        "genus": "Homo",
        "superkingdom": "Eukaryota",
        "phylum": "Chordata",
        "taxon_id": "9606",
    }


def test_process(tmp_path: Path) -> None:
    data_dir = tmp_path / "Data"
    data_dir.mkdir()
    (data_dir / "P12345.json").write_text(json.dumps(SAMPLE_ENTRY), encoding="utf-8")

    input_csv = tmp_path / "input.csv"
    input_csv.write_text("UniProt_id\nP12345\nP99999\n", encoding="utf-8")
    output_csv = tmp_path / "out.csv"

    uu.process(str(input_csv), str(output_csv), data_dir=str(data_dir))

    with output_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    assert rows[0]["UniProt_id"] == "P12345"
    assert rows[0]["names"] == "AP|Alt Protein|G1|GENE1|PK|Protein Kinase"
    assert rows[1]["UniProt_id"] == "P99999"
    assert rows[1]["names"] == ""
