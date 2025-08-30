from pathlib import Path
import sys
import pandas as pd
import pytest

sys.path.append(str(Path(__file__).resolve().parents[1]))

from library import chembl_library as cl
from get_target_data import read_ids

# Obtain a real ChEMBL ID from the provided targets.csv file
DATA_DIR = Path(__file__).resolve().parents[1]
SAMPLE_ID = read_ids(DATA_DIR / "targets.csv")[0]


class FakeResponse:
    def __init__(self, data):
        self._data = data
        self.status_code = 200

    def raise_for_status(self) -> None:
        return None

    def json(self):
        return self._data


SAMPLE_JSON = {
    "pref_name": "MAP3K14",
    "target_chembl_id": SAMPLE_ID,
    "target_components": {
        "target_component": {
            "component_description": "Mitogen-activated protein kinase kinase kinase 14",
            "component_id": 123,
            "relationship": "single",
            "target_component_synonyms": {
                "target_component_synonym": [
                    {"component_synonym": "MAP3K14", "syn_type": "GENE_SYMBOL"},
                    {"component_synonym": "NIK", "syn_type": "GENE_SYMBOL_OTHER"},
                    {"component_synonym": "2.7.11.25", "syn_type": "EC_NUMBER"},
                    {
                        "component_synonym": "Mitogen-activated protein kinase kinase kinase 14",
                        "syn_type": "UNIPROT",
                    },
                ]
            },
            "target_component_xrefs": {
                "target": [
                    {"xref_src_db": "UniProt", "xref_id": "Q99558"},
                    {
                        "xref_src_db": "HGNC",
                        "xref_id": "HGNC:6853",
                        "xref_name": "MAP3K14",
                    },
                ]
            },
        }
    },
}


SAMPLE_DF = pd.DataFrame(
    {
        "pref_name": ["MAP3K14"],
        "target_chembl_id": [SAMPLE_ID],
        "component_description": [
            "Mitogen-activated protein kinase kinase kinase 14",
        ],
        "component_id": ["123"],
        "relationship": ["single"],
        "gene": ["MAP3K14"],
        "uniprot_id": ["Q99558"],
        "mapping_uniprot_id": ["P12345"],
        "chembl_alternative_name": ["NIK"],
        "ec_code": ["2.7.11.25"],
        "hgnc_name": ["MAP3K14"],
        "hgnc_id": ["6853"],
    }
)


@pytest.fixture(autouse=True)
def _mock_mapping(monkeypatch) -> None:
    """Always return a fixed UniProt ID for mapping in tests."""

    monkeypatch.setattr(cl, "map_chembl_to_uniprot", lambda cid: "P12345")


def test_chunked_splits_list() -> None:
    assert list(cl._chunked([1, 2, 3, 4], 2)) == [[1, 2], [3, 4]]


def test_chunked_invalid_size() -> None:
    with pytest.raises(ValueError):
        list(cl._chunked([1], 0))


def test_get_target(monkeypatch) -> None:
    monkeypatch.setattr(cl._session, "get", lambda url, timeout=30: FakeResponse(SAMPLE_JSON))
    data = cl.get_target(SAMPLE_ID)
    assert data["uniprot_id"] == "Q99558"
    assert data["gene"] == "MAP3K14|NIK"
    assert data["mapping_uniprot_id"] == "P12345"


def test_get_targets(monkeypatch) -> None:
    bulk_json = {"targets": [SAMPLE_JSON]}
    monkeypatch.setattr(cl._session, "get", lambda url, timeout=30: FakeResponse(bulk_json))
    df = cl.get_targets([SAMPLE_ID])
    assert df.loc[0, "uniprot_id"] == "Q99558"
    assert df.loc[0, "mapping_uniprot_id"] == "P12345"
    assert df.shape[0] == 1


def test_extend_target(monkeypatch) -> None:
    monkeypatch.setattr(cl, "get_targets", lambda ids, chunk_size=50: SAMPLE_DF)
    input_df = pd.DataFrame({"task_chembl_id": [SAMPLE_ID]})
    out_df = cl.extend_target(input_df)
    assert out_df.loc[0, "chembl_pref_name"] == "MAP3K14"
    assert out_df.loc[0, "chembl_gene"] == "MAP3K14"
