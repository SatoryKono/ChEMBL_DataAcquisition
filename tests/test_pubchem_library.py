from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parents[1]))

from library import pubchem_library as pl


def test_url_encode() -> None:
    assert pl.url_encode("Caffeine 10%") == "Caffeine%2010%25"


def test_validate_cid() -> None:
    assert pl.validate_cid("") is None
    assert pl.validate_cid("0") is None
    assert pl.validate_cid("-1") is None
    assert pl.validate_cid("123") == "123"


def test_extract_cids() -> None:
    bindings = [
        {"cid": {"value": "123"}},
        {"cid": "http://rdf.ncbi.nlm.nih.gov/pubchem/compound/CID456"},
    ]
    assert pl._extract_cids(bindings) == ["123", "456"]


def test_get_cid(monkeypatch) -> None:
    sample = {"results": {"bindings": [{"cid": {"value": "123"}}, {"cid": {"value": "123"}}, {"cid": {"value": "456"}}]}}
    monkeypatch.setattr(pl, "make_request", lambda url, delay=3.0: sample)
    assert pl.get_cid("water") == "123|456"


def test_get_cid_from_smiles(monkeypatch) -> None:
    sample = {"IdentifierList": {"CID": [123, 123, 456]}}
    monkeypatch.setattr(pl, "make_request", lambda url, delay=3.0: sample)
    assert pl.get_cid_from_smiles("O") == "123|456"


def test_get_cid_from_inchi(monkeypatch) -> None:
    sample = {"IdentifierList": {"CID": [789]}}
    monkeypatch.setattr(pl, "make_request", lambda url, delay=3.0: sample)
    assert pl.get_cid_from_inchi("InChI=1S/H2O/h1H2") == "789"


def test_get_cid_from_inchikey(monkeypatch) -> None:
    sample = {"IdentifierList": {"CID": [321]}}
    monkeypatch.setattr(pl, "make_request", lambda url, delay=3.0: sample)
    assert pl.get_cid_from_inchikey("XLYOFNOQVPJJNP-UHFFFAOYSA-N") == "321"


def test_get_standard_name(monkeypatch) -> None:
    sample = {"InformationList": {"Information": [{"Title": "Water"}]}}
    monkeypatch.setattr(pl, "validate_cid", lambda cid: cid)
    monkeypatch.setattr(pl, "make_request", lambda url, delay=3.0: sample)
    assert pl.get_standard_name("123") == "Water"


def test_get_properties(monkeypatch) -> None:
    sample = {
        "PropertyTable": {
            "Properties": [
                {
                    "IUPACName": "Oxidane",
                    "MolecularFormula": "H2O",
                    "IsomericSMILES": "O",
                    "CanonicalSMILES": "O",
                    "InChI": "InChI=1S/H2O/h1H2",
                    "InChIKey": "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
                }
            ]
        }
    }
    monkeypatch.setattr(pl, "validate_cid", lambda cid: cid)
    monkeypatch.setattr(pl, "make_request", lambda url, delay=3.0: sample)
    props = pl.get_properties("123")
    assert props.MolecularFormula == "H2O"
    assert props.IUPACName == "Oxidane"


def test_process_compound(monkeypatch) -> None:
    monkeypatch.setattr(pl, "get_cid", lambda name: "123")
    monkeypatch.setattr(pl, "get_standard_name", lambda cid: "Water")
    monkeypatch.setattr(
        pl,
        "get_properties",
        lambda cid: pl.Properties("Oxidane", "H2O", "O", "O", "InChI", "Key"),
    )
    result = pl.process_compound("Water")
    assert result["CID"] == "123"
    assert result["Standard Name"] == "Water"
    assert result["MolecularFormula"] == "H2O"
