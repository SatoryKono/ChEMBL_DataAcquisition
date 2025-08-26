import json
import sys
from pathlib import Path

import responses

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from mylib import pubchem

DATA_DIR = Path(__file__).parent / "data"


@responses.activate
def test_get_cid(monkeypatch):
    monkeypatch.setattr(pubchem.time, "sleep", lambda *_: None)
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/rdf/query?graph=synonym&"
        "return=cid&format=json&name=aspirin"
    )
    with (DATA_DIR / "pubchem_getcid.json").open() as handle:
        responses.add(responses.GET, url, json=json.load(handle), status=200)
    assert pubchem.get_cid("aspirin") == "2244"


@responses.activate
def test_process_compound(monkeypatch):
    monkeypatch.setattr(pubchem.time, "sleep", lambda *_: None)
    base = "https://pubchem.ncbi.nlm.nih.gov/rest"

    with (DATA_DIR / "pubchem_getcid.json").open() as handle:
        responses.add(
            responses.GET,
            f"{base}/rdf/query?graph=synonym&return=cid&format=json&name=aspirin",
            json=json.load(handle),
            status=200,
        )
    with (DATA_DIR / "pubchem_standardname.json").open() as handle:
        responses.add(
            responses.GET,
            f"{base}/pug/compound/cid/2244/description/JSON",
            json=json.load(handle),
            status=200,
        )
    with (DATA_DIR / "pubchem_properties.json").open() as handle:
        responses.add(
            responses.GET,
            f"{base}/pug/compound/cid/2244/property/" "MolecularFormula,IUPACName,IsomericSMILES,CanonicalSMILES,InChI,InChIKey/JSON",
            json=json.load(handle),
            status=200,
        )
    result = pubchem.process_compound("aspirin")
    assert result["CID"] == "2244"
    assert result["Standard Name"] == "Aspirin"
    assert result["MolecularFormula"] == "C9H8O4"
