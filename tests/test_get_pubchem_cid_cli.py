from pathlib import Path
import sys
import argparse

sys.path.append(str(Path(__file__).resolve().parents[1]))

import get_pubchem_cid as gpc
from library import pubchem_library as pl


def test_run_smiles(monkeypatch, capsys) -> None:
    monkeypatch.setattr(pl, "get_cid_from_smiles", lambda s: "123")
    args = argparse.Namespace(smiles="O", inchi=None, inchikey=None, log_level="INFO")
    assert gpc.run(args) == 0
    assert capsys.readouterr().out.strip() == "123"


