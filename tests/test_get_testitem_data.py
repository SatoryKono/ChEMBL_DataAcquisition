from pathlib import Path
import sys

import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1]))

import get_testitem_data as gtd
from library import pubchem_library as pl


def test_add_pubchem_data(monkeypatch) -> None:
    df = pd.DataFrame({"molecule_structures.canonical_smiles": ["O", ""]})

    monkeypatch.setattr(pl, "get_cid_from_smiles", lambda s: "123" if s else None)
    monkeypatch.setattr(
        pl,
        "get_properties",
        lambda cid: pl.Properties(
            "Oxidane", "H2O", "O", "O", "InChI=1S/H2O/h1H2", "XLYOFNOQVPJJNP-UHFFFAOYSA-N"
        ),
    )

    result = gtd.add_pubchem_data(df)
    assert list(result["pubchem_cid"]) == ["123", ""]
    assert list(result["pubchem_iupac_name"]) == ["Oxidane", ""]

