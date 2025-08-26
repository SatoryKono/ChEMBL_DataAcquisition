import sys
from pathlib import Path

import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1] / "script"))
import chembl_lib  # type: ignore


def test_extend_target_merges(monkeypatch):
    df = pd.DataFrame({"task_chembl_id": ["A", "B"]})
    targets = pd.DataFrame({
        "target_chembl_id": ["A", "B"],
        "pref_name": ["nameA", "nameB"],
        "component_description": ["descA", "descB"],
        "component_id": [1, 2],
        "relationship": ["relA", "relB"],
        "gene": ["geneA", "geneB"],
        "chembl_alternative_name": ["altA", "altB"],
        "ec_code": ["ecA", "ecB"],
        "HGNC_name": ["hA", "hB"],
        "HGNC_id": ["1", "2"],
    })

    def fake_get_targets(ids):
        return targets[targets["target_chembl_id"].isin(ids)]

    monkeypatch.setattr(chembl_lib, "get_targets", fake_get_targets)
    merged = chembl_lib.extend_target(df)
    assert "chembl_pref_name" in merged.columns
    assert merged.loc[0, "chembl_pref_name"] == "nameA"
