import pandas as pd

from script import chembl_lib


def test_extend_target_uses_bulk(monkeypatch):
    calls: list[list[str]] = []

    def fake_get_targets(ids, chunk_size=50):  # noqa: D401 - test stub
        calls.append(list(ids))
        return pd.DataFrame(
            [
                {
                    "target_chembl_id": "CHEMBL1",
                    "pref_name": "name",
                    "component_description": "desc",
                    "component_id": 1,
                    "relationship": "rel",
                    "gene": "gene",
                    "chembl_alternative_name": "alt",
                    "ec_code": "ec",
                    "HGNC_name": "hgnc",
                    "HGNC_id": "123",
                }
            ]
        )

    monkeypatch.setattr(chembl_lib, "get_targets", fake_get_targets)

    df = pd.DataFrame({"task_chembl_id": ["CHEMBL1"]})
    result = chembl_lib.extend_target(df)

    assert calls == [["CHEMBL1"]]
    assert "chembl_pref_name" in result.columns
    assert result.loc[0, "chembl_pref_name"] == "name"
