import pandas as pd

import get_chembl_data
from script import chembl_lib as cc
from script import io_utils as io


def test_main_success(tmp_path, monkeypatch):
    ids_file = tmp_path / "ids.csv"
    ids_file.write_text("chembl_id\nCHEMBL1\n", encoding="utf8")

    monkeypatch.setattr(io, "read_ids", lambda *a, **k: ["CHEMBL1"])
    monkeypatch.setattr(
        cc,
        "get_targets",
        lambda ids: pd.DataFrame({"target_chembl_id": ids, "pref_name": ["n"]}),
    )

    out_file = tmp_path / "out.csv"
    argv = [
        "--type",
        "target",
        "--input",
        str(ids_file),
        "--output",
        str(out_file),
    ]

    result = get_chembl_data.main(argv)
    assert result == 0
    assert out_file.exists()


def test_main_read_ids_error(monkeypatch, tmp_path):
    ids_file = tmp_path / "ids.csv"

    def raise_error(*_a, **_k):
        raise ValueError("bad")

    monkeypatch.setattr(io, "read_ids", raise_error)
    argv = [
        "--type",
        "target",
        "--input",
        str(ids_file),
        "--output",
        str(tmp_path / "out.csv"),
    ]
    assert get_chembl_data.main(argv) == 1
