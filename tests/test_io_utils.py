from pathlib import Path

import pytest

from script.io_utils import read_ids

DATA_DIR = Path(__file__).parent / "data"


def test_read_ids(tmp_path):
    path = DATA_DIR / "ids.csv"
    assert read_ids(path) == ["CHEMBL1", "CHEMBL2"]


def test_read_ids_missing_column():
    path = DATA_DIR / "ids.csv"
    with pytest.raises(ValueError):
        read_ids(path, column="missing")


def test_read_ids_missing_file():
    with pytest.raises(FileNotFoundError):
        read_ids(Path("does_not_exist.csv"))


def test_read_ids_strips_whitespace(tmp_path):
    path = tmp_path / "ids.csv"
    path.write_text("chembl_id\n CHEMBL1 \nCHEMBL2\n", encoding="utf8")
    assert read_ids(path) == ["CHEMBL1", "CHEMBL2"]
