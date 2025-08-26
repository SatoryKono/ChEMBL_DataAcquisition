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
