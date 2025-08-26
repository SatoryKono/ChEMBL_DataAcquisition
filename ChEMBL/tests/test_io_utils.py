import sys
from pathlib import Path

import pytest

sys.path.append(str(Path(__file__).resolve().parents[1] / "script"))
import io_utils  # type: ignore


def test_read_ids_returns_filtered_list():
    src = Path(__file__).resolve().parent / "data" / "ids.csv"
    assert io_utils.read_ids(src) == ["CHEMBL1", "CHEMBL2"]


def test_read_ids_missing_column(tmp_path):
    bad = tmp_path / "bad.csv"
    bad.write_text("other\n1\n")
    with pytest.raises(ValueError):
        io_utils.read_ids(bad, column="chembl_id")
