from pathlib import Path
import sys

import pytest

sys.path.append(str(Path(__file__).resolve().parents[1]))

from get_target_data import read_ids


def test_read_ids(tmp_path: Path) -> None:
    csv_file = tmp_path / "ids.csv"
    csv_file.write_text("chembl_id\nCHEMBL1\n#N/A\n\nCHEMBL2\n", encoding="utf8")
    assert read_ids(csv_file) == ["CHEMBL1", "CHEMBL2"]


def test_read_ids_missing_column(tmp_path: Path) -> None:
    csv_file = tmp_path / "ids.csv"
    csv_file.write_text("other\n1\n2\n", encoding="utf8")
    with pytest.raises(ValueError):
        read_ids(csv_file, column="chembl_id")
