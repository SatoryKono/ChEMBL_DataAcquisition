from pathlib import Path
import sys
import json
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1]))

from library import uniprot_library as uu

DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "uniprot"


def test_iter_ids(tmp_path: Path) -> None:
    csv_file = tmp_path / "uids.csv"
    csv_file.write_text("uniprot_id\nQ99558\nQ99519\n", encoding="utf8")
    assert list(uu.iter_ids(csv_file)) == ["Q99558", "Q99519"]


def test_collect_info_real_json() -> None:
    info = uu.collect_info("Q99558", data_dir=str(DATA_DIR))
    assert info["uniprot_id"] == "Q99558"
    assert "MAP3K14" in info["names"]
    assert info["genus"] == "Homo"


def test_collect_info_downloads_missing(tmp_path: Path, monkeypatch) -> None:
    uid = "Q99558"
    sample = json.loads((DATA_DIR / f"{uid}.json").read_text(encoding="utf8"))

    def fake_fetch(uniprot_id: str) -> dict:
        assert uniprot_id == uid
        return sample

    monkeypatch.setattr(uu, "fetch_uniprot", fake_fetch)
    data_dir = tmp_path / "uniprot"
    info = uu.collect_info(uid, data_dir=str(data_dir))
    assert info["uniprot_id"] == uid
    assert (data_dir / f"{uid}.json").exists()


def test_process_writes_expected(tmp_path: Path) -> None:
    input_csv = tmp_path / "uids.csv"
    input_csv.write_text("uniprot_id\nQ99558\n", encoding="utf8")
    output_csv = tmp_path / "out.csv"
    uu.process(
        input_csv=str(input_csv),
        output_csv=str(output_csv),
        data_dir=str(DATA_DIR),
    )
    df = pd.read_csv(output_csv, dtype=str)
    assert df.loc[0, "uniprot_id"] == "Q99558"
    assert "Mitogen-activated protein kinase kinase kinase 14" in df.loc[0, "names"]
