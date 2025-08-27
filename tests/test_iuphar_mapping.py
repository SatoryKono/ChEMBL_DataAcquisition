from pathlib import Path
import sys

import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1]))

from library.iuphar_library import IUPHARData


def test_map_uniprot_file(tmp_path: Path) -> None:
    target_csv = tmp_path / "target.csv"
    family_csv = tmp_path / "family.csv"
    input_csv = tmp_path / "input.csv"
    output_csv = tmp_path / "output.csv"

    target_csv.write_text(
        "target_id,swissprot,hgnc_name,hgnc_id,gene_name,synonyms,family_id,target_name,type\n"
        "0001,Q12345,GeneX,1,GENE1,Syn1|Syn2,100,TargetX,Enzyme.Lyase\n",
        encoding="utf-8",
    )
    family_csv.write_text(
        "family_id,family_name,parent_family_id,target_id,type\n"
        "100,FamilyX,,0001,Enzyme.Lyase\n",
        encoding="utf-8",
    )
    input_csv.write_text(
        "uniprot_id;GuidetoPHARMACOLOGY;hgnc_name\n"
        "Q99999;0001;\n"
        "Q12345;;\n",
        encoding="utf-8",
    )

    data = IUPHARData.from_files(target_csv, family_csv)
    df = data.map_uniprot_file(input_csv, output_csv, sep=";")

    assert output_csv.exists()
    out_df = pd.read_csv(output_csv, sep=";", dtype=str)
    assert list(out_df["target_id"]) == ["0001", "0001"]
    assert list(df["target_id"]) == ["0001", "0001"]
