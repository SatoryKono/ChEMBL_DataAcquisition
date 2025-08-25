import csv
import json
import argparse
from pathlib import Path
from typing import List, Dict
import requests

API_URL = "https://rest.uniprot.org/uniprotkb/{id}.json"

def read_ids(csv_path: Path) -> List[str]:
    with csv_path.open(newline='') as f:
        # Attempt to read using DictReader to support column names
        reader = csv.DictReader(f)
        if reader.fieldnames and 'uniprot_id' in reader.fieldnames:
            return [row['uniprot_id'] for row in reader if row.get('uniprot_id')]
        # Fallback: use first column
        f.seek(0)
        reader2 = csv.reader(f)
        return [row[0] for row in reader2 if row]

def fetch_uniprot(uniprot_id: str) -> Dict:
    url = API_URL.format(id=uniprot_id)
    resp = requests.get(url)
    resp.raise_for_status()
    return resp.json()

def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch UniProt JSON data for IDs from CSV")
    parser.add_argument(
        "csv_file",
        type=Path,
        help="Path to CSV file with uniprot_id column or single column",
    )
    parser.add_argument(
        "output_dir",
        type=Path,
        help="Directory where individual UniProt JSON files will be saved",
    )
    args = parser.parse_args()

    ids = read_ids(args.csv_file)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for uid in ids:
        try:
            data = fetch_uniprot(uid)
        except Exception as exc:
            print(f"Failed to fetch {uid}: {exc}")
            continue

        out_file = args.output_dir / f"{uid}.json"
        with out_file.open("w") as out:
            json.dump(data, out, indent=2)

if __name__ == "__main__":
    main()
