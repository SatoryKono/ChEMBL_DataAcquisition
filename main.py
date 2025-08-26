"""Command-line interface for the PubChem client."""
from __future__ import annotations

import argparse
import csv
import json
import logging
from pathlib import Path
from typing import Iterable, List

from mylib.pubchem import process_compound


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Query PubChem for compound data")
    parser.add_argument(
        "--name",
        action="append",
        dest="names",
        help="Compound name. Can be used multiple times.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        help="Path to a text file with one compound name per line",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output file path",
    )
    parser.add_argument(
        "--format",
        choices=["json", "csv"],
        default="json",
        help="Output format",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        help="Logging level (e.g. INFO, DEBUG)",
    )
    return parser.parse_args()


def read_names(input_path: Path | None, cli_names: List[str] | None) -> List[str]:
    """Collect compound names from CLI or file."""
    names: List[str] = []
    if input_path:
        if not input_path.exists():
            raise FileNotFoundError(f"Input file not found: {input_path}")
        with input_path.open("r", encoding="utf-8") as handle:
            names.extend(line.strip() for line in handle if line.strip())
    if cli_names:
        names.extend(cli_names)
    if not names:
        raise ValueError("No compound names provided")
    return names


def write_output(records: Iterable[dict], path: Path, fmt: str) -> None:
    """Write *records* to *path* in given format."""
    if fmt == "json":
        with path.open("w", encoding="utf-8") as handle:
            json.dump(list(records), handle, ensure_ascii=False, indent=2)
    else:
        records = list(records)
        if records:
            headers = list(records[0].keys())
        else:
            headers = []
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=headers)
            writer.writeheader()
            writer.writerows(records)


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO))
    names = read_names(args.input, args.names)
    records = [process_compound(name) for name in names]
    write_output(records, args.output, args.format)


if __name__ == "__main__":
    main()
