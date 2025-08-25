"""Command line interface for ChEMBL data retrieval.

The script provides a small CLI for downloading target or assay
information from the ChEMBL REST API and storing the result as a CSV
file.  It demonstrates the use of :mod:`mylib.chembl_client`.

Example
-------
Fetch target records for two identifiers and write them to ``out.csv``::

    python main.py --type target --ids CHEMBL1,CHEMBL2 --output out.csv
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Iterable

from mylib import chembl_client as cc


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Fetch data from ChEMBL API")
    parser.add_argument(
        "--type",
        choices={"target", "assay"},
        required=True,
        help="Type of entity to retrieve",
    )
    parser.add_argument(
        "--ids",
        required=True,
        help="Comma separated list of ChEMBL identifiers",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Path to the output CSV file",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        help="Logging level (e.g. INFO, DEBUG)",
    )
    return parser.parse_args(argv)


def configure_logging(level: str) -> None:
    """Configure basic logging."""
    logging.basicConfig(level=getattr(logging, level.upper(), logging.INFO))


def main(argv: Iterable[str] | None = None) -> int:
    """Entry point for command line execution."""
    args = parse_args(argv)
    configure_logging(args.log_level)

    ids = args.ids.split(",")
    if args.type == "target":
        df = cc.get_targets(ids)
    else:
        df = cc.get_assays(ids)

    df.to_csv(args.output, index=False)
    logging.info("Wrote %d records to %s", len(df), args.output)
    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())

