"""Command line interface for ChEMBL data retrieval.

The script reads ChEMBL identifiers from a CSV file and downloads target,
assay, or document information from the ChEMBL REST API. Results are stored
in a CSV file. It demonstrates the use of :mod:`script` utilities.

Example
-------
Fetch target records listed in ``ids.csv`` and write them to ``out.csv``::

    python get_chembl_data.py --type target --input ids.csv --output out.csv
"""

from __future__ import annotations

import argparse
import logging
import sys

import argparse
import logging
from pathlib import Path

from pathlib import Path
from typing import Sequence

from library import chembl_library as cc
from library import uniprot_library as uu
from library.iuphar_library import IUPHARClassifier, IUPHARData



from __future__ import annotations

from pathlib import Path
import csv


def read_ids(
    path: str | Path,
    column: str = "chembl_id",
    sep: str = ",",
    encoding: str = "utf8",
) -> list[str]:
    """Read identifier values from a CSV file.

    Parameters
    ----------
    path : str or Path
        Path to the CSV file.
    column : str, optional
        Name of the column containing identifiers. Defaults to ``"chembl_id"``.
    sep : str, optional
        Field delimiter, by default a comma.
    encoding : str, optional
        File encoding, by default ``"utf8"``.

    Returns
    -------
    list of str
        Identifier values in the order they appear. Empty strings and ``"#N/A"``
        markers are discarded.

    Raises
    ------
    ValueError
        If ``column`` is not present in the input file.
    """
    try:
        with Path(path).open("r", encoding=encoding, newline="") as fh:
            reader = csv.DictReader(fh, delimiter=sep)
            if reader.fieldnames is None or column not in reader.fieldnames:
                raise ValueError(f"column '{column}' not found in {path}")
            ids: list[str] = []
            for row in reader:
                value = (row.get(column) or "").strip()
                if value and value != "#N/A":
                    ids.append(value)
            return ids
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"input file not found: {path}") from exc
    except csv.Error as exc:
        raise ValueError(f"malformed CSV in file: {path}: {exc}") from exc

def build_parser() -> argparse.ArgumentParser:
    """Create the CLI argument parser."""
    parser = argparse.ArgumentParser(description="IUPHAR data utilities")
    parser.add_argument(
        "--target-file",
        type=Path,
        default=Path("IUPHAR/_IUPHAR_target.csv"),
        help="Path to _IUPHAR_target.csv",
    )
    parser.add_argument(
        "--family-file",
        type=Path,
        default=Path("IUPHAR/_IUPHAR_family.csv"),
        help="Path to _IUPHAR_family.csv",
    )
    parser.add_argument("--target-id", help="Target ID for family lookup")
    parser.add_argument("--uniprot", help="UniProt accession for ID lookup")
    parser.add_argument(
        "--uniprot-file",
        type=Path,
        help="CSV file containing a uniprot_id column",
    )
    parser.add_argument(
        "--column",
        default="chembl_id",
        help="Name of the column with identifiers",
    )
    parser.add_argument(
        "--sep",
        default=",",
        help="CSV delimiter",
    )
    parser.add_argument(
        "--encoding",
        default="utf8",
        help="CSV file encoding",
    )

    parser.add_argument(
        "--output-file",
        type=Path,
        help="Output CSV file for UniProt mappings",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        help="Logging level (DEBUG, INFO, WARNING)",
    )
    return parser

def configure_logging(level: str) -> None:
    """Configure basic logging.

    Parameters
    ----------
    level:
        Logging verbosity such as ``"INFO"`` or ``"DEBUG"``.

    Returns
    -------
    None
        This function configures the root logger in-place and does not
        return a value.
    """
    logging.basicConfig(level=getattr(logging, level.upper(), logging.INFO))


def get_chembl_data(argv: Sequence[str] | None = None) -> int:
    """Entry point for command line execution.

    Parameters
    ----------
    argv:
        Optional sequence of command line arguments.

    Returns
    -------
    int
        Zero on success, non-zero on failure.
    """
    args = parse_args(argv)
    configure_logging(args.log_level)

    try:
        ids = read_ids(
            args.input, column=args.column, sep=args.sep, encoding=args.encoding
        )
    except (FileNotFoundError, ValueError) as exc:
        logging.error("Failed to read identifiers: %s", exc)
        return 1
    
        df = cc.get_targets(ids)
    

    try:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.output, index=False, encoding=args.encoding)
    except OSError as exc:
        logging.error("Failed to write output to %s: %s", args.output, exc)
        return 1

    logging.info("Wrote %d records to %s", len(df), args.output)
    return 0

def get_uniprot_data(argv: Sequence[str] | None = None) -> int:
    return uu.process(input_csv, output_csv)

if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
