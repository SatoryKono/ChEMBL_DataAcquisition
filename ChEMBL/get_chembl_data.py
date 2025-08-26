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
from pathlib import Path
from typing import Sequence

from script import chembl_lib as cc
from script import io_utils as io


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command line arguments.

    Parameters
    ----------
    argv:
        Optional sequence of argument strings. If ``None`` the arguments
        are read from ``sys.argv``.

    Returns
    -------
    argparse.Namespace
        Parsed command line options.
    """
    parser = argparse.ArgumentParser(description="Fetch data from ChEMBL API")
    parser.add_argument(
        "--type",
        choices={"target", "assay", "document"},
        required=True,
        help="Type of entity to retrieve",
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Path to CSV file containing identifiers",
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


def main(argv: Sequence[str] | None = None) -> int:
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
        ids = io.read_ids(
            args.input, column=args.column, sep=args.sep, encoding=args.encoding
        )
    except (FileNotFoundError, ValueError) as exc:
        logging.error("Failed to read identifiers: %s", exc)
        return 1
    if args.type == "target":
        df = cc.get_targets(ids)
    elif args.type == "assay":
        df = cc.get_assays(ids)
    else:
        df = cc.get_documents(ids)

    try:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.output, index=False, encoding=args.encoding)
    except OSError as exc:
        logging.error("Failed to write output to %s: %s", args.output, exc)
        return 1

    logging.info("Wrote %d records to %s", len(df), args.output)
    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
