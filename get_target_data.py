"""Command line interface for retrieving target data from external sources."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path
from typing import Sequence

from library import uniprot_library as uu
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
    """Create and return the top-level CLI argument parser.

    The command line interface is organised in sub-commands.  Only the
    ``uniprot`` sub-command is currently implemented and allows batch
    processing of UniProt identifiers contained in a CSV file.
    """

    parser = argparse.ArgumentParser(description="Target data utilities")
    parser.add_argument(
        "--log-level",
        default="INFO",
        help="Logging level (DEBUG, INFO, WARNING)",
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    uniprot = subparsers.add_parser(
        "uniprot", help="Extract information for UniProt accessions"
    )
    uniprot.add_argument(
        "input_csv",
        type=Path,
        help="CSV file containing a 'uniprot_id' column",
    )
    uniprot.add_argument(
        "output_csv",
        type=Path,
        help="Destination CSV file for the extracted information",
    )
    uniprot.add_argument(
        "--sep",
        default=",",
        help="CSV delimiter used for input and output files",
    )
    uniprot.add_argument(
        "--encoding",
        default="utf8",
        help="File encoding for input and output CSV files",
    )
    uniprot.add_argument(
        "--data-dir",
        default="uniprot",
        help="Directory containing '<uniprot_id>.json' files",
    )
    uniprot.set_defaults(func=run_uniprot)

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


def run_uniprot(args: argparse.Namespace) -> int:
    """Execute the ``uniprot`` sub-command.

    Parameters
    ----------
    args:
        Parsed command-line arguments specific to the ``uniprot`` sub-command.

    Returns
    -------
    int
        Zero on success, non-zero on failure.
    """

    uu.process(
        input_csv=str(args.input_csv),
        output_csv=str(args.output_csv),
        data_dir=str(args.data_dir),
        sep=args.sep,
        encoding=args.encoding,
    )
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point."""

    parser = build_parser()
    args = parser.parse_args(argv)
    configure_logging(args.log_level)
    if hasattr(args, "func"):
        return args.func(args)
    parser.print_help()
    return 1


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
