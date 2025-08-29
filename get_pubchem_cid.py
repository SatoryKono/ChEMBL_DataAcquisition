"""Command line utility to look up PubChem CIDs from structural identifiers."""

from __future__ import annotations

import argparse
import logging
from typing import Optional

from library import pubchem_library as pl


logger = logging.getLogger(__name__)


def build_parser() -> argparse.ArgumentParser:
    """Create and return the CLI argument parser."""

    parser = argparse.ArgumentParser(
        description="Lookup PubChem CIDs from SMILES, InChI or InChIKey"
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        help="Logging level (DEBUG, INFO, WARNING)",
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--smiles", help="SMILES string to query")
    group.add_argument("--inchi", help="InChI string to query")
    group.add_argument("--inchikey", help="InChIKey to query")
    return parser


def run(args: argparse.Namespace) -> int:
    """Execute the lookup based on parsed ``args``."""

    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO))
    cid: Optional[str]
    if args.smiles:
        cid = pl.get_cid_from_smiles(args.smiles)
    elif args.inchi:
        cid = pl.get_cid_from_inchi(args.inchi)
    else:
        cid = pl.get_cid_from_inchikey(args.inchikey)
    if cid:
        print(cid)
        return 0
    logger.error("CID not found")
    return 1


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    raise SystemExit(run(args))


if __name__ == "__main__":
    main()
