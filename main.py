"""Command line interface for IUPHAR data utilities."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from mylib.transforms import IUPHARClassifier, IUPHARData


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


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO))

    data = IUPHARData.from_files(args.target_file, args.family_file)
    classifier = IUPHARClassifier(data)

    if args.uniprot:
        rec = classifier.by_uniprot_id(args.uniprot)
        print(f"Target ID for UniProt {args.uniprot}: {rec.IUPHAR_target_id}")
        print(f"Family chain for target {rec.IUPHAR_target_id}: {rec.IUPHAR_tree}")
        print(f"Full ID path: {data.all_id(rec.IUPHAR_target_id)}")
        print(f"Full name path: {data.all_name(rec.IUPHAR_target_id)}")

    if args.uniprot_file and args.output_file:
        classifier.data.map_uniprot_file(args.uniprot_file, args.output_file)
        print(f"Results written to {args.output_file}")

    if args.target_id:
        fam_id = data.from_target_family_id(args.target_id)
        print(
            f"Family chain for target {args.target_id}: {data.family_chain(fam_id) if fam_id else []}"
        )
        print(f"Full ID path: {data.all_id(args.target_id)}")
        print(f"Full name path: {data.all_name(args.target_id)}")


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    main()
