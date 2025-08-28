

"""Command-line interface for the classification pipeline."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from library.classification_library import main_process, read_table, write_table


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Run classification pipeline")
    parser.add_argument("--input", required=True, help="Input CSV file path")
    parser.add_argument("--output", required=True, help="Output CSV file path")
    parser.add_argument("--sep", default=",", help="Field delimiter")
    parser.add_argument("--encoding", default="utf-8", help="File encoding")
    parser.add_argument(
        "--log-level", default="INFO", help="Logging level (e.g. INFO, DEBUG)"
    )
    return parser.parse_args()


def configure_logging(level: str) -> None:
    """Configure basic logging."""

    logging.basicConfig(level=getattr(logging, level.upper(), logging.INFO))


def main() -> None:
    """Entry point for CLI."""

    args = parse_args()
    configure_logging(args.log_level)

    input_path = Path(args.input)
    output_path = Path(args.output)

    logging.info("Reading input data from %s", input_path)
    df = read_table(str(input_path), sep=args.sep, encoding=args.encoding)
    logging.info("Running transformation pipeline")
    result = main_process(df)
    logging.info("Writing output to %s", output_path)
    write_table(result, str(output_path), sep=args.sep, encoding=args.encoding)


if __name__ == "__main__":
    main()