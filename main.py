"""CLI entry point for target table post-processing."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Sequence

from mylib import postprocess_file

logger = logging.getLogger(__name__)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Target table post-processing")
    parser.add_argument("input_csv", type=Path, help="Raw merged CSV file")
    parser.add_argument("output_csv", type=Path, help="Destination for cleaned CSV")
    parser.add_argument("--sep", default=",", help="CSV delimiter")
    parser.add_argument("--encoding", default="utf8", help="CSV encoding")
    parser.add_argument("--log-level", default="INFO", help="Logging level")
    return parser


def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=level.upper(), format="%(levelname)s:%(name)s:%(message)s"
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    configure_logging(args.log_level)
    try:
        postprocess_file(
            args.input_csv, args.output_csv, sep=args.sep, encoding=args.encoding
        )
        return 0
    except FileNotFoundError as exc:
        logger.error("%s", exc)
        return 1
    except OSError as exc:
        logger.error("failed to write output: %s", exc)
        return 1


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
