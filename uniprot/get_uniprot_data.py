#!/usr/bin/env python3
"""Read UniProt IDs from a CSV and extract names from JSON files.

For each UniProt ID listed in the input CSV, this script expects a
corresponding JSON file named ``<ID>.json`` in the ``Data`` directory.
It extracts all protein and gene names using :func:`extract_names` from
``uniprot_names`` and writes the results to an output CSV.

Usage:
    python uniprot_batch_names.py <input_csv> <output_csv>

The input CSV must contain a column named ``UniProt_id``.
"""

import csv
import json
import os
import sys
from typing import Iterable

from script import uniprot_utils as uu





def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print('Usage: python uniprot_batch_names.py <input_csv> <output_csv>', file=sys.stderr)
        return 1
    input_csv, output_csv = argv[1], argv[2]
    uu.process(input_csv, output_csv)
    return 0


if __name__ == '__main__':
    raise SystemExit(main(sys.argv))
