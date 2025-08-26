# ChEMBL_DataAcquisition

Utilities for working with the IUPHAR portion of the ChEMBL data set.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

The project ships both a reusable library and a small command line interface.

```bash
python main.py --uniprot Q11111 --target-id T1
```

Batch process UniProt accessions from a CSV file containing a ``uniprot_id``
column and write the mapping to another CSV:

```bash
python main.py \
    --target-file tests/data/target.csv \
    --family-file tests/data/family.csv \
    --uniprot-file tests/data/uniprot_input.csv \
    --output-file results.csv
```

## Development

Formatting and linting are recommended:

```bash
black mylib tests main.py
ruff mylib tests main.py
mypy mylib
```

Run the unit tests with `pytest`:

```bash
pytest
```
