# ChEMBL_DataAcquisition

Utilities for working with the IUPHAR portion of the ChEMBL data set.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

The project ships both a reusable library and a small command line interface.

```bash
python main.py \
    --target-file tests/data/target.csv \
    --family-file tests/data/family.csv \
    --uniprot Q11111
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

The resulting file includes the resolved ``target_id`` together with the
IUPHAR classification (class, subclass and family chain) for each UniProt
accession. If a UniProt lookup fails, the script optionally consults
``HGNC_name``, ``HGNC_id``, ``gene_name`` and ``synonyms`` columns (when
present) to resolve the identifier.

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
