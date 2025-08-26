# ChEMBL_DataAcquisition


Utilities for retrieving and processing chemical information from public
sources such as ChEMBL, IUPHAR, PubChem and UniProt.  Each subdirectory
contains a small library of helper functions and a command line script to
perform common data acquisition tasks.


## Installation

The project targets **Python 3.12**. Install the required dependencies:

```bash
pip install -r requirements.txt
```

Main dependencies

- pandas >= 2.0
- requests >= 2.31
- urllib3 >= 1.26
- responses >= 0.23 (tests)
- pytest >= 7.4 (tests)
- ruff >= 0.1 (lint)
- black >= 23.0 (formatting)
- mypy >= 1.8 (type checking)

## Usage

Each data source provides both a reusable library and a command line
interface. Examples:

### IUPHAR

```bash
python get_IUPHAR.py \
    --target-file tests/data/target.csv \
    --family-file tests/data/family.csv \
    --uniprot Q11111
```

Batch process UniProt accessions from a CSV file containing a
``uniprot_id`` column and write the mapping to another CSV:

```bash
python get_IUPHAR.py \
    --target-file tests/data/target.csv \
    --family-file tests/data/family.csv \
    --uniprot-file tests/data/uniprot_input.csv \
    --output-file results.csv
```

The resulting file includes the resolved ``target_id`` together with the
IUPHAR classification (class, subclass and family chain) for each
UniProt accession. If a UniProt lookup fails, the script optionally
consults ``hgnc_name``, ``hgnc_id``, ``gene_name`` and ``synonyms``
columns (when present) to resolve the identifier. When no target can be
found, an optional ``ec_number`` column is used to infer
``IUPHAR_type``, ``IUPHAR_class`` and ``IUPHAR_subclass``.

### UniProt

```bash
python get_uniprot_data.py tests/data/uniprot_input.csv results.csv
```

### PubChem

```bash
python get_PubChem.py --name "aspirin"
```

### ChEMBL

```bash
python get_chembl_data.py --target --input tests/data/target.csv  --output ChEMBL_target.csv 

```

## Development

Run quality checks and tests before submitting changes:

```bash
ruff check .
mypy uniprot/get_uniprot_data.py  # adjust paths as needed
pytest
```

 
