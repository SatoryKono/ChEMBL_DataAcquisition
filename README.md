# ChEMBL Data Acquisition

Utilities for downloading and integrating target information from
[ChEMBL](https://www.ebi.ac.uk/chembl/),
[UniProt](https://www.uniprot.org/) and the
[IUPHAR Guide to Pharmacology](https://www.guidetopharmacology.org/).

The repository exposes a single command line script `get_target_data.py`
which can query each data source individually or run the combined
pipeline to generate a unified CSV table.

## Installation

```bash
pip install pandas>=2.0 requests>=2.31
# Optional: type stubs for development
pip install pandas-stubs types-requests
```

## Usage

Each sub-command reads an input CSV and writes a new CSV with the
requested annotations. Delimiters and encodings can be customised with
`--sep` and `--encoding`.

Fetch ChEMBL targets for the identifiers in `targets.csv`:

```bash
python get_target_data.py chembl targets.csv chembl_results.csv
```

Parse UniProt JSON files in `uniprot/` and enrich the accessions listed
in `ids.csv`:

```bash
python get_target_data.py uniprot ids.csv uniprot_results.csv --data-dir uniprot
```

Map UniProt IDs to IUPHAR classifications:

```bash
python get_target_data.py iuphar uniprot_results.csv iuphar_results.csv \
    --target-csv data/_IUPHAR_target.csv \
    --family-csv data/_IUPHAR_family.csv
```

Run the full pipeline (ChEMBL → UniProt → IUPHAR) starting from a CSV
with a single `chembl_id` column:

```bash
python get_target_data.py all targets.csv merged.csv \
    --data-dir uniprot \
    --target-csv data/_IUPHAR_target.csv \
    --family-csv data/_IUPHAR_family.csv
```

Intermediate results can be retained with `--chembl-out`, `--uniprot-out`
and `--iuphar-out`.

## Development

Format and lint the code:

```bash
black .
ruff check .
```

Run the test-suite:

```bash
pytest -q
```

Static type checking (requires stub packages):

```bash
mypy --explicit-package-bases get_target_data.py library
```

## Example Data

A small `targets.csv` file is provided in the repository containing a
`chembl_id` column for smoke testing.

