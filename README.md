# ChEMBL Data Acquisition

Utilities for retrieving target and assay information from the
[ChEMBL](https://www.ebi.ac.uk/chembl/) REST API.

## Installation

The project requires Python 3.12 or later.  Install dependencies with:

```bash
pip install -r requirements.txt
```

## Usage

Fetch two targets and write the result to `targets.csv`:

```bash
python main.py --type target --ids CHEMBL25,CHEMBL26 --output targets.csv
```

Fetch a list of assays and write the result to `assays.csv`:

```bash
python main.py --type assay --ids CHEMBL1,CHEMBL2 --output assays.csv
```

## Development

Code style and static analysis:

```bash
ruff check .
black .
mypy .
```

Run the tests:

```bash
pytest
```

## Dependencies

- `requests >= 2.31`
- `pandas >= 2.0`
- `pytest >= 7.0` (for running tests)
- `ruff`, `black`, `mypy` (optional, for development)
