# ChEMBL_DataAcquisition

Utilities and scripts for retrieving chemical data.

## PubChem Client

A small client library and CLI for fetching compound information from the
[PubChem](https://pubchem.ncbi.nlm.nih.gov/) REST API.

### Installation

```bash
pip install -r requirements.txt
```

### Example

Fetch data for "aspirin" and save as JSON:

```bash
python main.py --name aspirin --output aspirin.json --format json
```

### Testing and quality tools

```bash
pytest -q
ruff .
mypy mylib
```
