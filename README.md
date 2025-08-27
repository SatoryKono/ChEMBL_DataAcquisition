# ChEMBL Data Acquisition

Utilities for downloading and integrating target information from
[ChEMBL](https://www.ebi.ac.uk/chembl/),
[UniProt](https://www.uniprot.org/) and the
[IUPHAR Guide to Pharmacology](https://www.guidetopharmacology.org/).
It also provides helpers for collecting publication metadata from
PubMed, Semantic Scholar, OpenAlex and CrossRef.


Five command line tools are available:


`get_target_data.py`
    Query biological data sources individually or run the combined
    pipeline to generate a unified CSV table.

`get_document_data.py`
    Retrieve document information from the services listed above.

`get_assay_data.py`
    Fetch assay information from the ChEMBL API for a list of assay IDs.


`get_activity_data.py`
    Fetch activity information from the ChEMBL API for a list of activity IDs.

`get_testitem_data.py`
    Fetch compound information from the ChEMBL API for a list of molecule IDs.


## Installation

Install the runtime dependencies:

```bash
pip install -r requirements.txt
```

Optional type stubs for development:

```bash
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

### Assay metadata

Retrieve assay information from the ChEMBL API for identifiers listed in
`assays.csv`:

```bash
python get_assay_data.py assays.csv assay_results.csv
```


### Activity metadata

Retrieve activity information from the ChEMBL API for identifiers listed in
`activities.csv`:

```bash
python get_activity_data.py activities.csv activity_results.csv
```

### Test item metadata

Retrieve compound information from the ChEMBL API for identifiers listed in
`molecules.csv`:

```bash
python get_testitem_data.py molecules.csv compound_results.csv
```


### Document metadata

Fetch PubMed, Semantic Scholar, OpenAlex and CrossRef records for PMIDs
listed in `pmids.csv`:

```bash
python get_document_data.py pubmed pmids.csv document_data.csv
```

Retrieve document information from the ChEMBL API for the identifiers in
`docs.csv`:

```bash
python get_document_data.py chembl docs.csv chembl_docs.csv
```

Run both pipelines and merge the outputs:

```bash
python get_document_data.py all docs.csv merged_docs.csv
```

get_target(chembl_target_id) / get_targets(ids, chunk_size=50) - fetch single/multiple targets, normalize fields (name, ChEMBL ID, HGNC, EC, etc.), robust to network errors.
get_assay(chembl_assay_id) / get_assays(ids, chunk_size=50) - load and normalize assay records.
get_activity(activity_id) / get_activities(ids, chunk_size=50) - load and normalize activity records.
get_testitem(ids, chunk_size=50) - load and normalize compound records.

get_document(chembl_document_id) / get_documents(ids, chunk_size=50) - publication metadata from ChEMBL.
extend_target(df, chembl_column="task_chembl_id", chunk_size=50) - join an input table with extended ChEMBL target information.
load_targets(path), load_families(path) - read CSV as strings, normalize headers, validate required columns.
IUPHARData.from_files(target_path, family_path) - container with target_df and family_df.
family_chain(start_id) - build parent_family_id chain.
target_id_by_uniprot / hgnc_name / hgnc_id / gene / name - map target IDs by identifiers and synonyms.

family_id_by_name, all_id, all_name - recover full ID/name paths.
map_uniprot_file(input_csv, output_csv, sep=",", encoding="utf-8") - read UniProt CSV, resolve to IUPHAR, add type/class/subclass, chain, full paths, write output.

ClassificationRecord - classification result structure.

Pipeline: process_compound(compound_name) - single call to retrieve all fields.

Load: fetch_uniprot(uniprot_id) - direct REST request.
extract_names(data) - protein and gene names.
extract_organism(data) - taxonomy (genus, superkingdom, phylum, taxon_id).
extract_keywords(data) - MF/CC keywords, EC, subcellular location, topology, membrane flags.

extract_ptm(data) - PTM feature flags.

extract_isoform(data) - isoform names/IDs/synonyms.
extract_crossrefs(data) - selected external DBs.
extract_activity(data) - catalytic reactions and EC numbers.

iter_ids(csv_path, ...) - yield UniProt IDs from CSV.

collect_info(uid, data_dir="uniprot") - gather all fields from local <uid>.json.
process(input_csv, output_csv, data_dir="uniprot", ...) - batch processing to CSV.
read_ids(path, column="chembl_id", ...) - safe reading of IDs from CSV.

build_parser() - defines CLI with subcommands:

CLI (get_target_data.py) - entry point.
test_iuphar_mapping.py - verifies map_uniprot_file maps UniProt ID to target ID and writes expected CSV.

test_read_ids.py - validates read_ids: filtering of empty/#N/A, error when column missing.
Project structure and purpose of modules
chembl_library.py

Role: utilities for retrieving targets, assays, and documents from ChEMBL, with wrappers for assembling results into pandas.DataFrame.

Key functions:

get_target(chembl_target_id) / get_targets(ids, chunk_size=50) — fetch single/multiple targets, normalize fields (name, ChEMBL ID, HGNC, EC, etc.), robust to network errors.

get_assay(chembl_assay_id) / get_assays(ids, chunk_size=50) — load and normalize assay records.

get_document(chembl_document_id) / get_documents(ids, chunk_size=50) — publication metadata from ChEMBL.

extend_target(df, chembl_column="task_chembl_id", chunk_size=50) — join an input table with extended ChEMBL target information.

iuphar_library.py

Role: core for working with preprocessed IUPHAR CSV: loading target/family tables, building classification chains, mapping UniProt ? IUPHAR, heuristic classification by EC/name.

Key components:

Data loading and validation:

load_targets(path), load_families(path) — read CSV as strings, normalize headers, validate required columns.

IUPHARData.from_files(target_path, family_path) — container with target_df and family_df.

Navigation and mapping:

family_chain(start_id) — build parent_family_id chain.

target_id_by_uniprot / hgnc_name / hgnc_id / gene / name — map target IDs by identifiers and synonyms.

family_id_by_name, all_id, all_name — recover full ID/name paths.

Batch mapping:

map_uniprot_file(input_csv, output_csv, sep=",", encoding="utf-8") — read UniProt CSV, resolve to IUPHAR, add type/class/subclass, chain, full paths, write output.

Classification:

ClassificationRecord — classification result structure.

IUPHARClassifier with methods by_target_id, by_uniprot_id, by_family_id, by_ec_number, by_name, and get(...).

pubchem_library.py

Role: client for the PubChem REST API: CID lookup, compound names, properties, and structured records.

Key functions:

Utilities: url_encode, make_request(url, delay=3.0), validate_cid.

CID lookup: get_cid(compound_name) (exact), get_all_cid(compound_name) (partial).

Metadata: get_standard_name(cid), get_properties(cid) ? Properties (IUPAC, formula, SMILES, InChI/Key).

Pipeline: process_compound(compound_name) — single call to retrieve all fields.

uniprot_library.py

Role: parsing UniProt JSON and batch processing local annotation files.

Key functions:

Load: fetch_uniprot(uniprot_id) — direct REST request.

Extract from entry:

extract_names(data) — protein and gene names.

extract_organism(data) — taxonomy (genus, superkingdom, phylum, taxon_id).

extract_keywords(data) — MF/CC keywords, EC, subcellular location, topology, membrane flags.

extract_ptm(data) — PTM feature flags.

extract_isoform(data) — isoform names/IDs/synonyms.

extract_crossrefs(data) — selected external DBs.

extract_activity(data) — catalytic reactions and EC numbers.

Streams:

iter_ids(csv_path, ...) — yield UniProt IDs from CSV.

collect_info(uid, data_dir="uniprot") — gather all fields from local <uid>.json.

process(input_csv, output_csv, data_dir="uniprot", ...) — batch processing to CSV.

get_target_data.py

Role: CLI wrapper around UniProt/ChEMBL/IUPHAR libraries.

Key parts:

read_ids(path, column="chembl_id", ...) — safe reading of IDs from CSV.

build_parser() — defines CLI with subcommands:

uniprot ? calls uniprot_library.process.

chembl ? reads IDs and calls chembl_library.get_targets.

iuphar ? uses IUPHARData.from_files and map_uniprot_file.

Runtime: configure_logging, run_uniprot, run_chembl, run_iuphar, main.

Module interaction

CLI (get_target_data.py) — entry point.

uniprot subcommand ? uniprot_library.process.

chembl subcommand ? chembl_library.get_targets.

iuphar subcommand ? IUPHARData + map_uniprot_file.

Tests

test_iuphar_mapping.py — verifies map_uniprot_file maps UniProt ID to target ID and writes expected CSV.

test_read_ids.py — validates read_ids: filtering of empty/#N/A, error when column missing.