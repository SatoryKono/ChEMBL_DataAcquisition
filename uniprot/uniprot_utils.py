"""Utility functions for extracting information from UniProt data.

This module consolidates the logic previously spread across
``uniprot_names.py`` and ``uniprot_batch_names.py`` into a reusable
library.  It provides helpers to parse UniProt JSON structures, gather
protein and gene names, organism taxonomy, and batch-process multiple
entries from a CSV file.

The most commonly used functions are:

``fetch_uniprot(uniprot_id)``
    Retrieve a UniProt JSON entry from the REST API given an accession ID.

``extract_names(data)``
    Parse a UniProt JSON object and return a set of all protein and gene
    names found in the entry.

``extract_organism(data)``
    Extract genus, superkingdom, phylum and taxon ID information from a
    UniProt JSON object.  Returns a dictionary with these fields.

``iter_ids(csv_path)``
    Read a CSV file containing a ``UniProt_id`` column and yield each ID.

``collect_info(uid, data_dir="Data")``
    Given a UniProt accession and directory containing ``<uid>.json``
    files, return a dictionary with the accession, all names, and
    organism taxonomy data.

``process(input_csv, output_csv, data_dir="Data")``
    Batch-process a CSV of UniProt IDs and write an output CSV with
    names and organism information for each ID.
"""

from __future__ import annotations

import csv
import json
import os
from typing import Any, Dict, Iterable, List, Set

import requests

API_URL = "https://rest.uniprot.org/uniprotkb/{id}.json"

__all__ = [
    "fetch_uniprot",
    "extract_names",
    "extract_organism",
    "iter_ids",
    "collect_info",
    "process",
]


def fetch_uniprot(uniprot_id: str) -> Dict[str, Any]:
    """Fetch a UniProt JSON record from the public REST API.

    Args:
        uniprot_id: UniProt accession identifier to retrieve.

    Returns:
        The JSON-decoded response as a dictionary.

    Raises:
        requests.HTTPError: If the HTTP request returned an unsuccessful
            status code.
    """
    url = API_URL.format(id=uniprot_id)
    resp = requests.get(url)
    resp.raise_for_status()
    return resp.json()


def _collect_name_fields(name_obj: Dict[str, Any]) -> Iterable[str]:
    """Yield all full and short names from a UniProt name object."""
    if not isinstance(name_obj, dict):
        return []
    names: List[str] = []
    full = name_obj.get("fullName")
    if isinstance(full, dict):
        value = full.get("value")
        if value:
            names.append(value)
    short = name_obj.get("shortName")
    if isinstance(short, dict):
        value = short.get("value")
        if value:
            names.append(value)
    elif isinstance(short, list):
        for item in short:
            if isinstance(item, dict):
                value = item.get("value")
                if value:
                    names.append(value)
    return names


def _extract_protein_names(desc: Dict[str, Any]) -> Set[str]:
    names: Set[str] = set()
    if not isinstance(desc, dict):
        return names
    rec = desc.get("recommendedName")
    if isinstance(rec, dict):
        names.update(_collect_name_fields(rec))
    for key in ("alternativeNames", "submissionNames"):
        items = desc.get(key) or []
        for item in items:
            names.update(_collect_name_fields(item))
    return names


def _extract_gene_names(entry: Dict[str, Any]) -> Set[str]:
    names: Set[str] = set()
    for gene in entry.get("genes", []):
        if not isinstance(gene, dict):
            continue
        main = gene.get("geneName")
        if isinstance(main, dict):
            value = main.get("value")
            if value:
                names.add(value)
        for syn in gene.get("synonyms", []):
            if isinstance(syn, dict):
                value = syn.get("value")
                if value:
                    names.add(value)
    return names


def extract_names(data: Any) -> Set[str]:
    """Return all protein and gene names found in ``data``.

    Args:
        data: A UniProt JSON structure, list of entries, or search results
            containing UniProt entries.

    Returns:
        A set of name strings aggregated from protein and gene sections.
    """
    names: Set[str] = set()
    if isinstance(data, dict) and "results" in data:
        entries = data["results"]
    elif isinstance(data, list):
        entries = data
    else:
        entries = [data]
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        names.update(_extract_protein_names(entry.get("proteinDescription", {})))
        names.update(_extract_gene_names(entry))
    return names


def extract_organism(data: Any) -> Dict[str, str]:
    """Return organism taxonomy information for the entry in ``data``.

    Args:
        data: A UniProt JSON structure, list of entries, or search results
            containing UniProt entries.

    Returns:
        A dictionary with keys ``genus``, ``superkingdom``, ``phylum`` and
        ``taxon_id``. Empty strings are returned when a field is missing.
    """
    result = {"genus": "", "superkingdom": "", "phylum": "", "taxon_id": ""}
    if isinstance(data, dict) and "results" in data:
        entries = data["results"]
    elif isinstance(data, list):
        entries = data
    else:
        entries = [data]
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        org = entry.get("organism", {})
        if not isinstance(org, dict):
            continue
        taxon_id = org.get("taxonId")
        if taxon_id is not None:
            result["taxon_id"] = str(taxon_id)
        lineage = org.get("lineage") or []
        if isinstance(lineage, list) and lineage:
            result["superkingdom"] = lineage[0]
            if len(lineage) >= 2:
                candidate = lineage[1]
                if (
                    isinstance(candidate, str)
                    and candidate.endswith("zoa")
                    and len(lineage) >= 3
                ):
                    result["phylum"] = lineage[2]
                else:
                    result["phylum"] = candidate
            result["genus"] = lineage[-1]
        sci_name = org.get("scientificName")
        if sci_name and not result["genus"]:
            result["genus"] = sci_name.split()[0]
        break
    return result


def iter_ids(csv_path: str) -> Iterable[str]:
    """Yield UniProt IDs from a CSV file with a ``UniProt_id`` column.

    Args:
        csv_path: Path to a CSV file containing a ``UniProt_id`` column.

    Yields:
        Each UniProt accession ID as a string.
    """
    with open(csv_path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if "UniProt_id" not in reader.fieldnames:
            raise ValueError("Input CSV must have a UniProt_id column")
        for row in reader:
            uid = row.get("UniProt_id")
            if uid:
                yield uid.strip()


def collect_info(uid: str, data_dir: str = "Data") -> Dict[str, str]:
    """Return names and organism data for ``uid``.

    Args:
        uid: UniProt accession identifier.
        data_dir: Directory containing ``<uid>.json`` files with UniProt data.

    Returns:
        A dictionary with keys ``UniProt_id``, ``names`` and organism
        taxonomy fields. Missing or invalid files leave fields empty.
    """
    json_path = os.path.join(data_dir, f"{uid}.json")
    result = {
        "UniProt_id": uid,
        "names": "",
        "genus": "",
        "superkingdom": "",
        "phylum": "",
        "taxon_id": "",
    }
    try:
        with open(json_path, "r", encoding="utf-8") as handle:
            data = json.load(handle)
        names = extract_names(data)
        org = extract_organism(data)
        result["names"] = "|".join(sorted(names))
        result.update(org)
    except FileNotFoundError:
        pass
    except json.JSONDecodeError:
        pass
    return result


def process(input_csv: str, output_csv: str, data_dir: str = "Data") -> None:
    """Read IDs from ``input_csv`` and write extracted data to ``output_csv``.

    Args:
        input_csv: Path to the CSV file listing UniProt IDs.
        output_csv: Destination path for the output CSV file.
        data_dir: Directory where JSON files for each ID are stored.

    Returns:
        None. The processed information is written to ``output_csv``.
    """
    rows = [collect_info(uid, data_dir) for uid in iter_ids(input_csv)]
    with open(output_csv, "w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "UniProt_id",
            "names",
            "genus",
            "superkingdom",
            "phylum",
            "taxon_id",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
