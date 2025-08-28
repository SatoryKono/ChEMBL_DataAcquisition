"""Data transformation utilities for classification data.

This module provides functions to read, transform, and write
classification tables.  It implements a pipeline similar to a Power
Query transformation, but using pandas.  The code is written to be
compatible with Python 3.12 and adheres to PEP 8 style guidelines.
"""

from __future__ import annotations

import logging
from itertools import chain
from typing import Iterable, List

import pandas as pd

logger = logging.getLogger(__name__)


def read_table(path: str, sep: str = ",", encoding: str = "utf-8") -> pd.DataFrame:
    """Read a CSV file into a :class:`pandas.DataFrame`.

    Parameters
    ----------
    path:
        Path to the CSV file.
    sep:
        Field delimiter, defaults to comma.
    encoding:
        File encoding, defaults to UTF-8.

    Returns
    -------
    pandas.DataFrame
        Parsed table.
    """

    logger.debug("Reading CSV file %s", path)
    return pd.read_csv(path, sep=sep, encoding=encoding)


def write_table(
    df: pd.DataFrame, path: str, sep: str = ",", encoding: str = "utf-8"
) -> None:
    """Write ``df`` to a CSV file.

    Parameters
    ----------
    df:
        Data to write.
    path:
        Output file path.
    sep:
        Field delimiter, defaults to comma.
    encoding:
        File encoding, defaults to UTF-8.
    """

    logger.debug("Writing CSV file %s", path)
    df.to_csv(path, sep=sep, encoding=encoding, index=False)


def to_lower_if_text(df: pd.DataFrame, columns: Iterable[str]) -> pd.DataFrame:
    """Convert specified text columns to lowercase.

    Non-string values are left unchanged.
    """

    result = df.copy()
    for col in columns:
        if col in result.columns:
            result[col] = result[col].apply(
                lambda x: x.lower() if isinstance(x, str) else x
            )
    return result


def split_and_clean(texts: Iterable[str]) -> List[str]:
    """Split pipe-delimited strings and return a list of unique tokens."""

    non_null = [t for t in texts if t is not None and str(t) != ""]
    split_lists = [str(t).split("|") for t in non_null]
    flat = list(chain.from_iterable(split_lists))
    trimmed = [t.strip() for t in flat]
    no_tag = [t.replace("synonyms=", "") for t in trimmed]
    non_empty = [t for t in no_tag if t]
    seen: set[str] = set()
    unique: List[str] = []
    for item in non_empty:
        if item not in seen:
            seen.add(item)
            unique.append(item)
    return unique


def validate_columns(df: pd.DataFrame, required: Iterable[str]) -> None:
    """Validate that required columns are present.

    Parameters
    ----------
    df:
        DataFrame to validate.
    required:
        Column names that must exist in ``df``.

    Raises
    ------
    ValueError
        If any required column is missing.
    """

    missing = [col for col in required if col not in df.columns]
    if missing:
        names = ", ".join(missing)
        raise ValueError(f"Missing required columns: {names}")


def get_table(source: pd.DataFrame) -> pd.DataFrame:
    """Replicate the initial Power Query transformations."""

    required = [
        "chembl_id",
        "target_id",
        "IUPHAR_family_id",
        "IUPHAR_type",
        "IUPHAR_class",
        "IUPHAR_subclass",
        "IUPHAR_chain",
        "full_id_path",
        "full_name_path",
        "gene",
        "component_description",
        "names_x",
        "chembl_alternative_name",
    ]
    validate_columns(source, required)

    lowered = to_lower_if_text(
        source,
        [
            "component_description",
            "pref_name",
            "gene",
            "chembl_alternative_name",
            "names_x",
            "cellular_component_x",
            "subcellular_location_x",
            "topology_x",
        ],
    )

    typed = lowered.copy()
    type_map = {
        "pref_name": "string",
        "chembl_id": "string",
        "component_description": "string",
        "component_id": "Int64",
        "relationship": "string",
        "gene": "string",
        "uniprot_id": "string",
        "chembl_alternative_name": "string",
        "names_x": "string",
        "target_id": "Int64",
        "IUPHAR_family_id": "Int64",
        "IUPHAR_type": "string",
        "IUPHAR_class": "string",
        "IUPHAR_subclass": "string",
        "IUPHAR_chain": "string",
        "full_id_path": "string",
        "full_name_path": "string",
    }
    for col, t in type_map.items():
        if col in typed.columns:
            typed[col] = typed[col].astype(t)

    typed["gene_name"] = typed["gene"].str.split("|").str[0]
    typed["_alternative_names"] = typed.apply(
        lambda r: split_and_clean(
            [
                r.get("gene"),
                r.get("component_description"),
                r.get("names_x"),
                r.get("chembl_alternative_name"),
            ]
        ),
        axis=1,
    )

    selected = typed[
        [
            "chembl_id",
            "target_id",
            "IUPHAR_family_id",
            "IUPHAR_type",
            "IUPHAR_class",
            "IUPHAR_subclass",
            "IUPHAR_chain",
            "full_id_path",
            "full_name_path",
            "gene_name",
            "_alternative_names",
        ]
    ]

    group_cols = [
        "chembl_id",
        "target_id",
        "IUPHAR_family_id",
        "IUPHAR_type",
        "IUPHAR_class",
        "IUPHAR_subclass",
        "IUPHAR_chain",
        "full_id_path",
        "full_name_path",
        "gene_name",
    ]

    def combine_unique(items: Iterable[List[str]]) -> List[str]:
        return list(dict.fromkeys(chain.from_iterable(items)))

    grouped = (
        selected.groupby(group_cols, dropna=False)["_alternative_names"].apply(combine_unique).reset_index()
    )
    expanded = grouped.explode("_alternative_names")
    filtered = expanded[expanded["_alternative_names"].notna() & (expanded["_alternative_names"] != "")]
    filtered = filtered.reset_index(drop=True)
    filtered["Index"] = range(len(filtered))
    return filtered


def build_base(source: pd.DataFrame) -> pd.DataFrame:
    """Helper replicating the ``BuildBase`` function."""

    done = source[
        (source["_alternative_names"].notna())
        & (source["_alternative_names"] != "")
        & (source["IUPHAR_type"].notna())
    ]
    filter_null_type = source[source["IUPHAR_type"].isna()]

    join_cols = [
        "chembl_id",
        "target_id",
        "IUPHAR_family_id",
        "IUPHAR_type",
        "IUPHAR_class",
        "IUPHAR_subclass",
        "IUPHAR_chain",
        "full_id_path",
        "full_name_path",
        "gene_name",
        "_alternative_names",
    ]
    joined = filter_null_type.merge(
        done[join_cols],
        on="_alternative_names",
        how="left",
        suffixes=("", "_right"),
    )

    joined = joined.drop(
        columns=[
            "target_id",
            "IUPHAR_family_id",
            "IUPHAR_type",
            "IUPHAR_class",
            "IUPHAR_subclass",
            "IUPHAR_chain",
            "full_id_path",
            "full_name_path",
        ]
    )

    joined = joined.rename(
        columns={
            "chembl_id_right": "chembl_id.1",
            "target_id_right": "target_id",
            "IUPHAR_family_id_right": "IUPHAR_family_id",
            "IUPHAR_type_right": "IUPHAR_type",
            "IUPHAR_class_right": "IUPHAR_class",
            "IUPHAR_subclass_right": "IUPHAR_subclass",
            "IUPHAR_chain_right": "IUPHAR_chain",
            "full_id_path_right": "full_id_path",
            "full_name_path_right": "full_name_path",
            "gene_name_right": "gene_name.1",
        }
    )

    joined["gene_match"] = joined["gene_name"] == joined["gene_name.1"]
    keep_with_type = joined.loc[joined["IUPHAR_type"].notna()].copy()
    keep_with_type["target_id"] = keep_with_type["target_id"].astype("string")
    dup_index = keep_with_type.copy()
    dup_index["Index_copy"] = dup_index["Index"]
    dup_index["IUPHAR_family_id_copy"] = dup_index["IUPHAR_family_id"]
    dup_index["IUPHAR_family_id_copy"] = dup_index["IUPHAR_family_id_copy"].astype("string")
    dup_index["Index_copy"] = dup_index["Index_copy"].astype("string")
    dup_index["Merged"] = dup_index["IUPHAR_family_id_copy"] + dup_index["Index_copy"]
    excluded = dup_index[~dup_index["Merged"].isin(["521825", "3271956"])]
    return excluded


def get_multyply(source: pd.DataFrame) -> pd.DataFrame:
    """Select whitelisted multiple matches."""

    base = build_base(source)
    group_cols = ["chembl_id", "gene_name", "Index"]
    counts = base.groupby(group_cols).size().reset_index(name="Count")
    multi = counts[counts["Count"] != 1]
    expanded = multi.merge(base, on=group_cols, how="left")
    keep = {"176838", "3231184", "49639", "5721806", "7381101", "765326", "8401880"}
    filtered = expanded[expanded["Merged"].isin(keep)]
    return filtered.drop(columns=["Count"])


def main_process(df: pd.DataFrame) -> pd.DataFrame:
    """Run the full transformation pipeline."""

    table1 = get_table(df)
    base = build_base(table1)
    group_cols = ["chembl_id", "gene_name", "Index"]
    counts = base.groupby(group_cols).size().reset_index(name="Count")
    singles = counts[counts["Count"] == 1]
    expanded = singles.merge(base, on=group_cols, how="left")
    combined = pd.concat([expanded, get_multyply(table1)], ignore_index=True)
    if "Count" in combined.columns:
        combined = combined.drop(columns=["Count"])
    return combined