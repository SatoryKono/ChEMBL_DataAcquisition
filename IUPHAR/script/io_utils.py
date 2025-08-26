"""Utility functions for reading IUPHAR CSV files."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd


EXPECTED_TARGET_COLUMNS: tuple[str, ...] = (
    "target_id",
    "swissprot",
    "hgnc_name",
    "hgnc_id",
    "gene_name",
    "synonyms",
    "family_id",
    "target_name",
)

EXPECTED_FAMILY_COLUMNS: tuple[str, ...] = (
    "family_id",
    "family_name",
    "parent_family_id",
    "target_id",
    "type",
)


def _validate_columns(df: pd.DataFrame, expected: Iterable[str]) -> None:
    """Validate that *df* contains the *expected* columns.

    Parameters
    ----------
    df:
        DataFrame to validate.
    expected:
        Iterable of expected column names.

    Raises
    ------
    ValueError
        If any expected column is missing.
    """
    missing = [c for c in expected if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns: {', '.join(missing)}")


def load_targets(path: str | Path, *, encoding: str = "utf-8") -> pd.DataFrame:
    """Load the ``_IUPHAR_target.csv`` file.

    The ``target_id`` and ``family_id`` columns in the official files are
    zero-padded strings. Reading them as numeric values would drop the
    padding and break lookups.  To avoid this, all columns are read as
    strings and missing values are replaced with empty strings.

    Parameters
    ----------
    path:
        Path to the target CSV file.
    encoding:
        File encoding. Defaults to UTF-8.

    Returns
    -------
    pandas.DataFrame
        Loaded target data with string identifiers.
    """
    df = pd.read_csv(path, encoding=encoding, dtype=str).fillna("")
    # Normalise legacy column names to the expected lowercase form
    df = df.rename(
        columns={
            "HGNC_NAME": "hgnc_name",
            "HGNC_name": "hgnc_name",
            "HGNC_ID": "hgnc_id",
            "HGNC_id": "hgnc_id",
        }
    )
    _validate_columns(df, EXPECTED_TARGET_COLUMNS)
    return df


def load_families(path: str | Path, *, encoding: str = "utf-8") -> pd.DataFrame:
    """Load the ``_IUPHAR_family.csv`` file.

    Similar to :func:`load_targets`, family identifiers are stored as
    zero-padded strings.  The file is therefore read with ``dtype=str`` and
    missing values are normalised to empty strings to retain formatting.

    Parameters
    ----------
    path:
        Path to the family CSV file.
    encoding:
        File encoding. Defaults to UTF-8.

    Returns
    -------
    pandas.DataFrame
        Loaded family data with string identifiers.
    """
    df = pd.read_csv(path, encoding=encoding, dtype=str).fillna("")
    _validate_columns(df, EXPECTED_FAMILY_COLUMNS)
    return df
