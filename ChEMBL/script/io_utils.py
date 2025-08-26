from __future__ import annotations

from pathlib import Path
from typing import List

import pandas as pd


def read_ids(
    path: str | Path,
    column: str = "chembl_id",
    sep: str = ",",
    encoding: str = "utf8",
) -> List[str]:
    """Read identifier values from a CSV file.

    Parameters
    ----------
    path : str or Path
        Path to the CSV file.
    column : str, optional
        Name of the column containing identifiers. Defaults to ``"chembl_id"``.
    sep : str, optional
        Field delimiter, by default a comma.
    encoding : str, optional
        File encoding, by default ``"utf8"``.

    Returns
    -------
    list of str
        Identifier values in the order they appear. Empty strings and ``"#N/A"``
        markers are discarded.

    Raises
    ------
    ValueError
        If ``column`` is not present in the input file.
    """
    df = pd.read_csv(path, sep=sep, encoding=encoding, dtype=str)
    if column not in df.columns:
        raise ValueError(f"column '{column}' not found in {path}")

    ids = df[column].dropna().astype(str)
    return [i for i in ids if i and i != "#N/A"]