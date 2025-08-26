from __future__ import annotations

from pathlib import Path

import pandas as pd


def read_ids(
    path: str | Path,
    column: str = "chembl_id",
    sep: str = ",",
    encoding: str = "utf8",
) -> list[str]:
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
    try:
        df = pd.read_csv(path, sep=sep, encoding=encoding, dtype=str)
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"input file not found: {path}") from exc
    except pd.errors.EmptyDataError as exc:
        raise ValueError(f"no data in file: {path}") from exc
    except pd.errors.ParserError as exc:
        raise ValueError(f"malformed CSV in file: {path}: {exc}") from exc

    if column not in df.columns:
        raise ValueError(f"column '{column}' not found in {path}")

    ids = df[column].dropna().astype(str)
    return [i for i in ids if i and i != "#N/A"]
