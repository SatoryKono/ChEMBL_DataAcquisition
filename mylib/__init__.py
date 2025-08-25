"""Utility package for accessing ChEMBL data."""

from .chembl_client import (
    get_target,
    get_targets,
    get_assay,
    get_assays,
    extend_target,
)

__all__ = [
    "get_target",
    "get_targets",
    "get_assay",
    "get_assays",
    "extend_target",
]

