"""Core transformation logic for working with IUPHAR data."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional

import logging

import pandas as pd
import requests

from .io_utils import load_families, load_targets

logger = logging.getLogger(__name__)


@dataclass
class IUPHARData:
    """Container for IUPHAR target and family data."""

    target_df: pd.DataFrame
    family_df: pd.DataFrame

    @classmethod
    def from_files(
        cls,
        target_path: str | Path,
        family_path: str | Path,
        *,
        encoding: str = "utf-8",
    ) -> "IUPHARData":
        """Load CSV files and return an :class:`IUPHARData` instance."""

        target_df = load_targets(target_path, encoding=encoding)
        family_df = load_families(family_path, encoding=encoding)
        return cls(target_df=target_df, family_df=family_df)

    # ------------------------------------------------------------------
    # Chain-related helpers
    # ------------------------------------------------------------------

    def family_chain(self, start_id: str) -> List[str]:
        """Return the family hierarchy starting from *start_id*.

        The chain is built by following ``parent_family_id`` links until a
        root family is reached.
        """

        chain: List[str] = []
        current = start_id
        while current:
            chain.append(current)
            parent_series = self.family_df.loc[
                self.family_df["family_id"] == current, "parent_family_id"
            ]
            if parent_series.empty:
                break
            parent = str(parent_series.iloc[0]) if pd.notna(parent_series.iloc[0]) else ""
            if not parent:
                break
            current = parent
        return chain

    def all_id(self, target_id: str) -> str:
        """Return the full family ID path for a target."""
        family_id = self.from_target_family_id(target_id)
        if not family_id:
            return ""
        chain = self.family_chain(family_id)
        return f"{target_id}#" + ">".join(chain)

    def all_name(self, target_id: str) -> str:
        """Return the full family name path for a target."""
        family_id = self.from_target_family_id(target_id)
        if not family_id:
            return ""
        chain = []
        for fam_id in self.family_chain(family_id):
            record = self.from_family_record(fam_id)
            if record is None:
                continue
            name = record.get("family_name")
            if name and name.lower() not in {"enzyme"}:
                chain.append(name)
        target_name = self.from_target_name(target_id)
        return f"{target_name}#" + ">".join(chain)

    # ------------------------------------------------------------------
    # Target ID mapping functions
    # ------------------------------------------------------------------

    def _select_target_ids(self, mask: pd.Series) -> List[str]:
        rows = self.target_df.loc[mask, "target_id"].dropna().astype(str)
        return list(rows.unique())

    def target_id_by_uniprot(self, uniprot_id: str) -> str:
        ids = self._select_target_ids(self.target_df["swissprot"].eq(uniprot_id))
        return ids[0] if ids else ""

    def target_id_by_hgnc_name(self, hgnc_name: str) -> str:
        if not hgnc_name:
            return ""
        ids = self._select_target_ids(self.target_df["HGNC_NAME"].eq(hgnc_name))
        return "|".join(ids) if ids else ""

    def target_id_by_hgnc_id(self, hgnc_id: str) -> str:
        ids = self._select_target_ids(self.target_df["HGNC_ID"].eq(hgnc_id))
        return "|".join(ids) if ids else ""

    def target_id_by_gene(self, gene_name: str) -> str:
        ids = self._select_target_ids(self.target_df["gene_name"].eq(gene_name))
        return "|".join(ids) if ids else ""

    def target_id_by_name(self, target_name: str) -> str:
        if not target_name:
            return ""
        mask = self.target_df["synonyms"].fillna("").str.contains(
            target_name, case=False, na=False
        )
        ids = self._select_target_ids(mask)
        return "|".join(ids) if ids else ""

    def target_ids_by_synonyms(self, synonyms: Iterable[str]) -> str:
        valid = [s for s in synonyms if s and len(s) > 3]
        ids: List[str] = []
        for syn in valid:
            mapped = self.target_id_by_name(syn)
            if mapped and "|" not in mapped:
                ids.append(mapped)
        unique = sorted(set(ids))
        return "|".join(unique) if unique else ""

    def target_id_from_row(self, row: pd.Series) -> str:
        uniprot = self.target_id_by_uniprot(row.get("uniprot_id", ""))
        hgnc_name = self.target_id_by_hgnc_name(row.get("HGNC_name", ""))
        hgnc_id = self.target_id_by_hgnc_id(row.get("HGNC_id", ""))
        name = ""
        synonyms = ""
        all_ids = [x for x in [uniprot, hgnc_name, hgnc_id, name, synonyms] if x]
        unique = sorted({i for part in all_ids for i in part.split("|") if i})
        if not unique:
            return "N/A"
        return "|".join(unique)

    def map_uniprot_file(
        self,
        input_path: str | Path,
        output_path: str | Path,
        *,
        encoding: str = "utf-8",
    ) -> pd.DataFrame:
        """Map UniProt accessions from *input_path* and write results.

        The input CSV must contain a column named ``uniprot_id``. For each
        accession the corresponding ``target_id`` is looked up along with the
        full family ID and name paths. The resulting table is written to
        ``output_path`` and also returned.

        Parameters
        ----------
        input_path:
            CSV file with a ``uniprot_id`` column.
        output_path:
            Where to write the resulting CSV file.
        encoding:
            Encoding used for both reading and writing. Defaults to UTF-8.

        Returns
        -------
        pandas.DataFrame
            DataFrame containing the mapping results.
        """

        df = pd.read_csv(input_path, dtype=str, encoding=encoding).fillna("")
        if "uniprot_id" not in df.columns:
            raise ValueError("Input file must contain 'uniprot_id' column")

        df["target_id"] = df["uniprot_id"].apply(self.target_id_by_uniprot)
        df["full_id_path"] = df["target_id"].apply(self.all_id)
        df["full_name_path"] = df["target_id"].apply(self.all_name)

        df.to_csv(output_path, index=False, encoding=encoding)
        logger.info("Wrote %d rows to %s", len(df), output_path)
        return df

    # ------------------------------------------------------------------
    # fromTargetID functions
    # ------------------------------------------------------------------

    def from_target_record(self, target_id: str) -> Optional[pd.Series]:
        rows = self.target_df[self.target_df["target_id"] == target_id]
        if rows.empty:
            return None
        return rows.iloc[0]

    def from_target_family_record(self, target_id: str) -> Optional[pd.Series]:
        mask = self.family_df["target_id"].fillna("").str.split("|").apply(
            lambda ids: target_id in ids
        )
        rows = self.family_df[mask]
        if rows.empty:
            return None
        return rows.iloc[0]

    def from_target_name(self, target_id: str) -> str:
        record = self.from_target_record(target_id)
        return str(record.get("target_name")) if record is not None else ""

    def from_target_type(self, target_id: str) -> str:
        record = self.from_target_record(target_id)
        if record is None:
            return ""
        value = record.get("type")
        return str(value) if pd.notna(value) else ""

    def from_target_synonyms(self, target_id: str) -> str:
        record = self.from_target_record(target_id)
        if record is None:
            return ""
        value = record.get("synonyms")
        return str(value) if pd.notna(value) else ""

    def from_target_family_id(self, target_id: str) -> str:
        record = self.from_target_record(target_id)
        if record is None:
            return ""
        value = record.get("family_id")
        return str(value) if pd.notna(value) else ""

    def from_target_parent_family(self, target_id: str) -> str:
        record = self.from_target_family_record(target_id)
        if record is None:
            return ""
        value = record.get("parent_family_id")
        return str(value) if pd.notna(value) else ""

    # ------------------------------------------------------------------
    # fromFamilyID functions
    # ------------------------------------------------------------------

    def from_family_type(self, family_id: str) -> str:
        record = self.from_family_record(family_id)
        if record is None:
            return ""
        value = record.get("type")
        return str(value) if pd.notna(value) else ""

    def from_family_record(self, family_id: str) -> Optional[pd.Series]:
        rows = self.family_df[self.family_df["family_id"] == family_id]
        if rows.empty:
            return None
        return rows.iloc[0]

    def from_family_parent(self, family_id: str) -> str:
        record = self.from_family_record(family_id)
        if record is None:
            return ""
        value = record.get("parent_family_id")
        return str(value) if pd.notna(value) else ""

    # ------------------------------------------------------------------
    # Web search
    # ------------------------------------------------------------------

    def websearch_gene_to_id(self, gene_name: str) -> dict:
        """Query the IUPHAR web API for *gene_name*.

        Returns a dictionary with the first hit or an empty dict on failure.
        """

        url = (
            "https://www.guidetopharmacology.org/services/targets/?geneSymbol="
            + gene_name
        )
        try:
            response = requests.get(url, timeout=10)
            response.raise_for_status()
            data = response.json()
            return data[0] if data else {}
        except requests.RequestException as exc:  # pragma: no cover - network errors
            logger.error("IUPHAR web request failed: %s", exc)
            return {}

    # ------------------------------------------------------------------
    # IUPHAR upload processing
    # ------------------------------------------------------------------

    def iuphar_upload(self) -> pd.DataFrame:
        """Reproduce the ``IUPHAR_upload`` transformation.

        The function downloads mapping files from the IUPHAR service and
        combines them with the local target and family tables.
        """

        uni_url = (
            "https://www.guidetopharmacology.org/DATA/GtP_to_UniProt_mapping.csv"
        )
        hgnc_url = (
            "https://www.guidetopharmacology.org/DATA/GtP_to_HGNC_mapping.csv"
        )

        uni_df = pd.read_csv(uni_url)
        hgnc_df = pd.read_csv(hgnc_url)
        hgnc_df = hgnc_df.rename(columns={"IUPHAR ID": "GtoPdb IUPHAR ID"})
        mapping = pd.merge(
            hgnc_df,
            uni_df,
            on="GtoPdb IUPHAR ID",
            how="outer",
            suffixes=("_hgnc", "_uni"),
        )
        mapping["Target id"] = mapping["GtoPdb IUPHAR ID"].fillna(mapping["IUPHAR ID"])

        joined = pd.merge(
            mapping,
            self.target_df,
            left_on="Target id",
            right_on="target_id",
            how="left",
        )
        joined = pd.merge(
            joined,
            self.family_df,
            left_on="family_id",
            right_on="family_id",
            how="left",
            suffixes=("", "_family"),
        )
        # Final cleanup of column names to match the M script
        result = joined.rename(
            columns={
                "UniProtKB ID": "chembl_swissprot",
                "HGNC ID": "chembl_HGNC_id",
                "IUPHAR Name": "chembl_name",
                "Family name": "family_name",
                "Family id": "family_id",
            }
        )
        return result
