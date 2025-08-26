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

    def family_id_by_name(self, target_name: str) -> str:
        """Return family identifiers associated with *target_name*.

        The function searches the synonym column for matches to ``target_name``
        and collects the corresponding family IDs. Multiple IDs are joined by
        a pipe character.
        """

        if not target_name:
            return ""
        target_ids = self.target_id_by_name(target_name).split("|")
        family_ids = [self.from_target_family_id(tid) for tid in target_ids if tid]
        unique = sorted({fid for fid in family_ids if fid})
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
        """Map UniProt IDs to classification data and write a CSV output.

        The input CSV must contain a ``uniprot_id`` column. Each accession is
        translated to the corresponding IUPHAR target identifier and a full
        classification record. The resulting table includes the target ID,
        class, subclass, and family chain in addition to the full ID and name
        paths.

        Parameters
        ----------
        input_path:
            CSV file with a ``uniprot_id`` column.
        output_path:
            Destination for the enriched CSV file.
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

        classifier = IUPHARClassifier(self)

        df["target_id"] = df["uniprot_id"].apply(self.target_id_by_uniprot)

        def _classify(tid: str) -> pd.Series:
            if not tid:
                return pd.Series(
                    {
                        "IUPHAR_family_id": "",
                        "IUPHAR_type": "",
                        "IUPHAR_class": "",
                        "IUPHAR_subclass": "",
                        "IUPHAR_chain": "",
                    }
                )
            record = classifier.by_target_id(tid)
            return pd.Series(
                {
                    "IUPHAR_family_id": record.IUPHAR_family_id,
                    "IUPHAR_type": record.IUPHAR_type,
                    "IUPHAR_class": record.IUPHAR_class,
                    "IUPHAR_subclass": record.IUPHAR_subclass,
                    "IUPHAR_chain": ">".join(record.IUPHAR_tree),
                }
            )

        class_df = df["target_id"].apply(_classify)
        df = pd.concat([df, class_df], axis=1)

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


@dataclass
class ClassificationRecord:
    """Represents a classification result derived from IUPHAR data."""

    IUPHAR_target_id: str = "N/A"
    IUPHAR_family_id: str = "N/A"
    IUPHAR_class: str = "Other Protein Target"
    IUPHAR_subclass: str = "Other Protein Target"
    IUPHAR_tree: List[str] = None  # set in __post_init__
    IUPHAR_type: str = "Other Protein Target.Other Protein Target"
    IUPHAR_name: str = "N/A"
    IUPHAR_ecNumber: List[str] = None
    STATUS: str = "N/A"

    def __post_init__(self) -> None:
        if self.IUPHAR_tree is None:
            self.IUPHAR_tree = ["0864-1", "0864"]
        if self.IUPHAR_ecNumber is None:
            self.IUPHAR_ecNumber = []


class IUPHARClassifier:
    """Utility class replicating Power Query classification logic."""

    def __init__(self, data: IUPHARData):
        self.data = data

    # ------------------------------------------------------------------
    # Validation helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _is_valid_parameter(parameter: Optional[str]) -> bool:
        return bool(parameter) and parameter not in {"N/A", "Other Protein Target"}

    @staticmethod
    def _is_valid_list(values: Iterable[str]) -> bool:
        lst = [v for v in values if v and v != "N/A"]
        if not lst:
            return False
        if lst[0] == "0864-1" and len(lst) == 2:
            return False
        return True

    # ------------------------------------------------------------------
    # Lookups
    # ------------------------------------------------------------------

    def _family_to_type(self, family_id: str) -> str:
        if self._is_valid_parameter(family_id):
            return self.data.from_family_type(family_id) or "N/A"
        return "N/A"

    def _family_to_chain(self, family_id: str) -> List[str]:
        if self._is_valid_parameter(family_id):
            return self.data.family_chain(family_id)
        return []

    def _target_record(self, target_id: str) -> Optional[pd.Series]:
        if self._is_valid_parameter(target_id):
            return self.data.from_target_record(target_id)
        return None

    # ------------------------------------------------------------------
    # Type helpers
    # ------------------------------------------------------------------

    def _target_to_type(self, target_id: str) -> str:
        record = self._target_record(target_id)
        if record is None:
            return "Other Protein Target.Other Protein Target"
        family_id = str(record.get("family_id", ""))
        type1 = str(record.get("type", ""))
        type2 = self._family_to_type(family_id)
        if self._is_valid_parameter(type1) and not self._is_valid_parameter(type2):
            type_val = type1
        elif self._is_valid_parameter(type2) and not self._is_valid_parameter(type1):
            type_val = type2
        elif not (
            self._is_valid_parameter(type1) or self._is_valid_parameter(type2)
        ):
            type_val = "Other Protein Target.Other Protein Target"
        elif type1 == type2:
            type_val = type1
        elif self._is_valid_parameter(type1):
            type_val = type1
        else:
            type_val = "N/A.N/A"
        return type_val or "Other Protein Target.Other Protein Target"

    @staticmethod
    def _name_to_type(target_name: str) -> str:
        if not IUPHARClassifier._is_valid_parameter(target_name):
            return "Other Protein Target.Other Protein Target"
        name = target_name.lower()
        if "kinase" in name:
            return "Enzyme.Transferase"
        if "oxidase" in name or "reductase" in name:
            return "Enzyme.Oxidoreductase"
        if "hydrolase" in name or "protease" in name or "phosphatases" in name:
            return "Enzyme.Hydrolase"
        if "atpase" in name:
            return "Transporter.N/A"
        if "solute carrier" in name:
            return "Transporter.SLC superfamily of solute carrier"
        if "transport" in name:
            return "Transporter.N/A"
        if "channel" in name:
            return "Ion channel.N/A"
        if "hormone" in name:
            return "Receptor.Nuclear hormone receptor"
        return "Other Protein Target.Other Protein Target"

    # Mapping for EC-number derived types to classification chain
    _CHAIN_MAP = {
        "Enzyme.Oxidoreductase": ["0690-1", "0690"],
        "Enzyme.Transferase": ["0690-2", "0690"],
        "Enzyme.Multifunctional": ["0690-3", "0690"],
        "Enzyme.Hydrolase": ["0690-4", "0690"],
        "Enzyme.Isomerase": ["0690-5", "0690"],
        "Enzyme.Lyase": ["0690-6", "0690"],
        "Enzyme.Ligase": ["0690-6", "0690"],
        "Receptor.Catalytic receptor": ["0862", "0688"],
        "Receptor.G protein-coupled receptor": ["0694", "0688"],
        "Receptor.Nuclear hormone receptor": ["0095", "0688"],
        "Transporter.ATP-binding cassette transporter family": ["0136", "0691"],
        "Transporter.F-type and V-type ATPase": ["0137", "0691"],
        "Transporter.P-type ATPase": ["0138", "0691"],
        "Transporter.SLC superfamily of solute carrier": ["0863", "0691"],
        "Ion channel.Ligand-gated ion channel": ["0697", "0689"],
        "Ion channel.Other ion channel": ["0861", "0689"],
        "Ion channel.Voltage-gated ion channel": ["0696", "0689"],
    }

    @staticmethod
    def _ec_number_to_type(ec_numbers: List[str]) -> str:
        if not IUPHARClassifier._is_valid_list(ec_numbers):
            return ""
        prefixes = {
            n.split(".")[0] for n in ec_numbers if "." in n and n.split(".")[0]
        }
        if not prefixes:
            return ""
        if len(prefixes) > 1:
            return "Enzyme.Multifunctional"
        code = prefixes.pop()
        mapping = {
            "1": "Enzyme.Oxidoreductase",
            "2": "Enzyme.Transferase",
            "3": "Enzyme.Hydrolase",
            "4": "Enzyme.Lyase",
            "5": "Enzyme.Isomerase",
            "6": "Enzyme.Ligase",
            "7": "Enzyme.Translocase",
        }
        return mapping.get(code, "")

    @classmethod
    def _ec_number_to_chain(cls, ec_numbers: List[str]) -> List[str]:
        target_type = cls._ec_number_to_type(ec_numbers)
        return cls._CHAIN_MAP.get(target_type, ["0864-1", "0864"])

    # ------------------------------------------------------------------
    # Public classification methods
    # ------------------------------------------------------------------

    def set_record(
        self,
        iuphar_target_id: str,
        iuphar_family_id: str,
        iuphar_name: str,
        status: Optional[str] = None,
        ec_numbers: Optional[List[str]] = None,
    ) -> ClassificationRecord:
        target_id = (
            iuphar_target_id if self._is_valid_parameter(iuphar_target_id) else "N/A"
        )
        family_id = iuphar_family_id
        if not self._is_valid_parameter(family_id) and self._is_valid_parameter(target_id):
            family_id = self.data.from_target_family_id(target_id) or "N/A"
        name = iuphar_name if self._is_valid_parameter(iuphar_name) else "N/A"
        numbers = ec_numbers if ec_numbers and self._is_valid_list(ec_numbers) else []
        if status is None:
            if self._is_valid_parameter(target_id):
                status = "target_id"
            elif self._is_valid_parameter(family_id):
                status = "family_id"
            elif self._is_valid_list(numbers):
                status = "ec_number"
            else:
                status = "N/A"

        if self._is_valid_parameter(target_id):
            type_val = self._target_to_type(target_id)
        elif self._is_valid_parameter(family_id):
            type_val = self._family_to_type(family_id)
        else:
            type_val = "Other Protein Target.Other Protein Target"

        parts = type_val.split(".")
        cls_part = parts[0] if parts else "Other Protein Target"
        sub_part = parts[1] if len(parts) > 1 else "Other Protein Target"

        tree = (
            self._family_to_chain(family_id)
            if self._is_valid_parameter(family_id)
            else ["0864-1", "0864"]
        )

        return ClassificationRecord(
            IUPHAR_target_id=target_id,
            IUPHAR_family_id=family_id if self._is_valid_parameter(family_id) else "N/A",
            IUPHAR_class=cls_part,
            IUPHAR_subclass=sub_part,
            IUPHAR_tree=tree,
            IUPHAR_type=type_val,
            IUPHAR_name=name,
            IUPHAR_ecNumber=numbers,
            STATUS=status,
        )

    def by_target_id(
        self, iuphar_target_id: str, optional_name: Optional[str] = None
    ) -> ClassificationRecord:
        if not self._is_valid_parameter(iuphar_target_id) or "|" in iuphar_target_id:
            return ClassificationRecord()
        family_id = self.data.from_target_family_id(iuphar_target_id)
        return self.set_record(iuphar_target_id, family_id, optional_name or "")

    def by_uniprot_id(self, uniprot_id: str) -> ClassificationRecord:
        """Classify a UniProt accession.

        The function resolves the accession to a target identifier and then
        delegates to :meth:`by_target_id`. If the accession cannot be mapped,
        a default :class:`ClassificationRecord` is returned.
        """

        target_id = self.data.target_id_by_uniprot(uniprot_id)
        if not target_id:
            return ClassificationRecord()
        return self.by_target_id(target_id)

    def by_family_id(
        self, iuphar_family_id: str, optional_name: Optional[str] = None
    ) -> ClassificationRecord:
        if not self._is_valid_parameter(iuphar_family_id) or "|" in iuphar_family_id:
            return ClassificationRecord()
        return self.set_record("N/A", iuphar_family_id, optional_name or "")

    def by_ec_number(
        self, iuphar_ec_number: str, optional_name: Optional[str] = None
    ) -> ClassificationRecord:
        numbers = (
            iuphar_ec_number.split("|")
            if "." in iuphar_ec_number or "|" in iuphar_ec_number
            else []
        )
        if not self._is_valid_list(numbers):
            return ClassificationRecord()
        type_val = self._ec_number_to_type(numbers) or "Other Protein Target.Other Protein Target"
        parts = type_val.split(".")
        cls_part = parts[0] if parts else "Other Protein Target"
        sub_part = parts[1] if len(parts) > 1 else "Other Protein Target"
        tree = self._ec_number_to_chain(numbers)
        name = optional_name if self._is_valid_parameter(optional_name) else "N/A"
        return ClassificationRecord(
            IUPHAR_target_id="N/A",
            IUPHAR_family_id="N/A",
            IUPHAR_type=type_val,
            IUPHAR_class=cls_part,
            IUPHAR_subclass=sub_part,
            IUPHAR_tree=tree,
            IUPHAR_name=name,
            IUPHAR_ecNumber=numbers,
            STATUS="ec_number",
        )

    def by_name(self, iuphar_name: str) -> ClassificationRecord:
        if not self._is_valid_parameter(iuphar_name):
            return ClassificationRecord()
        target_id = self.data.target_id_by_name(iuphar_name)
        family_id = self.data.family_id_by_name(iuphar_name)
        if self._is_valid_parameter(target_id):
            type_val = self._target_to_type(target_id)
        elif self._is_valid_parameter(family_id):
            type_val = self._family_to_type(family_id)
        else:
            type_val = self._name_to_type(iuphar_name)
        tree = (
            self._family_to_chain(family_id)
            if self._is_valid_parameter(family_id)
            else ["864-1", "864"]
        )
        parts = type_val.split(".")
        cls_part = parts[0] if parts else "Other Protein Target"
        sub_part = parts[1] if len(parts) > 1 else "Other Protein Target"
        return ClassificationRecord(
            IUPHAR_target_id=target_id if self._is_valid_parameter(target_id) else "N/A",
            IUPHAR_family_id=family_id if self._is_valid_parameter(family_id) else "N/A",
            IUPHAR_type=type_val,
            IUPHAR_class=cls_part,
            IUPHAR_subclass=sub_part,
            IUPHAR_tree=tree,
            IUPHAR_name=iuphar_name,
            STATUS="name",
        )

    def get(
        self,
        iuphar_target_id: str,
        iuphar_family_id: str,
        iuphar_ec_number: str,
        iuphar_name: str,
    ) -> ClassificationRecord:
        target_rec = (
            self.by_target_id(iuphar_target_id, iuphar_name)
            if self._is_valid_parameter(iuphar_target_id)
            else ClassificationRecord()
        )
        family_rec = (
            self.by_family_id(iuphar_family_id, iuphar_name)
            if self._is_valid_parameter(iuphar_family_id)
            else ClassificationRecord()
        )
        ec_rec = self.by_ec_number(iuphar_ec_number, iuphar_name)
        name_rec = self.by_name(iuphar_name)
        if target_rec.STATUS in {"target_id", "family_id"}:
            return target_rec
        if family_rec.STATUS == "family_id":
            return family_rec
        if ec_rec.STATUS == "ec_number":
            return ec_rec
        if name_rec.STATUS == "name":
            return name_rec
        return ClassificationRecord()

    # ------------------------------------------------------------------
    # Activity and initialisation helpers
    # ------------------------------------------------------------------

    def merge_activity(
        self, input_df: pd.DataFrame, activity_df: pd.DataFrame
    ) -> pd.DataFrame:
        merged = input_df.merge(
            activity_df, on="task_uniprot_id", how="left", suffixes=("", "_act")
        )
        merged["activity.ec_number"] = (
            merged.get("activity.ec_number", "").replace("EC-not-assigned", "").fillna("")
        )
        merged["ec_number"] = (
            merged["activity.ec_number"].fillna("")
            + "|" + merged.get("ec_number", "").fillna("")
        )
        merged["ec_number"] = (
            merged["ec_number"].str.split("|").apply(
                lambda x: "|".join(sorted({i.strip() for i in x if i.strip()}))
            )
        )
        return merged.drop(columns=["activity.ec_number"], errors="ignore")

    def init_class(
        self,
        input_df: pd.DataFrame,
        family_table: pd.DataFrame,
        activity_df: pd.DataFrame,
        db_name: str,
    ) -> pd.DataFrame:
        df = self.merge_activity(input_df, activity_df)

        def classify(row: pd.Series) -> pd.Series:
            record = self.get(
                row.get("guidetopharmacology_id", ""),
                row.get("guidetopharmacology_family", ""),
                row.get("ec_number", ""),
                row.get("chembl_component_description", ""),
            )
            return pd.Series(
                {
                    "guidetopharmacology_id": record.IUPHAR_target_id,
                    "guidetopharmacology_family": record.IUPHAR_family_id,
                    "IUPHAR_type": record.IUPHAR_type,
                    "IUPHAR_class": record.IUPHAR_class,
                    "IUPHAR_subclass": record.IUPHAR_subclass,
                    "IUPHAR_chain": ">".join(record.IUPHAR_tree),
                    "IUPHAR_name": record.IUPHAR_name,
                }
            )

        classified = df.apply(classify, axis=1)
        return pd.concat([df, classified], axis=1)

    def by_reference(
        self, target_df: pd.DataFrame, db_df: pd.DataFrame, db_column: str
    ) -> pd.DataFrame:
        if "task_uniprot_id" not in target_df.columns:
            raise ValueError("target_df must contain 'task_uniprot_id'")
        if db_column not in db_df.columns:
            raise ValueError(f"db_df must contain '{db_column}' column")
        merged = target_df.merge(db_df, on="task_uniprot_id", how="left")
        merged["guidetopharmacology_family"] = merged["guidetopharmacology_family"].fillna(
            merged.get("guidetopharmacology_family_y")
        )
        merged["guidetopharmacology_type"] = merged["guidetopharmacology_type"].fillna(
            merged.get("guidetopharmacology_type_y")
        )
        return merged.drop(
            columns=[
                col
                for col in [
                    "guidetopharmacology_family_y",
                    "guidetopharmacology_type_y",
                ]
                if col in merged.columns
            ]
        )

    def get_database(
        self, family_table: pd.DataFrame, name_table: pd.DataFrame, db_column: str
    ) -> pd.DataFrame:
        valid = family_table[family_table[db_column].notna() & (family_table[db_column] != "N/A")]
        merged = valid.merge(name_table, on="task_uniprot_id", how="left")
        records = []
        for key, group in merged.groupby(db_column):
            family = (
                group["guidetopharmacology_family"].dropna().mode().iloc[0]
                if not group["guidetopharmacology_family"].dropna().empty
                else ""
            )
            type_series = group["guidetopharmacology_id"].dropna().apply(
                self.data.from_target_type
            )
            type_val = type_series.mode().iloc[0] if not type_series.empty else ""
            records.append(
                {
                    db_column: key,
                    "guidetopharmacology_family": family,
                    "guidetopharmacology_type": type_val,
                }
            )
        return pd.DataFrame(records)

