"""Functions to map ChEMBL target IDs to UniProt accessions using the UniProt ID Mapping API."""

from __future__ import annotations

import json
import logging
import time
import urllib.parse
import urllib.request
from urllib.error import HTTPError
from typing import Optional, Callable, Any

API_BASE = "https://rest.uniprot.org/idmapping"

logger = logging.getLogger(__name__)


def map_chembl_to_uniprot(
    chembl_target_id: str,
    poll_interval: float = 0.5,
    timeout: float = 300.0,
    opener: Optional[Callable[..., Any]] = None,
) -> str:
    """Map a ChEMBL target identifier to a UniProt accession.

    Parameters
    ----------
    chembl_target_id:
        ChEMBL target identifier (e.g., ``"CHEMBL204"``).
    poll_interval:
        Seconds to wait between polling the UniProt API for job completion.
    timeout:
        Maximum number of seconds to wait for the mapping job to finish.
    opener:
        Optional callable with the same signature as :func:`urllib.request.urlopen`
        used to perform HTTP requests. Primarily intended for testing.

    Returns
    -------
    str
        UniProt accession corresponding to ``chembl_target_id``.

    Raises
    ------
    ValueError
        If the API reports failure, no UniProt ID is found, or a UniProt API
        request returns an HTTP error.
    TimeoutError
        If the mapping job does not complete within ``timeout`` seconds.
    URLError
        If a network-related error occurs.
    """

    if opener is None:
        opener = urllib.request.urlopen

    def _open_json(url: str, data: bytes | None = None) -> Any:
        """Open ``url`` and parse the JSON response."""
        try:
            with opener(url, data=data) as response:
                return json.load(response)
        except HTTPError as exc:  # pragma: no cover - network failure simulation
            body = ""
            if exc.fp is not None:
                try:
                    body = exc.fp.read().decode()
                except Exception:  # pragma: no cover - fallback if decode fails
                    body = ""
            raise ValueError(
                f"UniProt API request to {url} failed with status {exc.code}: {body or exc.reason}"
            ) from exc

    data = urllib.parse.urlencode(
        {"from": "ChEMBL", "to": "UniProtKB", "ids": chembl_target_id}
    ).encode()
    logger.debug("Submitting ID mapping job for %s", chembl_target_id)
    run_data = _open_json(f"{API_BASE}/run", data=data)
    job_id = run_data.get("jobId")
    if not job_id:
        raise ValueError("UniProt ID Mapping API did not return a job ID")

    status_url = f"{API_BASE}/status/{job_id}"
    start = time.time()
    result_data = {}
    while True:
        status_data = _open_json(status_url)
        if "results" in status_data:
            result_data = status_data
            status = "FINISHED"
        else:
            status = status_data.get("jobStatus") or status_data.get("status")
        logger.debug("Job %s status: %s", job_id, status)
        if status == "FINISHED":
            break
        if status == "FAILED":
            raise ValueError("UniProt ID mapping job failed")
        if time.time() - start > timeout:
            raise TimeoutError("UniProt ID mapping job timed out")
        time.sleep(poll_interval)

    if not result_data:
        result_url = f"{API_BASE}/uniprotkb/results/{job_id}?format=json"
        result_data = _open_json(result_url)

    results = result_data.get("results", [])
    if not results:
        raise ValueError(f"No UniProt ID found for {chembl_target_id}")

    first = results[0]
    to = first.get("to", {})
    accession = to.get("primaryAccession")
    if not accession:
        raise ValueError("Unexpected response format from UniProt ID mapping API")

    return accession


__all__ = ["map_chembl_to_uniprot"]
