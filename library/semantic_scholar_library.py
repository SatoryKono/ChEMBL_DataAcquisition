"""Utilities for retrieving metadata from the Semantic Scholar API."""

from __future__ import annotations

import json
import logging
import time
from typing import Any, Dict, Tuple

import requests

logger = logging.getLogger(__name__)

TIMEOUT = 10


def _do_request(
    session: requests.Session,
    url: str,
    sleep: float,
    expect_json: bool = True,
    retries: int = 2,
    **kwargs: Any,
) -> Tuple[Dict[str, Any] | str | None, str]:
    """Perform a GET request with retry and error handling."""
    for attempt in range(retries + 1):
        if attempt:
            time.sleep(sleep * attempt)
        try:
            resp = session.get(url, timeout=TIMEOUT, **kwargs)
        except requests.RequestException as exc:  # pragma: no cover - network errors
            if attempt >= retries:
                return None, str(exc)
            continue

        if resp.status_code in (429, 500, 502, 503, 504):
            if attempt >= retries:
                return None, f"HTTP {resp.status_code}: {resp.text[:100]}"
            continue
        if resp.status_code != 200:
            return None, f"HTTP {resp.status_code}: {resp.text[:100]}"

        if expect_json:
            try:
                return resp.json(), ""
            except ValueError as exc:
                return None, f"Invalid JSON: {exc}"
        return resp.text, ""
    return None, "Request failed"


def fetch_semantic_scholar(
    session: requests.Session, pmid: str, sleep: float
) -> Dict[str, str]:
    """Fetch metadata for a PMID from Semantic Scholar."""
    fields = "publicationTypes,externalIds,paperId,venue"
    headers = {"Accept": "application/json"}
    url = f"https://api.semanticscholar.org/graph/v1/paper/PMID:{pmid}"
    data, error = _do_request(
        session,
        url,
        sleep * 5,
        headers=headers,
        params={"fields": fields},
    )
    if error or not isinstance(data, dict):
        return {
            "scholar.PMID": pmid,
            "scholar.Venue": "",
            "scholar.PublicationTypes": "",
            "scholar.SemanticScholarId": "",
            "scholar.ExternalIds": "",
            "scholar.DOI": "",
            "scholar.Error": error or "Invalid response",
        }
    pubtypes = data.get("publicationTypes") or []
    external_ids = data.get("externalIds") or {}
    doi = external_ids.get("DOI") or ""
    return {
        "scholar.PMID": pmid,
        "scholar.Venue": data.get("venue", ""),
        "scholar.PublicationTypes": "; ".join(pubtypes) if pubtypes else "",
        "scholar.SemanticScholarId": data.get("paperId", ""),
        "scholar.ExternalIds": json.dumps(external_ids, ensure_ascii=False),
        "scholar.DOI": doi,
        "scholar.Error": "",
    }

