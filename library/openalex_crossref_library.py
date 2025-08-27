"""Utilities for retrieving metadata from OpenAlex and CrossRef APIs."""

from __future__ import annotations

import logging
import time
from typing import Any, Dict, List, Tuple
from urllib.parse import quote

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


def combine(items: List[str]) -> str:
    """Combine non-empty items into a pipe-separated string."""
    return "|".join([x for x in items if x])


def fetch_openalex(session: requests.Session, pmid: str, sleep: float) -> Dict[str, str]:
    """Fetch metadata for a PMID from the OpenAlex API."""
    url = f"https://api.openalex.org/works/pmid:{pmid}"
    data, error = _do_request(session, url, sleep)
    if error or not isinstance(data, dict):
        return {
            "OpenAlex.PublicationTypes": "",
            "OpenAlex.TypeCrossref": "",
            "OpenAlex.Genre": "",
            "OpenAlex.Id": "",
            "OpenAlex.Venue": "",
            "OpenAlex.MeshDescriptors": "",
            "OpenAlex.MeshQualifiers": "",
            "OpenAlex.Error": error or "Invalid response",
        }
    mesh_entries = data.get("mesh") or []
    descriptors: List[str] = []
    qualifiers: List[str] = []
    for entry in mesh_entries:
        d = entry.get("descriptor_name")
        if d:
            descriptors.append(d)
        for q in entry.get("qualifiers") or []:
            qn = q.get("qualifier_name")
            if qn:
                qualifiers.append(qn)
    return {
        "OpenAlex.PublicationTypes": data.get("type", ""),
        "OpenAlex.TypeCrossref": data.get("type_crossref", ""),
        "OpenAlex.Genre": data.get("genre", ""),
        "OpenAlex.Id": data.get("id", ""),
        "OpenAlex.Venue": data.get("host_venue", {}).get("display_name", ""),
        "OpenAlex.MeshDescriptors": combine(descriptors),
        "OpenAlex.MeshQualifiers": combine(qualifiers),
        "OpenAlex.Error": "",
    }


def fetch_crossref(session: requests.Session, doi: str, sleep: float) -> Dict[str, str]:
    """Fetch metadata for a DOI from the CrossRef API."""
    if not doi:
        return {
            "crossref.Type": "",
            "crossref.Subtype": "",
            "crossref.Title": "",
            "crossref.Subtitle": "",
            "crossref.Subject": "",
            "crossref.Error": "Missing DOI",
        }

    url = f"https://api.crossref.org/works/{quote(doi, safe='')}"
    data, error = _do_request(session, url, sleep)
    if error or not isinstance(data, dict):
        return {
            "crossref.Type": "",
            "crossref.Subtype": "",
            "crossref.Title": "",
            "crossref.Subtitle": "",
            "crossref.Subject": "",
            "crossref.Error": error or "Invalid response",
        }
    message = data.get("message", {})
    title = message.get("title") or [""]
    subtitle = message.get("subtitle") or [""]
    subject = "; ".join(message.get("subject") or [])
    return {
        "crossref.Type": message.get("type", ""),
        "crossref.Subtype": message.get("subtype", ""),
        "crossref.Title": title[0] if title else "",
        "crossref.Subtitle": subtitle[0] if subtitle else "",
        "crossref.Subject": subject,
        "crossref.Error": "",
    }

