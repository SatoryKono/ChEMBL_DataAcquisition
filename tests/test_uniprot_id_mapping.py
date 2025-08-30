import json
import json
from pathlib import Path

import pytest

import sys
sys.path.append(str(Path(__file__).resolve().parents[1]))

from library.uniprot_id_mapping import map_chembl_to_uniprot


class FakeResponse:
    def __init__(self, payload):
        self._payload = payload

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def read(self):
        return json.dumps(self._payload).encode()


def test_map_chembl_to_uniprot_success():
    responses = [
        {"jobId": "1"},
        {"jobStatus": "RUNNING"},
        {"results": [{"to": {"primaryAccession": "Q99558"}}]},
    ]
    def opener(url, data=None):
        return FakeResponse(responses.pop(0))

    accession = map_chembl_to_uniprot("CHEMBL123", opener=opener)
    assert accession == "Q99558"
