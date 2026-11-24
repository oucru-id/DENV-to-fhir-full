import os
from pathlib import Path
import pytest

REPO = Path(__file__).resolve().parent.parent


@pytest.fixture(autouse=True)
def fixed_timestamp(monkeypatch):
    monkeypatch.setenv("DENVFHIR_FIXED_TIMESTAMP", "2026-06-01T00:00:00+00:00")


@pytest.fixture
def repo_root():
    return REPO


@pytest.fixture
def panel_tsv():
    return REPO / "data" / "references" / "panel.tsv"


@pytest.fixture
def metadata_dir():
    return REPO / "tests" / "fixtures" / "metadata"
