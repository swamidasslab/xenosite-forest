"""Pytest/Hypothesis shared configuration."""

from __future__ import annotations

import os
from pathlib import Path

from hypothesis import settings
from hypothesis.database import DirectoryBasedExampleDatabase

_DB_PATH = Path(__file__).resolve().parents[1] / ".hypothesis" / "examples"
_DB_PATH.mkdir(parents=True, exist_ok=True)

settings.register_profile(
    "default",
    database=DirectoryBasedExampleDatabase(str(_DB_PATH)),
    print_blob=True,
)
settings.register_profile(
    "ci",
    database=DirectoryBasedExampleDatabase(str(_DB_PATH)),
    print_blob=True,
    deadline=None,
)
settings.load_profile(os.getenv("HYPOTHESIS_PROFILE", "default"))
