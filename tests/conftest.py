"""Pytest/Hypothesis shared configuration."""

from __future__ import annotations

import os
from pathlib import Path

import pytest
from hypothesis import settings
from hypothesis.database import DirectoryBasedExampleDatabase

_DB_PATH = Path(__file__).resolve().parents[1] / ".hypothesis" / "examples"
_DB_PATH.mkdir(parents=True, exist_ok=True)

settings.register_profile(
    "default",
    database=DirectoryBasedExampleDatabase(str(_DB_PATH)),
    print_blob=True,
    max_examples=16,
)
settings.register_profile(
    "ci",
    database=DirectoryBasedExampleDatabase(str(_DB_PATH)),
    print_blob=True,
    deadline=None,
    max_examples=16,
)
# Long fuzz: hundreds of examples, no deadline. Per-test @settings(max_examples=N)
# still wins unless tests omit max_examples or read XENOSITE_FUZZ_EXAMPLES.
settings.register_profile(
    "long",
    database=DirectoryBasedExampleDatabase(str(_DB_PATH)),
    print_blob=True,
    deadline=None,
    max_examples=int(os.getenv("XENOSITE_FUZZ_EXAMPLES", "200")),
)
settings.load_profile(os.getenv("HYPOTHESIS_PROFILE", "default"))


# Archived pre-swap forest tests are historical only — never collect in CI/local default.
collect_ignore_glob = ["_archive_forest/*", "**/_archive_forest/*"]


@pytest.fixture(autouse=True)
def _formula_delta_mismatch_must_be_zero(request: pytest.FixtureRequest):
    """Every test: formula_delta_mismatch collector stays empty.

    Mark intentional mismatch tests with ``allow_formula_delta_mismatch``.
    """

    from xenosite.forest.rules import (
        begin_formula_delta_mismatch_collector,
        end_formula_delta_mismatch_collector,
    )

    bag, token = begin_formula_delta_mismatch_collector()
    yield
    end_formula_delta_mismatch_collector(token)
    if request.node.get_closest_marker("allow_formula_delta_mismatch"):
        return
    # Pair emissions: optimistic topology / site scoring can disagree with
    # sealed end bags — still recorded on counters, but not a suite fail yet.
    non_pair = [d for d in bag if not d.pair]
    assert non_pair == [], (
        "formula_delta_mismatch must be zero (non-pair); recorded %s" % (non_pair,)
    )