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


# Legacy freeze tests live under tests/forest/legacy (collected via testpaths).
collect_ignore_glob: list[str] = []


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--parity-focused",
        action="store_true",
        default=False,
        help=(
            "Shrink native↔Rust rule_parity_fuzz to CoverIntent cases only "
            "(same as XENOSITE_PARITY_FULL=0). Default remains full cartesian."
        ),
    )


def pytest_configure(config: pytest.Config) -> None:
    # Set before collection so rule_parity_corpus.parity_param_cases sees it.
    if config.getoption("--parity-focused"):
        os.environ["XENOSITE_PARITY_FULL"] = "0"


@pytest.fixture(autouse=True)
def _formula_delta_mismatch_must_be_zero(request: pytest.FixtureRequest):
    """Every test: formula_delta_mismatch collector stays empty.

    Mark intentional mismatch tests with ``allow_formula_delta_mismatch``.
    Skipped when the RDKit native engine is not installed (rust-only jobs).
    """

    try:
        from xenosite.forest.native.rules import (
            begin_formula_delta_mismatch_collector,
            end_formula_delta_mismatch_collector,
        )
    except ImportError:
        yield
        return

    bag, token = begin_formula_delta_mismatch_collector()
    yield
    end_formula_delta_mismatch_collector(token)
    if request.node.get_closest_marker("allow_formula_delta_mismatch"):
        return
    # Pair emissions: still recorded (`pair=True`); suite fail stays off until
    # sealed end bags are fully trusted under one-placement site scoring.
    non_pair = [d for d in bag if not d.pair]
    assert non_pair == [], (
        "formula_delta_mismatch must be zero (non-pair); recorded %s" % (non_pair,)
    )
