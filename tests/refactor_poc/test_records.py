"""The records import, and pyright accepts that file."""

import subprocess
from pathlib import Path

from xenosite.refactor_poc.records import (
    Addition,
    AtomTrace,
    Forest,
    Formula,
    FragmentSplit,
    McsResult,
    Structure,
)
from xenosite.refactor_poc.rules import ReactionRule

_ROOT = Path(__file__).resolve().parents[2]
_RECORDS = "src/xenosite/refactor_poc/records.py"


class _StubRule(ReactionRule):
    """A ReactionRule with no chemistry, so the test does not import a SMARTS rule."""


def test_records_are_tuples_and_dicts():
    formula: Formula = {"counts": {"C": 2, "H": 6}, "charge": 0}
    assert formula["counts"] == {"C": 2, "H": 6}
    assert formula["charge"] == 0
    assert isinstance(formula, dict)

    rule: ReactionRule = _StubRule(name="Hydroxylation")
    addition = Addition(
        site=1,
        rules=(rule,),
        info={},
        effect={"adds": "O", "removes": "H"},
        name="Hydroxylation",
        phase1=None,
        depth=0,
    )
    assert addition.site == 1
    assert addition.name == "Hydroxylation"
    assert addition.depth == 0

    found = McsResult(embeddings=((0, 1),))
    assert found.embeddings == ((0, 1),)

    split = FragmentSplit(pieces=())
    assert split.pieces == ()

    structure: Structure = {"formula": formula, "csmi": "CC"}
    trace: AtomTrace = {
        "records": {"0": {"idx": [0], "depth": [0]}},
        "deletes": {},
        "transforms": [],
        "additions": {},
        "formula": formula,
        "delta_formula": {},
        "depth": 0,
        "last_tag": 0,
        "next_transform": 1,
    }
    forest: Forest = {"structure": structure, "atom_trace": trace}
    assert forest["structure"]["formula"]["charge"] == 0
    assert forest["atom_trace"]["next_transform"] == 1


def test_pyright_accepts_records():
    completed = subprocess.run(
        ["pyright", _RECORDS],
        cwd=_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
