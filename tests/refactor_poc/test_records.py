"""The records import, and pyright accepts that file."""

import subprocess
from pathlib import Path

from xenosite.refactor_poc.records import (
    Addition,
    Forest,
    Formula,
    FragmentSplit,
    McsResult,
    Structure,
)

_ROOT = Path(__file__).resolve().parents[2]
_RECORDS = "src/xenosite/refactor_poc/records.py"


def test_records_are_tuples_and_dicts():
    formula: Formula = {"counts": {"C": 2, "H": 6}, "charge": 0}
    assert formula["counts"] == {"C": 2, "H": 6}
    assert formula["charge"] == 0
    assert isinstance(formula, dict)

    addition = Addition(
        site=(1,),
        rules=("Hydroxylation",),
        info={},
        effect={},
        name="Hydroxylation",
        phase1=None,
        depth=0,
    )
    assert addition.site == (1,)
    assert addition.name == "Hydroxylation"
    assert addition.depth == 0

    found = McsResult(embeddings=((0, 1),))
    assert found.embeddings == ((0, 1),)

    split = FragmentSplit(pieces=())
    assert split.pieces == ()

    structure: Structure = {"formula": formula, "csmi": "CC"}
    forest: Forest = {"structure": structure, "atom_trace": {}}
    assert forest["structure"]["formula"]["charge"] == 0
    assert isinstance(forest["atom_trace"], dict)


def test_pyright_accepts_records():
    completed = subprocess.run(
        ["pyright", _RECORDS],
        cwd=_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
