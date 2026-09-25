"""The records import, and pyright accepts that file."""

import subprocess
from pathlib import Path

from xenosite.forest.records import (
    Addition,
    AtomRef,
    AtomTrace,
    Forest,
    Formula,
    FragmentSplit,
    McsResult,
    PatternInfo,
    Structure,
)
from xenosite.forest.rules import ReactionRule

_ROOT = Path(__file__).resolve().parents[2]
_RECORDS = "src/xenosite/forest/records.py"


class _StubRule(ReactionRule):
    """A ReactionRule with no chemistry, so the test does not import a SMARTS rule."""


def test_records_are_tuples_and_dicts():
    formula: Formula = {"counts": {"C": 2, "H": 6}, "charge": 0}
    assert formula["counts"] == {"C": 2, "H": 6}
    assert formula["charge"] == 0
    assert isinstance(formula, dict)

    rule: ReactionRule = _StubRule(name="Hydroxylation")
    pattern: PatternInfo = {"name": "h", "possibilities": ()}
    addition = Addition(
        site=1,
        rules=(rule,),
        info={},
        effect={"adds": "O", "removes": "H"},
        name="Hydroxylation",
        phase1=None,
        depth=0,
        pattern=pattern,
    )
    assert addition.site == 1
    assert addition.name == "Hydroxylation"
    assert addition.depth == 0
    assert addition.pattern is pattern
    assert addition.pattern is not None
    assert addition.pattern.get("name") == "h"
    # PatternInfo is linked, not a chain entry.
    assert addition.pattern not in addition.rules

    # FutureSite top-level leaf is AtomRef; bare int stays a Site.
    future = AtomRef(0, "O")
    assert future.idx == 0
    assert future.element == "O"
    assert future.depth == 0
    nested: tuple[int | AtomRef, ...] = (0, future)
    assert nested[0] == 0
    assert nested[1] is future

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
        "dedup_smi": ["CC"],
        "depth": 0,
        "last_tag": 0,
        "next_transform": 1,
    }
    forest: Forest = {"cache": structure, "atom_trace": trace}
    assert forest["cache"]["formula"]["charge"] == 0
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
