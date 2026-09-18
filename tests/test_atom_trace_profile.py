"""Opt-in AtomTracker.tags profiling bench (~1–3s).

Run::

    uv run pytest tests/test_atom_trace_profile.py -q -p no:xdist

Not collected under normal ``pytest -n auto`` unless this file is targeted
(or ``-m profile``). No hard wall-time gate — numbers go in LOG.md.
"""

from __future__ import annotations

import cProfile
import pstats
import time

import pytest
from rdkit import Chem

from xenosite.forest.phaseone import PhaseOneRS

pytestmark = pytest.mark.profile


def _expand(mols):
    nxt = []
    for m in mols:
        for _site, prods in PhaseOneRS.metabolize(Chem.Mol(m)):
            nxt.extend(p for p in prods if p is not None)
    return nxt


def test_benzene_phaseone_chained_metabolize_profile(capsys):
    """Benzene depth-1 then expand those products (depth-2 frontier)."""
    parents = _expand([Chem.MolFromSmiles("c1ccccc1")])
    assert parents, "depth-1 should emit products"

    pr = cProfile.Profile()
    t0 = time.perf_counter()
    pr.enable()
    children = _expand(parents)
    pr.disable()
    wall = time.perf_counter() - t0

    st = pstats.Stats(pr)
    watch = {
        "literal_eval": {"tt": 0.0, "ct": 0.0},
        "deepcopy": {"tt": 0.0, "ct": 0.0},
        "tags": {"tt": 0.0, "ct": 0.0},
    }
    for (_fn, _line, name), (_cc, _nc, tt, ct, _caller) in st.stats.items():
        if name in watch:
            watch[name]["tt"] += tt
            watch[name]["ct"] = max(watch[name]["ct"], ct)

    line = (
        "atom_trace_profile wall=%.3fs products=%d parents=%d "
        "literal_eval_cum=%.3fs deepcopy_cum=%.3fs tags_cum=%.3fs"
        % (
            wall,
            len(children),
            len(parents),
            watch["literal_eval"]["ct"],
            watch["deepcopy"]["ct"],
            watch["tags"]["ct"],
        )
    )
    # Soft bounds: should finish in a few seconds and do real work.
    assert len(children) > 50
    assert wall < 30.0
    with capsys.disabled():
        print(line)
