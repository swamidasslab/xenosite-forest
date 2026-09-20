"""Random phase-I walks, then ``find_path``.

Invariants, not a recipe replay: a hit's product is the target and is one
component; a miss is an empty result, not an exception. Terminal conjugates
are not expanded.

``canonical_emitted_sites`` is drawn with Hypothesis ``st.booleans()`` and
passed as kwargs on metabolize / find_path (no env / custom flip helpers).
"""

from __future__ import annotations

import os as _os

from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st
from rdkit import Chem

from xenosite.forest.find_path import PathCounters, find_path
from xenosite.forest.rules import (
    Acetylation,
    Dehydrogenation,
    Hydroxylation,
    QuinoneFormation,
)
from xenosite.forest.rulesets import RuleSet


def _fuzz_examples(default: int) -> int:
    raw = _os.environ.get("XENOSITE_FUZZ_EXAMPLES")
    if raw:
        return int(raw)
    return default

_CORPUS = ("CCO", "CC", "C=C", "c1ccccc1", "Oc1ccccc1", "CCN", "COc1ccccc1")
_MAX_HEAVY = 16


def _rules():
    return [Hydroxylation(), Dehydrogenation()]


def _ruleset():
    return RuleSet(_rules(), name="FuzzPhase1")


def _collect(mol, rules, seen: set[str], **site_kw):
    out = []
    for rule in rules:
        for products, _info in rule.metabolize(mol, **site_kw):
            for product in products:
                smi = product.xf.csmi
                if not smi or "." in smi or smi in seen:
                    continue
                out.append((rule, product, smi))
    return out


def _walk(draw, start: str, rules, n_steps: int, *, canonical_emitted_sites: bool):
    mol = Chem.MolFromSmiles(start)
    assume(mol is not None)
    assume(mol.GetNumHeavyAtoms() <= _MAX_HEAVY)
    path = [mol.xf.csmi]
    current = mol
    site_kw = {"canonical_emitted_sites": canonical_emitted_sites}
    for _step in range(n_steps):
        candidates = _collect(current, rules, set(path), **site_kw)
        assume(candidates)
        candidates.sort(key=lambda item: item[1].GetNumHeavyAtoms())
        pool = candidates[:8]
        _rule, product, smi = draw(st.sampled_from(pool))
        path.append(smi)
        current = product
        assume(current.GetNumHeavyAtoms() <= _MAX_HEAVY + 4)
    return path[0], path[-1]


def _search(
    reactant: str, target: str, *, max_nodes: int, canonical_emitted_sites: bool
):
    counters = PathCounters()
    hits = list(
        find_path(
            reactant,
            target,
            ruleset=_ruleset(),
            counters=counters,
            max_nodes=max_nodes,
            max_paths=1,
            canonical_emitted_sites=canonical_emitted_sites,
        )
    )
    return hits, counters


def _assert_hit_or_miss(hits, target: str):
    if not hits:
        return
    assert len(hits) == 1
    assert "." not in hits[0].smiles
    assert hits[0].smiles == target


@given(
    start=st.sampled_from(_CORPUS),
    data=st.data(),
    canonical_emitted_sites=st.booleans(),
)
@settings(
    max_examples=_fuzz_examples(6),
    deadline=20_000,
    derandomize=True,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.filter_too_much],
)
def test_fuzz_depth1_hit_matches_target(
    start: str, data, canonical_emitted_sites: bool
):
    reactant, target = _walk(
        data.draw,
        start,
        _rules(),
        1,
        canonical_emitted_sites=canonical_emitted_sites,
    )
    assume(reactant != target)
    hits, counters = _search(
        reactant,
        target,
        max_nodes=80,
        canonical_emitted_sites=canonical_emitted_sites,
    )
    assert counters.nodes <= 80
    # A one-step phase-I product is in the ruleset. A miss is a search bug.
    assert hits, (reactant, target, counters.nodes)
    _assert_hit_or_miss(hits, target)


@given(
    start=st.sampled_from(_CORPUS),
    data=st.data(),
    canonical_emitted_sites=st.booleans(),
)
@settings(
    max_examples=_fuzz_examples(4),
    deadline=20_000,
    derandomize=True,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.filter_too_much],
)
def test_fuzz_depth2_hit_or_honest_miss(
    start: str, data, canonical_emitted_sites: bool
):
    reactant, target = _walk(
        data.draw,
        start,
        _rules(),
        2,
        canonical_emitted_sites=canonical_emitted_sites,
    )
    assume(reactant != target)
    hits, counters = _search(
        reactant,
        target,
        max_nodes=120,
        canonical_emitted_sites=canonical_emitted_sites,
    )
    assert counters.nodes <= 120
    _assert_hit_or_miss(hits, target)


@given(canonical_emitted_sites=st.booleans())
@settings(max_examples=_fuzz_examples(4), deadline=10_000, derandomize=True)
def test_walk_candidates_skip_the_current_molecule(canonical_emitted_sites: bool):
    mol = Chem.MolFromSmiles("CCO")
    seen = {mol.xf.csmi}
    for _rule, _product, smi in _collect(
        mol,
        _rules(),
        seen,
        canonical_emitted_sites=canonical_emitted_sites,
    ):
        assert smi not in seen
        assert "." not in smi


@given(canonical_emitted_sites=st.booleans())
@settings(max_examples=_fuzz_examples(4), deadline=10_000, derandomize=True)
def test_max_nodes_stops(canonical_emitted_sites: bool):
    counters = PathCounters()
    hits = list(
        find_path(
            "c1ccccc1",
            "O=C1C=CC(=O)C=C1",
            ruleset=_ruleset(),
            counters=counters,
            max_nodes=3,
            max_paths=1,
            canonical_emitted_sites=canonical_emitted_sites,
        )
    )
    assert counters.nodes <= 3
    assert hits == [] or "." not in hits[0].smiles


def test_terminal_conjugate_is_not_expanded():
    phenol = Chem.MolFromSmiles("Oc1ccccc1")
    products, _info = next(Acetylation().metabolize(phenol))
    acetyl = products[0]
    assert acetyl.xf.is_terminal
    assert list(Hydroxylation().metabolize(acetyl)) == []
    assert list(QuinoneFormation().metabolize(acetyl)) == []
    counters = PathCounters()
    hits = list(
        find_path(
            acetyl,
            "O=C1C=CC(=O)C=C1",
            ruleset=RuleSet([Hydroxylation(), Dehydrogenation(), QuinoneFormation()]),
            counters=counters,
            max_nodes=5,
            max_paths=1,
        )
    )
    assert hits == []
    assert counters.nodes == 1
