"""Opt-in canonical lex-orbit site emission.

Default off: discovery order unchanged. Opt-in: ``site`` is the lex
representative; ``discovered_site`` holds pre-canonical indexes when they
differ; filters see discovery; unique-edit still collapses the orbit.
"""

from __future__ import annotations

import importlib.util

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from xenosite.forest.find_path import bfs, find_path
from xenosite.forest.graph_isomorphism import (
    all_site_pair_orbits_nauty,
    canonical_emitted_site,
    canonicalize_smarts_match,
    ensure_lexical_orbit_representatives,
    lexical_orbit_representatives,
    normalize_orbit_candidate,
)
from xenosite.forest.rdkitutil import MolFromSmiles
from xenosite.forest.rules import Dehydrogenation, Hydroxylation, QuinoneFormation

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("pynauty") is None,
    reason="pynauty required for lex-orbit representatives",
)


def _mol(smi: str):
    mol = MolFromSmiles(smi)
    assert mol is not None
    return mol


def _as_fs(site) -> frozenset[int]:
    if isinstance(site, frozenset):
        return site
    if isinstance(site, int):
        return frozenset({site})
    return frozenset(site)


def test_lexical_rep_stable_under_member_order():
    mol = _mol("c1ccccc1")
    families = all_site_pair_orbits_nauty(mol)
    ortho = families["atom_atom_unordered"][0]
    forward = lexical_orbit_representatives(
        [ortho], kind="atom_atom", ordered=False
    )
    reverse = lexical_orbit_representatives(
        [list(reversed(ortho))], kind="atom_atom", ordered=False
    )
    assert forward == reverse
    rep = min(
        normalize_orbit_candidate(m, kind="atom_atom", ordered=False)
        for m in ortho
    )
    for member in ortho:
        assert (
            canonical_emitted_site(
                member, forward, kind="atom_atom", ordered=False
            )
            == rep
        )


def test_canonicalize_smarts_match_maps_non_lex_to_rep():
    mol = _mol("c1ccccc1")
    families = all_site_pair_orbits_nauty(mol)
    ortho = families["atom_atom_unordered"][0]
    rep = min(ortho)
    non = next(m for m in ortho if m != rep)
    mapped = {1: non[0], 2: non[1]}
    out = canonicalize_smarts_match(
        mol, mapped, frozenset(non), unique_orbit="atom_atom"
    )
    assert out is not None
    emit_mapped, emit_site = out
    assert sorted(emit_site) == sorted(rep)
    assert sorted(emit_mapped.values()) == sorted(rep)


def test_canonicalize_smarts_match_singleton_atom_to_lex_rep():
    """One-atom sites remap onto the lex-smallest atom in the orbit."""

    mol = _mol("c1ccccc1")
    tables = ensure_lexical_orbit_representatives(mol)
    assert tables is not None
    # Benzene carbons share one orbit; lex rep is the min index.
    members = [a for a, rep in tables.atom.items() if rep == tables.atom[0]]
    assert len(members) == 6
    lex = min(members)
    non = max(members)
    assert non != lex
    mapped = {1: non}
    out = canonicalize_smarts_match(
        mol, mapped, frozenset({non}), unique_orbit="atom_atom"
    )
    assert out is not None
    emit_mapped, emit_site = out
    assert emit_site == frozenset({lex})
    assert emit_mapped[1] == lex


def test_default_off_no_discovered_site():
    mol = _mol("Oc1ccc(O)cc1")
    for _prod, info in Dehydrogenation().metabolize(mol):
        assert "discovered_site" not in info


def test_opt_in_emits_lex_site_with_discovered_escape_hatch():
    """Filter accepts only non-lex discoveries → emit site is lex + discovered_site."""

    mol = _mol("c1ccccc1")
    tables = ensure_lexical_orbit_representatives(mol)
    assert tables is not None

    filtered: list[frozenset[int]] = []

    def only_non_lex(_mol, site, _info):
        fs = _as_fs(site)
        filtered.append(fs)
        if len(fs) != 2:
            return True
        a, b = sorted(fs)
        try:
            rep = canonical_emitted_site(
                (a, b),
                tables.atom_atom_unordered,
                kind="atom_atom",
                ordered=False,
            )
        except KeyError:
            return True
        return (a, b) != rep

    products = list(
        QuinoneFormation().metabolize(
            mol,
            filter_sites=only_non_lex,
            canonical_emitted_sites=True,
        )
    )
    assert products
    assert filtered

    remapped = 0
    for prod, info in products:
        site_fs = _as_fs(info["site"])
        if "discovered_site" not in info:
            continue
        remapped += 1
        disc_fs = _as_fs(info["discovered_site"])
        assert disc_fs in filtered
        assert site_fs != disc_fs
        a, b = sorted(site_fs)
        assert canonical_emitted_site(
            (a, b),
            tables.atom_atom_unordered,
            kind="atom_atom",
            ordered=False,
        ) == (a, b)
        # Trace mirrors SiteInfo.
        tid = prod._forest["atom_trace"]["transforms"][-1]
        addition = prod._forest["atom_trace"]["additions"][tid]
        assert "discovered_site" in addition
        assert _as_fs(addition["site"]) == site_fs
        assert _as_fs(addition["discovered_site"]) == disc_fs

    assert remapped >= 1


def test_opt_in_does_not_double_emit_hq():
    mol = _mol("Oc1ccc(O)cc1")
    off = list(Dehydrogenation().metabolize(mol, canonical_emitted_sites=False))
    on = list(Dehydrogenation().metabolize(mol, canonical_emitted_sites=True))
    assert len(on) == len(off)
    assert {p.xf.csmi for p, _ in on} == {p.xf.csmi for p, _ in off}


def test_opt_in_singleton_emits_lex_atom_with_discovered_site():
    """Hydroxylation on benzene: filter non-lex → emit lex atom + discovered_site."""

    mol = _mol("c1ccccc1")
    tables = ensure_lexical_orbit_representatives(mol)
    assert tables is not None
    lex = min(a for a, rep in tables.atom.items() if rep == tables.atom[0])

    filtered: list[frozenset[int]] = []

    def only_non_lex(_mol, site, _info):
        fs = _as_fs(site)
        filtered.append(fs)
        if len(fs) != 1:
            return True
        atom = next(iter(fs))
        return atom != lex

    products = list(
        Hydroxylation().metabolize(
            mol,
            filter_sites=only_non_lex,
            canonical_emitted_sites=True,
        )
    )
    assert products
    assert filtered
    remapped = 0
    for prod, info in products:
        site_fs = _as_fs(info["site"])
        assert site_fs == frozenset({lex})
        if "discovered_site" not in info:
            continue
        remapped += 1
        disc_fs = _as_fs(info["discovered_site"])
        assert disc_fs != site_fs
        assert disc_fs in filtered
        tid = prod._forest["atom_trace"]["transforms"][-1]
        addition = prod._forest["atom_trace"]["additions"][tid]
        assert _as_fs(addition["site"]) == site_fs
        assert _as_fs(addition["discovered_site"]) == disc_fs
    assert remapped >= 1


def test_lex_reps_read_from_parent_after_product_cache_clear():
    """of_products / clear_structure wipe product cache; lex reps live on parent."""

    from xenosite.forest.rdkitutil import (
        copy_mol,
        restamp_product_forest_last_layer,
    )

    parent = _mol("c1ccccc1")
    parent = parent.xf.tracing._stamp()
    tables = ensure_lexical_orbit_representatives(parent)
    assert tables is not None
    assert parent._forest["cache"].get("lexical_orbit_representatives") is tables

    # Simulate a finished product whose structure cache was cleared.
    product = copy_mol(parent)
    product.xf.clear_structure()
    assert product._forest["cache"] == {}

    # Parent cache still intact.
    assert ensure_lexical_orbit_representatives(parent) is tables

    # Lookup via cleared child must pass parent= for the cache host.
    via_parent = ensure_lexical_orbit_representatives(product, parent=parent)
    assert via_parent is tables
    assert product._forest["cache"] == {}  # child still empty

    restamp_product_forest_last_layer(product, parent=parent)
    assert product._forest["cache"] == {}
    assert parent._forest["cache"].get("lexical_orbit_representatives") is tables


def test_canonical_emission_after_of_products_still_sets_discovered_site():
    """End-to-end: remap + restamp after of_products clear still works."""

    mol = _mol("c1ccccc1")
    tables = ensure_lexical_orbit_representatives(mol)
    assert tables is not None

    def only_non_lex(_mol, site, _info):
        fs = _as_fs(site)
        if len(fs) != 2:
            return True
        a, b = sorted(fs)
        try:
            rep = canonical_emitted_site(
                (a, b),
                tables.atom_atom_unordered,
                kind="atom_atom",
                ordered=False,
            )
        except KeyError:
            return True
        return (a, b) != rep

    products = list(
        QuinoneFormation().metabolize(
            mol,
            filter_sites=only_non_lex,
            canonical_emitted_sites=True,
        )
    )
    remapped = [info for _p, info in products if "discovered_site" in info]
    assert remapped
    for _p, info in products:
        if "discovered_site" not in info:
            continue
        # Product cache was cleared; addition still records the split.
        tid = _p._forest["atom_trace"]["transforms"][-1]
        addition = _p._forest["atom_trace"]["additions"][tid]
        assert "discovered_site" in addition
        assert _p._forest.get("cache", {}).get("lexical_orbit_representatives") is None


@given(canonical_emitted_sites=st.booleans())
@settings(max_examples=4, deadline=20_000, derandomize=True)
def test_bfs_forwards_canonical_emitted_sites(canonical_emitted_sites: bool):
    """``bfs`` splat: Hypothesis draws the flag; remapping only when True."""

    mol = _mol("c1ccccc1")
    tables = ensure_lexical_orbit_representatives(mol)
    assert tables is not None
    lex = min(a for a, rep in tables.atom.items() if rep == tables.atom[0])

    def only_non_lex(_mol, site, _info):
        fs = _as_fs(site)
        if len(fs) != 1:
            return True
        return next(iter(fs)) != lex

    products = list(
        bfs(
            mol,
            Hydroxylation(),
            depth=1,
            filter_sites=only_non_lex,
            canonical_emitted_sites=canonical_emitted_sites,
        )
    )
    assert products
    remapped = [(p, info) for p, info in products if "discovered_site" in info]
    if canonical_emitted_sites:
        assert remapped
        for _p, info in remapped:
            assert _as_fs(info["site"]) == frozenset({lex})
    else:
        assert remapped == []


@given(canonical_emitted_sites=st.booleans())
@settings(max_examples=4, deadline=20_000, derandomize=True)
def test_find_path_forwards_canonical_emitted_sites(canonical_emitted_sites: bool):
    """``find_path`` splat: Hypothesis draws the flag; either mode finds phenol."""

    hits = list(
        find_path(
            "c1ccccc1",
            "Oc1ccccc1",
            ruleset=Hydroxylation(),
            max_paths=1,
            max_nodes=50,
            canonical_emitted_sites=canonical_emitted_sites,
        )
    )
    assert hits
    assert hits[0].smiles == "Oc1ccccc1"

