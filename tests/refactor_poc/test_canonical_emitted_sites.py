"""Opt-in canonical lex-orbit site emission.

Default off: discovery order unchanged. Opt-in: ``site`` is the lex
representative; ``discovered_site`` holds pre-canonical indexes when they
differ; filters see discovery; unique-edit still collapses the orbit.
"""

from __future__ import annotations

import importlib.util

import pytest

from xenosite.refactor_poc.graph_isomorphism import (
    all_site_pair_orbits_nauty,
    canonical_emitted_site,
    canonicalize_smarts_match,
    ensure_lexical_orbit_representatives,
    lexical_orbit_representatives,
    normalize_orbit_candidate,
    set_canonical_emitted_sites,
)
from xenosite.refactor_poc.rdkitutil import MolFromSmiles
from xenosite.refactor_poc.rules import Dehydrogenation, QuinoneFormation

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


@pytest.fixture(autouse=True)
def _clear_canonical_override():
    set_canonical_emitted_sites(None)
    yield
    set_canonical_emitted_sites(None)


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
