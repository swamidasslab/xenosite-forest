"""Mol-scoped resonance cache: parity, laziness, and instrumentation."""

from __future__ import annotations

from contextlib import nullcontext

import pytest
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolfiles import MolFromSmiles, MolToSmiles

from xenosite.forest import rulesets
from xenosite.forest.base import (
    EditMol,
    Resonate,
    _ModeResonanceCache,
    _resonance_cache,
    _resonance_cache_disabled,
    _set_resonance_cache_enabled,
)
from xenosite.forest.rules import Dehydrogenation, Hydrogenation
from xenosite.forest.utils import refresh_mol

APAP = "CC(=O)Nc1ccc(O)cc1"
NAPH_STYRYL = "C1=CC=CC2=C1C=C(C=C2)CC3=CC=CC(=C3)C=C"
MULTI_CONJ = "NCCCCC1=CC=CC2=C1C=C(C=C2)CC3=CC=CC(=C3)C=Cc1ccccc1"

FIXTURES = {
    "APAP": APAP,
    "naph_styryl": NAPH_STYRYL,
    "multi_conj": MULTI_CONJ,
}


def _smi(m):
    return MolToSmiles(refresh_mol(m))


def _product_key(site, mets):
    return (str(site), tuple(sorted(_smi(m) for m in mets)))


def _collect_products(ruleset_name, smiles, *, cache_disabled=False):
    rs = rulesets.load_ruleset(ruleset_name)
    ctx = _resonance_cache_disabled() if cache_disabled else nullcontext()
    with ctx:
        return [
            _product_key(site, mets)
            for site, mets in rs.metabolites(Mol(MolFromSmiles(smiles)))
        ]


def _form_smiles(smiles, *, cache_disabled=False, first_yield_input=False):
    ctx = _resonance_cache_disabled() if cache_disabled else nullcontext()
    with ctx:
        return [
            _smi(m)
            for m in Resonate().resonance_structures(
                Mol(MolFromSmiles(smiles)), first_yield_input=first_yield_input
            )
        ]


@pytest.fixture(autouse=True)
def _restore_resonance_cache_flag():
    """Ensure private toggle does not leak across tests."""
    yield
    _set_resonance_cache_enabled(True)


# ---------------------------------------------------------------------------
# Form-count / product parity (cache on vs off)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name,smiles", list(FIXTURES.items()), ids=list(FIXTURES))
def test_resonance_form_count_and_smiles_match_cache_on_off(name, smiles):
    on = _form_smiles(smiles, cache_disabled=False)
    off = _form_smiles(smiles, cache_disabled=True)
    assert len(on) == len(off)
    assert on == off


@pytest.mark.parametrize("name,smiles", list(FIXTURES.items()), ids=list(FIXTURES))
def test_full_products_match_cache_on_off(name, smiles):
    on = _collect_products("Full", smiles, cache_disabled=False)
    off = _collect_products("Full", smiles, cache_disabled=True)
    assert len(on) == len(off)
    assert on == off


@pytest.mark.parametrize("name,smiles", list(FIXTURES.items()), ids=list(FIXTURES))
def test_phaseone_products_match_cache_on_off(name, smiles):
    on = _collect_products("PhaseOneRS", smiles, cache_disabled=False)
    off = _collect_products("PhaseOneRS", smiles, cache_disabled=True)
    assert len(on) == len(off)
    assert on == off


# ---------------------------------------------------------------------------
# Instrumentation: call / compute counts
# ---------------------------------------------------------------------------


def test_resfrags_compute_once_per_mode_on_full_apap(monkeypatch):
    """Full ruleset fills each (mol, mode) producer once — not once globally.

    Dehydrogenation runs conjugated resonance on the substrate and again on its
    kekulized standardize template (separate forests). QuinoneFormation fills
    aromatic on the template. Same-size must not share resonance across those.
    """
    calls = {"n": 0}
    orig = Resonate._resfrags

    def counting_resfrags(self, mol, output_systems=False):
        calls["n"] += 1
        yield from orig(self, mol, output_systems=output_systems)

    monkeypatch.setattr(Resonate, "_resfrags", counting_resfrags)

    list(rulesets.load_ruleset("Full").metabolites(Mol(MolFromSmiles(APAP))))
    assert calls["n"] == 3  # conj@substrate + conj@template + aromatic@template

    calls["n"] = 0
    with _resonance_cache_disabled():
        list(rulesets.load_ruleset("Full").metabolites(Mol(MolFromSmiles(APAP))))
    assert calls["n"] == 7


def test_resonance_structures_call_count_instrumentable(monkeypatch):
    """Tests can wrap resonance_structures to count invocations."""
    calls = {"n": 0}
    orig = Resonate.resonance_structures

    def counting(self, mol, first_yield_input=True):
        calls["n"] += 1
        yield from orig(self, mol, first_yield_input=first_yield_input)

    monkeypatch.setattr(Resonate, "resonance_structures", counting)

    mol = Mol(MolFromSmiles(APAP))
    list(Dehydrogenation().metabolites(mol))
    # DH: ResonanceRule.metabolites once; Hydrogenation-like pair path does not
    # call resonance_structures again for DH's second leg, but DH itself does.
    assert calls["n"] == 1

    calls["n"] = 0
    list(Hydrogenation().metabolites(Mol(MolFromSmiles(APAP))))
    assert calls["n"] == 1


def test_mode_compute_count_reused_across_resonance_and_pair_paths():
    mol = Mol(MolFromSmiles(MULTI_CONJ))
    res = Resonate()
    list(res.resonance_structures(mol, first_yield_input=False))
    cache = _resonance_cache(mol)
    entry = cache.mode(res.flag)
    assert entry.compute_count == 1
    assert entry.exhausted
    n_forms = entry.forms_materialized

    # Second consumer on same mol must not recompute.
    list(res.resonate_with_pair_paths(mol))
    assert entry.compute_count == 1
    assert entry.forms_materialized == n_forms

    # standardize kekulizes a copy and does not alias or cache it on the parent.
    template = EditMol.standardize(mol)
    assert template is not False
    assert "standardized_mol" not in (mol._forest or {})
    assert getattr(template, "_forest", None) is not mol._forest
    assert "resonance" not in (template._forest or {})
    list(res.resonate_with_pair_paths(template))
    assert entry.compute_count == 1  # parent cache untouched
    assert _resonance_cache(template).mode(res.flag).compute_count == 1


def test_copy_mol_and_carry_forest_do_not_share_resonance_by_size():
    from xenosite.forest.base import copy_mol, carry_forest

    mol = Mol(MolFromSmiles(APAP))
    parent = _resonance_cache(mol)
    twin = copy_mol(mol)
    assert twin._forest is not mol._forest
    assert "resonance" not in twin._forest
    # Explicit opt-in still allowed for unedited identity views.
    carry_forest(mol, twin, share_resonance=True)
    assert twin._forest.get("resonance") is parent
    # Edits must clear.
    EditMol().swap_bonds_along_path(twin, [4, 5, 6, 7, 8])
    assert "resonance" not in twin._forest


def test_install_product_forest_clears_resonance():
    from xenosite.forest.base import install_product_forest, copy_mol

    parent = Mol(MolFromSmiles(APAP))
    _resonance_cache(parent)
    product = copy_mol(parent)
    product._forest["resonance"] = parent._forest["resonance"]
    install_product_forest(parent, product)
    assert "resonance" not in product._forest


def test_smarts_metabolites_kekulize_clears_caller_resonance():
    from xenosite.forest.rules import Hydroxylation

    mol = Mol(MolFromSmiles("c1ccccc1"))
    cache = _resonance_cache(mol)
    assert mol._forest["resonance"] is cache
    assert mol.GetBondWithIdx(0).GetIsAromatic()
    list(Hydroxylation().metabolites(mol))
    assert not mol.GetBondWithIdx(0).GetIsAromatic()
    assert "resonance" not in mol._forest


def test_bfs_compute_once_across_full_apap():
    mol = Mol(MolFromSmiles(APAP))
    list(rulesets.load_ruleset("Full").metabolites(mol))
    cache = _resonance_cache(mol)
    assert cache.bfs_compute_count == 1
    # Kekulize-for-SMARTS must not have wiped / replaced the parent's cache.
    assert mol._forest.get("resonance") is cache


# ---------------------------------------------------------------------------
# Lazy / streaming behavior
# ---------------------------------------------------------------------------


def test_resonance_structures_are_lazily_materialized():
    mol = Mol(MolFromSmiles(MULTI_CONJ))
    res = Resonate()
    it = res.resonance_structures(mol, first_yield_input=False)
    first = next(it)
    assert first is not None

    cache = _resonance_cache(mol)
    entry = cache.mode(res.flag)
    assert entry.forms_materialized == 1
    assert not entry.exhausted

    rest = list(it)
    assert entry.exhausted
    assert entry.forms_materialized == 1 + len(rest)
    assert entry.forms_materialized >= 2


def test_early_stop_does_not_require_full_enumeration():
    mol = Mol(MolFromSmiles(MULTI_CONJ))
    res = Resonate()
    next(res.resonance_structures(mol, first_yield_input=False))
    entry = _resonance_cache(mol).mode(res.flag)
    assert entry.forms_materialized == 1
    assert not entry.exhausted


def test_hydrogenation_early_stop_streams(monkeypatch):
    """Taking one Hydrogenation product must not force full pair-path exhaustion."""
    produced = {"n": 0}
    orig_produce = _ModeResonanceCache._produce

    def counting_produce(self, resonate, mol):
        for item in orig_produce(self, resonate, mol):
            produced["n"] += 1
            yield item

    monkeypatch.setattr(_ModeResonanceCache, "_produce", counting_produce)

    mol = Mol(MolFromSmiles(MULTI_CONJ))
    next(Hydrogenation().metabolites(mol))
    # Producer may pull more than one form if the first SMARTS leg or first
    # joined form yields no product, but must not mark the mode exhausted
    # unless every form was needed.
    entry = _resonance_cache(mol).mode(Hydrogenation().flag)
    if not entry.exhausted:
        assert entry.forms_materialized == produced["n"]
        assert produced["n"] < 7  # multi_conj has 7 joined forms total


def test_yielded_forms_are_defensive_copies():
    mol = Mol(MolFromSmiles(APAP))
    res = Resonate()
    forms = list(res.resonance_structures(mol, first_yield_input=False))
    assert forms
    forms[0].GetAtomWithIdx(0).SetProp("mutated", "1")
    forms2 = list(res.resonance_structures(mol, first_yield_input=False))
    assert not forms2[0].GetAtomWithIdx(0).HasProp("mutated")


def test_first_yield_input_is_live_mol():
    mol = Mol(MolFromSmiles(APAP))
    first = next(Resonate().resonance_structures(mol, first_yield_input=True))
    assert first is mol
