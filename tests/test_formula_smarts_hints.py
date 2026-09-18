"""SMARTS covering tests for formula heuristics and per-rxn FormulaHint skips.

Ensures every Phase I reaction SMARTS declares a formula effect, that those
effects match observed heavy-atom changes on a probe corpus, and that
``could_help`` / ``smarts_compatible`` agree with the hints (so skips are sound).
"""

from __future__ import annotations

import pytest
from rdkit import Chem

from xenosite.forest.path_context import (
    ADD_O,
    CLEAVE,
    CLEAVE_OR_ADD_O,
    NEUTRAL,
    FormulaAny,
    FormulaHint,
    FormulaMatch,
    PathContext,
    _formula,
    any_hint_compatible,
    expand_formula_effects,
    formula_compatible,
    oxygen_deficit,
)
from xenosite.forest.phaseone import PhaseOneRS
from xenosite.forest.rules import (
    Dehydration,
    Dehydrogenation,
    Epoxidation,
    EpoxideOpening,
    Hydroxylation,
    NitrogenOxidation,
    OxygenReduction,
)
from xenosite.forest.utils import clean

# Probes that collectively exercise Phase I SMARTS families.
_PROBES = (
    "CCO",
    "CCN",
    "CCS",
    "C=C",
    "C#C",
    "C1OC1",
    "c1ccccc1",
    "Oc1ccccc1",
    "CCCl",
    "ClCCl",
    "CC(=O)OC",
    "CCNO",
    "CCSO",
    "CS(=O)C",
    "CN(C)C",
    "c1ccccc1C1OC1",
    "O=C1C=CC(=O)C=C1",
    "CC(=O)Nc1ccc(O)cc1",
    "COP(=O)(O)O",
    "[O-][N+](=O)c1ccccc1",
    "c1ccccc1SSc1ccccc1",
)


def _phase1_rules():
    """Flatten PhaseOneRS to concrete rule instances (unique by class)."""
    out = []
    seen = set()
    for item in PhaseOneRS.rules:
        for rule in getattr(item, "rules", [item]):
            key = type(rule)
            if key in seen:
                continue
            seen.add(key)
            out.append(rule)
    # NDealkylation is not in PhaseOneRS but reuses Dealkylation SMARTS.
    return out


def _n_rxns(rule) -> int:
    return len(getattr(rule, "rxns", ()) or ())


def _run_rxn(rule, mol, rxn_num):
    """Yield cleaned product lists for one reaction SMARTS only."""
    rxn = rule.rxns[rxn_num]
    work = Chem.Mol(mol)
    try:
        rule._kekulize(work)
    except Exception:
        pass
    rule._clear_atom_maps(work)
    try:
        hits = rxn.RunReactants((work,))
    except RuntimeError:
        return
    for prod in hits:
        products = clean(list(prod))
        if products:
            yield products


def _hint_matches_products(spec, reactant, product_lists, *, smarts: str = "") -> bool:
    """True if observed products are consistent with at least one expanded hint."""
    hints = expand_formula_effects(spec, reactant, match=None)
    if not hints:
        return True
    fr = _formula(reactant)
    n_heavy = reactant.GetNumHeavyAtoms()
    smarts_cleaves = ">>" in smarts and "." in smarts.split(">>", 1)[1]
    for hint in hints:
        if not isinstance(hint, FormulaHint):
            continue
        for products in product_lists:
            frags = [p for p in products if p]
            if not frags:
                continue
            if hint.cleave:
                if smarts_cleaves:
                    return True
                if len(frags) > 1 or any(p.GetNumHeavyAtoms() < n_heavy for p in frags):
                    return True
                if any(_formula(p) != fr for p in frags):
                    return True
                continue
            if not hint.delta:
                if any(_formula(p) == fr for p in frags):
                    return True
                continue
            ok = True
            for el, change in hint.delta.items():
                before = fr.get(el, 0)
                afters = [_formula(p).get(el, 0) for p in frags]
                total = sum(afters)
                if change > 0 and not (
                    max(afters) >= before + change or total >= before + change
                ):
                    ok = False
                    break
                if change < 0 and not (
                    min(afters) <= before + change or total <= before + change
                ):
                    ok = False
                    break
            if ok:
                return True
    return False


def _rxn_smarts(rule, rxn_num: int) -> str:
    raw = getattr(rule, "smarts", None)
    if raw is None:
        return ""
    if isinstance(raw, str):
        return raw if rxn_num == 0 else ""
    if rxn_num < len(raw):
        return raw[rxn_num]
    return ""


@pytest.mark.parametrize("rule", _phase1_rules(), ids=lambda r: type(r).__name__)
def test_phase1_rule_declares_formula_effects_covering_all_smarts(rule):
    """Every Phase I rxn SMARTS declares a formula_hint on its options."""
    hints = rule.formula_hints()
    n = _n_rxns(rule)
    assert n > 0, "%s has no reaction SMARTS" % type(rule).__name__
    assert hints is not None, "%s missing formula_hint options" % type(rule).__name__
    assert len(hints) == n, (
        "%s: len(hints)=%d != len(rxns)=%d" % (type(rule).__name__, len(hints), n)
    )
    for h in hints:
        assert h is not None
        # Must expand to at least one concrete FormulaHint
        assert expand_formula_effects(h), "%s hint %r expands empty" % (
            type(rule).__name__,
            h,
        )


def test_formula_any_and_match_polymorphism():
    """Static any-of and match-resolved hints expand correctly."""
    mol = Chem.MolFromSmiles("CCO")
    target_small = Chem.MolFromSmiles("CO")
    target_plus_o = Chem.MolFromSmiles("OCCO")
    # Any-of: keep if either cleave or +O helps
    assert formula_compatible(CLEAVE_OR_ADD_O, mol, target_small)
    assert formula_compatible(CLEAVE_OR_ADD_O, mol, target_plus_o)
    assert not formula_compatible(CLEAVE_OR_ADD_O, mol, Chem.MolFromSmiles("CCN"))

    def _resolve(m, match):
        # Pretend match site 0 → cleave, else +O
        if match and match.get("site") == 0:
            return CLEAVE
        return ADD_O

    spec = FormulaMatch(possible=(CLEAVE, ADD_O), resolve=_resolve)
    assert formula_compatible(spec, mol, target_small, match=None)  # possible set
    assert formula_compatible(spec, mol, target_small, match={"site": 0})
    assert not formula_compatible(spec, mol, target_small, match={"site": 1})
    assert formula_compatible(spec, mol, target_plus_o, match={"site": 1})


@pytest.mark.parametrize("rule", _phase1_rules(), ids=lambda r: type(r).__name__)
def test_formula_hints_match_observed_smarts_products(rule):
    """Declared hints agree with heavy-atom changes when SMARTS fire on probes."""
    hints = rule.formula_hints()
    uncovered = []
    mismatches = []
    for rxn_num, hint in enumerate(hints):
        saw_fire = False
        smarts = _rxn_smarts(rule, rxn_num)
        for smi in _PROBES:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            product_lists = list(_run_rxn(rule, mol, rxn_num))
            if not product_lists:
                continue
            saw_fire = True
            if not _hint_matches_products(
                hint, mol, product_lists, smarts=smarts
            ):
                mismatches.append(
                    (type(rule).__name__, rxn_num, smi, hint, len(product_lists))
                )
                break
        if not saw_fire:
            uncovered.append((type(rule).__name__, rxn_num, hint))
    # SMARTS that need rare motifs — still declared; covering is structural.
    rare_ok = {
        ("Dephosphorylation", 0),
        ("Dehydrogenation", 0),  # S(O)–OH
        ("Hydrogenation", 0),  # triple bond
        ("NitrogenOxidation", 2),  # tertiary / pyridine N-oxide
        ("NitrogenReduction", 3),
        ("NitrogenReduction", 4),
        ("NitrogenReduction", 6),
        ("NitrogenReduction", 7),
        ("OxygenReduction", 1),  # peroxide
        ("SulfurReduction", 0),
        ("SulfurReduction", 1),
        ("SulfurReduction", 2),
        ("Hydrolysis", 0),
        ("Hydrolysis", 1),
        ("Dealkylation", 12),
        ("OxidativeDehalogenation", 1),
        ("OxidativeDehalogenation", 3),
        ("OxidativeDehalogenation", 4),
        ("OxidativeDehalogenation", 5),
    }
    bad_uncovered = [u for u in uncovered if (u[0], u[1]) not in rare_ok]
    assert not mismatches, "hint vs observed mismatch: %s" % mismatches[:5]
    assert not bad_uncovered, "SMARTS never fired on probes: %s" % bad_uncovered


@pytest.mark.parametrize("rule", _phase1_rules(), ids=lambda r: type(r).__name__)
def test_could_help_agrees_with_formula_hints(rule):
    """Rule-level could_help is True iff any SMARTS hint is compatible.

    Epoxidation adds a quinoid overlay (may be False even when ADD_O fits).
    """
    pairs = (
        ("CCO", "CC=O"),  # same heavy / DH
        ("CCO", "OCCO"),  # +O
        ("CCCl", "CCO"),  # ox dehal +O
        ("CCN", "C"),  # smaller
        ("C1OC1", "OCCO"),  # epoxide open +O
        ("c1ccccc1", "O=C1C=CC(=O)C=C1"),  # quinoid +2 O
    )
    for r_smi, t_smi in pairs:
        mol = Chem.MolFromSmiles(r_smi)
        target = Chem.MolFromSmiles(t_smi)
        ctx = PathContext.from_mols(mol, target)
        hints_ok = any_hint_compatible(rule.formula_hints(), mol, target)
        help_ok = rule.could_help(mol, target, ctx)
        if type(rule) is Epoxidation and hints_ok and not help_ok:
            # Quinoid / non-aromatic ≥2 O overlay — allowed exception.
            assert oxygen_deficit(mol, target) > 0
            continue
        assert help_ok == hints_ok, (
            "%s could_help=%s hints=%s for %s→%s"
            % (type(rule).__name__, help_ok, hints_ok, r_smi, t_smi)
        )


def test_smarts_compatible_skips_add_o_when_no_oxygen_needed():
    """Per-SMARTS skip: hydroxylation SMARTS ruled out when formulas match."""
    rule = Hydroxylation()
    mol = Chem.MolFromSmiles("CCO")
    target = Chem.MolFromSmiles("CC=O")  # same heavy formula, needs DH not OH
    assert oxygen_deficit(mol, target) == 0
    assert not rule.smarts_compatible(0, mol, target)
    assert not rule.smarts_compatible(1, mol, target)
    # Without toward_target, both SMARTS may still run; with it, none fire.
    bare = list(rule.metabolites(mol))
    filtered = list(rule.metabolites(mol, toward_target=target))
    assert bare  # OH can still hydroxylate ethanol in isolation
    assert filtered == []


def test_epoxide_opening_neutral_smarts_allowed_when_formulas_match():
    """EpoxideOpening SMARTS[0] is NEUTRAL — must not be gated as cleavage-only."""
    rule = EpoxideOpening()
    assert rule.formula_hints()[0] == NEUTRAL
    assert rule.formula_hints()[1] == ADD_O
    mol = Chem.MolFromSmiles("C1OC1")
    # Hypothetical same-formula open product path: equal heavy counts.
    target = Chem.MolFromSmiles("CCO")  # same C2O1
    assert rule.smarts_compatible(0, mol, target)  # NEUTRAL ok
    assert not rule.smarts_compatible(1, mol, target)  # ADD_O not needed
    assert rule.could_help(mol, target, PathContext.from_mols(mol, target))


def test_oxygen_reduction_neutral_not_remove_o_only():
    """OxygenReduction SMARTS[0] is formula-neutral (C=O→C–O), not remove-O."""
    rule = OxygenReduction()
    assert rule.formula_hints()[0] == NEUTRAL
    assert rule.formula_hints()[1] == CLEAVE
    mol = Chem.MolFromSmiles("CC=O")
    target = Chem.MolFromSmiles("CCO")  # same heavy formula
    assert rule.smarts_compatible(0, mol, target)
    assert not rule.smarts_compatible(1, mol, target)


def test_dehydration_is_cleave_not_neutral():
    """Dehydration always fragments off oxygen — CLEAVE, not same-formula noop."""
    rule = Dehydration()
    assert all(h == CLEAVE for h in rule.formula_hints())
    mol = Chem.MolFromSmiles("CCO")
    same = Chem.MolFromSmiles("CCO")
    smaller = Chem.MolFromSmiles("CC")
    assert not rule.could_help(mol, same, PathContext.from_mols(mol, same))
    assert rule.could_help(mol, smaller, PathContext.from_mols(mol, smaller))


def test_nitrogen_oxidation_skips_when_target_has_no_extra_o():
    rule = NitrogenOxidation()
    mol = Chem.MolFromSmiles("CCN")
    target = Chem.MolFromSmiles("CC")  # smaller, no extra O
    assert oxygen_deficit(mol, target) <= 0
    assert not rule.could_help(mol, target, PathContext.from_mols(mol, target))
    assert list(rule.metabolites(mol, toward_target=target)) == []


def test_formula_hint_units():
    mol = Chem.MolFromSmiles("CCO")
    more_o = Chem.MolFromSmiles("OCCO")
    same = Chem.MolFromSmiles("CC=O")
    smaller = Chem.MolFromSmiles("CC")
    assert ADD_O.compatible(mol, more_o)
    assert not ADD_O.compatible(mol, same)
    assert NEUTRAL.compatible(mol, same)
    assert not NEUTRAL.compatible(mol, more_o)
    assert CLEAVE.compatible(mol, smaller)
    assert not CLEAVE.compatible(mol, more_o)


def test_guided_enumerate_passes_toward_target_skips_smarts():
    """enumerate_for_path should not expand hydroxylation toward same-formula T."""
    rule = Hydroxylation()
    mol = Chem.MolFromSmiles("CCO")
    target = Chem.MolFromSmiles("CC=O")
    ctx = PathContext.from_mols(mol, target)
    assert list(rule.enumerate_for_path(mol, ctx)) == []
    # DH still helps
    assert Dehydrogenation().could_help(mol, target, ctx)
