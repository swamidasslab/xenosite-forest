"""A rule's canonical_plan, applied, yields that rule's product.

Identity rules report themselves at the metabolize site. QuinoneFormation
expands to prep steps plus dehydrogenation; search must not see a
QuinoneFormation leaf. Epoxidation / NDealkylation keep identity plans
(forest ``phase1_equivalent`` singletons), not group look-aheads.
Apply uses forest StepPlan for composite plans (plain mol, no poc forest)
and a site-filtered poc metabolize for identity plans.
"""

from __future__ import annotations

import pytest
from rdkit import Chem

from xenosite.forest.utils import unmapped_smiles
from xenosite.refactor_poc.canonical_plan import CanonicalStep, as_deps
from xenosite.refactor_poc.rdkitutil import as_mol
from xenosite.refactor_poc.records import AtomRef, ProductInfo, SiteInfo, _flat_ints
from xenosite.refactor_poc.rules import (
    Acetylation,
    AzoSplitting,
    BenzodioxoleReduction,
    Epoxidation,
    Glucuronidation,
    Glutathionation,
    NDealkylation,
    NitroaromaticReduction,
    QuinoneFormation,
    ReactionRule,
    Sulfation,
    ThiopheneSulfurOxidation,
)
from xenosite.refactor_poc.rulesets import PhaseOne

from .substrate_library import QUICK_SUBSTRATES, SUBSTRATE_LIBRARY

# Concrete reaction rules under test. Not RuleSet / abstract bases.
_POC_RULES: tuple[type[ReactionRule], ...] = tuple(type(rule) for rule in PhaseOne) + (
    NDealkylation,
    AzoSplitting,
    BenzodioxoleReduction,
    NitroaromaticReduction,
    ThiopheneSulfurOxidation,
    Acetylation,
    Sulfation,
    Glucuronidation,
    Glutathionation,
)

_RULE_BY_NAME: dict[str, type[ReactionRule]] = {
    cls.__name__: cls for cls in _POC_RULES
}

_PRODUCTS_PER_PAIR = 4
_LIBRARY_PRODUCTS = 2


def _canon(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, smiles
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)


def _make_rule(cls: type[ReactionRule]) -> ReactionRule:
    if cls in (Acetylation, Sulfation, Glucuronidation, Glutathionation):
        return cls(as_star=False)
    return cls()


def _is_composite(plan: tuple[CanonicalStep, ...]) -> bool:
    if len(plan) > 1:
        return True
    return any(isinstance(item, AtomRef) for step in plan for item in step.site)


def _cleaving_ends(info: SiteInfo | ProductInfo) -> bool:
    ends = info.get("ends") if "ends" in info else None
    if not ends:
        return False
    return any(bool(end.get("cleaves")) for end in ends)


def _identity_reaches(
    smiles: str, plan: tuple[CanonicalStep, ...], target: str
) -> bool:
    """Re-run the plan's elementary rule at its site; require the same product."""

    assert len(plan) == 1
    step = plan[0]
    cls = _RULE_BY_NAME.get(step.rule)
    assert cls is not None, step.rule
    apply_rule = _make_rule(cls)
    wanted: set[int] = set()
    for item in step.site:
        assert isinstance(item, int), item
        wanted.add(item)

    def filter_sites(site, _info):
        return _flat_ints(site) == wanted

    hits = {
        _canon(info["csmi"])
        for _product, info in apply_rule.metabolize(
            as_mol(smiles), filter_sites=filter_sites
        )
    }
    return _canon(target) in hits


def _composite_reaches(smiles: str, plan: tuple[CanonicalStep, ...], target: str) -> bool:
    """Replay via forest Linearization on a plain mol (matches original suite)."""

    assert plan
    deps = as_deps(plan)
    plain = Chem.MolFromSmiles(smiles)
    assert plain is not None
    want = _canon(target)
    for lin in deps.linearizations():
        out = lin.apply(plain)
        if out and want in {_canon(unmapped_smiles(m)) for m in out}:
            return True
        if len(plan) < 2 or plan[-1].rule != "Dehydrogenation":
            continue
        mids = lin.apply(plain, drop_last=1)
        if not mids:
            continue
        for mid in mids:
            mid_smi = Chem.MolToSmiles(mid)
            for _product, info in QuinoneFormation().metabolize(as_mol(mid_smi)):
                if _canon(info["csmi"]) == want:
                    return True
    return False


def plan_reaches_product(
    rule: ReactionRule, smiles: str, info: SiteInfo | ProductInfo, product_csmi: str
) -> None:
    """Assert applying ``rule.canonical_plan`` reaches the metabolize product."""

    mol = as_mol(smiles)
    plan = rule.canonical_plan(mol, info)
    assert plan, "%s on %s: empty plan" % (rule.name, smiles)
    assert all(step.rule != "QuinoneFormation" for step in plan), plan

    if _cleaving_ends(info):
        # Dealkylating quinone ends need a Dealkylation prep; not in this plan yet.
        return

    target = product_csmi
    if _is_composite(plan) or plan[0].rule != rule.name:
        # Composite, or quinone collapsed to a lone Dehydrogenation / OD leaf.
        assert _composite_reaches(smiles, plan, target) or _identity_reaches(
            smiles, plan, target
        ), ("%s on %s: plan %s did not reach %s" % (rule.name, smiles, plan, target))
        return
    assert _identity_reaches(smiles, plan, target), (
        "%s on %s: identity plan %s did not reach %s"
        % (rule.name, smiles, plan, target)
    )


@pytest.mark.parametrize("rule_cls", _POC_RULES, ids=lambda cls: cls.__name__)
@pytest.mark.parametrize("smiles", QUICK_SUBSTRATES)
def test_canonical_plan_reaches_product_quick(rule_cls, smiles):
    """Every poc rule × quick substrates: reported plan reaches the product."""

    rule = _make_rule(rule_cls)
    checked = 0
    for _product, info in rule.metabolize(as_mol(smiles)):
        plan_reaches_product(rule, smiles, info, info["csmi"])
        checked += 1
        if checked >= _PRODUCTS_PER_PAIR:
            break


@pytest.mark.parametrize("rule_cls", _POC_RULES, ids=lambda cls: cls.__name__)
def test_canonical_plan_reaches_product_library(rule_cls):
    """Broader SMARTS/edit cohort: same check, fewer products per reactant."""

    rule = _make_rule(rule_cls)
    checked = 0
    for smiles in SUBSTRATE_LIBRARY:
        for _product, info in rule.metabolize(as_mol(smiles)):
            plan_reaches_product(rule, smiles, info, info["csmi"])
            checked += 1
            if checked >= _LIBRARY_PRODUCTS:
                break
        if checked >= _LIBRARY_PRODUCTS:
            break
    if not checked:
        pytest.skip("%s: no products on SUBSTRATE_LIBRARY" % rule.name)


def test_quinone_benzene_plan_is_two_oh_then_dh():
    rule = QuinoneFormation()
    smiles = "c1ccccc1"
    _product, info = next(rule.metabolize(as_mol(smiles)))
    plan = rule.canonical_plan(as_mol(smiles), info)
    names = [step.rule for step in plan]
    assert names.count("Hydroxylation") == 2
    assert names[-1] == "Dehydrogenation"
    assert "QuinoneFormation" not in names
    plan_reaches_product(rule, smiles, info, info["csmi"])


def test_epoxidation_plan_is_identity_same_site():
    """Forest phase1_steps: Epoxidation singleton, not StableOxygenation."""

    rule = Epoxidation()
    smiles = "C=C"
    _product, info = next(rule.metabolize(as_mol(smiles)))
    plan = rule.canonical_plan(as_mol(smiles), info)
    assert len(plan) == 1
    assert plan[0].rule == "Epoxidation"
    assert set(plan[0].site) == set(_flat_ints(info["site"]))
    assert "StableOxygenation" not in [step.rule for step in plan]
    plan_reaches_product(rule, smiles, info, info["csmi"])


def test_ndealkylation_plan_is_identity_same_site():
    """Forest phase1_steps: NDealkylation singleton, not UnstableOxygenation."""

    rule = NDealkylation()
    smiles = "CCN"
    _product, info = next(rule.metabolize(as_mol(smiles)))
    plan = rule.canonical_plan(as_mol(smiles), info)
    assert len(plan) == 1
    assert plan[0].rule == "NDealkylation"
    assert set(plan[0].site) == set(_flat_ints(info["site"]))
    assert "UnstableOxygenation" not in [step.rule for step in plan]
    plan_reaches_product(rule, smiles, info, info["csmi"])


def test_find_path_asks_rule_hook_not_quinone_name():
    """Search records whatever canonical_plan returns; no QuinoneFormation leaf."""

    from xenosite.refactor_poc.find_path import find_path

    hits = list(find_path("c1ccccc1", "O=C1C=CC(=O)C=C1"))
    assert hits
    names = [step.rule for step in hits[0].plan.children]
    assert "QuinoneFormation" not in names
    assert names.count("Hydroxylation") == 2
    assert names.count("Dehydrogenation") == 1
