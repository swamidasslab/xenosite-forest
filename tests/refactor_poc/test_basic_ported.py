"""Port of ``tests/test_basic.py`` chemical checks. Path finds use poc ``find_path``.

Epoxide / dehydrogenation path cases assert a non-empty plan. Dealkylation
site and sanitize checks keep the forest asserts. Nevirapine SMILES is the
constant from the forest suite (imported, not copied).
"""

from __future__ import annotations

import pytest
from rdkit.Chem.rdmolops import SanitizeMol

from test_basic import nevirapine
from xenosite.refactor_poc.find_path import find_path
from xenosite.refactor_poc.rdkit_api import MolFromSmiles
from xenosite.refactor_poc.rules import (
    Dealkylation,
    Dehydrogenation,
    Epoxidation,
    EpoxideOpening,
    Hydroxylation,
    _unique_csmi_key,
)
from xenosite.refactor_poc.rulesets import RuleSet

def test_epoxide_opening_aromatic():
    hits = list(
        find_path(
            "c1ccccc1",
            "C1=CC=CC(O)C1O",
            ruleset=RuleSet((Epoxidation, EpoxideOpening), name="EO"),
            max_nodes=200,
        )
    )
    assert hits


def test_epoxide_opening_kekulized():
    hits = list(
        find_path(
            "C1=CC=CC=C1",
            "C1=CC=CC(O)C1O",
            ruleset=RuleSet((Epoxidation, EpoxideOpening), name="EO"),
            max_nodes=200,
        )
    )
    assert hits


def test_propane_dehydrogenation_to_propene():
    hits = list(
        find_path(
            "CCC",
            "C=CC",
            ruleset=RuleSet((Dehydrogenation,), name="DH"),
            max_nodes=50,
        )
    )
    assert hits


def test_propane_dehydrogenation_followed_by_epoxidation():
    hits = list(
        find_path(
            "CCC",
            "C1OC1C",
            # RuleSet names may not contain '_'.
            ruleset=RuleSet((Dehydrogenation, Epoxidation), name="DHE"),
            max_nodes=200,
        )
    )
    assert hits


def test_hydroxyl_should_not_be_dealkylated():
    mol = MolFromSmiles("CCO")
    sites = [info["site"] for _product, info in Dealkylation().metabolize(mol)]
    assert frozenset([1, 2]) not in sites


def test_epoxide_opening2_rejects_invalid_target():
    """Forest ``C1=CC=CC1OC1`` does not parse; forest ``find_path`` treated
    ``end_mol=None`` as any-path. Poc must not soft-pass: invalid SMILES
    raises, and a depth-1 budget still misses the diol on a valid target.
    """

    rs = RuleSet((Epoxidation, EpoxideOpening), name="EO")
    with pytest.raises(ValueError, match="could not parse"):
        list(find_path("c1ccccc1", "C1=CC=CC1OC1", ruleset=rs, max_nodes=2))

    hits = list(
        find_path(
            "c1ccccc1",
            "C1=CC=CC(O)C1O",
            ruleset=rs,
            max_nodes=2,
            max_paths=3,
        )
    )
    assert hits == []


def test_epoxide_opening1_fail():
    """Depth-1 budget must not reach the diol from benzene (forest depth=1)."""

    hits = list(
        find_path(
            "c1ccccc1",
            "C1=CC=CC(O)C1O",
            ruleset=RuleSet((Epoxidation, EpoxideOpening), name="EO"),
            # Forest default depth=1. Poc hits the diol by node 3; keep a
            # one-hop ceiling so the empty-hits assert stays meaningful.
            max_nodes=2,
            max_paths=3,
        )
    )
    assert hits == []


def _dealk_sanitize(reactant: str) -> None:
    mol = MolFromSmiles(reactant)
    for _product, _info in Dealkylation().metabolize(mol):
        SanitizeMol(_product)


def test_matt_problem1():
    _dealk_sanitize("CC(C)CNS(=O)(=O)c1ccc(CCC(=O)Nc2ccc(Cl)cc2C)cc1")


def test_matt_problem2():
    _dealk_sanitize("Oc1c(C(=O)Nc2cccnc2)c(=O)n2CCc3cccc1c23")

def test_matt_problem4():
    rmol = MolFromSmiles(nevirapine)
    for _product, info in Dealkylation().metabolize(rmol):
        site = info["site"]
        if isinstance(site, int):
            continue
        atoms = tuple(site)
        if len(atoms) != 2:
            continue
        assert rmol.GetBondBetweenAtoms(*atoms), (
            "NO BOND BETWEEN ATOMS %d and %d" % tuple(a + 1 for a in atoms)
        )


def test_unique_metabolites():
    """``unique_csmi`` key is ``(rule, PatternInfo.name|SMARTS, csmi)``.

    Not forest's site+product cross-collapse: two patterns may share a site
    when their products differ; the same product must not repeat under one
    pattern token.
    """

    seen: set[tuple[str, str | None, str]] = set()
    rmol = MolFromSmiles(nevirapine)
    for product, info in Dealkylation().metabolize(rmol):
        key = _unique_csmi_key(info, product.xf.csmi)
        assert key not in seen
        seen.add(key)
    assert seen


def test_hydroxylation_partitions_h_count():
    """``h`` / ``h2`` do not both emit the same alcohol (ethane → one ``CCO``)."""

    products = list(Hydroxylation().metabolize(MolFromSmiles("CC")))
    assert [p.xf.csmi for p, _ in products] == ["CCO"]
    assert products[0][1]["pattern"].get("name") == "h2"

    aryl = list(Hydroxylation().metabolize(MolFromSmiles("c1ccccc1")))
    assert {p.xf.csmi for p, _ in aryl} == {"Oc1ccccc1"}
    assert all(info["pattern"].get("name") == "h" for _, info in aryl)
