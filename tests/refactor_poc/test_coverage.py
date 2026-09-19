"""One row per reactant, rule, and target. The fourth field names the SMARTS.

The goal is a row for every rule, every SMARTS, and every ``when``. This
file only seeds the rules the easy pairs already use.
"""

import pytest
from rdkit import Chem

from xenosite.refactor_poc.rules import (
    Acetylation,
    AzoSplitting,
    BenzodioxoleReduction,
    ConjugationRule,
    Dealkylation,
    Dehydrogenation,
    Hydroxylation,
    NDealkylation,
    NitroaromaticReduction,
    ThiopheneSulfurOxidation,
)


def _canon(smiles):
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, smiles
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)


def _fragments(rule, smiles, pattern):
    found = set()

    def filter_rules(rule, info):
        return info is pattern

    for product, _info in rule.metabolize(
        Chem.MolFromSmiles(smiles), filter_rules=filter_rules
    ):
        parsed = Chem.MolFromSmiles(Chem.MolToSmiles(product))
        assert parsed is not None
        for piece in Chem.GetMolFrags(parsed, asMols=True, sanitizeFrags=True):
            found.add(Chem.MolToSmiles(piece, canonical=True, isomericSmiles=False))
    return found


def _pattern(rule, smarts):
    groups = []
    for name in ("smarts", "endpoints"):
        group = getattr(rule, name, None)
        if group:
            groups.extend(group)
    for text, info in groups:
        if text == smarts:
            return info
    raise AssertionError(smarts)


@pytest.mark.parametrize(
    ("reactant", "rule_cls", "target", "smarts"),
    [
        ("c1ccccc1", Hydroxylation, "Oc1ccccc1", "[#6h:1]>>[*:1]O"),
        ("CCC", Hydroxylation, "CC(O)C", "[#6h2:1]>>[*:1]O"),
        (
            "COc1ccccc1",
            Dealkylation,
            "Oc1ccccc1",
            "[#6H3:1][#7,#8H0,#16:2]>>([*:2].[*:1]=O)",
        ),
        (
            "Oc1ccc(O)cc1",
            Dehydrogenation,
            "O=C1C=CC(=O)C=C1",
            "[#6:1]-[#8H:2]",
        ),
        (
            "CN(C)C",
            NDealkylation,
            "CNC",
            "[#6H3:1][#7:2]>>([*:2].[*:1]=O)",
        ),
        (
            "OC(=O)c1cc(/N=N/c2ccc(c(c2)C(=O)O)O)ccc1O",
            AzoSplitting,
            "Nc1ccc(O)c(C(=O)O)c1",
            "[#7:1]=[#7:2]>>[*:1].[*:2]",
        ),
        (
            "c1ccc2c(c1)OCO2",
            BenzodioxoleReduction,
            "Oc1ccccc1O",
            "[#6R:1]-[#8R:2]-[#6H2R:3]-[#8R:4]-[#6R:5]>>([*:1]-[*:2].[*:3].[*:4]-[*:5])",
        ),
        (
            "O=[N+]([O-])c1ccccc1",
            NitroaromaticReduction,
            "O=Nc1ccccc1",
            "[#8-1:1]-[#7+1:2]([#6R:4])=[#8:3]>>[*:1].[*:2]([*:4])=[*:3]",
        ),
        (
            "c1ccsc1",
            ThiopheneSulfurOxidation,
            "[O-][s+]1cccc1",
            "[#6:2]1=[#6:3][#6:4]=[#6:5][#16;v2,v4:1]1>>[*:2]1=[*:3][*:4]=[*:5][*&H0&+:1]1[O-]",
        ),
        (
            "Oc1ccccc1",
            ConjugationRule,
            "*Oc1ccccc1",
            "[#7,#8,#16;h:1]>>[*:1][#6](=[#8])[#6]",
        ),
        (
            "Oc1ccccc1",
            Acetylation,
            "CC(=O)Oc1ccccc1",
            "[#7,#8,#16;h:1]>>[*:1][#6](=[#8])[#6]",
        ),
    ],
)
def test_rule_emits_the_target_for_that_pattern(reactant, rule_cls, target, smarts):
    # Star collapse stays on the conjugation class. A row whose target is the
    # acetyl, not the star, asks for that conjugate.
    if getattr(rule_cls, "as_star", False) and "*" not in target:
        rule = rule_cls(as_star=False)
    else:
        rule = rule_cls()
    pattern = _pattern(rule, smarts)
    assert _canon(target) in _fragments(rule, reactant, pattern)
