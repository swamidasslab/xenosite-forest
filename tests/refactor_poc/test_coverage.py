"""One row per reactant, rule, and target. The fourth field names the SMARTS.

The goal is a row for every rule, every SMARTS, and every ``when``. This
file only seeds the rules the easy pairs already use.
"""

import pytest
from rdkit import Chem

from xenosite.refactor_poc.rules import Dealkylation, Dehydrogenation, Hydroxylation


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
    ],
)
def test_rule_emits_the_target_for_that_pattern(reactant, rule_cls, target, smarts):
    rule = rule_cls()
    pattern = _pattern(rule, smarts)
    assert _canon(target) in _fragments(rule, reactant, pattern)
