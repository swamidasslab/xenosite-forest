"""AtomTrace: 1-based atom numbers, 0-based depths.

# atom 1 = RDKit idx 0
"""

from __future__ import annotations

import inspect
import re

import pytest
from rdkit.Chem.rdmolfiles import MolFromSmiles, MolToSmiles

from xenosite.forest import AtomTrace, bfs, rules
from xenosite.forest.base import AtomTracker, ReactionRule, SmartsReactionRule
from xenosite.forest.trace import atom_no, rdkit_idx

from test_rules import examples


# Representative reactant for each of the 26 rule classes.
# Prefer the first example in test_rules when one exists.
_FALLBACK_REACTANTS = {
    "Acetylation": "C1=CC(=C(C=C1N)C(=O)O)O",
    "SulfurReduction": "c1ccccc1SSc1ccccc1",
    "NitroaromaticReduction": "[O-][N+](C1=CC2=C(C=C1)NC(CN=C2C3=CC=CC=C3Cl)=O)=O",
    "ThiopheneSulfurOxidation": "O=C(c1ccc(cc1)C(C(=O)O)C)c2sccc2",
}


def _rule_classes():
    return [
        cls
        for _, cls in inspect.getmembers(rules, inspect.isclass)
        if issubclass(cls, ReactionRule) and cls.__module__ == rules.__name__
    ]


def _reactant_for(rule_name: str) -> str:
    if rule_name in examples:
        return examples[rule_name][0][1]
    return _FALLBACK_REACTANTS[rule_name]


def test_atom_no_helpers():
    assert atom_no(0) == 1
    assert rdkit_idx(1) == 0


def test_reactant_only_identity_maps():
    # atom 1 = RDKit idx 0
    mol = MolFromSmiles("CCO")
    _, products = next(rules.Hydroxylation().metabolize(mol))
    assert products
    t = AtomTrace(mol)
    n_heavy = mol.GetNumHeavyAtoms()
    assert t.map() == {i: i for i in range(1, n_heavy + 1)}
    smi = MolToSmiles(mol, canonical=False)
    assert ":1" in smi and ":2" in smi and ":3" in smi
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            assert atom.GetAtomMapNum() == 0
        else:
            assert atom.GetAtomMapNum() == atom_no(atom.GetIdx())


def test_hydroxylation_cco_carbon1_order_and_maps():
    # atom 1 = RDKit idx 0  (terminal carbon of CCO)
    mol = MolFromSmiles("CCO")
    site, products = next(
        rules.Hydroxylation().metabolites_from_sites(mol, [frozenset({0})])
    )
    p = products[0]
    t = AtomTrace(p)

    smi = MolToSmiles(p, canonical=False)
    maps = [int(x) for x in re.findall(r":(\d+)", smi)]
    assert maps == [1, 2, 3]
    assert t.map() == {1: 1, 2: 3, 3: 4}
    assert t.follow(1) == (1, 1)
    assert t.origin(2) is None
    assert t.added() == frozenset({2})
    assert t.removed() == frozenset()
    assert t.depths == (0, 1)

    for atom in p.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        origin = t.origin(atom_no(atom.GetIdx()))
        if origin is None:
            assert atom.GetAtomMapNum() == 0
        else:
            assert origin == atom.GetAtomMapNum()


def test_do_not_tag_atoms_skips_tags_maps_reorder():
    mol = MolFromSmiles("CCO")
    _, products = next(
        rules.Hydroxylation().metabolize(mol, do_not_tag_atoms=True)
    )
    p = products[0]
    assert not p.HasProp(AtomTracker.tag_name)
    assert all(atom.GetAtomMapNum() == 0 for atom in p.GetAtoms())
    with pytest.raises(ValueError):
        AtomTrace(p)


class _CarbonToOxygen(SmartsReactionRule):
    smarts = ["[C:1]Cl>>[O:1].[Cl:2]"]


def test_element_change_oxygen_is_not_added():
    mol = MolFromSmiles("CCl")
    _, products = next(_CarbonToOxygen().metabolize(mol))
    oxygen = [p for p in products if p.GetAtomWithIdx(0).GetAtomicNum() == 8]
    assert oxygen
    t = AtomTrace(oxygen[0])
    assert t.origin(1) == 1
    assert 1 not in t.added()
    assert t.map()[1] == 1


def test_hydroxylation_react_atom_idx_not_in_added():
    mol = MolFromSmiles("CCO")
    _, products = next(rules.Hydroxylation().metabolize(mol))
    p = products[0]
    t = AtomTrace(p)
    for atom in p.GetAtoms():
        if not atom.HasProp("react_atom_idx"):
            continue
        n = atom_no(atom.GetIdx())
        assert n not in t.added()
        assert n in t.map().values()


def _first_product(rule_cls, smiles: str):
    mol = MolFromSmiles(smiles)
    assert mol is not None, smiles
    site, products = next(rule_cls().metabolize(mol))
    products = [p for p in products if p]
    assert products, f"{rule_cls.__name__} produced no mols from {smiles}"
    return mol, site, products[0]


@pytest.mark.parametrize("rule_cls", _rule_classes(), ids=lambda c: c.__name__)
def test_trace_every_rule(rule_cls):
    smiles = _reactant_for(rule_cls.__name__)
    reactant, _site, p = _first_product(rule_cls, smiles)
    t = AtomTrace(p)
    n_react = reactant.GetNumHeavyAtoms()
    reactant_atoms = set(range(1, n_react + 1))

    mapping = t.map()
    assert set(mapping.keys()).issubset(reactant_atoms)

    for i in mapping:
        followed = t.follow(i)
        assert followed[-1] == mapping[i]

    for i in t.added():
        assert t.origin(i) is None
        assert p.GetAtomWithIdx(rdkit_idx(i)).GetAtomMapNum() == 0

    assert t.removed().isdisjoint(mapping.keys())

    for atom in p.GetAtoms():
        if atom.HasProp("react_atom_idx"):
            assert atom_no(atom.GetIdx()) not in t.added()

    smi = MolToSmiles(p, canonical=False)
    maps = [int(x) for x in re.findall(r":(\d+)", smi)]
    assert maps == sorted(maps)


def test_multistep_bfs_maps_refer_to_original_carbons():
    result = list(
        bfs(
            ["CC", "OCCO"],
            "Full",
            depth=2,
            all_paths=False,
            outmols=True,
            phase1=True,
            ismi=True,
        )
    )
    assert result
    path_mols = result[0][2]
    p = path_mols[-1]
    t = AtomTrace(p)
    assert t.depths == tuple(range(len(path_mols)))
    assert set(t.map().keys()) == {1, 2}
    for atom in p.GetAtoms():
        if atom.GetAtomicNum() == 8:
            assert atom.GetAtomMapNum() == 0
            assert t.origin(atom_no(atom.GetIdx())) is None
        elif atom.GetAtomicNum() == 6:
            assert atom.GetAtomMapNum() in {1, 2}
