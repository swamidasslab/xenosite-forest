"""Shared helpers for refactor_poc correctness ports. SMILES stay in the source suites."""

from __future__ import annotations

from rdkit.Chem.rdmolops import RemoveStereochemistry

from xenosite.refactor_poc.find_path import find_path
from xenosite.refactor_poc.rdkit_api import (
    Mol,
    MolFromSmiles,
    MolToSmiles,
)
from xenosite.refactor_poc.rulesets import PhaseOne


def canon(smiles_or_mol: Mol | str | None) -> str:
    """Canonical non-isomeric SMILES; empty string if the input does not parse."""

    if smiles_or_mol is None:
        return ""
    if isinstance(smiles_or_mol, str):
        mol = MolFromSmiles(smiles_or_mol)
    else:
        mol = Mol(smiles_or_mol)
    if mol is None:
        return ""
    RemoveStereochemistry(mol)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return MolToSmiles(mol, canonical=True, isomericSmiles=False)


def product_smiles(rule, reactant: str) -> set[str]:
    """Canonical SMILES of every product mol a rule emits from ``reactant``.

    Production yields one connected mol per piece. A SMILES round-trip
    matches historical aromatic canons; it does not split dotted products.
    """

    mol = MolFromSmiles(reactant)
    assert mol is not None, reactant
    RemoveStereochemistry(mol)
    found: set[str] = set()
    for product, _info in rule.metabolize(mol):
        text = MolToSmiles(product)
        assert "." not in text, text
        parsed = MolFromSmiles(text)
        if parsed is None:
            continue
        found.add(canon(parsed))
    return found


def emits_product(rule, reactant: str, product: str) -> bool:
    """True when ``rule`` emits a single-component mol matching ``product``."""

    return canon(product) in product_smiles(rule, reactant)


def find_phaseone(reactant: str, product: str, *, max_nodes: int = 200, max_paths: int = 3):
    """``find_path`` over catalog ``PhaseOne``, not the tiny default Poc set."""

    return list(
        find_path(
            reactant,
            product,
            ruleset=PhaseOne,
            max_nodes=max_nodes,
            max_paths=max_paths,
        )
    )
