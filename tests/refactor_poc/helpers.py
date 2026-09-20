"""Shared helpers for refactor_poc correctness ports. SMILES stay in the source suites."""

from __future__ import annotations

from rdkit.Chem.rdmolops import RemoveStereochemistry

from xenosite.refactor_poc.rdkit_api import (
    GetMolFrags,
    Mol,
    MolFromSmiles,
    MolToSmiles,
)


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
    """Canonical SMILES of every fragment a rule emits from ``reactant``."""

    mol = MolFromSmiles(reactant)
    assert mol is not None, reactant
    RemoveStereochemistry(mol)
    found: set[str] = set()
    for product, _info in rule.metabolize(mol):
        text = MolToSmiles(product)
        parsed = MolFromSmiles(text)
        if parsed is None:
            continue
        for piece in GetMolFrags(parsed, asMols=True, sanitizeFrags=True):
            found.add(canon(piece))
    return found


def emits_product(rule, reactant: str, product: str) -> bool:
    """True when ``rule`` emits a molecule matching ``product``."""

    target = canon(product)
    mol = MolFromSmiles(reactant)
    assert mol is not None, reactant
    return any(canon(pred) == target for pred, _info in rule.metabolize(mol))
