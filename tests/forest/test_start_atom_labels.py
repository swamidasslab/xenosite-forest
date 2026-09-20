"""Input-model CX ``atomLabel`` on start atoms survives stamp and metabolize.

Without ``_forest["immutable"]["start_labels"]`` capture + restamp on the
stamp/check path, RunReactants / fragment finishing drop those props and they
stay gone.
"""

from __future__ import annotations

from rdkit import Chem

from xenosite.forest.rules import Hydroxylation


def _set_start_labels(mol, labels: dict[int, str]) -> None:
    for idx, text in labels.items():
        mol.GetAtomWithIdx(idx).SetProp("atomLabel", text)


def _labels_by_root(mol) -> dict[int, str]:
    tracing = mol.xf.tracing
    out: dict[int, str] = {}
    for atom in mol.GetAtoms():
        if not atom.HasProp("atomLabel"):
            continue
        root = tracing.atom_root(atom.GetIdx())
        if root is None:
            continue
        out[root] = atom.GetProp("atomLabel")
    return out


def test_stamp_keeps_input_atom_labels():
    mol = Chem.MolFromSmiles("CCO")
    _set_start_labels(mol, {0: "alpha", 2: "oxy"})
    stamped = mol.xf.tracing._stamp()
    assert stamped.GetAtomWithIdx(0).GetProp("atomLabel") == "alpha"
    assert stamped.GetAtomWithIdx(2).GetProp("atomLabel") == "oxy"
    assert dict(stamped._forest["immutable"]["start_labels"]) == {0: "alpha", 2: "oxy"}


def test_metabolize_preserves_start_atom_labels():
    mol = Chem.MolFromSmiles("CCO")
    _set_start_labels(mol, {0: "alpha", 1: "beta", 2: "oxy"})
    product, info = next(Hydroxylation().metabolize(mol))
    assert info["site"] == frozenset({0})
    preserved = _labels_by_root(product)
    assert preserved == {0: "alpha", 1: "beta", 2: "oxy"}
    # New oxygen has no depth-0 root and must not steal a start label.
    for atom in product.GetAtoms():
        if product.xf.tracing.atom_root(atom.GetIdx()) is None:
            assert not atom.HasProp("atomLabel")


def test_restamp_restores_cleared_start_atom_labels():
    mol = Chem.MolFromSmiles("CCO")
    _set_start_labels(mol, {0: "kept"})
    mol.xf.tracing._stamp()
    mol.GetAtomWithIdx(0).ClearProp("atomLabel")
    assert not mol.GetAtomWithIdx(0).HasProp("atomLabel")
    mol.xf.tracing._stamp()
    assert mol.GetAtomWithIdx(0).GetProp("atomLabel") == "kept"
