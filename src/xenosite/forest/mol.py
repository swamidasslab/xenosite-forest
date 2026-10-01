"""Product ``ForestMol`` wrapper with optional lazy RDKit interop."""

from __future__ import annotations

import ast
from typing import Any

from xenosite.forest._ext import load


def _lazy_rdkit() -> Any:
    try:
        from rdkit import Chem  # noqa: F401
    except ImportError as e:
        raise ImportError(
            "RDKit is required for ForestMol RDKit interop; "
            "install with the xenosite-forest[rdkit] extra"
        ) from e
    from rdkit import Chem

    return Chem


def _is_rdkit_mol(spec: Any) -> bool:
    Chem = _lazy_rdkit()
    return isinstance(spec, Chem.Mol)


def _smiles_atom_output_order(mol: Any) -> list[int]:
    """Parse RDKit ``_smilesAtomOutputOrder`` set by the last ``MolToSmiles``."""

    raw = mol.GetProp("_smilesAtomOutputOrder")
    order = ast.literal_eval(raw)
    return [int(i) for i in order]


def _require_identity_smiles_atom_order(mol: Any) -> None:
    """Raise if the last ``MolToSmiles`` did not emit atoms in ``GetIdx()`` order."""

    order = _smiles_atom_output_order(mol)
    expected = list(range(mol.GetNumAtoms()))
    if order != expected:
        raise ValueError(
            "RDKit MolToSmiles reordered heavy atoms; "
            "cannot build ForestMol in the input GetIdx() site frame "
            f"(output order {order}, expected {expected})"
        )


def _rdkit_to_smiles_preserving_atom_order(mol: Any) -> str:
    """Non-canonical SMILES whose write order matches ``GetIdx()``, or raise."""

    Chem = _lazy_rdkit()
    smiles = Chem.MolToSmiles(mol, canonical=False)
    _require_identity_smiles_atom_order(mol)
    return smiles


class ForestMol:
    """Delegates to the Rust ``ForestMol`` pyclass; adds RDKit helpers.

    **Atom-index frame:** metabolize / emission ``site`` indexes are in the
    **input** heavy-atom order. For a SMILES string that is the parse order of
    that string (display ``csmi`` may still be canonical spelling). For an
    RDKit ``Mol``, ingest writes non-canonical SMILES and checks
    ``_smilesAtomOutputOrder`` is identity — do **not** round-trip through
    canonical ``MolToSmiles`` (that reorders, e.g. ``c1ccccc1OC`` →
    ``COc1ccccc1``). ``to_rdkit()`` rebuilds from ``csmi`` and does **not**
    preserve indexes.
    """

    __slots__ = ("_inner",)

    def __init__(self, spec: str | Any) -> None:
        rust = load().ForestMol
        if isinstance(spec, str):
            self._inner = rust(spec)
            return
        if _is_rdkit_mol(spec):
            smiles = _rdkit_to_smiles_preserving_atom_order(spec)
            self._inner = rust(smiles)
            return
        if isinstance(spec, rust):
            self._inner = spec
            return
        raise TypeError(
            "ForestMol(spec): spec must be SMILES str, rdkit.Chem.Mol, or Rust ForestMol; "
            f"got {type(spec)!r}"
        )

    def __getattr__(self, name: str) -> Any:
        return getattr(self._inner, name)

    def __str__(self) -> str:
        return str(self._inner)

    def __repr__(self) -> str:
        return repr(self._inner)

    def _repr_html_(self) -> str:
        from . import notebook

        return notebook._forest_mol_html(self._inner)

    def to_rdkit(self) -> Any:
        """RDKit view from chematic CSMI (tags and input atom order are lost)."""

        Chem = _lazy_rdkit()
        mol = Chem.MolFromSmiles(self.csmi)
        if mol is None:
            raise ValueError(f"RDKit could not parse ForestMol CSMI: {self.csmi!r}")
        return mol

    @classmethod
    def from_rdkit(cls, mol: Any) -> ForestMol:
        """Build from RDKit ``Mol`` preserving heavy-atom ``GetIdx()`` order."""

        return cls(mol)

    @property
    def _rust(self) -> Any:
        """Rust pyclass payload (for APIs that require the extension type)."""

        return self._inner
