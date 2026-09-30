"""Product ``ForestMol`` wrapper with optional lazy RDKit interop."""

from __future__ import annotations

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


class ForestMol:
    """Delegates to the Rust ``ForestMol`` pyclass; adds RDKit helpers."""

    __slots__ = ("_inner",)

    def __init__(self, spec: str | Any) -> None:
        rust = load().ForestMol
        if isinstance(spec, str):
            self._inner = rust(spec)
            return
        if _is_rdkit_mol(spec):
            Chem = _lazy_rdkit()
            smiles = Chem.MolToSmiles(spec)
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
        """RDKit view from chematic CSMI (atom tags are not round-tripped)."""

        Chem = _lazy_rdkit()
        mol = Chem.MolFromSmiles(self.csmi)
        if mol is None:
            raise ValueError(f"RDKit could not parse ForestMol CSMI: {self.csmi!r}")
        return mol

    @classmethod
    def from_rdkit(cls, mol: Any) -> ForestMol:
        return cls(mol)

    @property
    def _rust(self) -> Any:
        """Rust pyclass payload (for APIs that require the extension type)."""

        return self._inner
