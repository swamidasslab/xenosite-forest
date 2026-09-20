#!/usr/bin/env python3
"""Render PERFORMANCE.md molecule SVGs with xenopict.

Reactants: xenopict circle marks on every atom + atom indices.
Products: MCS-aligned to the paired reactant (``Xenopict.align_to``).

Re-run from repo root (needs ``xenopict``; lab path via ``tool.uv.sources``):

  uv run python src/xenosite/forest/performance_assets/_render.py
"""

from __future__ import annotations

from pathlib import Path

from rdkit import Chem
from xenopict import Xenopict

OUT = Path(__file__).resolve().parent

# (asset_stem, smiles) for standalone reactants that are drawn once.
# Products list (product_stem, product_smiles, reactant_stem) so each product
# is oriented to its reactant; 1,4-NQ has two orientations.
REACTANTS: list[tuple[str, str]] = [
    ("eugenol_reactant", "COc1ccc(CC=C)cc1O"),
    ("dimethoxy_pea_reactant", "COc1ccc(CCN)cc1OC"),
    ("meoph_oh_reactant", "COc1ccc(O)cc1"),
    ("tba_reactant", "CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12"),
    ("meo2_naph_reactant", "COc1ccc2ccccc2c1"),
    ("meo1_naph_reactant", "COc1cccc2ccccc12"),
]

PRODUCTS: list[tuple[str, str, str]] = [
    ("eugenol_product", "O=C1C=CC(=O)C(CC=C)=C1", "eugenol_reactant"),
    ("dimethoxy_pea_product", "NCCc1ccc(O)c(O)c1", "dimethoxy_pea_reactant"),
    ("meoph_oh_product", "O=C1C=C(O)C(=O)C(O)=C1", "meoph_oh_reactant"),
    ("tba_product", "CC(C)(C)C#CC=CC=O", "tba_reactant"),
    ("nq12_product", "O=C1C(=O)c2ccccc2C=C1", "meo2_naph_reactant"),
    # Same 1,4-NQ SMILES, two alignments (one per reactant orientation).
    ("nq14_from_2meo_product", "O=C1C=CC(=O)c2ccccc12", "meo2_naph_reactant"),
    ("nq14_from_1meo_product", "O=C1C=CC(=O)c2ccccc12", "meo1_naph_reactant"),
]


def _mol(smi: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise ValueError(f"bad SMILES: {smi}")
    return mol


def draw_reactant(smi: str, path: Path) -> Xenopict:
    """Circle-mark every atom and show RDKit atom indices (via xenopict)."""

    mol = _mol(smi)
    pict = Xenopict(mol, add_atom_indices=True)
    pict.mark_atoms(list(range(mol.GetNumAtoms())))
    path.write_text(pict.to_svg())
    return pict


def draw_product_aligned(smi: str, template: Xenopict, path: Path) -> None:
    """Draw product MCS-aligned to the reactant template."""

    pict = Xenopict(_mol(smi))
    pict.align_to(template)
    path.write_text(pict.to_svg())


def main() -> None:
    templates: dict[str, Xenopict] = {}
    for stem, smi in REACTANTS:
        path = OUT / f"{stem}.svg"
        templates[stem] = draw_reactant(smi, path)
        print(f"wrote {path.name}  reactant  {smi}")

    for stem, smi, ref in PRODUCTS:
        path = OUT / f"{stem}.svg"
        draw_product_aligned(smi, templates[ref], path)
        print(f"wrote {path.name}  product→{ref}  {smi}")

    # Drop obsolete unaligned 1,4-NQ asset if present.
    old = OUT / "nq14_product.svg"
    if old.exists():
        old.unlink()
        print(f"removed {old.name} (replaced by nq14_from_*_product)")


if __name__ == "__main__":
    main()
