#!/usr/bin/env python3
"""Render PERFORMANCE.md molecule SVGs into this directory.

Re-run from repo root:
  uv run python src/xenosite/forest/performance_assets/_render.py
"""

from __future__ import annotations

from pathlib import Path

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

OUT = Path(__file__).resolve().parent

MOOLS: dict[str, str] = {
    "eugenol_reactant": "COc1ccc(CC=C)cc1O",
    "eugenol_product": "O=C1C=CC(=O)C(CC=C)=C1",
    "dimethoxy_pea_reactant": "COc1ccc(CCN)cc1OC",
    "dimethoxy_pea_product": "NCCc1ccc(O)c(O)c1",
    "meoph_oh_reactant": "COc1ccc(O)cc1",
    "meoph_oh_product": "O=C1C=C(O)C(=O)C(O)=C1",
    "tba_reactant": "CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12",
    "tba_product": "CC(C)(C)C#CC=CC=O",
    "meo2_naph_reactant": "COc1ccc2ccccc2c1",
    "nq12_product": "O=C1C(=O)c2ccccc2C=C1",
    "nq14_product": "O=C1C=CC(=O)c2ccccc12",
    "meo1_naph_reactant": "COc1cccc2ccccc12",
}


def draw_svg(smi: str, path: Path, size: tuple[int, int] = (280, 200)) -> None:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise ValueError(f"bad SMILES: {smi}")
    Chem.rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
    drawer.drawOptions().clearBackground = True
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    path.write_text(drawer.GetDrawingText())


def main() -> None:
    for name, smi in MOOLS.items():
        path = OUT / f"{name}.svg"
        draw_svg(smi, path)
        print(f"wrote {path.name}  {smi}")


if __name__ == "__main__":
    main()
