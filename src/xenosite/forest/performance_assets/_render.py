#!/usr/bin/env python3
"""Render PERFORMANCE.md molecule SVGs with xenopict.

Reactants: xenopict ``mark_atoms`` circles on **sites of metabolism only**
(no atom-index labels). Sites come from the live ``find_path`` plan's integer
``AtomRef.origin`` values on the reactant.

Products: MCS-aligned to the paired reactant via ``Xenopict.align_to``.

Re-run from repo root (needs ``xenopict`` via ``tool.uv.sources``):

  uv run python src/xenosite/forest/performance_assets/_render.py
"""

from __future__ import annotations

from pathlib import Path

from rdkit import Chem
from xenopict import Xenopict

from xenosite.forest.find_path import find_path
from xenosite.forest.rulesets import PhaseOne

OUT = Path(__file__).resolve().parent

# (reactant_stem, reactant_smi, product_stem, product_smi)
# product_smi is also the find_path target used to resolve SOM on the reactant.
# For the no-path 2-MeO→1,4-NQ case, SOM is taken from demethylation to 2-naphthol.
CASES: list[tuple[str, str, str, str, str | None]] = [
    # reactant_stem, reactant_smi, product_stem, product_smi, som_target_override
    (
        "eugenol_reactant",
        "COc1ccc(CC=C)cc1O",
        "eugenol_product",
        "O=C1C=CC(=O)C(CC=C)=C1",
        None,
    ),
    (
        "dimethoxy_pea_reactant",
        "COc1ccc(CCN)cc1OC",
        "dimethoxy_pea_product",
        "NCCc1ccc(O)c(O)c1",
        None,
    ),
    (
        "meoph_oh_reactant",
        "COc1ccc(O)cc1",
        "meoph_oh_product",
        "O=C1C=C(O)C(=O)C(O)=C1",
        None,
    ),
    (
        "tba_reactant",
        "CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12",
        "tba_product",
        "CC(C)(C)C#CC=CC=O",
        None,
    ),
    (
        "meo2_naph_reactant",
        "COc1ccc2ccccc2c1",
        "nq12_product",
        "O=C1C(=O)c2ccccc2C=C1",
        None,
    ),
    # No PhaseOne path to 1,4-NQ: mark demethylation SOM (path to 2-naphthol).
    (
        "meo2_naph_reactant_nopath",
        "COc1ccc2ccccc2c1",
        "nq14_from_2meo_product",
        "O=C1C=CC(=O)c2ccccc12",
        "Oc1ccc2ccccc2c1",
    ),
    (
        "meo1_naph_reactant",
        "COc1cccc2ccccc12",
        "nq14_from_1meo_product",
        "O=C1C=CC(=O)c2ccccc12",
        None,
    ),
]


def _mol(smi: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise ValueError(f"bad SMILES: {smi}")
    return mol


def reactant_som_atoms(plan) -> list[int]:
    """Integer AtomRef origins on the reactant (skip AddedRef / nested step refs)."""

    atoms: set[int] = set()
    for step in getattr(plan, "steps", ()) or ():
        site = getattr(step, "site", None) or ()
        for ref in site:
            origin = getattr(ref, "origin", None)
            if origin is None:
                continue
            try:
                atoms.add(int(origin))
            except (TypeError, ValueError):
                continue
    return sorted(atoms)


def som_for_pair(reactant_smi: str, target_smi: str) -> list[int]:
    hits = list(
        find_path(
            reactant_smi,
            target_smi,
            ruleset=PhaseOne,
            max_paths=1,
            max_nodes=800,
        )
    )
    if not hits:
        raise RuntimeError(f"no path for SOM: {reactant_smi} → {target_smi}")
    return reactant_som_atoms(hits[0].plan)


def draw_reactant(smi: str, som: list[int], path: Path) -> Xenopict:
    """Circle-mark SOM atoms only; no atom-index labels."""

    mol = _mol(smi)
    n = mol.GetNumAtoms()
    bad = [i for i in som if i < 0 or i >= n]
    if bad:
        raise ValueError(f"SOM out of range for {smi}: {bad} (n={n})")
    pict = Xenopict(mol, add_atom_indices=False)
    if som:
        pict.mark_atoms(som)
    path.write_text(pict.to_svg())
    return pict


def draw_product_aligned(smi: str, template: Xenopict, path: Path) -> None:
    pict = Xenopict(_mol(smi), add_atom_indices=False)
    pict.align_to(template)
    path.write_text(pict.to_svg())


def main() -> None:
    for r_stem, r_smi, p_stem, p_smi, som_target in CASES:
        target = som_target or p_smi
        som = som_for_pair(r_smi, target)
        r_path = OUT / f"{r_stem}.svg"
        template = draw_reactant(r_smi, som, r_path)
        print(f"wrote {r_path.name}  SOM={som}  {r_smi}")

        p_path = OUT / f"{p_stem}.svg"
        draw_product_aligned(p_smi, template, p_path)
        print(f"wrote {p_path.name}  aligned→{r_stem}  {p_smi}")

    # Drop obsolete shared 2-MeO reactant if we now use the nopath stem for 5b.
    # Keep meo2_naph_reactant.svg for the 1,2-NQ case.


if __name__ == "__main__":
    main()
