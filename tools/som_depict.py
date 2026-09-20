#!/usr/bin/env python3
"""Aligned, SoM-marked molecule depictions (out-of-library helper).

Uses xenopict (repo dependency) for circle-marked atoms and MCS ``align_to``.
Not a package export — live under ``tools/`` so PERFORMANCE / notes can reuse
one depiction path without reinventing it.

Examples (repo root)::

  # Single molecule, explicit SoM atoms
  uv run python tools/som_depict.py COc1ccc(O)cc1 --som 0,1,4,5 -o /tmp/meoph.svg

  # Reactant→product pair: resolve SoM via find_path, align product to reactant
  uv run python tools/som_depict.py --pair COc1ccc(O)cc1 O=C1C=C(O)C(=O)C(O)=C1 \\
      --out-dir docs/forest/performance_assets --stem meoph_oh

  # PERFORMANCE presets (mid-size + larger H2H cases)
  uv run python tools/som_depict.py --preset performance --out-dir docs/forest/performance_assets
  uv run python tools/som_depict.py --preset larger --out-dir docs/forest/performance_assets
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence

from rdkit import Chem
from xenopict import Xenopict

from xenosite.forest.find_path import find_path
from xenosite.forest.rulesets import PhaseOne

# ---------------------------------------------------------------------------
# PERFORMANCE / H2H case tables (stems match docs/forest/PERFORMANCE.md assets)
# ---------------------------------------------------------------------------

# (reactant_stem, reactant_smi, product_stem, product_smi, som_target_override)
# som_target_override: when the product is unreachable, find_path to this instead
# for reactant SoM (e.g. 2-MeO→1,4-NQ uses demethylation to 2-naphthol).
PerformanceCase = tuple[str, str, str, str, str | None]

PERFORMANCE_CASES: list[PerformanceCase] = [
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

LARGER_CASES: list[PerformanceCase] = [
    (
        "tbu_bis_nd_reactant",
        "CN(C)Cc1ccc(CN(C)Cc2ccc(C(C)(C)C)cc2)cc1",
        "tbu_bis_nd_product",
        "O=Cc1ccc(C=O)cc1",
        None,
    ),
    (
        "macrocycle_nd_reactant",
        "C1CCCCCCNC2CCCC(CC2)NCCCC1",
        "macrocycle_nd_product",
        "NC1CCCC(=O)CC1",
        None,
    ),
    (
        "tribenzyl_reactant",
        "N(Cc1ccccc1)(Cc1ccccc1)Cc1ccccc1",
        "tribenzyl_product",
        "O=Cc1ccccc1",
        None,
    ),
    (
        "triph_butyl_reactant",
        "c1ccccc1CCCCc2ccccc2CCCCc3ccccc3",
        "triph_butyl_product",
        "Oc1ccccc1CCCCc2ccccc2CCCCc3ccccc3",
        None,
    ),
    (
        "meo_diphenyl_reactant",
        "COc1ccc(Cc2ccc(OC)cc2)cc1",
        "meo_diphenyl_product",
        "Oc1ccc(Cc2ccc(O)cc2)cc1",
        None,
    ),
]

PRESETS: dict[str, list[PerformanceCase]] = {
    "performance": PERFORMANCE_CASES,
    "larger": LARGER_CASES,
    "all": PERFORMANCE_CASES + LARGER_CASES,
}


# ---------------------------------------------------------------------------
# Core
# ---------------------------------------------------------------------------


def mol_from_smiles(smi: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise ValueError(f"bad SMILES: {smi}")
    return mol


def reactant_som_atoms(plan) -> list[int]:
    """Integer ``AtomRef.origin`` values on the reactant (SoM sites).

    Skips ``AddedRef`` / nested step refs that lack a reactant origin.
    Bond / directed_bond sites contribute both endpoint origins.
    """

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


def som_via_find_path(
    reactant_smi: str,
    target_smi: str,
    *,
    max_nodes: int = 800,
) -> list[int]:
    """Resolve reactant SoM atoms from the first live ``find_path`` plan."""

    hits = list(
        find_path(
            reactant_smi,
            target_smi,
            ruleset=PhaseOne,
            max_paths=1,
            max_nodes=max_nodes,
        )
    )
    if not hits:
        raise RuntimeError(f"no path for SoM: {reactant_smi} → {target_smi}")
    return reactant_som_atoms(hits[0].plan)


def parse_som(text: str | None) -> list[int]:
    if not text:
        return []
    parts = [p.strip() for p in text.replace(";", ",").split(",") if p.strip()]
    return [int(p) for p in parts]


def draw_reactant(
    smi: str,
    som: Sequence[int],
    path: Path | None = None,
) -> Xenopict:
    """Circle-mark SoM atoms; no atom-index labels. Optionally write SVG."""

    mol = mol_from_smiles(smi)
    n = mol.GetNumAtoms()
    bad = [i for i in som if i < 0 or i >= n]
    if bad:
        raise ValueError(f"SoM out of range for {smi}: {bad} (n={n})")
    pict = Xenopict(mol, add_atom_indices=False)
    if som:
        pict.mark_atoms(list(som))
    if path is not None:
        path.write_text(pict.to_svg())
    return pict


def draw_product_aligned(
    smi: str,
    template: Xenopict,
    path: Path | None = None,
) -> Xenopict:
    """MCS-align product to ``template`` (usually the marked reactant)."""

    pict = Xenopict(mol_from_smiles(smi), add_atom_indices=False)
    pict.align_to(template)
    if path is not None:
        path.write_text(pict.to_svg())
    return pict


def depict_pair(
    reactant_smi: str,
    product_smi: str,
    *,
    out_dir: Path,
    stem: str | None = None,
    reactant_stem: str | None = None,
    product_stem: str | None = None,
    som: Sequence[int] | None = None,
    som_target: str | None = None,
    max_nodes: int = 800,
) -> tuple[Path, Path, list[int]]:
    """Write ``{stem}_reactant.svg`` / ``{stem}_product.svg`` (or explicit stems).

    SoM: use ``som`` if given; else ``find_path`` to ``som_target or product_smi``.
    Product is MCS-aligned to the reactant depiction.
    """

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    r_stem = reactant_stem or f"{stem or 'pair'}_reactant"
    p_stem = product_stem or f"{stem or 'pair'}_product"
    r_path = out_dir / f"{r_stem}.svg"
    p_path = out_dir / f"{p_stem}.svg"

    if som is None:
        target = som_target or product_smi
        som_atoms = som_via_find_path(reactant_smi, target, max_nodes=max_nodes)
    else:
        som_atoms = list(som)

    template = draw_reactant(reactant_smi, som_atoms, r_path)
    draw_product_aligned(product_smi, template, p_path)
    return r_path, p_path, som_atoms


def depict_case(case: PerformanceCase, out_dir: Path, *, max_nodes: int = 800) -> None:
    r_stem, r_smi, p_stem, p_smi, som_target = case
    r_path, p_path, som = depict_pair(
        r_smi,
        p_smi,
        out_dir=out_dir,
        reactant_stem=r_stem,
        product_stem=p_stem,
        som_target=som_target,
        max_nodes=max_nodes,
    )
    print(f"wrote {r_path.name}  SoM={som}  {r_smi}")
    print(f"wrote {p_path.name}  aligned→{r_stem}  {p_smi}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Aligned SoM-marked depictions (xenopict). "
            "Out-of-library helper — not part of the xenosite package API."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument(
        "smiles",
        nargs="?",
        help="Single-molecule SMILES (use with --som / -o)",
    )
    p.add_argument(
        "--pair",
        nargs=2,
        metavar=("REACTANT", "PRODUCT"),
        help="Reactant and product SMILES; resolve SoM via find_path + align",
    )
    p.add_argument(
        "--stem",
        help="Basename for --pair outputs: {stem}_reactant.svg / {stem}_product.svg",
    )
    p.add_argument(
        "--som",
        help="Comma-separated reactant atom indices to mark (skip find_path)",
    )
    p.add_argument(
        "--som-target",
        help="Alternate find_path target for SoM when product is unreachable",
    )
    p.add_argument(
        "-o",
        "--out",
        type=Path,
        help="Output SVG for single-molecule mode",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=Path("."),
        help="Directory for --pair / --preset outputs (default: .)",
    )
    p.add_argument(
        "--preset",
        choices=sorted(PRESETS),
        help="Render a built-in PERFORMANCE / larger H2H case set",
    )
    p.add_argument(
        "--max-nodes",
        type=int,
        default=800,
        help="find_path max_nodes when resolving SoM (default: 800)",
    )
    return p


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)

    if args.preset:
        for case in PRESETS[args.preset]:
            depict_case(case, args.out_dir, max_nodes=args.max_nodes)
        return 0

    if args.pair:
        reactant, product = args.pair
        stem = args.stem or "pair"
        r_path, p_path, som = depict_pair(
            reactant,
            product,
            out_dir=args.out_dir,
            stem=stem,
            som=parse_som(args.som) or None,
            som_target=args.som_target,
            max_nodes=args.max_nodes,
        )
        print(f"wrote {r_path}  SoM={som}")
        print(f"wrote {p_path}  aligned")
        return 0

    if args.smiles:
        if args.out is None:
            print("single-molecule mode needs -o/--out", file=sys.stderr)
            return 2
        som = parse_som(args.som)
        draw_reactant(args.smiles, som, args.out)
        print(f"wrote {args.out}  SoM={som}")
        return 0

    _build_parser().print_help()
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
