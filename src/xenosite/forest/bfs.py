"""BFS/DFS metabolite enumeration and the ``xenosite-forest`` CLI entry.

``bfs`` / ``dfs`` are re-exported from :mod:`xenosite.forest.find_path`.
The CLI finds a path when two SMILES are given, or enumerates one generation
when only a reactant is given.
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence

from .find_path import bfs, dfs, find_path
from .rdkitutil import canon_smiles
from .rulesets import PhaseOne


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="xenosite-forest",
        description="Find metabolic pathways or enumerate metabolites.",
    )
    parser.add_argument(
        "molecules",
        nargs="+",
        help="Reactant SMILES, or reactant and product SMILES.",
    )
    parser.add_argument(
        "-d",
        "--depth",
        default=1,
        type=int,
        help="Maximum enumeration depth (single-molecule mode; default: 1).",
    )
    parser.add_argument(
        "-m",
        "--max",
        dest="max_paths",
        type=int,
        default=1,
        help="Stop after this many path hits (two-molecule mode; default: 1).",
    )
    parser.add_argument(
        "--max-nodes",
        dest="max_nodes",
        type=int,
        default=800,
        help="Path-search node budget (two-molecule mode; default: 800).",
    )
    parser.add_argument(
        "--search",
        choices=("bfs", "dfs"),
        default="bfs",
        help="Enumeration order for single-molecule mode (default: bfs).",
    )

    args = parser.parse_args(argv)
    mols = list(args.molecules)

    if len(mols) >= 2:
        reactant, target = mols[0], mols[1]
        n = 0
        for outcome in find_path(
            reactant,
            target,
            ruleset=PhaseOne,
            max_nodes=args.max_nodes,
            max_paths=args.max_paths,
        ):
            print((outcome.smiles, outcome.plan))
            n += 1
            if args.max_paths is not None and n >= args.max_paths:
                break
        return

    search = dfs if args.search == "dfs" else bfs
    for product, info in search(mols[0], ruleset=PhaseOne, depth=args.depth):
        site = info.get("site") if isinstance(info, dict) else getattr(info, "site", None)
        print((canon_smiles(product), site))


if __name__ == "__main__":
    main()
