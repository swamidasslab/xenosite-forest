"""Find metabolic pathways between a reactant and an optional product."""

import argparse
import sys

from rdkit import rdBase

from .phaseone import PhaseOneRS
from .rulesets import RULESETS, load_ruleset
from .utils import load

RULESETS["PhaseOneRS"] = PhaseOneRS

rdBase.DisableLog("rdApp.*")


def _run_search(
    molstrings,
    ruleset="Full",
    termination_ruleset=None,
    max_paths=None,
    **kwargs,
):
    inputs = load(molstrings)

    if None in inputs:
        raise ValueError("Problem with input.")

    rules = load_ruleset(ruleset)

    if termination_ruleset is not None:
        termination_rulenames = load_ruleset(termination_ruleset).rulenames
    else:
        termination_rulenames = []

    def _iter():
        if len(inputs) > 2:
            for num, mol in enumerate(inputs):
                if not mol.HasProp("_Name"):
                    mol.SetProp("_Name", "Molecule%d" % num)
                yield from rules.find_path(
                    mol, termination_rulenames=termination_rulenames, **kwargs
                )
            return

        yield from rules.find_path(
            *inputs, termination_rulenames=termination_rulenames, **kwargs
        )

    n = 0
    for item in _iter():
        yield item
        n += 1
        if max_paths is not None and n >= max_paths:
            return


def bfs(
    molstrings,
    ruleset="Full",
    termination_ruleset=None,
    max_paths=None,
    shuffle_rng=None,
    **kwargs,
):
    """Run BFS to find paths linking reactants and an optional putative product.

    Args:
        molstrings: SMILES strings, RDKit mols, or a mix.
        ruleset: named ruleset, RuleSet instance, or list of those.
        termination_ruleset: optional ruleset whose rule names stop expansion.
        max_paths: optional cap on yielded pathways (``None`` = unlimited).
        shuffle_rng: optional ``random.Random`` to shuffle frontier / reaction /
            product order.
        **kwargs: forwarded to RuleSet.find_path (depth, phase1, all_paths,
            expand_star_conjugates, ...). Star adducts are not expanded further
            unless ``expand_star_conjugates=True``.
    """
    kwargs.setdefault("search", "bfs")
    if shuffle_rng is not None:
        kwargs["shuffle_rng"] = shuffle_rng
    yield from _run_search(
        molstrings,
        ruleset=ruleset,
        termination_ruleset=termination_ruleset,
        max_paths=max_paths,
        **kwargs,
    )


def dfs(
    molstrings,
    ruleset="Full",
    termination_ruleset=None,
    max_paths=None,
    shuffle_rng=None,
    **kwargs,
):
    """Depth-first pathway search (same arguments as ``bfs``).

    Yields deep paths earlier than BFS, so sampling a few two-step pathways
    does not require finishing the full breadth-first frontier. Pass
    ``shuffle_rng`` to randomize reaction / product order. Pass ``max_paths``
    to stop after that many yields (``None`` = unlimited).
    """
    kwargs["search"] = "dfs"
    if shuffle_rng is not None:
        kwargs["shuffle_rng"] = shuffle_rng
    yield from _run_search(
        molstrings,
        ruleset=ruleset,
        termination_ruleset=termination_ruleset,
        max_paths=max_paths,
        **kwargs,
    )


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="xenosite-forest",
        description="Find metabolic pathways between reactant and product SMILES.",
    )
    parser.add_argument(
        "molecules",
        nargs="+",
        help="Reactant SMILES, or reactant and product SMILES.",
    )
    parser.add_argument(
        "-a",
        "--all-paths",
        dest="all_paths",
        action="store_true",
        help="Output all valid paths, instead of just the first one.",
    )
    parser.add_argument(
        "-b",
        "--do-not-tag-atoms",
        dest="do_not_tag_atoms",
        action="store_true",
        help="Do not record atom-index tags on products.",
    )
    parser.add_argument(
        "-d",
        "--depth",
        default=1,
        type=int,
        help="Maximum search depth.",
    )
    parser.add_argument(
        "-e",
        "--phase1",
        action="store_true",
        help="Report sites as 0-based atom indexes.",
    )
    parser.add_argument(
        "-m",
        "--max",
        dest="max_paths",
        type=int,
        help="Stop after returning this many paths.",
    )
    parser.add_argument(
        "-r",
        "--ruleset",
        default="Full",
        help="Ruleset to use (default: Full).",
    )
    parser.add_argument(
        "-y",
        "--termination-ruleset",
        dest="termination_ruleset",
        help="Optional ruleset whose reactions terminate the search.",
    )
    parser.add_argument(
        "--expand-star-conjugates",
        dest="expand_star_conjugates",
        action="store_true",
        help="Allow further metabolism of star (*) conjugate adducts (default: off).",
    )
    parser.add_argument(
        "--search",
        choices=("bfs", "dfs"),
        default="bfs",
        help="Pathway search order (default: bfs).",
    )

    args = parser.parse_args(argv)
    kwargs = {
        "ruleset": args.ruleset,
        "termination_ruleset": args.termination_ruleset,
        "all_paths": args.all_paths,
        "do_not_tag_atoms": args.do_not_tag_atoms,
        "depth": args.depth,
        "phase1": args.phase1,
        "expand_star_conjugates": args.expand_star_conjugates,
        "search": args.search,
        "max_paths": args.max_paths,
    }

    search = dfs if args.search == "dfs" else bfs
    for smi, rules_and_sites, _mols in search(args.molecules, **kwargs):
        sys.stdout.write(str((smi, rules_and_sites)) + "\n")


if __name__ == "__main__":
    main()
