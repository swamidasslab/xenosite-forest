"""Find metabolic pathways between a reactant and an optional product."""

import argparse
import sys

from rdkit import rdBase

from .phaseone import PhaseOneRS
from .rulesets import RULESETS, load_ruleset
from .utils import load

RULESETS["PhaseOneRS"] = PhaseOneRS

rdBase.DisableLog("rdApp.*")


def bfs(molstrings, ruleset="Full", termination_ruleset=None, **kwargs):
    """Run BFS to find paths linking reactants and an optional putative product.

    Args:
        molstrings: SMILES strings, RDKit mols, or a mix.
        ruleset: named ruleset, RuleSet instance, or list of those.
        termination_ruleset: optional ruleset whose rule names stop expansion.
        **kwargs: forwarded to RuleSet.find_path (depth, phase1, all_paths, ...).
    """
    inputs = load(molstrings)

    if None in inputs:
        raise ValueError("Problem with input.")

    rules = load_ruleset(ruleset)

    if termination_ruleset is not None:
        termination_rulenames = load_ruleset(termination_ruleset).rulenames
    else:
        termination_rulenames = []

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
        help="Format sites as Phase I strings instead of frozensets.",
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

    args = parser.parse_args(argv)
    kwargs = {
        "ruleset": args.ruleset,
        "termination_ruleset": args.termination_ruleset,
        "all_paths": args.all_paths,
        "do_not_tag_atoms": args.do_not_tag_atoms,
        "depth": args.depth,
        "phase1": args.phase1,
    }

    for rxnnum, (smi, rules_and_sites, _mols) in enumerate(bfs(args.molecules, **kwargs)):
        if args.max_paths is not None and rxnnum >= args.max_paths:
            break
        sys.stdout.write(str((smi, rules_and_sites)) + "\n")


if __name__ == "__main__":
    main()
