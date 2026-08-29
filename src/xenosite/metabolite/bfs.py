"""Command line tool for performing BFS between a reactant and putative product."""

# Standard Library
import argparse
import sys


# Third Party
from .base import load
from .phaseone import PhaseOneRS
#from .report import AnnotatedReactions
from .rulesets import RULESETS, load_ruleset
from rdkit import Chem, rdBase

RULESETS['PhaseOneRS'] = PhaseOneRS

# Prevents spammy rdkit messages
rdBase.DisableLog('rdApp.*')


def bfs(molstrings, ruleset='Full', termination_ruleset=None, **kwargs):
    """Run bfs to find paths linking reactants and an optional putative product.

    Args:
        - outmols: output rdmols and sites in addition to string
    """

    inputs = load(molstrings)

    if None in inputs:
        raise ValueError('Problem with input.')

    RST = load_ruleset(ruleset)

    if termination_ruleset is not None:
        TRS = load_ruleset(termination_ruleset)
        termination_rulenames = TRS.rulenames
    else:
        termination_rulenames = []

    if len(inputs) > 2:
        for num, mol in enumerate(inputs):
            if not mol.HasProp('_Name'):
                mol.SetProp('_Name', 'Molecule%d' % num)
            for line in RST.find_path(
                    mol, termination_rulenames=termination_rulenames,
                    **kwargs):
                yield line

    else:

        for line in RST.find_path(
                *inputs, termination_rulenames=termination_rulenames,
                **kwargs):
            yield line


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='metfor')

    parser.add_argument(
        '-a',
        '--all_paths',
        action="store_true",
        help="Output all valid paths, instead of just the first one.")

    parser.add_argument('-b', '--do_not_tag_atoms', action='store_true')

    parser.add_argument(
        '-c', '--clean', action="store_true", help="do not circle SOMs")

    parser.add_argument(
        '-d', '--depth', default=1, type=int, help="The maximum search depth.")

    parser.add_argument(
        '-e',
        '--phase1',
        action='store_true',
        help=
        "Instead of outputting sites as frozensets, output Phase I site strings."
    )

    parser.add_argument(
        '-g',
        '--align',
        action='store_true',
        help='Align each reaction based on maximum common substructure.')

    parser.add_argument('-i', '--sites')

    parser.add_argument(
        '-l', '--limit_to_phase1_sites_in_sdf', action='store_true')

    parser.add_argument(
        '-m',
        '--max',
        type=int,
        help="Terminate after returning a specified number of metabolites")

    parser.add_argument(
        '-n',
        '--numbered',
        action="store_true",
        help="Intersperse atom index numbered images.")

    parser.add_argument(
        '-o', '--osdf', default='', help="File to which to save path(s)")

    parser.add_argument(
        '-r',
        '--ruleset',
        default='Full',
        help="Specifies the ruleset to be used.")

    parser.add_argument('-t', '--atom_paths', action="store_true")

    parser.add_argument('-u', '--metabolite_table', default='')

    parser.add_argument(
        '-w', '--row_num', help="add row numbers to html report")

    parser.add_argument(
        '-x', '--html', help='html file to which to write report.')

    parser.add_argument('-y', '--termination_ruleset')

    parser.add_argument(
        '-z',
        '--pdf_dir',
        help='Directory in which to save individual molecule PDFs.')

    # parser.add_argument(
    # '-p',
    # '--pdf_dir',
    # action="store_true",
    # help='Directory in which to save individual molecule PDFs.')

    ####################################################################################################

    opts, args = parser.parse_known_args()

    visualize = False
    if opts.pdf_dir or opts.html or opts.osdf:

        visualize = True
        #AR = AnnotatedReactions(**vars(opts))

    for rxnnum, (smi, rules_and_sites, mols) in enumerate(
            bfs(args, **vars(opts))):

        if opts.max and rxnnum >= opts.max:
            break

        if visualize:
            AR.add(
                mols,
                rules_and_sites=rules_and_sites,
                rxnnum=rxnnum,
                **vars(opts))

        sys.stdout.write(str((smi, rules_and_sites)) + '\n')

    if visualize:
        AR(**vars(opts))
