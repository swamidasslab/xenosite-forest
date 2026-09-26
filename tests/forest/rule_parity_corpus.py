"""Parity corpus: mol + intended rule/pattern/when correspondences.

Each :class:`ParityEntry` is ``(smiles, covers)``. ``covers`` is a tuple of
:class:`CoverIntent` (rule, pattern, possibility index, when key) — a list of
correspondences when one mol is the designated cover for several inventory
possibilities.

Meta-tests verify every intent and inventory completeness.
Parametric parity / CSMI suites use :func:`parity_param_cases`: full
rule×mol cartesian by default; focused CoverIntent-only under
``XENOSITE_PARITY_FULL=0`` / ``pytest --parity-focused`` (CI may disable).
"""

from __future__ import annotations

import os
from typing import NamedTuple


class CoverIntent(NamedTuple):
    """One inventory possibility this mol is intended to cover."""

    rule: str
    pattern: str
    poss_i: int
    when: tuple[int | None, int | None, int | None] | None


class ParityEntry(NamedTuple):
    """Corpus row: substrate SMILES + correspondence list."""

    smiles: str
    covers: tuple[CoverIntent, ...]


PARITY_CORPUS: tuple[ParityEntry, ...] = (
    ParityEntry(
        'CCN',
        (
            CoverIntent('Acetylation', 'acetyl', 0, (1, 7, None)),
            CoverIntent('ConjugationRule', 'acetyl', 0, (1, 7, None)),
            CoverIntent('Dealkylation', 'methylene_carboxylic', 0, (2, 7, None)),
            CoverIntent('Dealkylation', 'methylene_carbonyl', 0, (2, 7, None)),
            CoverIntent('Dealkylation', 'methylene_alcohol', 0, (2, 7, None)),
            CoverIntent('Dehydrogenation', 'amine', 0, (2, 7, 2)),
            CoverIntent('Dehydrogenation', 'amine_end', 0, (2, 7, 2)),
            CoverIntent('NDealkylation', 'methylene_carboxylic', 0, None),
            CoverIntent('NDealkylation', 'methylene_carbonyl', 0, None),
            CoverIntent('NDealkylation', 'methylene_alcohol', 0, None),
            CoverIntent('NitrogenOxidation', 'hydroxylamine', 1, (1, 7, 2)),
            CoverIntent('NitrogenOxidation', 'nitroso', 0, None),
        ),
    ),
    ParityEntry(
        'CCO',
        (
            CoverIntent('Acetylation', 'acetyl', 1, (1, 8, None)),
            CoverIntent('ConjugationRule', 'acetyl', 1, (1, 8, None)),
            CoverIntent('Dealkylation', 'cc_alcohol', 2, (1, 6, 3)),
            CoverIntent('Dealkylation', 'cc_carbonyl', 2, (1, 6, 3)),
            CoverIntent('Dehydration', 'alcohol', 0, (1, 6, None)),
            CoverIntent('Dehydration', 'beta_elimination', 0, None),
            CoverIntent('Dehydrogenation', 'alcohol', 0, None),
            CoverIntent('Dehydrogenation', 'alkyl', 1, (2, 6, 2)),
            CoverIntent('Dehydrogenation', 'phenol_end', 0, None),
            CoverIntent('Dehydrogenation', 'methide_end', 1, (2, 6, 2)),
            CoverIntent('Glucuronidation', 'alcohol', 0, (2, 6, None)),
            CoverIntent('Hydrogenation', 'path_end', 0, None),
            CoverIntent('Hydroxylation', 'h2', 0, (1, 6, 2)),
            CoverIntent('Hydroxylation', 'h2', 1, (1, 6, 3)),
            CoverIntent('Sulfation', 'alcohol', 0, (1, 6, None)),
        ),
    ),
    ParityEntry(
        'CCS',
        (
            CoverIntent('Acetylation', 'acetyl', 2, (1, 16, None)),
            CoverIntent('ConjugationRule', 'acetyl', 2, (1, 16, None)),
            CoverIntent('Dealkylation', 'methylene_carboxylic', 2, (2, 16, None)),
            CoverIntent('Dealkylation', 'methylene_carbonyl', 2, (2, 16, None)),
            CoverIntent('Dealkylation', 'methylene_alcohol', 2, (2, 16, None)),
            CoverIntent('Glutathionation', 'thiol', 0, None),
            CoverIntent('SulfurOxidation', 'zwitterion', 0, None),
            CoverIntent('SulfurOxidation', 'hydroxy', 0, None),
            CoverIntent('SulfurOxidation', 'oxo', 0, None),
            CoverIntent('SulfurReduction', 'thioether', 0, (2, 6, None)),
        ),
    ),
    ParityEntry(
        'OC(=O)c1cc(/N=N/c2ccc(c(c2)C(=O)O)O)ccc1O',
        (
            CoverIntent('AzoSplitting', 'azo', 0, None),
        ),
    ),
    ParityEntry(
        'c1ccc2c(c1)OCO2',
        (
            CoverIntent('BenzodioxoleReduction', 'dioxole_methylene', 0, None),
        ),
    ),
    ParityEntry(
        'CN(C)C',
        (
            CoverIntent('Dealkylation', 'methyl_carboxylic', 0, (2, 7, None)),
            CoverIntent('Dealkylation', 'methyl_carbonyl', 0, (2, 7, None)),
            CoverIntent('Dealkylation', 'methyl_alcohol', 0, (2, 7, None)),
            CoverIntent('NDealkylation', 'methyl_carboxylic', 0, None),
            CoverIntent('NDealkylation', 'methyl_carbonyl', 0, None),
            CoverIntent('NDealkylation', 'methyl_alcohol', 0, None),
            CoverIntent('NitrogenOxidation', 'n_oxide', 0, None),
        ),
    ),
    ParityEntry(
        'CC(=O)OC',
        (
            CoverIntent('Dealkylation', 'methyl_carboxylic', 1, (2, 8, None)),
            CoverIntent('Dealkylation', 'methyl_carbonyl', 1, (2, 8, None)),
            CoverIntent('Dealkylation', 'methyl_alcohol', 1, (2, 8, None)),
            CoverIntent('Dealkylation', 'quaternary_alcohol', 1, (2, 8, None)),
            CoverIntent('Dehydration', 'carbonyl', 0, (1, 6, None)),
            CoverIntent('Dehydrogenation', 'methide_end', 0, (2, 6, 3)),
            CoverIntent('Hydrolysis', 'add_water', 1, (3, 8, None)),
            CoverIntent('Hydrolysis', 'cleave', 1, (3, 8, None)),
            CoverIntent('OxygenReduction', 'carbonyl', 0, (2, 6, None)),
        ),
    ),
    ParityEntry(
        'CS(=O)C',
        (
            CoverIntent('Dealkylation', 'methyl_carboxylic', 2, (2, 16, None)),
            CoverIntent('Dealkylation', 'methyl_carbonyl', 2, (2, 16, None)),
            CoverIntent('Dealkylation', 'methyl_alcohol', 2, (2, 16, None)),
            CoverIntent('SulfurReduction', 'sulfoxide', 0, None),
        ),
    ),
    ParityEntry(
        'C1OC1',
        (
            CoverIntent('Dealkylation', 'methylene_carboxylic', 1, (2, 8, None)),
            CoverIntent('Dealkylation', 'methylene_carbonyl', 1, (2, 8, None)),
            CoverIntent('Dealkylation', 'methylene_alcohol', 1, (2, 8, None)),
            CoverIntent('Dealkylation', 'cc_alcohol', 1, (1, 6, 2)),
            CoverIntent('Dealkylation', 'cc_carbonyl', 1, (1, 6, 2)),
            CoverIntent('EpoxideOpening', 'rearrange', 0, None),
            CoverIntent('EpoxideOpening', 'hydrate', 0, None),
            CoverIntent('Glutathionation', 'epoxide_ch2', 0, None),
            CoverIntent('QuinoneFormation', 'single_to_double', 4, (2, 6, 2)),
            CoverIntent('QuinoneFormation', 'dealkylate', 1, (2, 8, None)),
        ),
    ),
    ParityEntry(
        'c1ccc2[nH]ccc2c1',
        (
            CoverIntent('Dealkylation', 'methine_carbonyl', 0, (2, 7, None)),
            CoverIntent('Dealkylation', 'methine_alcohol', 0, (2, 7, None)),
            CoverIntent('NDealkylation', 'methine_carbonyl', 0, None),
            CoverIntent('NDealkylation', 'methine_alcohol', 0, None),
        ),
    ),
    ParityEntry(
        'c1ccccc1C1OC1',
        (
            CoverIntent('Dealkylation', 'methine_carbonyl', 1, (2, 8, None)),
            CoverIntent('Dealkylation', 'methine_alcohol', 1, (2, 8, None)),
            CoverIntent('Dehydrogenation', 'methide_end', 2, (2, 6, 1)),
            CoverIntent('Glutathionation', 'epoxide_ch', 0, None),
            CoverIntent('QuinoneFormation', 'single_to_double', 5, (2, 6, 1)),
        ),
    ),
    ParityEntry(
        'c1ccsc1',
        (
            CoverIntent('Dealkylation', 'methine_carbonyl', 2, (2, 16, None)),
            CoverIntent('Dealkylation', 'methine_alcohol', 2, (2, 16, None)),
            CoverIntent('ThiopheneSulfurOxidation', 'thiophene_s_oxide', 0, None),
        ),
    ),
    ParityEntry(
        'CC(=O)Nc1ccc(O)cc1',
        (
            CoverIntent('Dealkylation', 'quaternary_alcohol', 0, (2, 7, None)),
            CoverIntent('Hydrolysis', 'add_water', 0, (3, 7, None)),
            CoverIntent('Hydrolysis', 'cleave', 0, (3, 7, None)),
            CoverIntent('NDealkylation', 'quaternary_alcohol', 0, None),
            CoverIntent('QuinoneFormation', 'single_to_double', 2, (2, 7, 1)),
            CoverIntent('QuinoneFormation', 'dealkylate', 0, (2, 7, None)),
        ),
    ),
    ParityEntry(
        'c1ccccc1SSc1ccccc1',
        (
            CoverIntent('Dealkylation', 'quaternary_alcohol', 2, (2, 16, None)),
            CoverIntent('SulfurReduction', 'disulfide', 0, None),
        ),
    ),
    ParityEntry(
        'Oc1ccccc1',
        (
            CoverIntent('Dealkylation', 'cc_quaternary_alcohol', 0, None),
            CoverIntent('QuinoneFormation', 'single_to_double', 0, (2, 8, None)),
        ),
    ),
    ParityEntry(
        'c1ccccc1',
        (
            CoverIntent('Dealkylation', 'cc_alcohol', 0, (1, 6, 1)),
            CoverIntent('Dealkylation', 'cc_carbonyl', 0, (1, 6, 1)),
            CoverIntent('QuinoneFormation', 'add_carbonyl_o', 0, None),
        ),
    ),
    ParityEntry(
        'OCN',
        (
            CoverIntent('Dealkylation', 'hemiaminal', 0, (2, 7, None)),
            CoverIntent('NDealkylation', 'hemiaminal', 0, None),
        ),
    ),
    ParityEntry(
        'OCOC',
        (
            CoverIntent('Dealkylation', 'hemiaminal', 1, (2, 8, None)),
        ),
    ),
    ParityEntry(
        'OCS',
        (
            CoverIntent('Dealkylation', 'hemiaminal', 2, (2, 16, None)),
        ),
    ),
    ParityEntry(
        'CCNO',
        (
            CoverIntent('Dehydration', 'alcohol', 1, (1, 7, None)),
            CoverIntent('Dehydrogenation', 'amine', 1, (2, 7, 1)),
            CoverIntent('Dehydrogenation', 'amine_end', 1, (2, 7, 1)),
            CoverIntent('NitrogenOxidation', 'hydroxylamine', 0, (1, 7, 1)),
            CoverIntent('NitrogenReduction', 'hydroxylamine', 0, None),
        ),
    ),
    ParityEntry(
        '[O-][N+](=O)c1ccccc1',
        (
            CoverIntent('Dehydration', 'carbonyl', 1, (1, 7, None)),
            CoverIntent('NitroaromaticReduction', 'nitro_charged', 0, None),
            CoverIntent('NitroaromaticReduction', 'nitro_neutral', 0, None),
            CoverIntent('NitrogenReduction', 'nitro_charged', 0, None),
            CoverIntent('NitrogenReduction', 'nitro_anion', 0, None),
            CoverIntent('NitrogenReduction', 'nitro_neutral', 0, None),
            CoverIntent('NitrogenReduction', 'nitro_to_amine', 0, None),
            CoverIntent('NitrogenReduction', 'nitro_both', 0, None),
            CoverIntent('NitrogenReduction', 'nitro_both_any', 0, None),
            CoverIntent('OxygenReduction', 'carbonyl', 1, (2, 7, None)),
            CoverIntent('QuinoneFormation', 'iminium', 0, None),
        ),
    ),
    ParityEntry(
        'CS(=O)O',
        (
            CoverIntent('Dehydrogenation', 'sulfoxide', 0, None),
        ),
    ),
    ParityEntry(
        'CC=O',
        (
            CoverIntent('Dehydrogenation', 'alkyl', 0, (2, 6, 3)),
            CoverIntent('Glutathionation', 'carbonyl', 0, None),
        ),
    ),
    ParityEntry(
        'CC(C)Cc1ccc(C(C)C(=O)O)cc1',
        (
            CoverIntent('Dehydrogenation', 'alkyl', 2, (2, 6, 1)),
        ),
    ),
    ParityEntry(
        'COP(=O)(O)O',
        (
            CoverIntent('Dephosphorylation', 'phosphate_ester', 0, (2, 15, None)),
        ),
    ),
    ParityEntry(
        'C=C',
        (
            CoverIntent('Epoxidation', 'epoxide', 0, (2, 6, None)),
            CoverIntent('Glutathionation', 'alkene', 0, None),
            CoverIntent('Hydrogenation', 'alkene', 0, None),
        ),
    ),
    ParityEntry(
        'O=C=Nc1ccccc1',
        (
            CoverIntent('Epoxidation', 'epoxide', 1, (2, 7, None)),
            CoverIntent('Glutathionation', 'isocyanate', 0, (3, 8, None)),
        ),
    ),
    ParityEntry(
        'c1ccccc1C(=O)O',
        (
            CoverIntent('Glucuronidation', 'carboxylate', 0, None),
        ),
    ),
    ParityEntry(
        'CC1(C)OC1',
        (
            CoverIntent('Glutathionation', 'epoxide_c', 0, None),
        ),
    ),
    ParityEntry(
        'CF',
        (
            CoverIntent('Glutathionation', 'halide', 0, (2, 9, None)),
            CoverIntent('OxidativeDehalogenation', 'alcohol', 0, (1, 9, None)),
            CoverIntent('ReductiveDehalogenation', 'cleave', 0, (1, 9, None)),
        ),
    ),
    ParityEntry(
        'CCCl',
        (
            CoverIntent('Glutathionation', 'halide', 1, (2, 17, None)),
            CoverIntent('OxidativeDehalogenation', 'alcohol', 1, (1, 17, None)),
            CoverIntent('OxidativeDehalogenation', 'carboxylic', 1, (1, 17, None)),
            CoverIntent('ReductiveDehalogenation', 'cleave', 1, (1, 17, None)),
            CoverIntent('ReductiveDehalogenation', 'alkene', 1, (1, 17, None)),
        ),
    ),
    ParityEntry(
        'BrCc1ccccc1',
        (
            CoverIntent('Glutathionation', 'halide', 2, (2, 35, None)),
            CoverIntent('OxidativeDehalogenation', 'alcohol', 2, (1, 35, None)),
            CoverIntent('OxidativeDehalogenation', 'carboxylic', 2, (1, 35, None)),
            CoverIntent('ReductiveDehalogenation', 'cleave', 2, (1, 35, None)),
            CoverIntent('ReductiveDehalogenation', 'alkene', 2, (1, 35, None)),
        ),
    ),
    ParityEntry(
        'CI',
        (
            CoverIntent('Glutathionation', 'halide', 3, (2, 53, None)),
            CoverIntent('OxidativeDehalogenation', 'alcohol', 3, (1, 53, None)),
            CoverIntent('ReductiveDehalogenation', 'cleave', 3, (1, 53, None)),
        ),
    ),
    ParityEntry(
        'O=C1C=CC(=O)C=C1',
        (
            CoverIntent('Glutathionation', 'michael', 0, (4, 8, None)),
        ),
    ),
    ParityEntry(
        'N=CC=Cc1ccccc1',
        (
            CoverIntent('Glutathionation', 'michael', 1, (4, 7, None)),
        ),
    ),
    ParityEntry(
        'CC1CN1',
        (
            CoverIntent('Glutathionation', 'aziridine_ch', 0, None),
        ),
    ),
    ParityEntry(
        'c1ccccc1N1CC1',
        (
            CoverIntent('Glutathionation', 'aziridine_ch2', 0, None),
        ),
    ),
    ParityEntry(
        'CC1(C)NC1',
        (
            CoverIntent('Glutathionation', 'aziridine_c', 0, None),
        ),
    ),
    ParityEntry(
        'COS(=O)(=O)C',
        (
            CoverIntent('Glutathionation', 'mesylate', 0, None),
        ),
    ),
    ParityEntry(
        'S=C=Nc1ccccc1',
        (
            CoverIntent('Glutathionation', 'isocyanate', 1, (3, 16, None)),
        ),
    ),
    ParityEntry(
        'C#C',
        (
            CoverIntent('Hydrogenation', 'alkyne', 0, None),
            CoverIntent('Hydroxylation', 'h', 0, (1, 6, 1)),
        ),
    ),
    ParityEntry(
        'CC(=O)SC',
        (
            CoverIntent('Hydrolysis', 'add_water', 2, (3, 16, None)),
            CoverIntent('Hydrolysis', 'cleave', 2, (3, 16, None)),
        ),
    ),
    ParityEntry(
        'O=Nc1ccccc1',
        (
            CoverIntent('NitrogenReduction', 'nitroso', 0, None),
        ),
    ),
    ParityEntry(
        '[At]CC=C',
        (
            CoverIntent('OxidativeDehalogenation', 'alcohol', 4, (1, 85, None)),
            CoverIntent('OxidativeDehalogenation', 'carboxylic', 4, (1, 85, None)),
            CoverIntent('OxidativeDehalogenation', 'rearrange', 4, (1, 85, None)),
            CoverIntent('ReductiveDehalogenation', 'cleave', 4, (1, 85, None)),
            CoverIntent('ReductiveDehalogenation', 'alkene', 4, (1, 85, None)),
        ),
    ),
    ParityEntry(
        'CC(F)C',
        (
            CoverIntent('OxidativeDehalogenation', 'carbonyl', 0, (1, 9, None)),
        ),
    ),
    ParityEntry(
        'O=C(NCC(O)c1ccc([N+](=O)[O-])cc1)C(Cl)Cl',
        (
            CoverIntent('OxidativeDehalogenation', 'carbonyl', 1, (1, 17, None)),
        ),
    ),
    ParityEntry(
        'CC(Br)C',
        (
            CoverIntent('OxidativeDehalogenation', 'carbonyl', 2, (1, 35, None)),
        ),
    ),
    ParityEntry(
        'CC(I)C',
        (
            CoverIntent('OxidativeDehalogenation', 'carbonyl', 3, (1, 53, None)),
        ),
    ),
    ParityEntry(
        '[At]C(C)C',
        (
            CoverIntent('OxidativeDehalogenation', 'carbonyl', 4, (1, 85, None)),
        ),
    ),
    ParityEntry(
        'CCF',
        (
            CoverIntent('OxidativeDehalogenation', 'carboxylic', 0, (1, 9, None)),
            CoverIntent('ReductiveDehalogenation', 'alkene', 0, (1, 9, None)),
        ),
    ),
    ParityEntry(
        'CCI',
        (
            CoverIntent('OxidativeDehalogenation', 'carboxylic', 3, (1, 53, None)),
            CoverIntent('ReductiveDehalogenation', 'alkene', 3, (1, 53, None)),
        ),
    ),
    ParityEntry(
        'c1ccccc1F',
        (
            CoverIntent('OxidativeDehalogenation', 'rearrange', 0, (1, 9, None)),
            CoverIntent('QuinoneFormation', 'replace_halogen', 0, (2, 9, None)),
        ),
    ),
    ParityEntry(
        'Clc1ccc(O)cc1',
        (
            CoverIntent('OxidativeDehalogenation', 'rearrange', 1, (1, 17, None)),
            CoverIntent('QuinoneFormation', 'replace_halogen', 1, (2, 17, None)),
        ),
    ),
    ParityEntry(
        'BrCC=C',
        (
            CoverIntent('OxidativeDehalogenation', 'rearrange', 2, (1, 35, None)),
        ),
    ),
    ParityEntry(
        'ICC=C',
        (
            CoverIntent('OxidativeDehalogenation', 'rearrange', 3, (1, 53, None)),
        ),
    ),
    ParityEntry(
        'FC(F)F',
        (
            CoverIntent('OxidativeDehalogenation', 'gem_carboxylic', 0, (1, 9, None)),
            CoverIntent('OxidativeDehalogenation', 'gem_hydrate', 0, (1, 9, None)),
        ),
    ),
    ParityEntry(
        'ClCCl',
        (
            CoverIntent('OxidativeDehalogenation', 'gem_carboxylic', 1, (1, 17, None)),
            CoverIntent('OxidativeDehalogenation', 'gem_hydrate', 1, (1, 17, None)),
        ),
    ),
    ParityEntry(
        'BrCBr',
        (
            CoverIntent('OxidativeDehalogenation', 'gem_carboxylic', 2, (1, 35, None)),
            CoverIntent('OxidativeDehalogenation', 'gem_hydrate', 2, (1, 35, None)),
        ),
    ),
    ParityEntry(
        'ICI',
        (
            CoverIntent('OxidativeDehalogenation', 'gem_carboxylic', 3, (1, 53, None)),
            CoverIntent('OxidativeDehalogenation', 'gem_hydrate', 3, (1, 53, None)),
        ),
    ),
    ParityEntry(
        '[At]C[At]',
        (
            CoverIntent('OxidativeDehalogenation', 'gem_carboxylic', 4, (1, 85, None)),
            CoverIntent('OxidativeDehalogenation', 'gem_hydrate', 4, (1, 85, None)),
        ),
    ),
    ParityEntry(
        'COO',
        (
            CoverIntent('OxygenReduction', 'peroxide', 0, None),
        ),
    ),
    ParityEntry(
        'Nc1ccc(O)cc1',
        (
            CoverIntent('QuinoneFormation', 'single_to_double', 1, (2, 7, 2)),
        ),
    ),
    ParityEntry(
        'Cc1ccc(O)cc1',
        (
            CoverIntent('QuinoneFormation', 'single_to_double', 3, (2, 6, 3)),
        ),
    ),
    ParityEntry(
        'c1ccccc1Br',
        (
            CoverIntent('QuinoneFormation', 'replace_halogen', 2, (2, 35, None)),
        ),
    ),
    ParityEntry(
        'c1ccccc1I',
        (
            CoverIntent('QuinoneFormation', 'replace_halogen', 3, (2, 53, None)),
        ),
    ),
    ParityEntry(
        'C1=CC2OC2C=C1',
        (
            CoverIntent('Sulfation', 'epoxide_methyl_sulfone', 0, None),
        ),
    ),
    ParityEntry(
        'CCSO',
        (
            CoverIntent('SulfurReduction', 'thioether', 1, (2, 8, None)),
        ),
    ),
    ParityEntry('Oc1ccc(O)cc1', ()),
    ParityEntry('COc1ccc(O)cc1', ()),
    ParityEntry('Clc1ccccc1', ()),
    ParityEntry('CN(C)c1ccccc1', ()),
    ParityEntry('COc1ccccc1', ()),
    ParityEntry('c1ccc2ccccc2c1', ()),
    ParityEntry('C1=CC=CC2=C1C=C(C=C2)CC3=CC=CC(=C3)C=C', ()),
    ParityEntry('c1ccccc1CCCCc2ccccc2CCCCc3ccccc3', ()),
    ParityEntry('c1ccccc1C1CO1', ()),
    ParityEntry('ClCc1ccccc1', ()),
    ParityEntry('Nc1ccccc1', ()),
    ParityEntry('N=C=Nc1ccccc1', ()),
    ParityEntry('CC(=O)Nc1ccccc1', ()),
    ParityEntry('C=Cc1ccccc1', ()),
    ParityEntry('CC(=O)Oc1ccccc1C(=O)O', ()),
    ParityEntry('CN(C)CCOC(c1ccccc1)c1ccccc1', ()),
    ParityEntry('Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1', ()),
    ParityEntry('CSC', ()),
    ParityEntry('OC(C)N(C)C', ()),
    ParityEntry('NC(O)C', ()),
    ParityEntry('OCSC', ()),
    ParityEntry('O=CC=Cc1ccccc1', ()),
    ParityEntry('CC=CC=N', ()),
    ParityEntry('CC1CN1C', ()),
    ParityEntry('c1ccccc1C1CN1', ()),
    ParityEntry('ClC1CN1', ()),
    ParityEntry('CCN=O', ()),
    ParityEntry('COOC', ()),
    ParityEntry('CS(O)=O', ()),
    ParityEntry('CCl', ()),
    ParityEntry('CBr', ()),
    ParityEntry('CCBr', ()),
    ParityEntry('CC(Cl)C', ()),
    ParityEntry('ClC(Cl)Cl', ()),
    ParityEntry('IC(I)C', ()),
    ParityEntry('IC(Cl)C', ()),
    ParityEntry('ClC(I)Cl', ()),
    ParityEntry('ClCC=C', ()),
    ParityEntry('CC[At]', ()),
    ParityEntry('[At]CC', ()),
    ParityEntry('[At]Cc1ccccc1', ()),
    ParityEntry('[At]C(Cl)[At]', ()),
    ParityEntry('ClC([At])Cl', ()),
    ParityEntry('[At]C([At])C', ()),
    ParityEntry('[At]CI', ()),
    ParityEntry('c1ccccc1Cl', ()),
    ParityEntry('Fc1ccc(O)cc1', ()),
    ParityEntry('Brc1ccc(O)cc1', ()),
    ParityEntry('Ic1ccc(O)cc1', ()),
    ParityEntry('ICc1ccccc1', ()),
    ParityEntry('FCc1ccccc1', ()),
    ParityEntry('CNc1ccccc1', ()),
    ParityEntry('CSc1ccccc1', ()),
    ParityEntry('CCOc1ccccc1', ()),
    ParityEntry('CCNc1ccccc1', ()),
    ParityEntry('CCSc1ccccc1', ()),
    ParityEntry('CS', ()),
    ParityEntry('c1ccccc1S', ()),
    ParityEntry('CCN(C)C', ()),
    ParityEntry('CC(=O)NC', ()),
    ParityEntry('CC(=S)OC', ()),
    ParityEntry('CC(=S)NC', ()),
    ParityEntry('CNc1ccc(O)cc1', ()),
    ParityEntry('CN(C)c1ccc(O)cc1', ()),
    ParityEntry('CCCC', ()),
    ParityEntry('Oc1ccccc1O', ()),
    ParityEntry('Oc1c(O)cccc1', ()),
    ParityEntry('Oc1cc(O)ccc1', ()),
    ParityEntry('Nc1ccccc1N', ()),
    ParityEntry('Nc1ccc(N)cc1', ()),
    ParityEntry('Nc1c(N)cccc1', ()),
    ParityEntry('COc1ccccc1OC', ()),
    ParityEntry('COc1ccc(OC)c(OC)c1', ()),
    ParityEntry('Clc1ccccc1Cl', ()),
    ParityEntry('Brc1ccc(Br)cc1', ()),
    ParityEntry('CCc1ccccc1CC', ()),
    ParityEntry('CCc1c(C)cccc1', ()),
    ParityEntry('Cc1c(C)cccc1', ()),
    ParityEntry('CC(=C)c1ccccc1C(=C)C', ()),
    ParityEntry('O=CC=O', ()),
    ParityEntry('O=CC=CC=O', ()),
    ParityEntry('c1ccc(N(C)c2ccccc2)cc1', ()),
    ParityEntry('c1ccc2c(c1)Nc1ccccc1C2', ()),
    ParityEntry('c1ccc2c(c1)Nc1ccccc1O2', ()),
    ParityEntry('[nH]1cccc1', ()),
    ParityEntry('c1ccncc1', ()),
)


def parity_fuzz_mols() -> tuple[str, ...]:
    """Unique SMILES in ``PARITY_CORPUS`` (order preserved)."""

    return tuple(dict.fromkeys(e.smiles for e in PARITY_CORPUS))


# Back-compat alias for suites that still iterate mols alone.
PARITY_FUZZ_MOLS: tuple[str, ...] = parity_fuzz_mols()


def _filter_rules(
    rules: tuple[str, ...] | list[str] | None,
) -> set[str] | None:
    if rules is None:
        return None
    return set(rules)


def parity_cover_cases(
    rules: tuple[str, ...] | list[str] | None = None,
) -> tuple[tuple[str, str, str, int, object], ...]:
    """Every declared cover: ``(rule, smiles, pattern, poss_i, when)``.

    Parametrize meta / focused chemistry tests over this for guaranteed
    inventory coverage (one case per possibility when corpus is complete).
    Optional ``rules`` restricts to those rule names.
    """

    allowed = _filter_rules(rules)
    out: list[tuple[str, str, str, int, object]] = []
    for entry in PARITY_CORPUS:
        for cover in entry.covers:
            if allowed is not None and cover.rule not in allowed:
                continue
            out.append(
                (cover.rule, entry.smiles, cover.pattern, cover.poss_i, cover.when)
            )
    return tuple(out)


def parity_rule_mol_cases(
    rules: tuple[str, ...] | list[str] | None = None,
) -> tuple[tuple[str, str], ...]:
    """Focused ``(rule_name, smiles)`` from declared covers (deduped).

    Default mode for parametric parity / CSMI: every inventory possibility
    has ≥1 designated case; not the full rule×mol cartesian.
    """

    allowed = _filter_rules(rules)
    seen: dict[tuple[str, str], None] = {}
    for entry in PARITY_CORPUS:
        for cover in entry.covers:
            if allowed is not None and cover.rule not in allowed:
                continue
            seen.setdefault((cover.rule, entry.smiles), None)
    return tuple(seen)


def parity_rule_mol_cases_full(
    rules: tuple[str, ...] | list[str] | None = None,
) -> tuple[tuple[str, str], ...]:
    """Cartesian expansion: every selected rule × every corpus mol.

    Used when ``parity_full_enabled()`` — slower full sweep beyond designated
    covers. When ``rules`` is omitted, uses every rule that appears in at
    least one :class:`CoverIntent`.
    """

    if rules is None:
        rule_names = tuple(
            dict.fromkeys(c.rule for e in PARITY_CORPUS for c in e.covers)
        )
    else:
        rule_names = tuple(dict.fromkeys(rules))
    mols = parity_fuzz_mols()
    return tuple((rule, smiles) for rule in rule_names for smiles in mols)


def parity_full_enabled() -> bool:
    """True for full rule×mol coverage (default).

    Full cartesian is the local default so inventory coverage is guaranteed.
    Opt into focused CoverIntent-only with ``XENOSITE_PARITY_FULL=0`` /
    ``false`` / ``focused`` or ``pytest --parity-focused`` (CI may disable).
    """

    raw = os.environ.get("XENOSITE_PARITY_FULL", "1").strip().lower()
    if raw in {"0", "false", "no", "off", "focused"}:
        return False
    return True


def parity_param_cases(
    rules: tuple[str, ...] | list[str] | None = None,
) -> tuple[tuple[str, str], ...]:
    """Cases for parametric parity / CSMI: focused, or full when enabled."""

    if parity_full_enabled():
        return parity_rule_mol_cases_full(rules)
    return parity_rule_mol_cases(rules)
