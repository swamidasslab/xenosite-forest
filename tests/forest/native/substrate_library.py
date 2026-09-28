"""Shared substrate SMILES for forest coverage and canonical-plan checks.

Collected from forest fuzz / SMARTS-coverage probes (phase1 steps, formula
hints, bfs crashers, conjugation probes). Import this module; do not copy
the list into other test files.
"""

from __future__ import annotations

# Phase-I / quinone / conjugation probes plus a few drug-like crashers.
# Designed to hit the SMARTS families the forest rules declare.
SUBSTRATE_LIBRARY: tuple[str, ...] = (
    # Formula / Phase I SMARTS probes
    "CCO",
    "CCN",
    "CCS",
    "C=C",
    "C#C",
    "C1OC1",
    "c1ccccc1",
    "Oc1ccccc1",
    "CCCl",
    "ClCCl",
    "CC(=O)OC",
    "CCNO",
    "CCSO",
    "CS(=O)C",
    "CN(C)C",
    "c1ccccc1C1OC1",
    "O=C1C=CC(=O)C=C1",
    "CC(=O)Nc1ccc(O)cc1",
    "COP(=O)(O)O",
    "[O-][N+](=O)c1ccccc1",
    "c1ccccc1SSc1ccccc1",
    # Quinone / aromatic coverage
    "Oc1ccc(O)cc1",
    "Nc1ccc(O)cc1",
    "Clc1ccc(O)cc1",
    "COc1ccc(O)cc1",
    "Cc1ccc(O)cc1",
    "Clc1ccccc1",
    "CN(C)c1ccccc1",
    "COc1ccccc1",
    "c1ccc2ccccc2c1",
    "C1=CC=CC2=C1C=C(C=C2)CC3=CC=CC(=C3)C=C",
    "c1ccccc1CCCCc2ccccc2CCCCc3ccccc3",
    # Conjugation / epoxide / halide / Michael
    "c1ccccc1C(=O)O",
    "c1ccccc1C1CO1",
    "ClCc1ccccc1",
    "BrCc1ccccc1",
    "CCS",
    "C=C",
    "CC=O",
    "c1ccccc1N1CC1",
    "COS(=O)(=O)C",
    "O=C=Nc1ccccc1",
    "Nc1ccccc1",
    "N=C=Nc1ccccc1",
    "CC(=O)Nc1ccccc1",
    # Reductions / heteroatoms / drugs (bfs corpus)
    "C=Cc1ccccc1",
    "CC(=O)Oc1ccccc1C(=O)O",
    "c1ccsc1",
    "c1ccc2[nH]ccc2c1",
    "CN(C)CCOC(c1ccccc1)c1ccccc1",
    "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
    "O=C(NCC(O)c1ccc([N+](=O)[O-])cc1)C(Cl)Cl",
    "Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1",
    # SMARTS cohorts under-represented above
    "CSC",
    "c1ccc2c(c1)OCO2",
    "OC(=O)c1cc(/N=N/c2ccc(c(c2)C(=O)O)O)ccc1O",
)

# Small slice for @pytest.mark.parametrize. Every entry is in SUBSTRATE_LIBRARY.
_QUICK = frozenset(
    {
        "CCO",
        "CCN",
        "CCS",
        "C=C",
        "C#C",
        "C1OC1",
        "c1ccccc1",
        "Oc1ccccc1",
        "Oc1ccc(O)cc1",
        "CCCl",
        "COc1ccccc1",
        "CSC",
        "CN(C)C",
        "CC(=O)OC",
        "COP(=O)(O)O",
        "[O-][N+](=O)c1ccccc1",
        "Clc1ccccc1",
        "Nc1ccccc1",
    }
)
QUICK_SUBSTRATES: tuple[str, ...] = tuple(
    smiles for smiles in SUBSTRATE_LIBRARY if smiles in _QUICK
)
