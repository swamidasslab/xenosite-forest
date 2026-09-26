"""Molecule corpus for Rust↔RDKit leaf-rule product parity fuzz.

Large, diverse SMILES chosen so every PatternInfo possibility — including
every ``when`` arm — has at least one covering mol. Completeness is owned by
:mod:`tests.forest.test_rule_parity_corpus` (meta-test first; grow this list
until that test is green; no xfails for missing cover).

Sources: :data:`substrate_library.SUBSTRATE_LIBRARY`, the rare-OR extras used
by pattern-info coverage, plus quaternary epoxide/aziridine carbons that hit
``Glutathionation`` ``epoxide_c`` / ``aziridine_c``.
"""

from __future__ import annotations

from .substrate_library import SUBSTRATE_LIBRARY

# Rare OR / heteroatom / halide / At branches beyond SUBSTRATE_LIBRARY.
_PATTERN_EXTRAS: tuple[str, ...] = (
    "OCN",
    "OCOC",
    "OCS",
    "OC(C)N(C)C",
    "NC(O)C",
    "OCSC",
    "O=CC=Cc1ccccc1",
    "N=CC=Cc1ccccc1",
    "CC=CC=N",
    "CC1CN1",
    "CC1CN1C",
    "c1ccccc1C1CN1",
    "ClC1CN1",
    "O=Nc1ccccc1",
    "CCN=O",
    "COO",
    "COOC",
    "CS(=O)O",
    "CS(O)=O",
    "C1=CC2OC2C=C1",
    "CF",
    "CCl",
    "CBr",
    "CI",
    "CCF",
    "CCBr",
    "CCI",
    "CC(F)C",
    "CC(Cl)C",
    "CC(Br)C",
    "CC(I)C",
    "FC(F)F",
    "BrCBr",
    "ICI",
    "ClC(Cl)Cl",
    "IC(I)C",
    "IC(Cl)C",
    "ClC(I)Cl",
    "ClCC=C",
    "BrCC=C",
    "ICC=C",
    "[At]CC=C",
    "CC[At]",
    "[At]CC",
    "[At]C(C)C",
    "[At]Cc1ccccc1",
    "[At]C[At]",
    "[At]C(Cl)[At]",
    "ClC([At])Cl",
    "[At]C([At])C",
    "[At]CI",
    "c1ccccc1F",
    "c1ccccc1Cl",
    "c1ccccc1Br",
    "c1ccccc1I",
    "Fc1ccc(O)cc1",
    "Brc1ccc(O)cc1",
    "Ic1ccc(O)cc1",
    "ICc1ccccc1",
    "FCc1ccccc1",
    "CNc1ccccc1",
    "CSc1ccccc1",
    "CCOc1ccccc1",
    "CCNc1ccccc1",
    "CCSc1ccccc1",
    "CS",
    "c1ccccc1S",
    "CCN(C)C",
    "CC(=O)NC",
    "CC(=O)SC",
    "CC(=S)OC",
    "CC(=S)NC",
    "CNc1ccc(O)cc1",
    "CN(C)c1ccc(O)cc1",
    "CCCC",
    "S=C=Nc1ccccc1",
    # Quaternary ring C (H0 + non-H substituent) for GSH epoxide_c / aziridine_c.
    "CC1(C)OC1",
    "CC1(C)NC1",
)

PARITY_FUZZ_MOLS: tuple[str, ...] = tuple(
    dict.fromkeys((*SUBSTRATE_LIBRARY, *_PATTERN_EXTRAS))
)
