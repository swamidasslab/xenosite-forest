"""Molecule corpus for Rust↔RDKit leaf-rule product parity fuzz.

Large, diverse SMILES chosen so every PatternInfo possibility — including
every ``when`` arm — has at least one covering mol. Completeness is owned by
:mod:`tests.forest.test_rule_parity_corpus` (meta-test first; grow this list
until that test is green; no xfails for missing cover).

Also includes ResonancePair **close-end / identical-partner** substrates
(``_PAIR_CLOSE_EXTRAS``): ortho catechols, crowded ethers, adjacent alkyls,
neighboring path_ends. Those geometries often need special handling; Rust
is usually more correct, Python sometimes ignores them.

Sources: :data:`substrate_library.SUBSTRATE_LIBRARY`, the rare-OR extras used
by pattern-info coverage, quaternary epoxide/aziridine carbons that hit
``Glutathionation`` ``epoxide_c`` / ``aziridine_c``, and pair-close extras.
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

# ResonancePair ends that land close (graph distance ≤ 2) and/or share the
# same partner role on both ends (e.g. phenol×phenol). These often need
# special unique-edit / path handling; Rust is usually stricter/correct —
# Python sometimes ignores the close/identical-partner case. Owned by
# ``test_parity_fuzz_mols_cover_close_pair_ends``.
_PAIR_CLOSE_EXTRAS: tuple[str, ...] = (
    "Oc1ccccc1O",  # catechol ortho — identical partner O, short path
    "Oc1c(O)cccc1",  # catechol ortho (explicit)
    "Oc1cc(O)ccc1",  # resorcinol meta
    "Nc1ccccc1N",  # ortho diamine — identical partner N
    "Nc1ccc(N)cc1",  # para diamine
    "Nc1c(N)cccc1",
    "COc1ccccc1OC",  # ortho dimethoxy — QF dealkylate ends close
    "COc1ccc(OC)c(OC)c1",  # crowded methoxy
    "Clc1ccccc1Cl",  # ortho dihalo QF
    "Brc1ccc(Br)cc1",
    "CCc1ccccc1CC",  # ortho diethyl — DH alkyl ends close
    "CCc1c(C)cccc1",  # adjacent ethyl/methyl
    "Cc1c(C)cccc1",  # ortho xylene methide/alkyl
    "CC(=C)c1ccccc1C(=C)C",  # ortho bis-methide
    "O=CC=O",  # glyoxal — H path_end neighbors
    "O=CC=CC=O",
    # Bridging heteroatom partners (C14 exclusive_partner): both sides can
    # match; sharing that N/O partner atom must be refused by effect data.
    "c1ccc(N(C)c2ccccc2)cc1",  # N-methyl diphenylamine — shared N
    "c1ccc2c(c1)Nc1ccccc1C2",  # dihydroacridine — bridging NH
    "c1ccc2c(c1)Nc1ccccc1O2",  # phenoxazine — bridging N (+ O)
    "[nH]1cccc1",  # pyrrole — aromatic heteroatom `#` / specialize `n`
)

PARITY_FUZZ_MOLS: tuple[str, ...] = tuple(
    dict.fromkeys((*SUBSTRATE_LIBRARY, *_PATTERN_EXTRAS, *_PAIR_CLOSE_EXTRAS))
)
