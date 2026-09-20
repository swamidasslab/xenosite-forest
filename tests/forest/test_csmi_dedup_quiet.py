"""SiteDeduplicationWarning only for true unique-edit misses.

Isomorphic cleavage siblings from one site (two anilines from azobenzene;
two methanols from ethanol C–C dealkylation) and shared leaving groups across
distinct organics must not warn — unique-edit already did its job.
"""

from __future__ import annotations

import warnings

import pytest
from rdkit import Chem

from xenosite.forest.rules import (
    AzoSplitting,
    SiteDeduplicationWarning,
    Dealkylation,
    Dehydration,
    Dehydrogenation,
    NitrogenReduction,
    OxidativeDehalogenation,
    ReductiveDehalogenation,
)


@pytest.mark.parametrize(
    "rule_cls,smiles",
    [
        (AzoSplitting, "c1ccc(N=Nc2ccccc2)cc1"),
        (Dealkylation, "CCO"),
        (Dealkylation, "c1ccc2c(c1)OCO2"),
        # Distinct epoxide carbons → same diol CSMI (product iso, not orbit miss)
        (Dealkylation, "C1=CC2OC2c2ccc3nc4ccc5ccccc5c4cc3c21"),
        (NitrogenReduction, "[O-][N+](=O)c1ccccc1"),
        (OxidativeDehalogenation, "ClC(Cl)Cl"),
        (Dehydration, "CC(O)C(O)C"),
        (Dehydrogenation, "CC(O)C(O)C"),
        # Distinct C–Br sites → same alkene CSMI after reductive elimination
        (
            ReductiveDehalogenation,
            "Br/C/1=C/CCC(Br)C(Br)CCC(C(CC1)Br)Br",
        ),
    ],
)
def test_no_csmi_on_cleavage_or_symmetry(rule_cls: type, smiles: str) -> None:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", SiteDeduplicationWarning)
        list(rule_cls().metabolize(Chem.MolFromSmiles(smiles)))
    csmi = [w for w in caught if issubclass(w.category, SiteDeduplicationWarning)]
    assert not csmi, f"{rule_cls.__name__} {smiles}: {[str(w.message) for w in csmi]}"
