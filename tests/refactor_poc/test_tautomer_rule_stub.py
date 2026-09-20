"""TautomerRule stub raises NotImplementedError."""

from __future__ import annotations

import pytest
from rdkit import Chem

from xenosite.refactor_poc.rules import TautomerRule


def test_tautomer_rule_raises():
    with pytest.raises(NotImplementedError, match="design stub"):
        list(TautomerRule().metabolites(Chem.MolFromSmiles("O=C1CCCCC1")))
