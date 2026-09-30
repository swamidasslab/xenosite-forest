"""Rust↔native leaf product parity — same corpus cartesian as C18 baseline.

Port of ``tests/forest/test_rule_parity_fuzz.py`` from tag
``py-rust-parity-c18-baseline`` / ``a814b83``. Product identity is RDKit CSMI
(C1). Site-bag equality is omitted until ForestMol exposes ranks again; C18
asserted both.

Conjugation leaves compare **star** products (native default ``as_star=True``,
Rust always collapses). Full-adduct native (``as_star=False``) is inventory-only.

C18's 29 ``ProductParityXfail`` rows stay ``xfail(strict=True)``. New gaps
outside that set are hard fails — that is the HEAD vs baseline delta.
"""

from __future__ import annotations

from typing import Any

import pytest
from rdkit import Chem
from rdkit.Chem.rdmolops import RemoveStereochemistry

from xenosite.forest import available, leaf_rule, load
from xenosite.forest.native import rules as native_rules
from xenosite.forest.native.rdkit_api import MolFromSmiles, MolToSmiles
from xenosite.forest.native.rules import ReactionRule

from .rule_parity_corpus import (
    C18_OPEN_PRODUCT_XFAIL_COUNT,
    parity_full_enabled,
    parity_param_cases,
    parity_rule_mol_cases,
    parity_rule_mol_cases_full,
    product_parity_xfail_cases,
    product_parity_xfail_reason,
)
from .rule_parity_pairs import paired_rule_names

# Product parity while chasing C18: Effect/formula soft-check stays off in Rust
# materialize; do not fail the suite on recorded formula_delta_mismatch.
# SiteDeduplicationWarning stays an error suite-wide; this fuzz cares about
# product sets, not unique-edit completeness on native.
pytestmark = [
    pytest.mark.skipif(
        not available(),
        reason="xenosite.forest Rust extension required",
    ),
    pytest.mark.allow_formula_delta_mismatch,
    pytest.mark.filterwarnings("ignore::xenosite.forest.native.rules.SiteDeduplicationWarning"),
]

def _rdkit_csmi(smiles_or_mol: Any) -> str | None:
    if smiles_or_mol is None:
        return None
    if isinstance(smiles_or_mol, str):
        mol = MolFromSmiles(smiles_or_mol)
    else:
        mol = Chem.Mol(smiles_or_mol)
    if mol is None:
        return None
    RemoveStereochemistry(mol)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return MolToSmiles(mol, canonical=True, isomericSmiles=False)


def _native_products(rule: ReactionRule, smiles: str) -> set[str]:
    mol = MolFromSmiles(smiles)
    assert mol is not None, smiles
    found: set[str] = set()
    for products, _info in rule.metabolize(mol):
        for product in products:
            csmi = _rdkit_csmi(product)
            if csmi is None:
                continue
            found.add(csmi)
    return found


def _rust_products(rule_name: str, smiles: str) -> set[str]:
    ext = load()
    mol = ext.ForestMol(smiles)
    found: set[str] = set()
    for emission in leaf_rule(rule_name).metabolize(mol):
        for csmi_raw in emission.product_csmis():
            # product_csmis may be CXSMILES (`*OCC |$SO3;;$|`); compare graph only.
            bare = str(csmi_raw).split("|", 1)[0].strip()
            for piece in bare.split("."):
                csmi = _rdkit_csmi(piece)
                if csmi is None:
                    continue
                found.add(csmi)
    return found


def _assert_product_parity(rule_name: str, smiles: str) -> None:
    cls = getattr(native_rules, rule_name)
    # Conjugation target is star products (*OCC / CX labels). Rust always
    # collapses; native defaults as_star=True — do not force as_star=False
    # (that path is inventory-only via instantiate_rule).
    rule = cls()
    native = _native_products(rule, smiles)
    rust = _rust_products(rule_name, smiles)
    only_native = native - rust
    only_rust = rust - native
    assert not only_native and not only_rust, (
        f"{rule_name} on {smiles}: product-set drift\n"
        f"  only native ({len(only_native)}): {sorted(only_native)[:20]}\n"
        f"  only rust   ({len(only_rust)}): {sorted(only_rust)[:20]}"
    )


def _parity_params() -> list:
    if not available():
        return [
            pytest.param(
                "Hydroxylation",
                "CCO",
                marks=pytest.mark.skip(reason="xenosite.forest Rust extension required"),
            )
        ]
    paired = paired_rule_names()
    out = []
    for rule_name, smiles in parity_param_cases(paired):
        reason = product_parity_xfail_reason(rule_name, smiles)
        if reason:
            out.append(
                pytest.param(
                    rule_name,
                    smiles,
                    marks=pytest.mark.xfail(reason=reason, strict=True),
                    id=f"{rule_name}-{smiles}",
                )
            )
        else:
            out.append(pytest.param(rule_name, smiles, id=f"{rule_name}-{smiles}"))
    return out


_CASES = _parity_params()


@pytest.mark.parametrize("rule_name,smiles", _CASES)
def test_leaf_rule_product_parity(rule_name: str, smiles: str) -> None:
    _assert_product_parity(rule_name, smiles)


def test_parity_param_mode_defaults_to_full() -> None:
    paired = paired_rule_names()
    focused = parity_rule_mol_cases(paired)
    full = parity_rule_mol_cases_full(paired)
    if parity_full_enabled():
        assert len(_CASES) == len(full)
        assert len(full) > len(focused)
    else:
        assert len(_CASES) == len(focused)
        assert len(focused) < len(full)


def test_c18_product_xfail_count_locked() -> None:
    cases = product_parity_xfail_cases()
    assert len(cases) == C18_OPEN_PRODUCT_XFAIL_COUNT
