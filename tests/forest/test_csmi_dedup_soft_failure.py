"""Parametric CSMI-dedup compliance (C11).

Rules declare ``unique_csmi_compliant`` (class data). Yield-layer CSMI dedup
runs only when that flag is True. Goal: every leaf compliant — dedup from
sites / patterns / ``when`` alone.

This file is the **only** suite that xfails on non-compliant rules. Other
parity / chemistry tests run them normally.

- **Compliant:** hard-fail on duplicate ``(rule, pattern, emission CSMI)``
  keys or ``SiteDeduplicationWarning`` (unique-edit miss).
- **Non-compliant:** xfail when a hit is present; fail if the corpus no
  longer hits (then mark ``unique_csmi_compliant = True``).

Default cases: full leaf×corpus cartesian via :func:`parity_param_cases`.
Focused CoverIntent-only with ``XENOSITE_PARITY_FULL=0`` /
``pytest --parity-focused`` (CI may disable for speed).
"""

from __future__ import annotations

from collections import Counter
import warnings

import pytest

from xenosite.forest.rdkit_api import MolFromSmiles
from xenosite.forest.rules import (
    ReactionRule,
    SiteDeduplicationWarning,
    _unique_csmi_key,
)

from .pattern_info_inventory import instantiate_rule
from .rule_parity_corpus import (
    PARITY_FUZZ_MOLS,
    parity_full_enabled,
    parity_param_cases,
    parity_rule_mol_cases,
    parity_rule_mol_cases_full,
)
from .rule_parity_pairs import python_leaf_classes, python_parity_exception

# Formula-delta suite gate is orthogonal; this sweep only asserts CSMI soft fails.
pytestmark = pytest.mark.allow_formula_delta_mismatch


def _leaf_names() -> list[str]:
    out: list[str] = []
    for name, cls in sorted(python_leaf_classes().items()):
        if python_parity_exception(cls):
            continue
        out.append(name)
    return out


_LEAVES = _leaf_names()
_CASES = (
    parity_param_cases(_LEAVES)
    if _LEAVES
    else (("Hydroxylation", "CCO"),)
)


def _noncompliant_names() -> list[str]:
    out: list[str] = []
    for name in _LEAVES:
        cls = python_leaf_classes()[name]
        if not bool(getattr(cls, "unique_csmi_compliant", True)):
            out.append(name)
    return out


_NONCOMPLIANT = _noncompliant_names()


def _csmi_dup_hits(rule: ReactionRule, smiles: str) -> tuple[int, int]:
    """Return ``(duplicate_key_count, SiteDeduplicationWarning count)``.

    Counts under ``unique_csmi=False`` so detection does not depend on whether
    the rule's yield layer is armed (compliance gate).
    """

    mol = MolFromSmiles(smiles)
    assert mol is not None, smiles
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", SiteDeduplicationWarning)
        rows = list(rule.metabolize(mol, unique_csmi=False))
    n_warn = sum(
        1 for w in caught if issubclass(w.category, SiteDeduplicationWarning)
    )
    keys: list[tuple[str, str | None, frozenset[str]]] = []
    for products, info in rows:
        fragment: list[str] = []
        unstable = False
        for product in products:
            dedup = product.xf.tracing.dedup_smi
            if dedup is None:
                unstable = True
                break
            fragment.append(dedup)
        if unstable:
            continue
        keys.append(_unique_csmi_key(info, frozenset(fragment)))
    n_dup = sum(1 for _key, n in Counter(keys).items() if n > 1)
    return n_dup, n_warn


@pytest.mark.parametrize("rule_name,smiles", _CASES)
def test_unique_csmi_compliance(rule_name: str, smiles: str) -> None:
    """Compliant: no CSMI dups. Non-compliant: xfail only here when dups remain."""

    if not _LEAVES:
        pytest.skip("no Python leaf rules")

    cls = python_leaf_classes()[rule_name]
    rule = instantiate_rule(cls)
    compliant = bool(getattr(rule, "unique_csmi_compliant", True))
    n_dup, n_warn = _csmi_dup_hits(rule, smiles)
    has_hit = n_dup > 0 or n_warn > 0

    if not compliant:
        if has_hit:
            pytest.xfail(
                f"{rule_name} unique_csmi_compliant=False: CSMI dup / "
                f"SiteDeduplicationWarning on {smiles!r} "
                f"(dups={n_dup}, warns={n_warn}); fix unique-edit then mark compliant"
            )
        return

    if has_hit:
        raise AssertionError(
            f"compliant rule {rule_name} on {smiles!r}: "
            f"CSMI duplicate keys={n_dup}, SiteDeduplicationWarning={n_warn}. "
            f"Dedup must come from sites/patterns/whens (C11), or set "
            f"unique_csmi_compliant=False until fixed."
        )


@pytest.mark.parametrize("rule_name", _NONCOMPLIANT or ["_none_"])
def test_noncompliant_still_hits_somewhere(rule_name: str) -> None:
    """Non-compliant leaves must still show a corpus hit (else mark compliant)."""

    if rule_name == "_none_":
        pytest.skip("all leaves unique_csmi_compliant")
    cls = python_leaf_classes()[rule_name]
    rule = instantiate_rule(cls)
    assert not rule.unique_csmi_compliant
    # Sweep the full corpus: a non-compliant flag must still be justified by a
    # hit somewhere, not only on designated CoverIntent mols for this rule.
    for smiles in PARITY_FUZZ_MOLS:
        n_dup, n_warn = _csmi_dup_hits(rule, smiles)
        if n_dup > 0 or n_warn > 0:
            return
    raise AssertionError(
        f"{rule_name} is unique_csmi_compliant=False but PARITY_FUZZ_MOLS has no "
        f"CSMI dup / SiteDeduplicationWarning — set unique_csmi_compliant=True"
    )


def test_csmi_param_mode_defaults_to_full() -> None:
    """Default collection uses full cartesian unless focused is requested."""

    if not _LEAVES:
        pytest.skip("no Python leaf rules")
    focused = parity_rule_mol_cases(_LEAVES)
    full = parity_rule_mol_cases_full(_LEAVES)
    if parity_full_enabled():
        assert len(_CASES) == len(full)
        assert len(full) > len(focused)
    else:
        assert len(_CASES) == len(focused)
        assert len(focused) < len(full)
