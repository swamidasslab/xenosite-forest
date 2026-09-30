"""Native RDKit vs Rust chematic product-set parity on the shared coverage pool.

SoT substrates: ``tests/data/coverage_substrates.txt`` via
:mod:`.substrate_library` (same file Rust catalog scans). Leaves are looked up
by catalog name (``leaf_rule`` / native class) — not ``xf:`` short codes.

Compares RDKit-canonical fragment SMILES sets per leaf × substrate. A product
only one door emits is a failure unless ``docs/forest/DIVERGENCES.md`` records
both sets and the reason — do not soft-skip. SMIRKS edits must not introduce
drift unless they progress toward correct chemistry.
"""

from __future__ import annotations

import pytest
from rdkit import Chem

from xenosite.forest import available
from xenosite.forest.native import rules as native_rules

from .helpers import canon, rust_leaf_products_canon
from .rule_parity_corpus import product_parity_xfail_reason
from .substrate_library import COVERAGE_CANDIDATES

pytestmark = pytest.mark.filterwarnings(
    "ignore::xenosite.forest.native.rules.SiteDeduplicationWarning",
)

# Native PhaseOne leaves (+ NDealkylation) that also exist as Rust
# ``LEAF_CTORS``. EpoxideHydration is Rust-only (composite); conjugation
# adducts stay on the adduct xfail. Lookup is by catalog name, not ``xf:``.
_PARITY_LEAVES: tuple[str, ...] = (
    "Hydroxylation",
    "Epoxidation",
    "SulfurOxidation",
    "NitrogenOxidation",
    "Dehydrogenation",
    "QuinoneFormation",
    "Dephosphorylation",
    "EpoxideOpening",
    "Hydrolysis",
    "Dehydration",
    "Hydrogenation",
    "NitrogenReduction",
    "OxygenReduction",
    "ReductiveDehalogenation",
    "SulfurReduction",
    "Dealkylation",
    "OxidativeDehalogenation",
    "NDealkylation",
)


def _native_products(rule_cls: type, smiles: str) -> set[str] | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    found: set[str] = set()
    for products, _info in rule_cls().metabolize(mol):
        for product in products:
            text = Chem.MolToSmiles(product)
            if "." in text:
                for frag in text.split("."):
                    c = canon(frag)
                    if c:
                        found.add(c)
                continue
            c = canon(product)
            if c:
                found.add(c)
    return found


def _rust_products(leaf: str, smiles: str) -> set[str] | None:
    return rust_leaf_products_canon(leaf, smiles)


@pytest.fixture(scope="module", autouse=True)
def _require_rust_extension() -> None:
    if not available():
        pytest.fail(
            "xenosite.forest._rust extension required for native↔Rust coverage parity"
        )


def _coverage_parity_params() -> list:
    out: list = []
    for leaf in _PARITY_LEAVES:
        for smiles in COVERAGE_CANDIDATES:
            reason = product_parity_xfail_reason(leaf, smiles)
            if reason:
                out.append(
                    pytest.param(
                        leaf,
                        smiles,
                        marks=pytest.mark.xfail(reason=reason, strict=True),
                        id=f"{smiles}-{leaf}",
                    )
                )
            else:
                out.append(
                    pytest.param(leaf, smiles, id=f"{smiles}-{leaf}"),
                )
    return out


@pytest.mark.parametrize("leaf,smiles", _coverage_parity_params())
def test_native_rust_product_sets_match_on_coverage(leaf: str, smiles: str) -> None:
    rule_cls = getattr(native_rules, leaf)
    native = _native_products(rule_cls, smiles)
    rust = _rust_products(leaf, smiles)
    if native is None and rust is None:
        return
    assert native is not None, f"{leaf} on {smiles}: native could not parse substrate"
    assert rust is not None, f"{leaf} on {smiles}: Rust could not parse substrate"
    only_native = native - rust
    only_rust = rust - native
    assert not only_native and not only_rust, (
        f"{leaf} on {smiles}: product-set drift\n"
        f"  only native ({len(only_native)}): {sorted(only_native)[:20]}\n"
        f"  only rust   ({len(only_rust)}): {sorted(only_rust)[:20]}\n"
        f"  (record intentional divergence in docs/forest/DIVERGENCES.md; "
        f"do not soft-skip)"
    )
