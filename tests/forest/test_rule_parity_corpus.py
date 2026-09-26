"""Meta-test: :data:`PARITY_FUZZ_MOLS` covers every pattern / every ``when``.

Write and keep this green *before* relying on the mol set for Rust↔RDKit
parity fuzz. Failures list each uncovered possibility so the corpus can be
grown deliberately (no silent skips, no xfails for missing cover).

Also requires ResonancePair **close-end** and **identical-partner** coverage
(ortho catechols, crowded ethers, …) — geometries Python often mishandles
relative to Rust.
"""

from __future__ import annotations

from rdkit.Chem import rdmolops

from xenosite.forest.rdkit_api import MolFromSmiles
from xenosite.forest.rules import (
    Dehydrogenation,
    Hydrogenation,
    QuinoneFormation,
    ResonancePairRule,
)

from .pattern_cover import uncovered_possibilities
from .pattern_info_inventory import (
    iter_pattern_possibilities,
    when_key,
)
from .rule_parity_corpus import PARITY_FUZZ_MOLS


def test_parity_fuzz_mols_cover_every_pattern_and_when():
    """Every inventory possibility has ≥1 covering mol in ``PARITY_FUZZ_MOLS``."""

    rows = list(iter_pattern_possibilities())
    gaps = uncovered_possibilities(PARITY_FUZZ_MOLS)
    if not gaps:
        return

    lines = [
        f"PARITY_FUZZ_MOLS ({len(PARITY_FUZZ_MOLS)} mols) covers "
        f"{len(rows) - len(gaps)}/{len(rows)} possibilities; "
        f"{len(gaps)} uncovered:"
    ]
    for row, _info in gaps[:40]:
        lines.append(
            f"  {row.rule_cls.__name__}/{row.pattern_name or '?'}#"
            f"{row.poss_i} when={when_key(row.when)} "
            f"group={row.group} smarts={row.smarts!r}"
        )
    if len(gaps) > 40:
        lines.append(f"  … and {len(gaps) - 40} more")
    raise AssertionError("\n".join(lines))


def test_parity_fuzz_mols_cover_every_when_branch_explicitly():
    """When-bearing possibilities are a non-empty subset; all of them covered."""

    when_rows = [r for r in iter_pattern_possibilities() if r.when is not None]
    assert when_rows, "expected when-branches in inventory"
    gaps = [
        row
        for row, _ in uncovered_possibilities(PARITY_FUZZ_MOLS)
        if row.when is not None
    ]
    assert not gaps, (
        "when-branches uncovered in PARITY_FUZZ_MOLS:\n"
        + "\n".join(
            f"  {r.rule_cls.__name__}/{r.pattern_name or '?'}#"
            f"{r.poss_i} when={when_key(r.when)}"
            for r in gaps[:40]
        )
    )


def test_parity_fuzz_mols_inventory_is_nonempty_target():
    """Sanity: the pattern inventory we must cover is non-empty."""

    rows = list(iter_pattern_possibilities())
    assert rows, "pattern inventory is empty"
    assert len(PARITY_FUZZ_MOLS) >= 50, (
        f"expected a diverse corpus, got {len(PARITY_FUZZ_MOLS)} mols"
    )


def _pair_end_distance(mol, info) -> int | None:
    ends = info.get("end_atoms")
    if not ends or len(ends) != 2:
        site = info.get("discovered_site", info.get("site"))
        if isinstance(site, (set, frozenset)) and len(site) == 2:
            ends = tuple(site)
        else:
            return None
    a, b = int(ends[0]), int(ends[1])
    path = rdmolops.GetShortestPath(mol, a, b)
    if not path:
        return None
    return len(path) - 1


def _identical_partner_role(info) -> bool:
    """True when both pair ends declare the same non-empty ``partner`` string."""

    ends = info.get("ends")
    if not ends or len(ends) != 2:
        return False
    p0 = ends[0].get("partner") or ""
    p1 = ends[1].get("partner") or ""
    return bool(p0) and p0 == p1


def _close_pair_hits(
    rule: ResonancePairRule, smiles_list: tuple[str, ...]
) -> tuple[list[str], list[str]]:
    """Return ``(close_mols, identical_partner_and_close_mols)``."""

    close: list[str] = []
    identical_close: list[str] = []
    for smiles in smiles_list:
        mol = MolFromSmiles(smiles)
        if mol is None:
            continue
        saw_close = False
        saw_identical_close = False
        for _products, info in rule.metabolize(mol):
            if "ends" not in info:
                continue
            dist = _pair_end_distance(mol, info)
            if dist is None or dist > 2:
                continue
            saw_close = True
            if _identical_partner_role(info):
                saw_identical_close = True
        if saw_close:
            close.append(smiles)
        if saw_identical_close:
            identical_close.append(smiles)
    return close, identical_close


def test_parity_fuzz_mols_cover_close_pair_ends():
    """Corpus must include close ResonancePair ends (and identical-partner closes).

    Graph distance ≤ 2 between end atoms. Identical ``partner`` on both ends
    (e.g. phenol×phenol ortho) is a separate required subclass — Rust usually
    handles these; Python sometimes does not.
    """

    required: list[tuple[type[ResonancePairRule], bool]] = [
        (Dehydrogenation, True),  # need identical-partner close (catechol)
        (QuinoneFormation, True),
        (Hydrogenation, False),  # path_end has no partner field; close only
    ]
    gaps: list[str] = []
    for cls, need_identical in required:
        rule = cls()
        close, identical_close = _close_pair_hits(rule, PARITY_FUZZ_MOLS)
        if not close:
            gaps.append(
                f"{cls.__name__}: no pair emission with end distance ≤ 2 "
                f"in PARITY_FUZZ_MOLS — add ortho/adjacent substrates"
            )
        if need_identical and not identical_close:
            gaps.append(
                f"{cls.__name__}: no close pair with identical end "
                f"``partner`` (e.g. Oc1ccccc1O phenol×phenol) in corpus"
            )
    assert not gaps, "close / identical-partner pair cover gaps:\n" + "\n".join(
        f"  {g}" for g in gaps
    )
