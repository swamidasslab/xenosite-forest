"""Meta-test: :data:`PARITY_CORPUS` covers every pattern / every ``when``.

Write and keep this green *before* relying on the corpus for Rust↔RDKit
parity. Failures list each uncovered possibility so the corpus can be
grown deliberately (no silent skips, no xfails for missing cover).

Each :class:`~.rule_parity_corpus.ParityEntry` carries a list of
:class:`~.rule_parity_corpus.CoverIntent` correspondences (rule / pattern /
poss_i / when). Meta-tests verify intents hold and are inventory-complete.
Parametric suites use those intents via :func:`~.rule_parity_corpus.parity_param_cases`
(focused) or full cartesian under ``--parity-full`` / ``XENOSITE_PARITY_FULL``.

Also requires ResonancePair **close-end** and **identical-partner** coverage
(ortho catechols, crowded ethers, …) — geometries Python often mishandles
relative to Rust.

Mapped reactant ``[#Z]`` (aromaticable Z) must hit **aliphatic and aromatic**
in the corpus whenever both states are reachable — data indicator on the
SMARTS, not rule-name branches (equivalent ``#`` expand coverage).
"""

from __future__ import annotations

import pytest
from rdkit.Chem import rdmolops

from xenosite.forest.rdkit_api import MolFromSmiles
from xenosite.forest.rules import (
    Dehydrogenation,
    Hydrogenation,
    QuinoneFormation,
    ResonancePairRule,
)

from .pattern_cover import (
    find_cover,
    parse_mols,
    pattern_info_for,
    uncovered_possibilities,
)
from .pattern_info_inventory import (
    iter_pattern_possibilities,
    possibility_key,
    when_key,
)
from .rule_parity_corpus import (
    PARITY_CORPUS,
    PARITY_FUZZ_MOLS,
    parity_cover_cases,
)


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


def test_parity_fuzz_mols_cover_aliphatic_and_aromatic_smarts_branches():
    """Reactant ``#N`` maps must hit aliphatic and aromatic when both are reachable.

    Data indicator (not rule names): a mapped ``[#Z:…]`` with aromaticable Z
    (C/N/O/P/S) expands to aliphatic **and** aromatic organic spellings under
    chematic / ``specialize_smirks_for_maps``. Corpus coverage must include
    both aromaticity states whenever either state is reachable for that SMARTS
    (corpus ∪ small aliphatic/aromatic probes). Grow ``PARITY_FUZZ_MOLS`` when
    a branch is missing — do not special-case rule classes.
    """

    import re

    from xenosite.forest.rules import _kekule_forms

    from .pattern_info_inventory import iter_pattern_possibilities

    aromaticable = {6, 7, 8, 15, 16}
    bracket = re.compile(r"\[#(\d+)([^\]]*)\]")
    # Reachability probes only — gaps are fixed by growing PARITY_FUZZ_MOLS.
    aliphatic_probe = (
        "CCO",
        "CCN",
        "CS",
        "CC",
        "C=C",
        "CC=O",
        "CC=CC=O",
        "CCl",
        "NCC",
        "OCC",
        "SCC",
    )
    aromatic_probe = (
        "[nH]1cccc1",
        "c1ccccc1",
        "o1cccc1",
        "s1cccc1",
        "c1ccncc1",
        "c1ccoc1",
        "c1ccsc1",
        "c1ccc2[nH]ccc2c1",
        "n1ccccc1",
        "Oc1ccccc1",
        "Nc1ccccc1",
        "Sc1ccccc1",
    )

    def mapped_hash_atoms(smarts: str) -> list[tuple[int, int]]:
        reactant = smarts.split(">>", 1)[0]
        out: list[tuple[int, int]] = []
        for m in bracket.finditer(reactant):
            z = int(m.group(1))
            if z not in aromaticable:
                continue
            map_m = re.search(r":(\d+)", m.group(2))
            if not map_m:
                continue
            out.append((int(map_m.group(1)), z))
        return out

    def collect_hits(
        reactant: str, mapno: int, z: int, smiles_list: tuple[str, ...]
    ) -> tuple[set[str], set[str]]:
        ali: set[str] = set()
        aro: set[str] = set()
        for smiles in smiles_list:
            mol = MolFromSmiles(smiles)
            if mol is None:
                continue
            for work in _kekule_forms(mol):
                for mapped in work.xf.smarts_matches(reactant):
                    if mapno not in mapped:
                        continue
                    atom = work.GetAtomWithIdx(mapped[mapno])
                    if atom.GetAtomicNum() != z:
                        continue
                    if atom.GetIsAromatic():
                        aro.add(smiles)
                    else:
                        ali.add(smiles)
        return ali, aro

    seen_patterns: set[tuple[str, str, str]] = set()
    gaps: list[str] = []
    both = 0
    for row in iter_pattern_possibilities():
        key = (row.rule_cls.__name__, row.pattern_name or "?", row.smarts)
        if key in seen_patterns:
            continue
        seen_patterns.add(key)
        reactant = row.smarts.split(">>", 1)[0]
        for mapno, z in mapped_hash_atoms(row.smarts):
            ali, aro = collect_hits(reactant, mapno, z, PARITY_FUZZ_MOLS)
            p_ali, _ = collect_hits(reactant, mapno, z, aliphatic_probe)
            _, p_aro = collect_hits(reactant, mapno, z, aromatic_probe)
            need_ali = bool(ali or p_ali)
            need_aro = bool(aro or p_aro)
            if need_ali and need_aro:
                both += 1
                missing: list[str] = []
                if not ali:
                    missing.append("aliphatic")
                if not aro:
                    missing.append("aromatic")
                if missing:
                    gaps.append(
                        f"{row.rule_cls.__name__}/{row.pattern_name or '?'} "
                        f"map={mapno} z={z} missing {','.join(missing)} "
                        f"(# expand needs both; smarts={reactant!r})"
                    )

    assert both > 0, (
        "expected some reactant # maps reachable in both aromaticity states"
    )
    assert not gaps, (
        "SMARTS aliphatic/aromatic branch gaps — add covering mols to "
        "PARITY_FUZZ_MOLS (data indicator: mapped [#Z] on aromaticable Z):\n"
        + "\n".join(f"  {g}" for g in gaps[:40])
    )


def test_parity_fuzz_mols_inventory_is_nonempty_target():
    """Sanity: the pattern inventory we must cover is non-empty."""

    rows = list(iter_pattern_possibilities())
    assert rows, "pattern inventory is empty"
    assert len(PARITY_FUZZ_MOLS) >= 50, (
        f"expected a diverse corpus, got {len(PARITY_FUZZ_MOLS)} mols"
    )
    assert PARITY_CORPUS, "PARITY_CORPUS is empty"
    assert parity_cover_cases(), "expected CoverIntent correspondences"


def test_parity_cover_intents_complete_inventory():
    """Every inventory possibility appears exactly once as a CoverIntent key."""

    inv = {
        possibility_key(r.rule_cls, r.pattern_name or "?", r.poss_i, r.when)
        for r in iter_pattern_possibilities()
    }
    covers = {
        (rule, pattern, poss_i, when)
        for rule, _smiles, pattern, poss_i, when in parity_cover_cases()
    }
    missing = sorted(inv - covers)
    extra = sorted(covers - inv)
    lines: list[str] = []
    if missing:
        lines.append(f"{len(missing)} inventory keys lack CoverIntent:")
        lines.extend(f"  {m}" for m in missing[:40])
    if extra:
        lines.append(f"{len(extra)} CoverIntent keys not in inventory:")
        lines.extend(f"  {e}" for e in extra[:40])
    assert not lines, "\n".join(lines)


@pytest.mark.parametrize(
    "rule,smiles,pattern,poss_i,when",
    parity_cover_cases() or [("Hydroxylation", "CC", "h2", 0, None)],
)
def test_parity_cover_intent_holds(
    rule: str, smiles: str, pattern: str, poss_i: int, when: object
) -> None:
    """Declared CoverIntent actually selects that possibility on that mol."""

    if not parity_cover_cases():
        pytest.skip("no cover intents")
    rows = [
        r
        for r in iter_pattern_possibilities()
        if r.rule_cls.__name__ == rule
        and (r.pattern_name or "?") == pattern
        and r.poss_i == poss_i
        and when_key(r.when) == when
    ]
    assert rows, (
        f"no inventory row for CoverIntent "
        f"{rule}/{pattern}#{poss_i} when={when!r}"
    )
    row = rows[0]
    info = pattern_info_for(row)
    reactant = row.smarts.split(">>", 1)[0]
    hit = find_cover(info, reactant, poss_i, parse_mols([smiles]))
    assert hit is not None, (
        f"CoverIntent failed: {rule}/{pattern}#{poss_i} when={when!r} "
        f"on {smiles!r} (smarts={reactant!r})"
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
