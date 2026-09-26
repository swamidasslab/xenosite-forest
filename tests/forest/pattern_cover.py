"""Find covering mols for PatternInfo possibilities (including every ``when``).

Shared by pattern-info coverage and the Rust↔RDKit parity fuzz corpus
meta-test. A mol *covers* possibility ``poss_i`` when a SMARTS hit on that
mol makes ``resolve_effect`` select that branch (earlier ``when`` arms miss).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

from xenosite.forest.rdkit_api import Mol, MolFromSmiles
from xenosite.forest.records import PatternInfo, When
from xenosite.forest.rules import _kekule_forms, _site_indexes, _when_matches

from .pattern_info_inventory import (
    PatternPossibility,
    instantiate_rule,
    iter_pattern_possibilities,
    patterns_on,
)


def atom_satisfies(mol: Mol, mapped: Mapping[int, int], when: When | None) -> bool:
    if when is None:
        return True
    return _when_matches(mol, mapped, when)


def selects_possibility(
    context: Mol,
    mapped: Mapping[int, int],
    info: PatternInfo,
    poss_i: int,
) -> bool:
    """True when ``resolve_effect(context, ...)`` would pick ``possibilities[poss_i]``.

    Production resolves against the aromatic context mol, not the kekulé copy.
    """

    possibilities = info.get("possibilities") or ()
    if poss_i < 0 or poss_i >= len(possibilities):
        return False
    target_when = possibilities[poss_i].get("when")

    def hits(when: When | None) -> bool:
        return atom_satisfies(context, mapped, when)

    if target_when is not None and not hits(target_when):
        return False
    for earlier in possibilities[:poss_i]:
        if hits(earlier.get("when")):
            return False
    return True


def find_cover(
    info: PatternInfo,
    reactant: str,
    poss_i: int,
    mols: Sequence[tuple[str, Mol]],
) -> tuple[str, Mol, Mol, dict[int, int]] | None:
    """First ``(smiles, context, kekule, mapped)`` that selects ``poss_i``."""

    for smiles, mol in mols:
        for work in _kekule_forms(mol):
            for mapped in work.xf.smarts_matches(reactant):
                if not selects_possibility(mol, mapped, info, poss_i):
                    continue
                if not _site_indexes(mapped, info):
                    continue
                return smiles, mol, work, dict(mapped)
    return None


def parse_mols(smiles_list: Sequence[str]) -> list[tuple[str, Mol]]:
    out: list[tuple[str, Mol]] = []
    for smiles in smiles_list:
        mol = MolFromSmiles(smiles)
        if mol is not None:
            out.append((smiles, mol))
    return out


def pattern_info_for(row: PatternPossibility) -> PatternInfo:
    """Live PatternInfo object for ``row`` (same object resolve_effect sees)."""

    rule = instantiate_rule(row.rule_cls)
    for group, smarts, pattern in patterns_on(rule):
        if (
            group == row.group
            and smarts == row.smarts
            and (pattern.get("name") or "") == row.pattern_name
        ):
            return pattern
    raise AssertionError(
        f"no PatternInfo for {row.rule_cls.__name__}/"
        f"{row.pattern_name!r} group={row.group} smarts={row.smarts!r}"
    )


def uncovered_possibilities(
    smiles_list: Sequence[str],
) -> list[tuple[PatternPossibility, PatternInfo]]:
    """Every inventory possibility with no covering mol in ``smiles_list``."""

    mols = parse_mols(smiles_list)
    gaps: list[tuple[PatternPossibility, PatternInfo]] = []
    for row in iter_pattern_possibilities():
        info = pattern_info_for(row)
        reactant = row.smarts.split(">>", 1)[0]
        if find_cover(info, reactant, row.poss_i, mols) is None:
            gaps.append((row, info))
    return gaps


def covering_smiles_map(
    smiles_list: Sequence[str],
) -> dict[tuple[str, str, int], str]:
    """``(rule, pattern_name, poss_i) → first covering smiles`` for hits only."""

    mols = parse_mols(smiles_list)
    out: dict[tuple[str, str, int], str] = {}
    for row in iter_pattern_possibilities():
        info = pattern_info_for(row)
        reactant = row.smarts.split(">>", 1)[0]
        cover = find_cover(info, reactant, row.poss_i, mols)
        if cover is not None:
            key = (row.rule_cls.__name__, row.pattern_name or "?", row.poss_i)
            out[key] = cover[0]
    return out
