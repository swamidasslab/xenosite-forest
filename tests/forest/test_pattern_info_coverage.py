"""Parameterized PatternInfo validity + ``when``-branch resolution coverage.

One case per possibility on every catalog PatternInfo. A covering mol is
chosen so ``resolve_effect`` selects that possibility (including every
``when`` branch). Assertions stay on PatternInfo shape and resolution, not
full product chemistry.

Inventory comes from :mod:`.pattern_info_inventory` so this file and the
completeness meta-test cannot drift. Known unreachable branches are
explicit xfails — do not drop them.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import pytest
from rdkit import Chem

from xenosite.forest.rdkit_api import Mol, MolFromSmiles
from xenosite.forest.records import Effect, PatternInfo, When
from xenosite.forest.rules import (
    ReactionRule,
    _kekule_forms,
    _site_indexes,
    _site_map_aromatic,
    _span,
    _when_matches,
    resolve_effect,
)

from .pattern_info_inventory import (
    PATTERNLESS_REACTION_RULE_BASES,
    discover_reaction_rule_classes,
    instantiate_rule,
    iter_pattern_possibilities,
    pattern_possibility_keys,
    patterns_on,
    possibility_key,
    when_key,
)
from .substrate_library import SUBSTRATE_LIBRARY

# Extra substrates beyond SUBSTRATE_LIBRARY that hit rare OR branches
# (At / I geminal dihalides, hemiaminals, aziridines, arene oxide, etc.).
_PATTERN_SUBSTRATES: tuple[str, ...] = (
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
    "CCCl",
    "CCBr",
    "CCI",
    "CC(F)C",
    "CC(Cl)C",
    "CC(Br)C",
    "CC(I)C",
    "FC(F)F",
    "ClCCl",
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
    "Clc1ccc(O)cc1",
    "Brc1ccc(O)cc1",
    "Ic1ccc(O)cc1",
    "ClCc1ccccc1",
    "BrCc1ccccc1",
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
)

_CANDIDATES: tuple[str, ...] = tuple(
    dict.fromkeys((*SUBSTRATE_LIBRARY, *_PATTERN_SUBSTRATES))
)

# (rule, pattern_name, poss_index) -> reason. Unreachable with current SMARTS.
_UNREACHABLE: dict[tuple[str, str, int], str] = {}


def _atom_satisfies(mol: Mol, mapped: Mapping[int, int], when: When | None) -> bool:
    if when is None:
        return True
    return _when_matches(mol, mapped, when)


def _selects_possibility(
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
        return _atom_satisfies(context, mapped, when)

    if target_when is not None and not hits(target_when):
        return False
    for earlier in possibilities[:poss_i]:
        if hits(earlier.get("when")):
            return False
    return True


def _find_cover(
    info: PatternInfo,
    reactant: str,
    poss_i: int,
    mols: Sequence[tuple[str, Mol]],
) -> tuple[str, Mol, Mol, dict[int, int]] | None:
    for smiles, mol in mols:
        for work in _kekule_forms(mol):
            for mapped in work.xf.smarts_matches(reactant):
                if not _selects_possibility(mol, mapped, info, poss_i):
                    continue
                if not _site_indexes(mapped, info):
                    continue
                return smiles, mol, work, dict(mapped)
    return None


def _assert_pattern_valid(smarts: str, info: PatternInfo) -> None:
    possibilities = info.get("possibilities")
    assert possibilities, "PatternInfo needs at least one possibility"
    assert "span" in info, "PatternInfo needs span"
    assert info["span"] == _span(possibilities), "span must match collapsed possibilities"
    reactant = smarts.split(">>", 1)[0]
    query = Chem.MolFromSmarts(reactant)
    assert query is not None, f"reactant SMARTS does not parse: {reactant}"
    for poss in possibilities:
        when = poss.get("when")
        if when is None:
            continue
        assert "map" in when, when
        assert "z" in when or "h" in when, when


def _assert_resolution(
    context: Mol,
    mapped: Mapping[int, int],
    info: PatternInfo,
    poss_i: int,
) -> Effect:
    possibilities = info.get("possibilities") or ()
    chosen = possibilities[poss_i]
    effect = resolve_effect(context, mapped, info)
    assert "when" not in effect

    site_map = info.get("site_map", 1)
    if isinstance(site_map, int):
        assert site_map in mapped
        atom = context.GetAtomWithIdx(mapped[site_map])
        assert effect.get("symbol") == atom.GetSymbol()
        assert effect.get("h") == atom.GetTotalNumHs()
        assert effect.get("site_aromatic") == atom.GetIsAromatic()
    else:
        assert all(m in mapped for m in site_map)

    assert _site_indexes(mapped, info)

    # Declared outcome fields from the chosen possibility survive resolution.
    # breaks_ring / partner / partner_h may be filled by resolve_effect, so
    # they are not compared to the declared possibility here.
    # dearomatizes on PatternInfo is capability: resolve_effect narrows it to
    # whether a site_map atom is aromatic (same as ResonancePair merge).
    for key in (
        "adds",
        "removes",
        "cleaves",
        "leave_count",
        "methide",
        "exclusive_partner",
        "needs",
    ):
        if key in chosen:
            assert effect.get(key) == chosen[key], (key, effect.get(key), chosen[key])

    if chosen.get("dearomatizes"):
        assert effect.get("dearomatizes") == _site_map_aromatic(
            context, mapped, site_map
        ), (
            "dearomatizes",
            effect.get("dearomatizes"),
            _site_map_aromatic(context, mapped, site_map),
        )
    elif "dearomatizes" in chosen:
        assert effect.get("dearomatizes") is False

    when = chosen.get("when")
    if when is not None:
        assert _when_matches(context, mapped, when)
        idx = mapped.get(when.get("map", -1))
        assert idx is not None
        atom = context.GetAtomWithIdx(idx)
        if "z" in when:
            assert atom.GetAtomicNum() == when["z"]
        if "h" in when:
            assert atom.GetTotalNumHs() == when["h"]

    return effect


def _parsed_mols() -> list[tuple[str, Mol]]:
    out: list[tuple[str, Mol]] = []
    for smiles in _CANDIDATES:
        mol = MolFromSmiles(smiles)
        if mol is not None:
            out.append((smiles, mol))
    return out


def _build_cases():
    mols = _parsed_mols()
    cases = []
    when_total = 0
    for row in iter_pattern_possibilities():
        if row.when is not None:
            when_total += 1
        reactant = row.smarts.split(">>", 1)[0]
        # Re-resolve PatternInfo from the live rule so cover search sees
        # the same object the parameterized test will assert on.
        rule = instantiate_rule(row.rule_cls)
        info = None
        for group, smarts, pattern in patterns_on(rule):
            if (
                group == row.group
                and smarts == row.smarts
                and (pattern.get("name") or "") == row.pattern_name
            ):
                info = pattern
                break
        assert info is not None, row
        cover = _find_cover(info, reactant, row.poss_i, mols)
        gap = _UNREACHABLE.get((row.rule_cls.__name__, row.pattern_name, row.poss_i))
        case_id = (
            f"{row.rule_cls.__name__}/{row.pattern_name or '?'}#{row.poss_i}"
            f"/when={when_key(row.when)}"
        )
        if gap is not None:
            cases.append(
                pytest.param(
                    row.rule_cls,
                    row.group,
                    row.smarts,
                    row.pattern_name,
                    row.poss_i,
                    row.when,
                    None,
                    id=case_id,
                    marks=pytest.mark.xfail(reason=gap, strict=True),
                )
            )
        elif cover is None:
            cases.append(
                pytest.param(
                    row.rule_cls,
                    row.group,
                    row.smarts,
                    row.pattern_name,
                    row.poss_i,
                    row.when,
                    None,
                    id=case_id,
                    marks=pytest.mark.xfail(
                        reason=(
                            "no covering mol in substrate library + "
                            "_PATTERN_SUBSTRATES; add a substrate or "
                            "record an _UNREACHABLE gap"
                        ),
                        strict=True,
                    ),
                )
            )
        else:
            smiles, _mol, _work, _mapped = cover
            cases.append(
                pytest.param(
                    row.rule_cls,
                    row.group,
                    row.smarts,
                    row.pattern_name,
                    row.poss_i,
                    row.when,
                    smiles,
                    id=case_id,
                )
            )
    return cases, when_total


_CASES, _WHEN_BRANCH_COUNT = _build_cases()


@pytest.mark.parametrize(
    "rule_cls, group, smarts, pattern_name, poss_i, when, smiles",
    _CASES,
)
def test_pattern_info_possibility(
    rule_cls: type[ReactionRule],
    group: str,
    smarts: str,
    pattern_name: str,
    poss_i: int,
    when: When | None,
    smiles: str | None,
):
    rule = instantiate_rule(rule_cls)
    info = None
    for g, text, pattern in patterns_on(rule):
        if g == group and text == smarts and pattern.get("name", "") == pattern_name:
            info = pattern
            break
    assert info is not None, (rule_cls.__name__, group, smarts, pattern_name)
    _assert_pattern_valid(smarts, info)

    assert smiles is not None, "expected a covering smiles for this possibility"
    mol = MolFromSmiles(smiles)
    assert mol is not None, smiles
    reactant = smarts.split(">>", 1)[0]

    hit = None
    for work in _kekule_forms(mol):
        for mapped in work.xf.smarts_matches(reactant):
            if _selects_possibility(mol, mapped, info, poss_i):
                hit = dict(mapped)
                break
        if hit is not None:
            break
    assert hit is not None, (
        f"no match selecting poss[{poss_i}] when={when} on {smiles} "
        f"for {rule_cls.__name__}/{pattern_name}"
    )
    # resolve_effect reads aromatic flags from the caller's mol (context).
    _assert_resolution(mol, hit, info, poss_i)


def test_when_branch_inventory_matches_params():
    """Every ``when`` on every PatternInfo has a parameterized row."""

    when_params = sum(1 for param in _CASES if param.values[5] is not None)
    assert when_params == _WHEN_BRANCH_COUNT
    live = sum(1 for row in iter_pattern_possibilities() if row.when is not None)
    assert when_params == live


def test_coverage_rows_match_rule_derived_possibilities():
    """Coverage params == :func:`iter_pattern_possibilities` (no silent gaps)."""

    expected = pattern_possibility_keys()
    observed = {
        possibility_key(param.values[0], param.values[3], param.values[4], param.values[5])
        for param in _CASES
    }
    missing = sorted(expected - observed)
    extra = sorted(observed - expected)
    assert not missing and not extra, (
        f"pattern-info coverage drifted from rule inventory; "
        f"missing={missing[:20]}{'…' if len(missing) > 20 else ''} "
        f"extra={extra[:20]}{'…' if len(extra) > 20 else ''}"
    )


def test_patternless_whitelist_is_exact():
    """Whitelist entries subclass ReactionRule, exist, and stay patternless."""

    discovered = set(discover_reaction_rule_classes())
    for cls in PATTERNLESS_REACTION_RULE_BASES:
        assert issubclass(cls, ReactionRule), cls
        assert cls in discovered, (
            f"whitelist entry {cls.__name__} is not a discovered ReactionRule "
            f"subclass under xenosite.forest"
        )

    patterned_whitelist = sorted(
        cls.__name__
        for cls in PATTERNLESS_REACTION_RULE_BASES
        if list(patterns_on(instantiate_rule(cls)))
    )
    assert not patterned_whitelist, (
        "whitelist entries must be patternless; move coverage to the "
        f"inventory instead: {patterned_whitelist}"
    )


def test_every_concrete_rule_with_patterns_is_inventoried():
    """Non-whitelisted rules that declare patterns appear in the inventory."""

    keys = pattern_possibility_keys()
    missing: list[tuple[str, str, int, object]] = []
    for cls in discover_reaction_rule_classes():
        if cls in PATTERNLESS_REACTION_RULE_BASES:
            continue
        rule = instantiate_rule(cls)
        for _group, _smarts, info in patterns_on(rule):
            name = info.get("name") or ""
            for poss_i, poss in enumerate(info.get("possibilities") or ()):
                key = possibility_key(cls, name, poss_i, poss.get("when"))
                if key not in keys:
                    missing.append(key)
    assert not missing, (
        f"rule-derived possibilities missing from inventory: {missing[:20]}"
        f"{'…' if len(missing) > 20 else ''}"
    )
