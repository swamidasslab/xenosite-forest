"""Meta-test: :data:`PARITY_FUZZ_MOLS` covers every pattern / every ``when``.

Write and keep this green *before* relying on the mol set for Rust↔RDKit
parity fuzz. Failures list each uncovered possibility so the corpus can be
grown deliberately (no silent skips, no xfails for missing cover).
"""

from __future__ import annotations

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
