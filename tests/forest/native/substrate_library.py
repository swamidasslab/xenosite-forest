"""Shared substrate SMILES for forest coverage and canonical-plan checks.

**SoT:** ``tests/data/coverage_substrates.txt`` (shared with Rust
``include_str!``). Expand the ``[pattern]`` section when a PatternInfo / When
arm stays mute — do **not** bury probes here or grow ``_example_substrates``.

Historical sources (phase1 SMARTS, formula hints, conjugation, rare OR arms)
were collected into that file once; see the file header.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

_DATA = Path(__file__).resolve().parents[2] / "data" / "coverage_substrates.txt"


@lru_cache(maxsize=1)
def _sections() -> tuple[tuple[str, ...], tuple[str, ...]]:
    library: list[str] = []
    pattern: list[str] = []
    section: str | None = None
    text = _DATA.read_text(encoding="utf-8")
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line in ("[library]", "[pattern]"):
            section = line[1:-1]
            continue
        if section == "library":
            library.append(line)
        elif section == "pattern":
            pattern.append(line)
        else:
            raise ValueError(
                f"coverage_substrates.txt: {line!r} outside [library]/[pattern] "
                f"(section={section!r})"
            )
    if not library:
        raise ValueError("coverage_substrates.txt: empty [library]")
    if not pattern:
        raise ValueError("coverage_substrates.txt: empty [pattern]")
    return tuple(library), tuple(pattern)


def _unique_union(*parts: tuple[str, ...]) -> tuple[str, ...]:
    return tuple(dict.fromkeys(s for part in parts for s in part))


SUBSTRATE_LIBRARY: tuple[str, ...] = _sections()[0]
PATTERN_SUBSTRATES: tuple[str, ...] = _sections()[1]
COVERAGE_CANDIDATES: tuple[str, ...] = _unique_union(SUBSTRATE_LIBRARY, PATTERN_SUBSTRATES)

# Small slice for @pytest.mark.parametrize. Every entry is in SUBSTRATE_LIBRARY.
_QUICK = frozenset(
    {
        "CCO",
        "CCN",
        "CCS",
        "C=C",
        "C#C",
        "C1OC1",
        "c1ccccc1",
        "Oc1ccccc1",
        "Oc1ccc(O)cc1",
        "CCCl",
        "COc1ccccc1",
        "CSC",
        "CN(C)C",
        "CC(=O)OC",
        "COP(=O)(O)O",
        "[O-][N+](=O)c1ccccc1",
        "Clc1ccccc1",
        "Nc1ccccc1",
    }
)
QUICK_SUBSTRATES: tuple[str, ...] = tuple(
    smiles for smiles in SUBSTRATE_LIBRARY if smiles in _QUICK
)
