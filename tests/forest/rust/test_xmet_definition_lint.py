"""Lint mapped XMET chemist definitions for Forest mentions (sibling xmet.yaml).

When a sibling ``xmet.yaml`` with ``4000xxx`` chemist homes is present (prefer
``xenosite-tagger/xenosite-xmet``), every SSSOM subject that is an ``xmet:``
chemist home must not mention Forest / Metabolic Forest / ``forest.*`` in its
concept block (preferred_label, definition, synonyms). Forest-map alias spine
terms may name Forest and are out of scope (not ``4000xxx``).

Known offenders stay on ``KNOWN_FOREST_MENTION_SUBJECTS`` until upstream
``xmet.yaml`` patches land. New hits fail; cleared known IDs must be dropped
from the allowlist. Semantic mismatch notes live in LOG.md / TODO.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from xenosite.forest import available

pytestmark = pytest.mark.skipif(not available(), reason="Rust extension not built")

REPO = Path(__file__).resolve().parents[3]
_XMET_CANDIDATES = (
    REPO.parent / "xenosite-tagger" / "xenosite-xmet" / "data" / "ontology" / "xmet.yaml",
    REPO.parent / "xenosite-xmet" / "data" / "ontology" / "xmet.yaml",
)
XMET_YAML = next((p for p in _XMET_CANDIDATES if p.is_file()), _XMET_CANDIDATES[0])
FOREST_MENTION = re.compile(
    r"(?i)\b(?:metabolic\s+)?forest\b|forest\.|forest-map|rainbow\s+ruleset"
)
ID_RE = re.compile(r"^(\s*)-\s+id:\s*(xmet:\S+)\s*$")

# Known Forest mentions in mapped chemist homes (tagger xmet.yaml). Cleared
# only when upstream strips Forest/Rainbow-identity language from these concepts.
KNOWN_FOREST_MENTION_SUBJECTS = frozenset(
    {
        "xmet:4000046",
        "xmet:4000077",
        "xmet:4000113",
        "xmet:4000152",
        "xmet:4000291",
        "xmet:4000316",
        "xmet:4000317",
        "xmet:4000327",
        "xmet:4000328",
        "xmet:4000329",
        "xmet:4000330",
        "xmet:4000331",
        "xmet:4000332",
        "xmet:4000333",
        "xmet:4000334",
        "xmet:4000335",
        "xmet:4000336",
        "xmet:4000337",
        "xmet:4000344",
        "xmet:4000345",
        "xmet:4000346",
        "xmet:4000349",
        "xmet:4000350",
        "xmet:4000351",
        "xmet:4000381",
        "xmet:4000382",
        "xmet:4000383",
        "xmet:4000392",
        "xmet:4000406",
    }
)


def _load_concept_blocks(path: Path) -> dict[str, str]:
    """Map ``xmet:…`` → raw concept-block text (until next ``- id:``)."""

    out: dict[str, str] = {}
    current: str | None = None
    buf: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        m = ID_RE.match(line)
        if m:
            if current is not None:
                out[current] = "\n".join(buf)
            current = m.group(2)
            buf = [line]
            continue
        if current is not None:
            buf.append(line)
    if current is not None:
        out[current] = "\n".join(buf)
    return out


def _is_chemist_home(subject_id: str) -> bool:
    local = subject_id.removeprefix("xmet:")
    return local.isdigit() and local.startswith("4")


def _label_from_block(block: str) -> str:
    for line in block.splitlines():
        if "preferred_label:" in line:
            return line.split("preferred_label:", 1)[1].strip().strip("'\"")
    return ""


@pytest.mark.skipif(not XMET_YAML.is_file(), reason="sibling xenosite-xmet not checked out")
def test_chemist_definitions_omit_forest_mentions():
    from xenosite.forest import forest_xmet_sssom

    blocks = _load_concept_blocks(XMET_YAML)
    text = forest_xmet_sssom()
    lines = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    header = lines[0].split("\t")
    si = header.index("subject_id")
    new_failures: list[str] = []
    known_hits: list[str] = []
    seen: set[str] = set()
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) <= si:
            continue
        subject = parts[si]
        if not subject.startswith("xmet:") or not _is_chemist_home(subject):
            continue
        if subject in seen:
            continue
        seen.add(subject)
        block = blocks.get(subject)
        if block is None:
            continue
        if not FOREST_MENTION.search(block):
            continue
        label = _label_from_block(block)
        msg = f"{subject} ({label})"
        if subject in KNOWN_FOREST_MENTION_SUBJECTS:
            known_hits.append(subject)
        else:
            new_failures.append(msg)
    stale = sorted(KNOWN_FOREST_MENTION_SUBJECTS - set(known_hits))
    assert not new_failures, (
        "new chemist XMET concepts mention Forest (propose xmet.yaml patches upstream):\n"
        + "\n".join(new_failures)
    )
    assert not stale, (
        "KNOWN_FOREST_MENTION_SUBJECTS stale (defs cleaned upstream — drop from allowlist):\n"
        + "\n".join(stale)
    )
