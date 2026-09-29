"""Lint mapped XMET chemist definitions for Forest mentions (sibling xmet.yaml).

When a sibling ``xmet.yaml`` with ``4000xxx`` chemist homes is present (prefer
``xenosite-tagger/xenosite-xmet``), every SSSOM subject that is an ``xmet:``
chemist home must not mention Forest / Metabolic Forest / ``forest.*`` in its
concept block (preferred_label, definition, synonyms). Forest-map alias spine
terms may name Forest and are out of scope (not ``4000xxx``).

Today many mapped homes still mention Forest in tagger ``xmet.yaml``. That is
an **xfail** (visible, strict) until upstream prose is cleaned — not an
allowlist skip. See ``.cursor/rules/never-skip-tests.mdc``.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from xenosite.forest import available

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


def _mapped_chemist_forest_mentions() -> list[str]:
    from xenosite.forest import forest_xmet_sssom

    blocks = _load_concept_blocks(XMET_YAML)
    text = forest_xmet_sssom()
    lines = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    header = lines[0].split("\t")
    si = header.index("subject_id")
    failures: list[str] = []
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
        if FOREST_MENTION.search(block):
            failures.append(f"{subject} ({_label_from_block(block)})")
    return failures


@pytest.mark.xfail(
    reason=(
        "tagger xmet.yaml: mapped 4000xxx chemist homes still mention "
        "Forest/Rainbow identity; clean upstream then drop this xfail"
    ),
    strict=True,
)
def test_chemist_definitions_omit_forest_mentions():
    if not available():
        pytest.fail("Rust extension not built; rebuild before running this lint")
    if not XMET_YAML.is_file():
        pytest.fail(f"sibling xmet.yaml required for this lint: {XMET_YAML}")
    failures = _mapped_chemist_forest_mentions()
    assert not failures, (
        "chemist XMET concepts mention Forest (fix xmet.yaml upstream):\n"
        + "\n".join(failures)
    )
