#!/usr/bin/env python3
"""Offline harvest of Forest SMIRKS/SMARTS into XRM candidate JSONL.

Does **not** import ``xenosite.forest`` (avoids RDKit / namer coupling). Parses
``src/xenosite/forest/rules.py`` as text for ``Smirks("…")`` literals and nearby
pattern ``name=`` / rule class context.

Uses:
  - Seed reactant SMARTS for XRM assignment rules (review before promote).
  - Validation corpus: assert XRM structural rules cover Forest pattern families.

Output: ``data/candidates/forest-smarts.jsonl``
"""

from __future__ import annotations

import ast
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
RULES_PY = REPO / "src/xenosite/forest/rules.py"
OUT = Path(__file__).resolve().parents[1] / "data/candidates/forest-smarts.jsonl"

SMIRKS_RE = re.compile(r'Smirks\(\s*"([^"]+)"\s*\)')
# Pattern name near a smirks often looks like name="methyl_alcohol" or "name": "…"
NAME_NEAR_RE = re.compile(r'name\s*=\s*"([A-Za-z0-9_/\-]+)"')
CLASS_RE = re.compile(r"^class\s+(\w+)\s*\(")


def reactant_of(smirks: str) -> str:
    return smirks.split(">>", 1)[0] if ">>" in smirks else smirks


def product_of(smirks: str) -> str | None:
    if ">>" not in smirks:
        return None
    return smirks.split(">>", 1)[1]


def main() -> None:
    text = RULES_PY.read_text()
    lines = text.splitlines()
    current_class = None
    rows = []
    seen = set()

    for i, line in enumerate(lines):
        mclass = CLASS_RE.match(line)
        if mclass:
            current_class = mclass.group(1)
        for m in SMIRKS_RE.finditer(line):
            smirks = m.group(1)
            # Look nearby (±8 lines) for a pattern name.
            window = "\n".join(lines[max(0, i - 8) : i + 9])
            names = NAME_NEAR_RE.findall(window)
            pname = names[-1] if names else None
            key = (current_class, pname, smirks)
            if key in seen:
                continue
            seen.add(key)
            reactant = reactant_of(smirks)
            product = product_of(smirks)
            forest_tag = None
            if current_class and pname:
                forest_tag = f"forest.pattern:{current_class}/{pname}"
            elif current_class:
                forest_tag = f"forest.rule:{current_class}"
            rows.append(
                {
                    "candidate_id": f"cand:forest-smarts:{len(rows):04d}",
                    "source": "forest.rules.py",
                    "rule_class": current_class,
                    "pattern_name": pname,
                    "smirks": smirks,
                    "reactant_smarts": reactant,
                    "product_smarts": product,
                    "forest_tag": forest_tag,
                    "suggested_spines": [
                        "chemical transformation",
                        "rule provenance",
                        "Metabolic Forest map",
                    ],
                    "notes": (
                        "Harvested statically from Forest rules.py. Review before "
                        "promoting into xenobiotic.jsonl. Prefer as validation "
                        "coverage checks; namer must not import forest."
                    ),
                }
            )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w") as f:
        for r in rows:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(f"wrote {OUT} ({len(rows)} smirks from {RULES_PY.relative_to(REPO)})")

    # Coverage report vs existing XRM assignment reactant SMARTS
    asg = Path(__file__).resolve().parents[1] / "data/assignments/xenobiotic.jsonl"
    xrm_smarts = set()
    for line in asg.read_text().splitlines():
        if not line.startswith("{"):
            continue
        o = json.loads(line)
        if o.get("reactant_smarts"):
            xrm_smarts.add(o["reactant_smarts"])
    forest_reactants = {r["reactant_smarts"] for r in rows}
    overlap = xrm_smarts & forest_reactants
    print(
        f"xrm reactant_smarts={len(xrm_smarts)} "
        f"forest_unique_reactants={len(forest_reactants)} "
        f"exact_overlap={len(overlap)}"
    )
    print(
        "note: low exact overlap is expected (XRM uses simpler chemist SMARTS; "
        "Forest SMIRKS are denser). Use forest-smarts.jsonl for validation pairs."
    )
    # Refresh coverage report used by Rust validation tests.
    import subprocess
    import sys

    subprocess.check_call(
        [sys.executable, str(Path(__file__).with_name("validate_forest_coverage.py"))]
    )


if __name__ == "__main__":
    main()
