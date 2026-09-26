#!/usr/bin/env python3
"""Validate XRM coverage of harvested Forest SMIRKS (offline; no forest import).

Reads:
  - data/candidates/forest-smarts.jsonl
  - data/mappings/xrm-forest.sssom.tsv
  - data/assignments/xenobiotic.jsonl

Writes:
  - data/candidates/forest-coverage.json  (report for CI / review)

Exit 0 always when report is written; non-zero only on missing harvest file.
The Rust test asserts policy thresholds on the report.
"""

from __future__ import annotations

import json
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HARVEST = ROOT / "data/candidates/forest-smarts.jsonl"
SSSOM = ROOT / "data/mappings/xrm-forest.sssom.tsv"
ASSIGN = ROOT / "data/assignments/xenobiotic.jsonl"
OUT = ROOT / "data/candidates/forest-coverage.json"


def main() -> None:
    if not HARVEST.exists():
        raise SystemExit(f"missing {HARVEST}; run harvest_forest_smarts.py first")

    rows = [
        json.loads(line)
        for line in HARVEST.read_text().splitlines()
        if line.strip()
    ]

    # SSSOM object ids (opaque forest.* CURIEs)
    sssom_objects: set[str] = set()
    sssom_by_object: dict[str, list[str]] = defaultdict(list)
    for line in SSSOM.read_text().splitlines():
        if not line or line.startswith("#") or line.startswith("subject"):
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        subj, _pred, obj = parts[0], parts[1], parts[2]
        sssom_objects.add(obj)
        sssom_by_object[obj].append(subj)

    # Assignment tags_any forest.* strings
    assign_tags: set[str] = set()
    for line in ASSIGN.read_text().splitlines():
        if not line.startswith("{"):
            continue
        o = json.loads(line)
        for t in o.get("tags_any") or []:
            if t.startswith("forest."):
                assign_tags.add(t)

    by_class: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        cls = r.get("rule_class") or "unknown"
        by_class[cls].append(r)

    class_report = []
    mapped_classes = 0
    taggable_classes = 0
    for cls in sorted(by_class):
        rule_curie = f"forest.rule:{cls}"
        patterns = sorted(
            {
                r["forest_tag"]
                for r in by_class[cls]
                if r.get("forest_tag", "").startswith("forest.pattern:")
            }
        )
        pattern_mapped = [p for p in patterns if p in sssom_objects]
        rule_mapped = rule_curie in sssom_objects
        rule_tagged = rule_curie in assign_tags or any(
            t.startswith(f"forest.pattern:{cls}/") for t in assign_tags
        )
        if rule_mapped or pattern_mapped:
            mapped_classes += 1
        if rule_tagged or any(p in assign_tags for p in patterns):
            taggable_classes += 1
        class_report.append(
            {
                "rule_class": cls,
                "n_smirks": len(by_class[cls]),
                "forest_rule": rule_curie,
                "sssom_rule_mapped": rule_mapped,
                "sssom_subjects": sssom_by_object.get(rule_curie, []),
                "patterns": patterns,
                "patterns_sssom_mapped": pattern_mapped,
                "assignment_taggable": rule_tagged
                or any(p in assign_tags for p in patterns),
            }
        )

    # Pattern-level gaps (harvested pattern tags with no SSSOM row)
    all_pattern_tags = sorted(
        {
            r["forest_tag"]
            for r in rows
            if (r.get("forest_tag") or "").startswith("forest.pattern:")
        }
    )
    unmapped_patterns = [p for p in all_pattern_tags if p not in sssom_objects]

    report = {
        "source_harvest": str(HARVEST.relative_to(ROOT)),
        "n_smirks": len(rows),
        "n_rule_classes": len(by_class),
        "n_pattern_tags": len(all_pattern_tags),
        "classes_with_sssom": mapped_classes,
        "classes_assignment_taggable": taggable_classes,
        "unmapped_pattern_count": len(unmapped_patterns),
        "unmapped_patterns_sample": unmapped_patterns[:30],
        "classes": class_report,
        "policy": {
            "namer_imports_forest": False,
            "harvest_is_text_parse_of_rules_py": True,
            "use": (
                "Validate Forest→XRM SSSOM/tag coverage; promote reviewed SMARTS "
                "into assignments; do not link forest into the namer."
            ),
        },
    }
    OUT.write_text(json.dumps(report, indent=2) + "\n")
    print(
        f"wrote {OUT}: classes={len(by_class)} "
        f"sssom_covered={mapped_classes} taggable={taggable_classes} "
        f"unmapped_patterns={len(unmapped_patterns)}"
    )


if __name__ == "__main__":
    main()
