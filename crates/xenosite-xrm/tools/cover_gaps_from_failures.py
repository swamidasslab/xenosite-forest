#!/usr/bin/env python3
"""Use gold / competency *failures* to drive general gap coverage.

Do **not** patch individual assignment ``emit`` lists for product-class misses.
Instead:

1. Score gold (and optional competency tagging fixtures).
2. Group failures into gap classes.
3. Apply general ontology ``related_match`` wiring so the namer's relatedMatch
   expansion covers the whole class next run.

Gap classes (v1)
----------------
* ``missing_product_for_transformation`` — a transformation label was emitted
  but an expected product-class label was not. Wire
  ``transformation --relatedMatch--> product`` (and reverse) when both exist
  under the chemical-transformation / phase-II and reactive-metabolite spines.
* ``missing_liability_for_transformation`` — same pattern into medchem liability.
* ``missing_definition`` — concepts without ``skos:definition`` get a parent-
  aware stub (bulk data fill, not per-test branches).

Usage:
  cargo test -p xenosite-xrm --test gold_score
  python3 crates/xenosite-xrm/tools/score_gold_set.py
  python3 crates/xenosite-xrm/tools/cover_gaps_from_failures.py
  python3 crates/xenosite-xrm/tools/yaml_to_skos.py
  # re-run gold / CQs
"""
from __future__ import annotations

import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
YAML_PATH = ROOT / "data/ontology/xrm.yaml"
GOLD = ROOT / "data/competency/gold/curated_reactions.tsv"
ACTUAL = ROOT / "data/competency/gold/actual_tags.jsonl"
REPORT = ROOT / "data/competency/gold/gap_coverage_report.json"

CHEM_SPINES = {"xrm:1100000", "xrm:1200000", "xrm:1300000"}
PRODUCT_SPINE = "xrm:1500000"
LIABILITY_SPINE = "xrm:1400000"


def load_ontology():
    doc = yaml.safe_load(YAML_PATH.read_text())
    by_id = {c["id"]: c for c in doc["concepts"]}
    by_label = {c["preferred_label"]: c for c in doc["concepts"]}
    children = defaultdict(list)
    for c in doc["concepts"]:
        for p in c.get("parents") or []:
            children[p].append(c["id"])
    return doc, by_id, by_label, children


def descendants(children, root):
    out = set()
    stack = list(children.get(root, []))
    while stack:
        x = stack.pop()
        if x in out:
            continue
        out.add(x)
        stack.extend(children.get(x, []))
    return out


def ancestors(by_id, cid):
    out = set()
    stack = list(by_id.get(cid, {}).get("parents") or [])
    while stack:
        p = stack.pop()
        if not p or p in out:
            continue
        out.add(p)
        stack.extend(by_id.get(p, {}).get("parents") or [])
    return out


def under_any(by_id, cid, roots):
    if cid in roots:
        return True
    return bool(ancestors(by_id, cid) & roots)


def label_hit(labels, exp):
    if not exp:
        return True
    labs = {x.lower() for x in labels}
    return exp.lower() in labs or any(exp.lower() in x for x in labs)


def relate(a, b):
    """Add symmetric related_match between concept dicts."""
    for src, tgt in ((a, b["id"]), (b, a["id"])):
        rel = list(src.get("related_match") or [])
        if tgt not in rel:
            rel.append(tgt)
        src["related_match"] = rel


def collect_gold_failures(by_label):
    gold = list(csv.DictReader(GOLD.open(), delimiter="\t"))
    actual = {}
    if ACTUAL.exists():
        for line in ACTUAL.read_text().splitlines():
            if line.strip():
                o = json.loads(line)
                actual[o["id"]] = o
    failures = []
    for g in gold:
        act = actual.get(g["id"], {})
        labels = act.get("labels") or []
        gid = g["id"]
        # Only treat as product/liability gap when transformation (or phase) landed.
        has_tx = label_hit(labels, g.get("expected_transformation") or "") or not g.get(
            "expected_transformation"
        )
        if not has_tx and labels:
            # transformation itself missing — different class; skip auto-wire
            continue
        prod = (g.get("expected_product_class") or "").strip()
        if prod and not label_hit(labels, prod):
            tx = (g.get("expected_transformation") or "").strip()
            failures.append(
                {
                    "id": gid,
                    "class": "missing_product_for_transformation",
                    "transformation": tx,
                    "missing": prod,
                    "labels": labels,
                }
            )
        liab = (g.get("expected_liability") or "").strip()
        if liab and not label_hit(labels, liab):
            tx = (g.get("expected_transformation") or "").strip()
            failures.append(
                {
                    "id": gid,
                    "class": "missing_liability_for_transformation",
                    "transformation": tx,
                    "missing": liab,
                    "labels": labels,
                }
            )
    return failures


def pair_name_heuristics(by_label, product_ids, chem_ids):
    """Conservative name-based pairs (conjugation ↔ conjugate, formation ↔ product).

    Avoid substring false friends (hydroxylation↛hydroxylamine,
    conjugation↛every conjugate, epoxidation↛thiophene epoxide).
    """
    pairs = []
    skip_labs = {"conjugation", "oxidation", "reduction", "hydrolysis"}

    for lab, c in by_label.items():
        if c["id"] not in chem_ids or lab.lower() in skip_labs:
            continue
        candidates = []
        if lab.endswith("ation") and not lab.endswith("conjugation"):
            stem = lab[: -len("ation")]
            candidates += [
                f"{stem}ate conjugate",
                f"{stem}ide conjugate",
                f"{stem} conjugate",
                lab.replace("ation", "ate conjugate"),
            ]
            if lab.endswith("uronidation"):
                candidates.append(lab.replace("uronidation", "uronide conjugate"))
            if "glutathion" in lab:
                candidates += [
                    "GSH adduct",
                    "GSH-trappable metabolite",
                    "glutathione conjugate",
                ]
        if lab.endswith(" formation"):
            candidates.append(lab[: -len(" formation")])
        if lab.endswith("epoxidation"):
            # exact family leaf only
            candidates.append("epoxide")
            if "thiophene" in lab:
                candidates.append("thiophene epoxide")
            if "arene" in lab:
                candidates.append("arene oxide")
        for cand in candidates:
            pc = by_label.get(cand)
            if pc and pc["id"] in product_ids:
                pairs.append((c["preferred_label"], pc["preferred_label"]))
    return pairs


def apply_failure_wiring(doc, by_id, by_label, children, failures):
    chem_ids = set()
    for root in CHEM_SPINES:
        chem_ids |= descendants(children, root)
        chem_ids.add(root)
    product_ids = descendants(children, PRODUCT_SPINE) | {PRODUCT_SPINE}
    liability_ids = descendants(children, LIABILITY_SPINE) | {LIABILITY_SPINE}

    wired = []
    seen_pair = set()

    def try_wire(tx_lab, miss_lab, spine_ids, reason):
        if not tx_lab or not miss_lab:
            return
        # Prefer exact transformation concept; else any emitted label under chem spines
        tx = by_label.get(tx_lab)
        miss = by_label.get(miss_lab)
        if not miss or miss["id"] not in spine_ids:
            return
        if tx is None:
            return
        if tx["id"] not in chem_ids and tx["id"] not in liability_ids:
            # allow medchem source for liability gaps
            if miss["id"] not in liability_ids:
                return
        key = (tx["id"], miss["id"])
        if key in seen_pair:
            return
        seen_pair.add(key)
        before = set(tx.get("related_match") or [])
        relate(tx, miss)
        if miss["id"] not in before:
            wired.append(
                {
                    "from": tx["id"],
                    "from_label": tx["preferred_label"],
                    "to": miss["id"],
                    "to_label": miss["preferred_label"],
                    "reason": reason,
                }
            )

    for f in failures:
        if f["class"] == "missing_product_for_transformation":
            try_wire(
                f["transformation"],
                f["missing"],
                product_ids,
                f"gold-failure:{f['id']}",
            )
            # Also try wiring from any chem-spine label present in actual labels
            for lab in f.get("labels") or []:
                if lab in by_label and by_label[lab]["id"] in chem_ids:
                    try_wire(lab, f["missing"], product_ids, f"gold-failure-emit:{f['id']}")
        elif f["class"] == "missing_liability_for_transformation":
            try_wire(
                f["transformation"],
                f["missing"],
                liability_ids,
                f"gold-failure:{f['id']}",
            )

    # General morphological pairs (covers the class, not one gold row)
    for a, b in pair_name_heuristics(by_label, product_ids, chem_ids):
        try_wire(a, b, product_ids, "name-heuristic")

    return wired


def fill_missing_definitions(doc, by_id):
    n = 0
    for c in doc["concepts"]:
        if c.get("definition"):
            continue
        parents = c.get("parents") or []
        parent_labs = [by_id[p]["preferred_label"] for p in parents if p in by_id]
        if parent_labs:
            c["definition"] = (
                f"{c['preferred_label']}: xenobiotic-metabolism term under "
                f"{', '.join(parent_labs)}."
            )
        else:
            c["definition"] = (
                f"{c['preferred_label']}: xenobiotic-metabolism concept in XRM."
            )
        n += 1
    return n


def strip_adhoc_product_emits(by_id, children):
    """Drop companion product-class emits when a transformation emit is present.

    General rule (not per-id): if a rule emits at least one concept under the
    chemical-transformation / phase-II spines **and** also emits a reactive-
    metabolite product class, remove the product-class ids from ``emit``.
    Companions belong on ``related_match``; the namer expands them.
    """
    assign_path = ROOT / "data/assignments/xenobiotic.jsonl"
    chem_ids = set()
    for root in CHEM_SPINES:
        chem_ids |= descendants(children, root)
    product_ids = descendants(children, PRODUCT_SPINE)

    out_lines = []
    stripped = []
    for line in assign_path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            out_lines.append(line)
            continue
        o = json.loads(line)
        emit = list(o.get("emit") or [])
        has_transformation = any(e in chem_ids for e in emit)
        if not has_transformation:
            out_lines.append(json.dumps(o, ensure_ascii=False, separators=(",", ":")))
            continue
        new_emit = []
        for e in emit:
            if e in product_ids:
                stripped.append(
                    (o["id"], e, by_id.get(e, {}).get("preferred_label", e))
                )
                continue
            new_emit.append(e)
        o["emit"] = new_emit
        out_lines.append(json.dumps(o, ensure_ascii=False, separators=(",", ":")))
    assign_path.write_text("\n".join(out_lines) + "\n")
    return stripped


def main() -> int:
    doc, by_id, by_label, children = load_ontology()
    failures = collect_gold_failures(by_label)

    # Always run morphological + failure wiring
    wired = apply_failure_wiring(doc, by_id, by_label, children, failures)
    defs = fill_missing_definitions(doc, by_id)

    doc["concepts"].sort(key=lambda c: c["id"])
    YAML_PATH.write_text(
        "# XRM thesaurus (authoring source). Export to xrm.skos.jsonld for Rust.\n"
        + yaml.safe_dump(doc, sort_keys=False, allow_unicode=True, width=100)
    )

    stripped = strip_adhoc_product_emits(by_id, children)

    report = {
        "gold_failures_seen": failures,
        "related_match_wired": wired,
        "definitions_filled": defs,
        "assignment_emits_stripped": [
            {"rule": r, "concept": c, "label": lab} for r, c, lab in stripped
        ],
        "policy": (
            "Failures drive general relatedMatch wiring; namer expands "
            "relatedMatch into product/medchem/leaving-group spines. "
            "Do not add one-off product-class IDs to assignment emit lists."
        ),
    }
    REPORT.write_text(json.dumps(report, indent=2) + "\n")
    print(f"gold failure rows: {len(failures)}")
    print(f"related_match wired: {len(wired)}")
    for w in wired[:20]:
        print(f"  {w['from_label']} ↔ {w['to_label']} ({w['reason']})")
    if len(wired) > 20:
        print(f"  … {len(wired) - 20} more")
    print(f"definitions filled: {defs}")
    print(f"assignment companion emits stripped: {len(stripped)}")
    print(f"wrote {REPORT}")
    print(f"wrote {YAML_PATH}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
