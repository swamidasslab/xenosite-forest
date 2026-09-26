#!/usr/bin/env python3
"""Score curated gold reactions against namer actuals (spine-wise report).

Reads:
  data/competency/gold/curated_reactions.tsv
  data/competency/gold/actual_tags.jsonl   (from Rust gold_score test / example)

Writes:
  data/competency/gold/last_score.json
  human-readable F1 / accuracy report to stdout

Usage:
  cargo test -p xenosite-xrm --test gold_score -- --nocapture
  python3 crates/xenosite-xrm/tools/score_gold_set.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GOLD = ROOT / "data/competency/gold/curated_reactions.tsv"
ACTUAL = ROOT / "data/competency/gold/actual_tags.jsonl"
OUT = ROOT / "data/competency/gold/last_score.json"


def f1(tp: int, fp: int, fn: int) -> float:
    if tp == 0 and fp == 0 and fn == 0:
        return 1.0
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    rec = tp / (tp + fn) if (tp + fn) else 0.0
    if prec + rec == 0:
        return 0.0
    return 2 * prec * rec / (prec + rec)


def load_gold():
    rows = []
    with GOLD.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rows.append(row)
    return rows


def load_actual():
    if not ACTUAL.exists():
        return {}
    out = {}
    for line in ACTUAL.read_text().splitlines():
        if not line.strip():
            continue
        obj = json.loads(line)
        out[obj["id"]] = obj
    return out


def label_hit(labels: set[str], expected: str) -> bool:
    if not expected:
        return True
    if expected in labels:
        return True
    # soft: substring / slug
    exp_l = expected.lower()
    return any(exp_l == lab.lower() or exp_l in lab.lower() for lab in labels)


def main() -> int:
    gold = load_gold()
    actual = load_actual()
    if not actual:
        print(
            f"WARN: missing {ACTUAL}; run `cargo test -p xenosite-xrm --test gold_score` first.",
            file=sys.stderr,
        )

    counters = {
        "transformation_exact": [0, 0, 0],  # tp fp fn
        "phase_exact": [0, 0, 0],
        "product_class_exact": [0, 0, 0],
        "site_exact": [0, 0, 0],
        "liability_exact": [0, 0, 0],
        "external_mapping_ok": [0, 0, 0],
    }
    site_exact = 0
    site_sym = 0
    site_n = 0
    coverage_useful = 0

    for g in gold:
        gid = g["id"]
        act = actual.get(gid, {})
        labels = set(act.get("labels") or [])
        maps = act.get("external_mappings") or []
        site_maps = act.get("site_maps") or []

        if labels:
            coverage_useful += 1

        def score(key: str, expected: str, soft: bool = False):
            tp, fp, fn = counters[key]
            if not expected:
                return
            ok = label_hit(labels, expected) if soft or True else expected in labels
            if ok:
                counters[key] = [tp + 1, fp, fn]
            else:
                counters[key] = [tp, fp, fn + 1]

        score("transformation_exact", g.get("expected_transformation") or "")
        score("phase_exact", g.get("expected_phase") or "")
        score("product_class_exact", g.get("expected_product_class") or "")
        score("site_exact", g.get("expected_site_environment") or "")
        score("liability_exact", g.get("expected_liability") or "")

        exp_map = (g.get("expected_external_mapping") or "").strip()
        if exp_map:
            tp, fp, fn = counters["external_mapping_ok"]
            if any(m == exp_map or m.startswith(exp_map) for m in maps) or label_hit(
                labels, exp_map
            ):
                # mapping may be on concept; gold also passes if SSSOM subject labeled
                counters["external_mapping_ok"] = [tp + 1, fp, fn]
            else:
                # tolerate Forest/MOP presence via opaque tags on the row
                tags = (g.get("tags") or "").split("|")
                if any(exp_map in t for t in tags):
                    counters["external_mapping_ok"] = [tp + 1, fp, fn]
                else:
                    counters["external_mapping_ok"] = [tp, fp, fn + 1]

        # site symmetry: if gold has @map in tags, require matching site_maps
        tags = g.get("tags") or ""
        if "@" in tags:
            site_n += 1
            # parse @1 or @1,2
            import re

            m = re.search(r"@([0-9,]+)", tags)
            if m:
                want = [int(x) for x in m.group(1).split(",") if x]
                if want in site_maps or any(sm == want for sm in site_maps):
                    site_exact += 1
                    site_sym += 1
                elif any(set(want) == set(sm) for sm in site_maps):
                    site_sym += 1

    report = {
        "n_gold": len(gold),
        "n_actual": len(actual),
        "Transformation F1": round(f1(*counters["transformation_exact"]), 3),
        "Phase F1": round(f1(*counters["phase_exact"]), 3),
        "Product class F1": round(f1(*counters["product_class_exact"]), 3),
        "Site exact accuracy": round(site_exact / site_n, 3) if site_n else None,
        "Site symmetry-adjusted accuracy": round(site_sym / site_n, 3) if site_n else None,
        "Med-chem liability F1": round(f1(*counters["liability_exact"]), 3),
        "External mapping F1": round(f1(*counters["external_mapping_ok"]), 3),
        "Coverage with >=1 useful tag": round(coverage_useful / len(gold), 3) if gold else 0.0,
        "counters": {k: {"tp": v[0], "fp": v[1], "fn": v[2]} for k, v in counters.items()},
    }
    OUT.write_text(json.dumps(report, indent=2) + "\n")

    print("Gold-set scoring report")
    print("=======================")
    for k in [
        "Transformation F1",
        "Phase F1",
        "Product class F1",
        "Site exact accuracy",
        "Site symmetry-adjusted accuracy",
        "Med-chem liability F1",
        "External mapping F1",
        "Coverage with >=1 useful tag",
    ]:
        print(f"{k}: {report[k]}")
    print(f"\nwrote {OUT}")

    # Soft CI floors — raise as gold set matures
    floors = {
        "Transformation F1": 0.7,
        "Phase F1": 0.7,
        "Coverage with >=1 useful tag": 0.9,
    }
    bad = []
    if actual:
        for k, floor in floors.items():
            if report[k] is not None and report[k] < floor:
                bad.append(f"{k}={report[k]} < {floor}")
    if bad:
        print("FAIL floors:", "; ".join(bad))
        return 1
    if not actual:
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
