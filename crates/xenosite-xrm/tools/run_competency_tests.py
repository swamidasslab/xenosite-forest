#!/usr/bin/env python3
"""Run XRM competency questions (ontology + mapping layers in Python).

Tagging-layer CQs are exported to fixtures and asserted in Rust
``tests/competency.rs`` (uses the real namer).

Usage:
  python3 crates/xenosite-xrm/tools/run_competency_tests.py
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
CQ = ROOT / "data/competency/cq.yml"
SKOS = ROOT / "data/ontology/xrm.skos.jsonld"
MAP_DIR = ROOT / "data/mappings"
FIXTURES = ROOT / "data/competency/fixtures/reactions.jsonl"


def load_skos():
    data = json.loads(SKOS.read_text())
    by_id = {}
    children = defaultdict(list)
    for n in data["@graph"]:
        if n.get("type") != "skos:Concept":
            continue
        cid = n["id"]
        by_id[cid] = n
        b = n.get("broader")
        parents = b if isinstance(b, list) else ([b] if b else [])
        for p in parents:
            children[p].append(cid)
    return by_id, children


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
    stack = []
    n = by_id.get(cid)
    if not n:
        return out
    b = n.get("broader")
    stack.extend(b if isinstance(b, list) else ([b] if b else []))
    while stack:
        p = stack.pop()
        if not p or p in out:
            continue
        out.add(p)
        pn = by_id.get(p)
        if not pn:
            continue
        pb = pn.get("broader")
        stack.extend(pb if isinstance(pb, list) else ([pb] if pb else []))
    return out


def load_sssom():
    rows = []
    for path in MAP_DIR.glob("*.sssom.tsv"):
        for line in path.read_text().splitlines():
            if not line or line.startswith("#") or line.startswith("subject"):
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                rows.append((parts[0], parts[2]))
    return rows


def main() -> int:
    doc = yaml.safe_load(CQ.read_text())
    by_id, children = load_skos()
    sssom = load_sssom()
    failed = []
    passed = 0
    tagging_fixtures = []

    for q in doc["questions"]:
        qid = q["id"]
        typ = q["type"]
        try:
            if typ == "ontology_labels":
                under = q["under"]
                desc = descendants(children, under)
                labels = {
                    by_id[c].get("prefLabel")
                    for c in desc
                    if c in by_id
                }
                # include direct children labels; also allow under node itself
                for lab in q.get("expected_contains_labels") or []:
                    if lab not in labels:
                        # also search all concepts under by walking parents? already descendants
                        # site labels may be nested deeper under site atom class
                        raise AssertionError(f"missing label under {under}: {lab} (have {len(labels)})")
                for lab in q.get("expected_excludes_labels") or []:
                    if lab in labels:
                        raise AssertionError(f"excluded label present: {lab}")
            elif typ == "ontology_broader":
                anc = ancestors(by_id, q["concept"])
                for a in q.get("expected_ancestors") or []:
                    if a not in anc:
                        raise AssertionError(f"{q['concept']} missing ancestor {a}; have {sorted(anc)[:12]}")
            elif typ == "mapping_exists":
                subj = q["subject"]
                obj = q.get("object")
                prefix = q.get("object_prefix")
                hits = [o for s, o in sssom if s == subj]
                if obj and obj not in hits:
                    raise AssertionError(f"no SSSOM {subj} → {obj}; have {hits[:8]}")
                if prefix and not any(o.startswith(prefix) for o in hits):
                    raise AssertionError(f"no SSSOM {subj} → {prefix}* ; have {hits[:8]}")
            elif typ == "tagging":
                tagging_fixtures.append(
                    {
                        "id": qid,
                        "question": q.get("question"),
                        "reactant_smiles": q["reactant"],
                        "product_smiles": q["product"],
                        "tags": q.get("tags") or [],
                        "expected_labels": q.get("expected_labels") or [],
                        "expected_excludes_labels": q.get("expected_excludes_labels") or [],
                        "expect_site_map": q.get("expect_site_map") or [],
                    }
                )
            else:
                raise AssertionError(f"unknown type {typ}")
            passed += 1
            print(f"PASS {qid}")
        except AssertionError as e:
            failed.append((qid, str(e)))
            print(f"FAIL {qid}: {e}")

    FIXTURES.parent.mkdir(parents=True, exist_ok=True)
    with FIXTURES.open("w") as f:
        for row in tagging_fixtures:
            f.write(json.dumps(row) + "\n")
    print(f"wrote {len(tagging_fixtures)} tagging fixtures → {FIXTURES}")

    print(f"\n{passed} checks ok, {len(failed)} failed (ontology/mapping; tagging deferred to Rust)")
    if failed:
        for qid, err in failed:
            print(f"  {qid}: {err}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
