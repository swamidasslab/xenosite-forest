#!/usr/bin/env python3
"""Run XRM competency questions (ontology / mapping / SPARQL / SHACL).

Tagging-layer CQs are exported to fixtures and asserted in Rust
``tests/competency.rs`` (uses the real namer).

Usage:
  python3 crates/xenosite-xrm/tools/run_competency_tests.py
  python3 crates/xenosite-xrm/tools/run_competency_tests.py --refresh-ttl
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
CQ = ROOT / "data/competency/cq.yml"
SKOS_JSON = ROOT / "data/ontology/xrm.skos.jsonld"
SKOS_TTL = ROOT / "data/ontology/xrm.skos.ttl"
SHACL = ROOT / "data/ontology/xrm.shacl.ttl"
MAP_DIR = ROOT / "data/mappings"
FIXTURES = ROOT / "data/competency/fixtures/reactions.jsonl"
EXPECTED_TAGS = ROOT / "data/competency/fixtures/expected_tags.jsonl"
SPARQL_DIR = ROOT / "data/competency"


def load_skos():
    data = json.loads(SKOS_JSON.read_text())
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


def iri_to_curie(term) -> str:
    s = str(term)
    prefix = "https://xenosite.org/ontology/xrm#"
    if s.startswith(prefix):
        return "xrm:" + s[len(prefix) :]
    if s.startswith("xrm:"):
        return s
    m = re.search(r"xrm[:#]([0-9A-Za-z_.:/-]+)$", s)
    if m:
        return "xrm:" + m.group(1).split("/")[-1]
    return s


def run_sparql(query: str | None = None, query_path: Path | None = None) -> set[str]:
    from rdflib import Graph

    if not SKOS_TTL.exists():
        raise AssertionError(f"missing {SKOS_TTL}; run tools/jsonld_to_ttl.py")
    g = Graph()
    g.parse(SKOS_TTL, format="turtle")
    if query_path is not None:
        q = query_path.read_text()
    else:
        q = query or ""
    out = set()
    for row in g.query(q):
        # first projected variable is the term id
        out.add(iri_to_curie(row[0]))
    return out


def run_shacl(data_path: Path, shapes_path: Path, allow_warnings: bool) -> None:
    from rdflib import Graph
    from pyshacl import validate

    data = Graph()
    data.parse(data_path, format="turtle")
    shapes = Graph()
    shapes.parse(shapes_path, format="turtle")
    ok, _, text = validate(
        data,
        shacl_graph=shapes,
        inference="rdfs",
        abort_on_first=False,
        allow_infos=True,
        allow_warnings=allow_warnings,
    )
    if not ok:
        raise AssertionError(f"SHACL failed:\n{text[:2500]}")


def refresh_artifacts() -> None:
    import subprocess

    subprocess.check_call([sys.executable, str(ROOT / "tools/jsonld_to_ttl.py")])
    subprocess.check_call([sys.executable, str(ROOT / "tools/export_mapping_views.py")])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh-ttl", action="store_true", help="Regenerate TTL + mapping views first")
    args = ap.parse_args()
    if args.refresh_ttl or not SKOS_TTL.exists():
        refresh_artifacts()

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
                labels = {by_id[c].get("prefLabel") for c in desc if c in by_id}
                for lab in q.get("expected_contains_labels") or []:
                    if lab not in labels:
                        raise AssertionError(
                            f"missing label under {under}: {lab} (have {len(labels)})"
                        )
                for lab in q.get("expected_excludes_labels") or []:
                    if lab in labels:
                        raise AssertionError(f"excluded label present: {lab}")
            elif typ == "ontology_broader":
                anc = ancestors(by_id, q["concept"])
                for a in q.get("expected_ancestors") or []:
                    if a not in anc:
                        raise AssertionError(
                            f"{q['concept']} missing ancestor {a}; have {sorted(anc)[:12]}"
                        )
            elif typ == "mapping_exists":
                subj = q["subject"]
                obj = q.get("object")
                prefix = q.get("object_prefix")
                hits = [o for s, o in sssom if s == subj]
                if obj and obj not in hits:
                    raise AssertionError(f"no SSSOM {subj} → {obj}; have {hits[:8]}")
                if prefix and not any(o.startswith(prefix) for o in hits):
                    raise AssertionError(f"no SSSOM {subj} → {prefix}* ; have {hits[:8]}")
            elif typ == "sparql":
                if q.get("query_text"):
                    got = run_sparql(query=q["query_text"])
                else:
                    got = run_sparql(query_path=SPARQL_DIR / q["query"])
                for need in q.get("expected_contains") or []:
                    if need not in got:
                        raise AssertionError(
                            f"SPARQL missing {need}; got {sorted(got)[:20]} (n={len(got)})"
                        )
                for bad in q.get("expected_excludes") or []:
                    if bad in got:
                        raise AssertionError(f"SPARQL unexpectedly contains {bad}")
                if q.get("min_results") is not None and len(got) < q["min_results"]:
                    raise AssertionError(f"SPARQL min_results {q['min_results']} > {len(got)}")
            elif typ == "shacl":
                data = ROOT / "data/competency" / q.get("data", "../ontology/xrm.skos.ttl")
                shapes = ROOT / "data/competency" / q.get("shapes", "../ontology/xrm.shacl.ttl")
                # resolve .. paths cleanly
                data = (SPARQL_DIR / q["data"]).resolve() if "data" in q else SKOS_TTL
                shapes = (SPARQL_DIR / q["shapes"]).resolve() if "shapes" in q else SHACL
                run_shacl(data, shapes, allow_warnings=bool(q.get("allow_warnings", True)))
            elif typ == "definition_coverage":
                total = len(by_id)
                with_def = sum(1 for n in by_id.values() if n.get("definition"))
                frac = with_def / total if total else 0.0
                need = float(q.get("min_fraction", 0.8))
                if frac < need:
                    raise AssertionError(
                        f"definition coverage {frac:.3f} < {need:.3f} ({with_def}/{total})"
                    )
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
                        "expect_site_label_substring": q.get("expect_site_label_substring"),
                        "expected_site_kind": q.get("expected_site_kind"),
                        "expected_site_environment": q.get("expected_site_environment") or [],
                    }
                )
            else:
                raise AssertionError(f"unknown type {typ}")
            passed += 1
            print(f"PASS {qid}")
        except Exception as e:
            failed.append((qid, str(e)))
            print(f"FAIL {qid}: {e}")

    FIXTURES.parent.mkdir(parents=True, exist_ok=True)
    with FIXTURES.open("w") as f:
        for row in tagging_fixtures:
            f.write(json.dumps(row) + "\n")
    with EXPECTED_TAGS.open("w") as f:
        for row in tagging_fixtures:
            f.write(
                json.dumps(
                    {
                        "id": row["id"],
                        "expected_tags": row["expected_labels"],
                        "expected_excludes": row["expected_excludes_labels"],
                        "expect_site_map": row["expect_site_map"],
                    }
                )
                + "\n"
            )
    print(f"wrote {len(tagging_fixtures)} tagging fixtures → {FIXTURES}")
    print(f"wrote expected tags → {EXPECTED_TAGS}")

    print(
        f"\n{passed} checks ok, {len(failed)} failed "
        f"(ontology/mapping/sparql/shacl; tagging deferred to Rust)"
    )
    if failed:
        for qid, err in failed:
            print(f"  {qid}: {err}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
