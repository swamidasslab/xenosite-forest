#!/usr/bin/env python3
"""Harvest XRM candidate terms / synonyms / examples from public sources.

Stdlib only. See README.md in this directory.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from sources import REGISTRY  # noqa: E402
from sources.base import merge_candidates  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--sources",
        default="seed,chebi,pubchem,kegg,rhea,go,reactome",
        help="Comma-separated source keys from: " + ",".join(sorted(REGISTRY)),
    )
    p.add_argument(
        "--out",
        type=Path,
        default=ROOT.parents[1] / "data" / "candidates" / "round-001.jsonl",
        help="Output JSONL path",
    )
    p.add_argument("--chebi-limit", type=int, default=12)
    p.add_argument("--pubchem-limit", type=int, default=14)
    p.add_argument("--kegg-limit", type=int, default=60)
    p.add_argument("--rhea-limit", type=int, default=40)
    p.add_argument("--go-limit", type=int, default=30)
    p.add_argument("--no-merge", action="store_true", help="Do not merge by pref_label")
    args = p.parse_args(argv)

    cfg = {
        "chebi_limit": args.chebi_limit,
        "pubchem_limit": args.pubchem_limit,
        "kegg_limit": args.kegg_limit,
        "rhea_limit": args.rhea_limit,
        "go_limit": args.go_limit,
    }

    wanted = [s.strip() for s in args.sources.split(",") if s.strip()]
    unknown = [s for s in wanted if s not in REGISTRY]
    if unknown:
        print("unknown sources:", ", ".join(unknown), file=sys.stderr)
        print("known:", ", ".join(sorted(REGISTRY)), file=sys.stderr)
        return 2

    rows = []
    for key in wanted:
        print(f"harvesting {key}…", file=sys.stderr)
        try:
            part = REGISTRY[key](cfg)
        except Exception as e:  # noqa: BLE001
            print(f"  FAILED {key}: {e}", file=sys.stderr)
            continue
        print(f"  {len(part)} records", file=sys.stderr)
        rows.extend(part)

    if not args.no_merge:
        rows = merge_candidates(rows)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8") as fh:
        for row in rows:
            fh.write(json.dumps(row, ensure_ascii=False) + "\n")
    print(f"wrote {len(rows)} candidates → {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
