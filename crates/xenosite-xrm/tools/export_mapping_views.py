#!/usr/bin/env python3
"""Derive compact mapping views (mesh/kegg/rxno_mop) from SSSOM TSVs."""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAP = ROOT / "data/mappings"
OUT = ROOT / "data/mappings/views"


def rows(path: Path):
    for line in path.read_text().splitlines():
        if not line or line.startswith("#") or line.startswith("subject"):
            continue
        parts = line.split("\t")
        if len(parts) >= 3:
            yield parts[0], parts[1], parts[2], parts[4] if len(parts) > 4 else ""


def write_view(name: str, predicate_filter, object_prefixes: tuple[str, ...]) -> int:
    out_path = OUT / name
    lines = ["subject_id\tpredicate_id\tobject_id\tsubject_label"]
    n = 0
    for path in sorted(MAP.glob("*.sssom.tsv")):
        for s, p, o, lab in rows(path):
            if object_prefixes and not any(o.startswith(pref) for pref in object_prefixes):
                continue
            if predicate_filter and p not in predicate_filter:
                continue
            lines.append(f"{s}\t{p}\t{o}\t{lab}")
            n += 1
    out_path.write_text("\n".join(lines) + "\n")
    return n


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    counts = {
        "mesh.tsv": write_view("mesh.tsv", None, ("mesh:",)),
        "kegg.tsv": write_view(
            "kegg.tsv",
            None,
            ("kegg:", "KEGG:", "kegg.reaction:", "kegg.pathway:", "kegg.rclass:", "kegg.enzyme:"),
        ),
        "rxno_mop.tsv": write_view("rxno_mop.tsv", None, ("rxno:", "mop:", "RXNO:", "MOP:")),
        "go_chebi_rhea.tsv": write_view(
            "go_chebi_rhea.tsv", None, ("GO:", "CHEBI:", "chebi:", "RHEA:", "rhea:")
        ),
    }
    for name, n in counts.items():
        print(f"wrote {OUT / name} ({n} rows)")


if __name__ == "__main__":
    main()
