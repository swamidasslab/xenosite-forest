"""Rhea reactions for xenobiotic conjugation / oxidation keywords."""

from __future__ import annotations

import urllib.parse

from .base import candidate, http_get_text

QUERIES = [
    "glucuron",
    "hydroxylat",
    "epoxid",
    "glutathion",
    "demethyl",
    "sulfat*transferase",
]


def harvest(cfg=None):
    limit = int((cfg or {}).get("rhea_limit", 40))
    per_q = max(5, limit // max(1, len(QUERIES)))
    rows = []
    for q in QUERIES:
        params = urllib.parse.urlencode(
            {
                "query": q,
                "columns": "rhea-id,equation,ec,chebi-id,go,reaction-xref(KEGG),reaction-xref(Reactome)",
                "format": "tsv",
                "limit": str(per_q),
            }
        )
        url = f"https://www.rhea-db.org/rhea/?{params}"
        try:
            text = http_get_text(url)
        except Exception as e:  # noqa: BLE001
            rows.append(
                candidate(
                    candidate_id=f"cand:rhea:error:{q}",
                    pref_label=f"rhea fetch error {q}",
                    notes=str(e),
                    sources=[{"system": "rhea", "id": q, "url": url}],
                )
            )
            continue
        lines = [ln for ln in text.splitlines() if ln.strip()]
        if not lines:
            continue
        header = lines[0].split("\t")
        for ln in lines[1:]:
            cols = ln.split("\t")
            rec = {header[i]: (cols[i] if i < len(cols) else "") for i in range(len(header))}
            rid = rec.get("Reaction identifier") or rec.get("rhea-id") or ""
            eq = rec.get("Equation") or rec.get("equation") or ""
            if not rid:
                continue
            pref = eq.split(" = ")[0].strip()[:80] if eq else rid
            syns = [eq, rid]
            if rec.get("EC number") or rec.get("ec"):
                syns.append(rec.get("EC number") or rec.get("ec"))
            sources = [
                {
                    "system": "rhea",
                    "id": rid if rid.startswith("RHEA:") else f"RHEA:{rid}",
                    "url": f"https://www.rhea-db.org/rhea/{rid.replace('RHEA:', '')}",
                }
            ]
            for system, key in (
                ("kegg", "Cross-reference (KEGG)"),
                ("reactome", "Cross-reference (Reactome)"),
                ("go", "Gene Ontology"),
            ):
                val = rec.get(key) or rec.get(key.lower()) or ""
                if val:
                    sources.append({"system": system, "id": val.split(";")[0], "url": ""})
            rows.append(
                candidate(
                    candidate_id=f"cand:rhea:{rid.replace(':', '_')}",
                    pref_label=pref or rid,
                    synonyms=[s for s in syns if s],
                    definition=eq or None,
                    sources=sources,
                    examples=[
                        {
                            "reactant_smiles": None,
                            "product_smiles": None,
                            "label": eq,
                            "refs": [{"system": "rhea", "id": rid}],
                        }
                    ],
                    suggested_spines=["chemist reaction type"],
                    notes=f"Rhea query={q}; map ChEBI participants to SMILES for goldens",
                )
            )
            if len(rows) >= limit:
                return rows
    return rows
