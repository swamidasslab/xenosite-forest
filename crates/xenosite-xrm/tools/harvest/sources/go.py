"""Gene Ontology xenobiotic metabolism process terms via QuickGO."""

from __future__ import annotations

import json
import urllib.parse
import urllib.request

from .base import candidate

# QuickGO search for xenobiotic-related biological processes.
QUERIES = [
    "xenobiotic metabolic process",
    "drug metabolic process",
    "glucuronidation",
    "glutathione conjugation",
    "epoxygenase P450 pathway",
]


def harvest(cfg=None):
    limit = int((cfg or {}).get("go_limit", 40))
    rows = []
    for q in QUERIES:
        params = urllib.parse.urlencode({"query": q, "limit": "20", "page": "1"})
        url = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/search?{params}"
        req_headers = {
            "Accept": "application/json",
            "User-Agent": "xenosite-xrm-harvest/0.1",
        }
        try:
            request = urllib.request.Request(url, headers=req_headers)
            with urllib.request.urlopen(request, timeout=45) as resp:
                data = json.loads(resp.read().decode())
        except Exception as e:  # noqa: BLE001
            rows.append(
                candidate(
                    candidate_id=f"cand:go:error:{q.replace(' ', '_')}",
                    pref_label=f"go fetch error {q}",
                    notes=str(e),
                    sources=[{"system": "go", "id": q, "url": url}],
                )
            )
            continue
        results = data.get("results") or data.get("response", {}).get("docs") or []
        for hit in results:
            go_id = hit.get("id") or hit.get("goId") or ""
            name = hit.get("name") or hit.get("label") or go_id
            definition = None
            if isinstance(hit.get("definition"), dict):
                definition = hit["definition"].get("text")
            elif isinstance(hit.get("definition"), str):
                definition = hit["definition"]
            rows.append(
                candidate(
                    candidate_id=f"cand:go:{go_id.replace(':', '_')}",
                    pref_label=name,
                    synonyms=[go_id, q],
                    definition=definition,
                    sources=[
                        {
                            "system": "go",
                            "id": go_id,
                            "url": f"https://www.ebi.ac.uk/QuickGO/term/{go_id}",
                        }
                    ],
                    suggested_spines=["chemist reaction type", "bioactivation"],
                    notes="GO process term; use as broader/relatedMatch, not enzyme prefLabel",
                )
            )
            if len(rows) >= limit:
                return rows
    return rows
