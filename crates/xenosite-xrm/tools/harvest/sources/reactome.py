"""Reactome pathway names for phase II / xenobiotic conjugation."""

from __future__ import annotations

import json
import urllib.request

from .base import candidate

# Stable pathway ids from Reactome content service.
PATHWAYS = [
    "R-HSA-156588",  # Glucuronidation
    "R-HSA-156580",  # Phase II - Conjugation of compounds
    "R-HSA-211859",  # Biological oxidations
    "R-HSA-211945",  # Phase I - Functionalization of compounds
    "R-HSA-156584",  # Cytosolic sulfonation of small molecules
    "R-HSA-174403",  # Glutathione conjugation
]


def harvest(cfg=None):
    rows = []
    for pid in PATHWAYS:
        url = f"https://reactome.org/ContentService/data/query/{pid}"
        try:
            req = urllib.request.Request(
                url,
                headers={
                    "Accept": "application/json",
                    "User-Agent": "xenosite-xrm-harvest/0.1",
                },
            )
            with urllib.request.urlopen(req, timeout=45) as resp:
                d = json.loads(resp.read().decode())
        except Exception as e:  # noqa: BLE001
            rows.append(
                candidate(
                    candidate_id=f"cand:reactome:error:{pid}",
                    pref_label=f"reactome fetch error {pid}",
                    notes=str(e),
                    sources=[{"system": "reactome", "id": pid, "url": url}],
                )
            )
            continue
        name = d.get("displayName") or d.get("name") or pid
        syns = [pid]
        if d.get("stId"):
            syns.append(d["stId"])
        definition = d.get("summation")
        if isinstance(definition, list) and definition:
            # summation entries often have 'text'
            texts = []
            for s in definition:
                if isinstance(s, dict) and s.get("text"):
                    texts.append(s["text"])
            definition = " ".join(texts)[:500] if texts else None
        rows.append(
            candidate(
                candidate_id=f"cand:reactome:{pid}",
                pref_label=name,
                synonyms=syns,
                definition=definition if isinstance(definition, str) else None,
                sources=[
                    {
                        "system": "reactome",
                        "id": pid,
                        "url": f"https://reactome.org/content/detail/{pid}",
                    }
                ],
                suggested_spines=["chemist reaction type", "phase II"],
                notes="Reactome pathway; good for phase-level and conjugation synonyms",
            )
        )
    return rows
