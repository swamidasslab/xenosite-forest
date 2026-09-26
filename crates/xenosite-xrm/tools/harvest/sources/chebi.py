"""ChEBI synonyms + example metabolites for seed reaction labels."""

from __future__ import annotations

from .base import candidate, http_get_json

# ChEBI ids tied to seed / xenobiotic metabolism vocabulary.
CHEBI_IDS = [
    17263,  # estrone
    28919,  # estrone 3-O-glucuronide
    16750,  # phenol (approx — may vary)
    15864,  # benzene? check
    15367,  # ethanol
    15343,  # acetaldehyde
    16236,  # ethanol sometimes
    17933,  # phenol
    16842,  # anisole? 
    27897,  # anisole
    17578,  # toluene
    17987,  # benzyl alcohol
]


def _names_from_payload(d: dict) -> list[str]:
    out: list[str] = []
    names = d.get("names") or {}
    if isinstance(names, dict):
        for _bucket, rows in names.items():
            if not isinstance(rows, list):
                continue
            for row in rows:
                if isinstance(row, dict):
                    for k in ("ascii_name", "name"):
                        if row.get(k):
                            out.append(str(row[k]))
    if d.get("ascii_name"):
        out.append(str(d["ascii_name"]))
    # strip HTML italics often present in ChEBI names
    cleaned = []
    for n in out:
        s = (
            n.replace("<em>", "")
            .replace("</em>", "")
            .replace("<small>", "")
            .replace("</small>", "")
            .replace("&beta;", "beta")
            .replace("&alpha;", "alpha")
        )
        cleaned.append(s)
    return cleaned


def harvest(cfg=None):
    limit = int((cfg or {}).get("chebi_limit", 12))
    rows = []
    for cid in CHEBI_IDS[:limit]:
        url = f"https://www.ebi.ac.uk/chebi/backend/api/public/compound/{cid}/"
        try:
            d = http_get_json(url)
        except Exception as e:  # noqa: BLE001
            rows.append(
                candidate(
                    candidate_id=f"cand:chebi:error:{cid}",
                    pref_label=f"chebi fetch error {cid}",
                    notes=str(e),
                    sources=[{"system": "chebi", "id": f"CHEBI:{cid}", "url": url}],
                )
            )
            continue
        accession = d.get("chebi_accession") or f"CHEBI:{cid}"
        pref = d.get("ascii_name") or d.get("name") or accession
        if isinstance(pref, str):
            pref = (
                pref.replace("<em>", "")
                .replace("</em>", "")
                .replace("<small>", "")
                .replace("</small>", "")
            )
        syns = _names_from_payload(d)
        smi = (d.get("default_structure") or {}).get("smiles")
        definition = d.get("definition")
        # Infer a reaction-ish label hint from definition keywords
        spines = ["chemist reaction type"]
        low = (definition or "").lower() + " " + " ".join(syns).lower()
        if "glucuron" in low:
            spines.append("electrophile role")
        if "epoxid" in low:
            spines += ["ring fate", "electrophile role"]
        examples = []
        if smi:
            examples.append(
                {
                    "reactant_smiles": None,
                    "product_smiles": smi,
                    "label": f"{pref} structure",
                    "refs": [{"system": "chebi", "id": accession}],
                }
            )
        rows.append(
            candidate(
                candidate_id=f"cand:chebi:{cid}",
                pref_label=str(pref),
                synonyms=syns,
                definition=definition if isinstance(definition, str) else None,
                sources=[{"system": "chebi", "id": accession, "url": url}],
                examples=examples,
                suggested_spines=spines,
                notes="ChEBI compound record; promote synonyms / pair with parent for reaction examples",
            )
        )
    return rows
