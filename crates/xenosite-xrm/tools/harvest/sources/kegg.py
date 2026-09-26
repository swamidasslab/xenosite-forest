"""KEGG reaction names matching xenobiotic / phase I–II keywords."""

from __future__ import annotations

import re

from .base import candidate, http_get_text

KEYWORDS = [
    "hydroxylat",
    "demethyl",
    "dealkyl",
    "epoxid",
    "glucuron",
    "sulfat",  # sulfation / sulfate conjugation noise filtered below
    "glutathion",
    "acetylat",
    "dehalogen",
    "nitroreduct",
    "quinone",
]


def _label_from_name(name: str) -> str:
    # Prefer the short enzyme-style clause before ';' when present.
    head = name.split(";")[0].strip()
    # Strip EC-ish prefixes
    head = re.sub(r"^[A-Za-z0-9\-]+:", "", head).strip()
    return head[:120] if head else name[:120]


def harvest(cfg=None):
    limit = int((cfg or {}).get("kegg_limit", 80))
    text = http_get_text("https://rest.kegg.jp/list/reaction")
    rows = []
    for ln in text.splitlines():
        if "\t" not in ln:
            continue
        rid, name = ln.split("\t", 1)
        low = name.lower()
        if not any(k in low for k in KEYWORDS):
            continue
        # Drop obvious non-xenobiotic sulfate-reduction noise somewhat
        if "sulfat" in low and "conjug" not in low and "transferase" not in low:
            if "hydroxylat" not in low and "glucuron" not in low:
                continue
        pref = _label_from_name(name)
        syns = [name.strip(), rid]
        spines = ["chemist reaction type"]
        if "glucuron" in low:
            spines += ["electrophile role"]
        if "epoxid" in low:
            spines += ["ring fate", "electrophile role"]
        if "demethyl" in low or "dealkyl" in low:
            spines += ["oxygenation outcome", "metabolite cardinality"]
        if "hydroxylat" in low:
            spines += ["redox polarity", "formula-delta class"]
        rows.append(
            candidate(
                candidate_id=f"cand:kegg:{rid}",
                pref_label=pref,
                synonyms=syns,
                definition=name.strip(),
                sources=[
                    {
                        "system": "kegg",
                        "id": rid,
                        "url": f"https://www.kegg.jp/entry/{rid}",
                    }
                ],
                examples=[
                    {
                        "reactant_smiles": None,
                        "product_smiles": None,
                        "label": name.strip(),
                        "refs": [{"system": "kegg", "id": rid}],
                    }
                ],
                suggested_spines=spines,
                notes="KEGG reaction name; resolve COMPOUND SMILES via EQUATION for goldens",
            )
        )
        if len(rows) >= limit:
            break
    return rows
