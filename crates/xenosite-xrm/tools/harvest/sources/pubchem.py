"""PubChem CID → SMILES/title enrichment for seed examples."""

from __future__ import annotations

import re

from .base import candidate, http_get_json

# CIDs used in seed examples and common xenobiotic pairs.
CIDS = [
    241,  # benzene
    996,  # phenol
    1140,  # toluene
    244,  # benzyl alcohol
    7519,  # anisole
    702,  # ethanol
    177,  # acetaldehyde
    6325,  # ethene
    6354,  # oxirane
    1983,  # acetaminophen
    243,  # benzoic acid
    464,  # hippuric acid
    6324,  # ethane
    297,  # methane? may be wrong
    887,  # methanol
]


def harvest(cfg=None):
    limit = int((cfg or {}).get("pubchem_limit", len(CIDS)))
    rows = []
    for cid in CIDS[:limit]:
        url = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
            f"{cid}/property/ConnectivitySMILES,Title/JSON"
        )
        try:
            data = http_get_json(url)
            prop = data["PropertyTable"]["Properties"][0]
        except Exception as e:  # noqa: BLE001
            rows.append(
                candidate(
                    candidate_id=f"cand:pubchem:error:{cid}",
                    pref_label=f"pubchem fetch error {cid}",
                    notes=str(e),
                    sources=[{"system": "pubchem", "id": f"CID:{cid}", "url": url}],
                )
            )
            continue
        title = prop.get("Title") or f"CID:{cid}"
        smi = prop.get("ConnectivitySMILES") or prop.get("CanonicalSMILES")
        # Suggest a tokenized synonym from the title
        syn = [title]
        slug = re.sub(r"[^a-z0-9]+", " ", title.lower()).strip()
        if slug and slug != title.lower():
            syn.append(slug)
        rows.append(
            candidate(
                candidate_id=f"cand:pubchem:{cid}",
                pref_label=title,
                synonyms=syn,
                sources=[
                    {
                        "system": "pubchem",
                        "id": f"CID:{cid}",
                        "url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}",
                    }
                ],
                examples=[
                    {
                        "reactant_smiles": smi,
                        "product_smiles": None,
                        "label": title,
                        "refs": [{"system": "pubchem", "id": f"CID:{cid}"}],
                    }
                ]
                if smi
                else [],
                suggested_spines=["chemist reaction type"],
                notes="PubChem compound; pair reactant/product CIDs in seed or review",
            )
        )
    return rows
