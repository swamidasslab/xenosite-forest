#!/usr/bin/env python3
"""Fill under-band spines with well-motivated leaves; expand localization templates.

Does not trim high spines. Extra clear terms remain welcome.
"""

from __future__ import annotations

import re
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
YAML_PATH = ROOT / "data/ontology/xrm.yaml"

NEXT = {
    "reactive": 1500200,
    "bio": 2100200,
    "prov": 1900200,
    "loc": 1600200,
    "product_status": 1800200,
}


def next_id(bucket: str) -> str:
    n = NEXT[bucket]
    NEXT[bucket] = n + 1
    return f"xrm:{n:07d}"


def main():
    doc = yaml.safe_load(YAML_PATH.read_text())
    existing = {c["preferred_label"].lower() for c in doc["concepts"]}
    ids = {c["id"] for c in doc["concepts"]}

    def add(pref, parent, definition, synonyms=None, bucket="reactive"):
        if pref.lower() in existing:
            return
        cid = next_id(bucket)
        while cid in ids:
            cid = next_id(bucket)
        entry = {
            "id": cid,
            "preferred_label": pref,
            "definition": definition,
            "parents": [parent] if isinstance(parent, str) else list(parent),
        }
        if synonyms:
            entry["synonyms"] = synonyms
        doc["concepts"].append(entry)
        existing.add(pref.lower())
        ids.add(cid)

    # Product / metabolite class (need ~20 more to reach 70)
    for pref, definition, syn in [
        ("hydroquinone metabolite", "Hydroquinone product class.", []),
        ("aminophenol metabolite", "Aminophenol product class.", []),
        ("nitrosoarenes", "Nitrosoarene metabolite class.", ["nitrosoarene"]),
        ("nitrenium ion", "Nitrenium reactive intermediate.", []),
        ("carbocation intermediate", "Carbocation metabolic intermediate.", []),
        ("benzylic radical", "Benzylic radical intermediate.", []),
        ("thiyl radical", "Thiyl radical intermediate.", []),
        ("sulfenic acid", "Sulfenic acid intermediate/product.", []),
        ("sulfinic acid metabolite", "Sulfinic acid metabolite.", []),
        ("sulfonic acid metabolite", "Sulfonic acid metabolite.", []),
        ("disulfide metabolite", "Disulfide metabolite product.", []),
        ("lactam metabolite", "Lactam ring metabolite.", []),
        ("lactone metabolite", "Lactone ring metabolite.", []),
        ("carboxylic acid glucuronide", "Acyl-linked glucuronide of a carboxylic acid.", []),
        ("ether glucuronide", "O-glucuronide of an alcohol/phenol.", []),
        ("N-glucuronide", "N-linked glucuronide product.", []),
        ("quaternary ammonium glucuronide", "Quaternary N-glucuronide product.", []),
        ("sulfate ester metabolite", "Sulfate ester conjugate product.", []),
        ("cysteine adduct", "Cysteine conjugate / adduct product.", []),
        ("N-acetylcysteine adduct", "NAC / mercapturate-path adduct.", ["NAC adduct"]),
        ("glycine conjugate", "Glycine-conjugated metabolite.", []),
        ("taurine conjugate", "Taurine-conjugated metabolite.", []),
        ("CoA conjugate", "Coenzyme A conjugate intermediate/product.", []),
        ("methylated metabolite", "Methyl-conjugated metabolite product.", []),
        ("acetylated metabolite", "Acetyl-conjugated metabolite product.", []),
        ("dealkylated amine", "Amine product after N-dealkylation.", []),
        ("dealkylated phenol", "Phenol product after O-dealkylation.", []),
        ("nor-metabolite", "N- or O-dealkylated nor-metabolite.", []),
        ("carboxylic acid from alcohol", "Carboxylic acid from alcohol oxidation cascade.", []),
        ("ketone metabolite", "Ketone metabolite product class.", []),
    ]:
        add(pref, "xrm:1500000", definition, syn or None, "reactive")

    # Biological context (need ~5–45 more; aim mid-band)
    for pref, definition, syn in [
        ("CYP3A context", "CYP3A family context (orthogonal).", ["CYP3A"]),
        ("CYP2D6 context", "CYP2D6 isoform context (orthogonal).", []),
        ("CYP2C9 context", "CYP2C9 isoform context (orthogonal).", []),
        ("CYP2C19 context", "CYP2C19 isoform context (orthogonal).", []),
        ("CYP1A2 context", "CYP1A2 isoform context (orthogonal).", []),
        ("CYP2E1 context", "CYP2E1 isoform context (orthogonal).", []),
        ("UGT1A1 context", "UGT1A1 isoform context (orthogonal).", []),
        ("UGT2B7 context", "UGT2B7 isoform context (orthogonal).", []),
        ("SULT1A1 context", "SULT1A1 isoform context (orthogonal).", []),
        ("GST-P context", "GST-P family context (orthogonal).", []),
        ("NAT2 context", "NAT2 isoform context (orthogonal).", []),
        ("FMO3 context", "FMO3 isoform context (orthogonal).", []),
        ("skin tissue context", "Skin tissue framing.", []),
        ("placenta tissue context", "Placenta tissue framing.", []),
        ("adrenal tissue context", "Adrenal tissue framing.", []),
        ("rabbit species context", "Rabbit species framing.", []),
        ("minipig species context", "Minipig species framing.", []),
        ("serum matrix context", "Serum matrix framing.", []),
        ("whole blood matrix context", "Whole-blood matrix framing.", []),
        ("tissue homogenate matrix", "Tissue homogenate incubation matrix.", []),
        ("recombinant enzyme matrix", "Recombinant enzyme incubation matrix.", []),
        ("supersome matrix context", "Insect-cell supersome matrix.", []),
        ("lysosome compartment", "Lysosomal compartment framing.", []),
        ("peroxisome compartment", "Peroxisome compartment framing.", []),
        ("nucleus compartment", "Nuclear compartment framing.", []),
        ("extracellular compartment", "Extracellular compartment framing.", []),
        ("male sex context", "Male biological sex framing.", []),
        ("female sex context", "Female biological sex framing.", []),
        ("pediatric age context", "Pediatric age framing.", []),
        ("geriatric age context", "Geriatric age framing.", []),
        ("fed state context", "Fed-state ADME framing.", []),
        ("fasted state context", "Fasted-state ADME framing.", []),
    ]:
        add(pref, "xrm:2100000", definition, syn or None, "bio")

    # Rule / model provenance (need ~8+)
    for pref, definition, syn in [
        ("repo rule derived", "Tag derived from a repository structure rule.", ["forest rule derived"]),
        ("SMARTS pattern match", "Tag from a matched SMARTS pattern.", []),
        ("formula delta match", "Tag from elemental formula delta.", []),
        ("caller CURIE match", "Tag from an opaque caller CURIE/tag.", []),
        ("Rainbow SO class map", "Mapped to Rainbow stable oxygenation class.", []),
        ("Rainbow UO class map", "Mapped to Rainbow unstable oxygenation class.", []),
        ("Rainbow DH class map", "Mapped to Rainbow dehydrogenation class.", []),
        ("Rainbow HD class map", "Mapped to Rainbow hydrolysis class.", []),
        ("Rainbow RD class map", "Mapped to Rainbow reduction class.", []),
        ("XenoNet edge score", "Provenance from a XenoNet/Forest path edge score.", []),
        ("XenoNet path score", "Provenance from a multi-step path score.", []),
        ("legacy SOM model output", "Mapped from a legacy site-of-metabolism model.", []),
        ("legacy reactivity model output", "Mapped from a legacy reactivity model.", []),
        ("exclusion SMARTS fired", "Negative/exclusion SMARTS suppressed a tag.", []),
        ("priority tie-break", "Tag chosen by priority among competing rules.", []),
        ("manual curator override", "Curator overrode automatic tagging.", []),
        ("harvest promotion", "Term/synonym promoted from offline harvest.", []),
        ("atom-map required match", "Assignment required mapped atoms to fire.", []),
        ("site aromaticity gate", "Assignment gated on site aromaticity.", []),
        ("multi-site split", "Emission split across distinct localized sites.", []),
    ]:
        add(pref, "xrm:1900000", definition, syn or None, "prov")

    # Localization templates (need ~15–30; not per-atom enumerations)
    for pref, definition, syn in [
        ("atom site template", "Template localizing a reaction to an atom site.", ["atom_site"]),
        ("bond site template", "Template localizing a reaction to a bond site.", ["bond_site"]),
        ("ring site template", "Template localizing a reaction to a ring site.", ["ring_site"]),
        ("attachment atom template O", "Template for O-attachment localization.", []),
        ("attachment atom template N", "Template for N-attachment localization.", []),
        ("attachment atom template S", "Template for S-attachment localization.", []),
        ("attachment atom template C", "Template for C-attachment localization.", []),
        ("attachment atom template acyl", "Template for acyl-attachment localization.", []),
        ("canonical site name", "Canonical generated site-localized display name.", ["canonical_site_name"]),
        ("map-number site label", "Site label using atom-map numbers (e.g. @1).", []),
        ("element-index site label", "Site label using element + index (e.g. C7).", []),
        ("motif-localized name template", "Motif-aware localized name (e.g. para on anisole).", []),
        ("product-class name template", "Product-class display name (e.g. phenol formation).", []),
        ("medchem tag name template", "Med-chem liability display phrase.", []),
        ("regio descriptor template", "Ortho/meta/para (or similar) regio descriptor.", []),
        ("stereo descriptor template", "Stereo descriptor in a localized name.", []),
        ("multi-site list template", "Template listing multiple localized sites.", []),
        ("bond endpoints template", "Template naming both ends of a bond site.", []),
        ("ring index template", "Template naming a ring by index/system.", []),
        ("heteroatom symbol template", "Template inserting N/O/S/halogen symbol.", []),
        ("chain position template", "Template for ω / ω-1 / terminal position names.", []),
        ("benzylic position template", "Template for benzylic position naming.", []),
        ("allylic position template", "Template for allylic position naming.", []),
        ("aniline nitrogen template", "Template for aniline N localization.", []),
        ("phenol oxygen template", "Template for phenolic O localization.", []),
        ("tertiary amine nitrogen template", "Template for tertiary amine N localization.", []),
        ("site environment phrase template", "Phrase combining site environment + reaction.", []),
        ("underspecified site template", "Display when site localization is ambiguous.", []),
    ]:
        add(pref, "xrm:1600000", definition, syn or None, "loc")

    # Cross-spine related_match for pathway logic examples
    by_label = {c["preferred_label"]: c for c in doc["concepts"]}

    def relate(a, b):
        ca, cb = by_label.get(a), by_label.get(b)
        if not ca or not cb:
            return
        for src, tgt in ((ca, cb["id"]), (cb, ca["id"])):
            rel = list(src.get("related_match") or [])
            if tgt not in rel:
                rel.append(tgt)
            src["related_match"] = rel

    relate("aromatic hydroxylation", "quinone-imine formation")
    relate("aromatic hydroxylation", "quinone formation")
    relate("N-dealkylation", "aldehyde metabolite")
    relate("N-dealkylation", "aldehyde forming")
    relate("epoxidation", "epoxide GSH conjugation")
    relate("acyl glucuronidation", "acyl glucuronide risk")
    relate("phenol", "aromatic soft spot")
    relate("tertiary amine", "amine soft spot")

    # Slug synonyms for new leaves
    for c in doc["concepts"]:
        if c["id"] not in ids:
            continue
        # only freshly style: ensure slug for leaves under target spines
        parents = c.get("parents") or []
        if not any(
            p in ("xrm:1500000", "xrm:1600000", "xrm:1900000", "xrm:2100000") for p in parents
        ):
            continue
        slug = re.sub(r"[^a-z0-9]+", "_", c["preferred_label"].lower()).strip("_")
        syns = list(c.get("synonyms") or [])
        if slug and slug not in {s.lower() for s in syns} and slug != c["preferred_label"].lower():
            syns.append(slug)
            c["synonyms"] = syns

    doc["concepts"].sort(key=lambda c: c["id"])
    YAML_PATH.write_text(
        "# XRM thesaurus (authoring source). Export to xrm.skos.jsonld for Rust.\n"
        "# Site-localized names (e.g. C4 hydroxylation) are NOT concepts; generate from templates.\n"
        + yaml.safe_dump(doc, sort_keys=False, allow_unicode=True, width=100)
    )
    print(f"concepts={len(doc['concepts'])}")


if __name__ == "__main__":
    main()
