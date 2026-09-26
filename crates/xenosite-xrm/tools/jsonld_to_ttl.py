#!/usr/bin/env python3
"""Export data/ontology/xrm.skos.jsonld → xrm.skos.ttl (+ optional xmet attachment props).

Also materializes ``xmet:hasAttachmentAtomType`` when a concept ``skos:relatedMatch``
points at an attachment-atom site concept — so SPARQL competency queries can use
the operational property from ANNOTATION.md without inventing SKOS parents.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
JSONLD = ROOT / "data/ontology/xrm.skos.jsonld"
TTL = ROOT / "data/ontology/xrm.skos.ttl"

XRM = "https://xenosite.org/ontology/xrm#"
XMET = "https://xenosite.org/ontology/xmet#"
SKOS = "http://www.w3.org/2004/02/skos/core#"

ATTACHMENT = {
    "xrm:1600135": "oxygen",  # attachment atom O
    "xrm:1600136": "nitrogen",
    "xrm:1600137": "sulfur",
    "xrm:1600138": "carbon",
    "xrm:1600139": "acyl",
}


def curie_to_iri(curie: str) -> str:
    if curie.startswith("xrm:"):
        return XRM + curie.split(":", 1)[1]
    if curie.startswith("http://") or curie.startswith("https://"):
        return curie
    # opaque / external CURIEs kept as string literals in match objects when needed
    return curie


def ttl_iri(curie: str) -> str:
    if curie.startswith("xrm:"):
        return f"xrm:{curie.split(':', 1)[1]}"
    if curie.startswith("http://") or curie.startswith("https://"):
        return f"<{curie}>"
    # non-xrm CURIEs → angle-bracket synthetic IRIs under xrm:ext/
    safe = curie.replace(":", "/")
    return f"<https://xenosite.org/ontology/xrm/ext/{safe}>"


def esc(s: str) -> str:
    return s.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n")


def as_list(v):
    if v is None:
        return []
    if isinstance(v, list):
        return v
    return [v]


def main() -> None:
    data = json.loads(JSONLD.read_text())
    lines = [
        "@prefix skos: <http://www.w3.org/2004/02/skos/core#> .",
        "@prefix xrm: <https://xenosite.org/ontology/xrm#> .",
        "@prefix xmet: <https://xenosite.org/ontology/xmet#> .",
        "@prefix dct: <http://purl.org/dc/terms/> .",
        "",
    ]
    for n in data["@graph"]:
        nid = n["id"]
        subj = ttl_iri(nid)
        typ = n.get("type", "")
        if typ == "skos:ConceptScheme":
            lines.append(f"{subj} a skos:ConceptScheme ;")
            if n.get("prefLabel"):
                lines.append(f'  skos:prefLabel "{esc(n["prefLabel"])}" ;')
            if n.get("definition"):
                lines.append(f'  skos:definition "{esc(n["definition"])}" ;')
            lines[-1] = lines[-1].rstrip(" ;") + " ."
            lines.append("")
            continue
        if typ != "skos:Concept":
            continue
        stmts = ["a skos:Concept"]
        if n.get("prefLabel"):
            stmts.append(f'skos:prefLabel "{esc(n["prefLabel"])}"')
        if n.get("definition"):
            stmts.append(f'skos:definition "{esc(n["definition"])}"')
        if n.get("inScheme"):
            stmts.append(f"skos:inScheme {ttl_iri(n['inScheme'])}")
        for lab in as_list(n.get("altLabel")):
            stmts.append(f'skos:altLabel "{esc(lab)}"')
        for pred in ("broader", "narrower", "exactMatch", "closeMatch", "broadMatch", "relatedMatch"):
            for obj in as_list(n.get(pred)):
                stmts.append(f"skos:{pred} {ttl_iri(obj)}")
                if pred == "relatedMatch" and obj in ATTACHMENT:
                    stmts.append(f"xmet:hasAttachmentAtomType {ttl_iri(obj)}")
        lines.append(f"{subj} {stmts[0]} ;")
        for s in stmts[1:-1]:
            lines.append(f"  {s} ;")
        if len(stmts) > 1:
            lines.append(f"  {stmts[-1]} .")
        else:
            lines[-1] = lines[-1].rstrip(" ;") + " ."
        lines.append("")

    TTL.write_text("\n".join(lines) + "\n")
    print(f"wrote {TTL}")


if __name__ == "__main__":
    main()
