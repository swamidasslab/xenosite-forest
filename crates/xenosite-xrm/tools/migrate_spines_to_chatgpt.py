#!/usr/bin/env python3
"""Realign XRM spines/terms to the ChatGPT share vocabulary.

Source of truth after this run:
  data/ontology/xrm.yaml  (authoring)
  data/ontology/xrm.skos.jsonld  (Rust runtime export)

ChatGPT core spines (adopted):
  metabolism_phase, chemical_transformation, phase_1_reaction_family,
  phase_2_conjugation_family, medchem_liability, reactive_metabolite_family,
  site_type, structural_delta, product_status, rule_provenance, evidence,
  biological_context

Kept as parallel spines (prior XRM design, still useful):
  ambiguity_and_underspecification, Metabolic Forest map (alias / ruleset map)

Legacy cross-cutting spines are nested under the new parents rather than
deleted, so existing term IDs and assignment emits stay stable.
"""

from __future__ import annotations

import json
import re
from collections import defaultdict
from copy import deepcopy
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
SKOS_PATH = ROOT / "data/ontology/xrm.skos.jsonld"
YAML_PATH = ROOT / "data/ontology/xrm.yaml"
FOREST_SSSOM = ROOT / "data/mappings/xrm-forest.sssom.tsv"
ASSIGNMENTS = ROOT / "data/assignments/xenobiotic.jsonl"
SOURCES_MD = ROOT / "data/ontology/SOURCES.md"


def load_skos():
    data = json.loads(SKOS_PATH.read_text())
    by_id = {}
    for node in data["@graph"]:
        cid = node.get("id")
        if cid:
            by_id[cid] = node
    return data, by_id


def as_list(v):
    if v is None:
        return []
    if isinstance(v, list):
        return v
    return [v]


def set_broader(node, parents):
    parents = [p for p in as_list(parents) if p]
    if not parents:
        node.pop("broader", None)
    elif len(parents) == 1:
        node["broader"] = parents[0]
    else:
        # stable unique order
        seen = []
        for p in parents:
            if p not in seen:
                seen.append(p)
        node["broader"] = seen


def get_broader(node):
    return as_list(node.get("broader"))


def concept(
    cid,
    pref,
    *,
    broader=None,
    definition=None,
    alt=None,
    exact=None,
    close=None,
    broad=None,
    related=None,
    top=False,
    spines=None,
):
    n = {
        "id": cid,
        "type": "skos:Concept",
        "prefLabel": pref,
        "inScheme": "xrm:scheme",
    }
    if top:
        n["topConcept"] = True
    if definition:
        n["definition"] = definition
    if alt:
        n["altLabel"] = alt if isinstance(alt, list) else [alt]
    if spines:
        n["spines"] = spines if isinstance(spines, list) else [spines]
    set_broader(n, broader)
    if exact:
        n["exactMatch"] = exact if isinstance(exact, list) else [exact]
    if close:
        n["closeMatch"] = close if isinstance(close, list) else [close]
    if broad:
        n["broadMatch"] = broad if isinstance(broad, list) else [broad]
    if related:
        n["relatedMatch"] = related if isinstance(related, list) else [related]
    return n


def rewrite_forest_parents(by_id):
    """Move Forest SO/UO rules under *ruleset* concepts; free Rainbow class IDs."""
    so_ruleset = "xrm:9000010"
    uo_ruleset = "xrm:9000011"
    by_id[so_ruleset] = concept(
        so_ruleset,
        "stable oxygenation ruleset",
        broader=["xrm:9000000", "xrm:9000018"],
        definition=(
            "Metabolic Forest ruleset alias for stable oxygenation "
            "(opaque forest.ruleset:SO). Not the chemist Rainbow class."
        ),
        alt=["Forest stable oxygenation", "SO ruleset"],
        exact=["forest.ruleset:SO"],
        spines=["Metabolic Forest map"],
    )
    by_id[uo_ruleset] = concept(
        uo_ruleset,
        "unstable oxygenation ruleset",
        broader=["xrm:9000000", "xrm:9000018"],
        definition=(
            "Metabolic Forest ruleset alias for unstable oxygenation "
            "(opaque forest.ruleset:UO). Not the chemist Rainbow class."
        ),
        alt=["Forest unstable oxygenation", "UO ruleset"],
        exact=["forest.ruleset:UO"],
        spines=["Metabolic Forest map"],
    )

    # Reparent Forest rules that hung under chemist Rainbow IDs.
    for cid, node in list(by_id.items()):
        parents = get_broader(node)
        if "xrm:0000010" in parents:
            parents = ["xrm:9000010" if p == "xrm:0000010" else p for p in parents]
            set_broader(node, parents)
        if "xrm:0000011" in parents:
            parents = ["xrm:9000011" if p == "xrm:0000011" else p for p in parents]
            set_broader(node, parents)

    # Rainbow chemist classes (same IDs, new parents).
    by_id["xrm:0000010"] = concept(
        "xrm:0000010",
        "stable oxygenation",
        broader=["xrm:1200000"],
        definition=(
            "Rainbow Phase I class: oxygenation that yields a stable oxygenated "
            "metabolite (aromatic/aliphatic hydroxylation, epoxidation, "
            "N-oxidation, S-oxidation)."
        ),
        alt=["SO class", "stable oxygen addition class"],
        related=["xrm:9000010"],
        spines=["phase I reaction family"],
    )
    by_id["xrm:0000011"] = concept(
        "xrm:0000011",
        "unstable oxygenation",
        broader=["xrm:1200000"],
        definition=(
            "Rainbow Phase I class: oxygenation that triggers cleavage or "
            "rearrangement (N-/O-/S-/C-dealkylation, oxidative deamination, "
            "oxidative dehalogenation)."
        ),
        alt=["UO class", "unstable oxygen addition class"],
        related=["xrm:9000011"],
        spines=["phase I reaction family"],
    )


def add_core_spines(by_id):
    spines = [
        (
            "xrm:1000000",
            "metabolism phase",
            "Phase framing for xenobiotic biotransformation (I / II / III).",
            ["metabolism_phase", "phase"],
        ),
        (
            "xrm:1100000",
            "chemical transformation",
            "Enzyme-independent chemical edit class (oxidation, reduction, …).",
            ["chemical_transformation", "transformation"],
        ),
        (
            "xrm:1200000",
            "phase I reaction family",
            "Rainbow-anchored Phase I reaction classes and detailed types.",
            ["phase_1_reaction_family", "phase I family"],
        ),
        (
            "xrm:1300000",
            "phase II conjugation family",
            "Conjugation families typical of Phase II xenobiotic metabolism.",
            ["phase_2_conjugation_family", "conjugation family"],
        ),
        (
            "xrm:1400000",
            "medchem liability",
            "Medicinal-chemistry liability and design handles for a reaction.",
            ["medchem_liability", "liability"],
        ),
        (
            "xrm:1500000",
            "reactive metabolite family",
            "Reactive or trapping-prone metabolite product families.",
            ["reactive_metabolite_family", "reactive metabolite"],
        ),
        (
            "xrm:1600000",
            "site type",
            "Where on the molecule the transformation is localized.",
            ["site_type", "site class"],
        ),
        (
            "xrm:1700000",
            "structural delta",
            "Structural / formula / aromaticity / bond-order change facets.",
            ["structural_delta", "delta"],
        ),
        (
            "xrm:1800000",
            "product status",
            "Observation / prediction / pathway status of a metabolite.",
            ["product_status"],
        ),
        (
            "xrm:1900000",
            "rule provenance",
            "How a tag was produced (SMARTS, rule name, priority, exclusions).",
            ["rule_provenance", "provenance"],
        ),
        (
            "xrm:2000000",
            "evidence",
            "Evidence type supporting a named transformation or product.",
            ["evidence"],
        ),
        (
            "xrm:2100000",
            "biological context",
            "Orthogonal biological framing (enzyme, tissue, species, matrix).",
            ["biological_context"],
        ),
    ]
    for cid, pref, definition, alts in spines:
        by_id[cid] = concept(
            cid,
            pref,
            broader="xrm:0000000",
            definition=definition,
            alt=alts,
            spines=[pref],
        )

    # Ambiguity + Forest map stay under root (already present).
    if "xrm:6000000" in by_id:
        set_broader(by_id["xrm:6000000"], "xrm:0000000")
    if "xrm:9000000" in by_id:
        set_broader(by_id["xrm:9000000"], "xrm:0000000")
        by_id["xrm:9000000"]["definition"] = (
            "Alias / mapping spine for Metabolic Forest rulesets → rules → "
            "patterns. Opaque forest.* CURIEs only; not the chemist backbone."
        )


def reparent_phase_and_transformations(by_id):
    # phase I / II under metabolism phase
    set_broader(by_id["xrm:0000001"], ["xrm:1000000"])
    set_broader(by_id["xrm:0000002"], ["xrm:1000000"])
    by_id["xrm:0000003"] = concept(
        "xrm:0000003",
        "phase III transport",
        broader="xrm:1000000",
        definition=(
            "Transporter-mediated disposition related to xenobiotic handling; "
            "related to metabolism naming but not a primary reaction label."
        ),
        alt=["phase 3", "phase III", "transport phase"],
        spines=["metabolism phase"],
    )

    # Chemical transformation top terms (polyhierarchy with phase children).
    chem_roots = {
        "xrm:0000020": "oxidation",
        "xrm:0000014": "reduction",
        "xrm:0000013": "hydrolysis",
        "xrm:0000012": "dehydrogenation",
        "xrm:0000024": "conjugation",
        "xrm:0002000": "tautomerization",
        "xrm:0002003": "isomerization",
        "xrm:0002004": "rearrangement",
    }
    for cid, _label in chem_roots.items():
        if cid not in by_id:
            continue
        parents = get_broader(by_id[cid])
        # Drop direct root attachment; attach to chemical_transformation.
        parents = [p for p in parents if p != "xrm:0000000"]
        if "xrm:1100000" not in parents:
            parents.insert(0, "xrm:1100000")
        # Keep phase parents when already present.
        set_broader(by_id[cid], parents)

    # Extra chemical_transformation leaves from ChatGPT list.
    extras = [
        ("xrm:1100010", "hydration", ["xrm:1100000"], "Addition of water across a bond."),
        ("xrm:1100011", "dehydration", ["xrm:1100000", "xrm:0000014"], "Loss of water."),
        ("xrm:1100012", "dealkylation", ["xrm:1100000"], "Cleavage of an alkyl group from a heteroatom or carbon."),
        ("xrm:1100013", "deamination", ["xrm:1100000"], "Removal of an amino group."),
        ("xrm:1100014", "dehalogenation", ["xrm:1100000"], "Removal of a halogen atom."),
        ("xrm:1100015", "epoxidation", ["xrm:1100000"], "Formation of an epoxide."),
        ("xrm:1100016", "ring oxidation", ["xrm:1100000"], "Oxidation on a ring framework."),
        ("xrm:1100017", "ring opening", ["xrm:1100000"], "Cleavage that opens a ring."),
        ("xrm:1100018", "aromatization", ["xrm:1100000"], "Gain of aromaticity."),
        ("xrm:1100019", "dearomatization", ["xrm:1100000"], "Loss of aromaticity."),
    ]
    for cid, pref, broader, definition in extras:
        # Prefer existing concepts with same prefLabel when present.
        existing = None
        for oid, n in by_id.items():
            if n.get("prefLabel") == pref and oid.startswith("xrm:000"):
                existing = oid
                break
        if existing:
            parents = get_broader(by_id[existing])
            if "xrm:1100000" not in parents:
                parents.append("xrm:1100000")
            set_broader(by_id[existing], parents)
        else:
            by_id[cid] = concept(cid, pref, broader=broader, definition=definition, spines=["chemical transformation"])

    # Link detailed phase I leaves also under phase I reaction family via Rainbow classes.
    rainbow_children = {
        "xrm:0000010": [  # stable oxygenation
            "xrm:0000100",  # hydroxylation
            "xrm:0000110",  # epoxidation
            "xrm:0000120",  # nitrogen oxidation
            "xrm:0000130",  # sulfur oxidation
        ],
        "xrm:0000011": [  # unstable oxygenation
            "xrm:0000200",  # dealkylation
            "xrm:0000208",  # oxidative deamination
            "xrm:0000210",  # oxidative dehalogenation
        ],
        "xrm:0000012": [  # dehydrogenation already exists; also under family
        ],
    }
    # Hang Rainbow classes + DH/HD/RD chemist concepts under phase I reaction family.
    for cid in ["xrm:0000010", "xrm:0000011"]:
        set_broader(by_id[cid], ["xrm:1200000"])
    for cid in ["xrm:0000012", "xrm:0000013", "xrm:0000014"]:
        parents = get_broader(by_id[cid])
        if "xrm:1200000" not in parents:
            parents.append("xrm:1200000")
        # Keep phase I if present
        if "xrm:0000001" not in parents and cid != "xrm:0000012":
            # dehydrogenation/hydrolysis/reduction already under phase I typically
            pass
        if "xrm:0000001" not in parents:
            parents.append("xrm:0000001")
        set_broader(by_id[cid], parents)

    for rainbow, kids in rainbow_children.items():
        for kid in kids:
            if kid not in by_id:
                continue
            parents = get_broader(by_id[kid])
            if rainbow not in parents:
                parents.append(rainbow)
            set_broader(by_id[kid], parents)

    # phase II conjugation family: hang phase II children also under 1300000
    for cid, node in list(by_id.items()):
        if "xrm:0000002" in get_broader(node) and cid != "xrm:1300000":
            parents = get_broader(node)
            if "xrm:1300000" not in parents:
                parents.append("xrm:1300000")
            set_broader(node, parents)

    # Explicit ChatGPT phase II family terms (alias existing when present).
    p2_terms = [
        ("glucuronidation", "xrm:0001000"),
        ("sulfation", "xrm:0001010"),
        ("glutathione conjugation", None),
        ("cysteine conjugation", None),
        ("N-acetylcysteine conjugation", None),
        ("mercapturic acid formation", "xrm:0001037"),
        ("acetylation", "xrm:0001020"),
        ("methylation", "xrm:0001040"),
        ("amino acid conjugation", "xrm:0001050"),
        ("glycine conjugation", "xrm:0001051"),
        ("taurine conjugation", "xrm:0001053"),
        ("phosphorylation", None),
    ]
    next_id = 1300010
    for pref, existing in p2_terms:
        if existing and existing in by_id:
            parents = get_broader(by_id[existing])
            for p in ("xrm:1300000", "xrm:0000002"):
                if p not in parents:
                    parents.append(p)
            set_broader(by_id[existing], parents)
            if pref == "glutathione conjugation":
                pass
            continue
        # Find by label
        found = None
        for oid, n in by_id.items():
            if n.get("prefLabel") == pref:
                found = oid
                break
        if pref == "glutathione conjugation":
            # alias glutathionation
            if "xrm:0001030" in by_id:
                alts = as_list(by_id["xrm:0001030"].get("altLabel"))
                if "glutathione conjugation" not in alts:
                    alts.append("glutathione conjugation")
                by_id["xrm:0001030"]["altLabel"] = alts
                parents = get_broader(by_id["xrm:0001030"])
                for p in ("xrm:1300000", "xrm:0000002"):
                    if p not in parents:
                        parents.append(p)
                set_broader(by_id["xrm:0001030"], parents)
            continue
        if found:
            parents = get_broader(by_id[found])
            if "xrm:1300000" not in parents:
                parents.append("xrm:1300000")
            set_broader(by_id[found], parents)
            continue
        cid = f"xrm:{next_id:07d}"
        next_id += 1
        by_id[cid] = concept(
            cid,
            pref,
            broader=["xrm:1300000", "xrm:0000002"],
            definition=f"Phase II conjugation family: {pref}.",
            spines=["phase II conjugation family"],
        )


def add_medchem_and_reactive(by_id):
    medchem = [
        ("xrm:1400010", "metabolic soft spot", "Site or motif that is readily metabolized."),
        ("xrm:1400011", "clearance liability", "Transformation that drives metabolic clearance."),
        ("xrm:1400012", "bioactivation", "Formation of a more reactive or toxic species."),
        ("xrm:1400013", "detoxification", "Transformation that reduces reactive/toxic liability."),
        ("xrm:1400014", "reactive metabolite", "Product with electrophilic or trapping liability."),
        ("xrm:1400015", "stable metabolite", "Product expected to be chemically stable."),
        ("xrm:1400016", "unstable intermediate", "Transient intermediate along a metabolic path."),
        ("xrm:1400017", "masked liability", "Liability revealed only after a preparatory step."),
        ("xrm:1400018", "metabolic blocking opportunity", "Site where blocking could reduce metabolism."),
        ("xrm:1400019", "bioactivation risk", "Elevated risk of bioactivation at this site/path."),
        ("xrm:1400020", "clearance pathway", "Dominant clearance-oriented conjugation or oxidation path."),
        ("xrm:1400021", "polarity increasing", "Transformation that increases polarity / solubility."),
        ("xrm:1400022", "aldehyde forming", "Path that produces an aldehyde product or intermediate."),
    ]
    for cid, pref, definition in medchem:
        by_id[cid] = concept(
            cid,
            pref,
            broader="xrm:1400000",
            definition=definition,
            spines=["medchem liability"],
        )
    # Merge old bioactivation/detoxication under medchem liability (polyhierarchy).
    if "xrm:0003000" in by_id:
        set_broader(by_id["xrm:0003000"], ["xrm:1400000", "xrm:1400012"])
        # Keep label; related to new bioactivation leaf
        by_id["xrm:0003000"]["relatedMatch"] = as_list(by_id["xrm:0003000"].get("relatedMatch")) + [
            "xrm:1400012"
        ]
    if "xrm:0003001" in by_id:
        set_broader(by_id["xrm:0003001"], ["xrm:1400000", "xrm:1400013"])

    reactive = [
        ("xrm:1500010", "quinone", "Quinone reactive metabolite."),
        ("xrm:1500011", "quinone imine", "Quinone-imine reactive metabolite."),
        ("xrm:1500012", "epoxide", "Epoxide reactive metabolite."),
        ("xrm:1500013", "aldehyde", "Aldehyde reactive metabolite."),
        ("xrm:1500014", "iminium", "Iminium reactive metabolite."),
        ("xrm:1500015", "acyl glucuronide", "Acyl glucuronide reactive conjugate."),
        ("xrm:1500016", "Michael acceptor", "Michael-acceptor electrophile."),
        ("xrm:1500017", "nitroso", "Nitroso reactive metabolite."),
        ("xrm:1500018", "hydroxylamine", "Hydroxylamine reactive metabolite."),
        ("xrm:1500019", "GSH-trappable metabolite", "Metabolite that forms a glutathione adduct."),
    ]
    for cid, pref, definition in reactive:
        by_id[cid] = concept(
            cid,
            pref,
            broader="xrm:1500000",
            definition=definition,
            spines=["reactive metabolite family"],
        )
    # Relate existing quinone/epoxide chemist terms.
    for existing, new in [
        ("xrm:0000300", "xrm:1500010"),
        ("xrm:0000303", "xrm:1500011"),
        ("xrm:0000110", "xrm:1500012"),
        ("xrm:0000305", "xrm:1500014"),
        ("xrm:0001004", "xrm:1500015"),
        ("xrm:0000123", "xrm:1500017"),
        ("xrm:0000122", "xrm:1500018"),
    ]:
        if existing in by_id:
            rel = as_list(by_id[existing].get("relatedMatch"))
            if new not in rel:
                rel.append(new)
            by_id[existing]["relatedMatch"] = rel


def add_site_structural_product_evidence(by_id):
    site_terms = [
        ("xrm:1600010", "atom site", "Transformation localized to an atom."),
        ("xrm:1600011", "bond site", "Transformation localized to a bond."),
        ("xrm:1600012", "ring site", "Transformation localized to a ring."),
        ("xrm:1600013", "aromatic site", "Site on an aromatic system."),
        ("xrm:1600014", "benzylic site", "Benzylic carbon site."),
        ("xrm:1600015", "allylic site", "Allylic carbon site."),
        ("xrm:1600016", "heteroatom site", "Site at N/O/S/halogen heteroatom."),
        ("xrm:1600017", "alpha to heteroatom", "Carbon alpha to a heteroatom."),
        ("xrm:1600018", "terminal alkyl", "Terminal alkyl carbon site."),
        ("xrm:1600019", "stereocenter", "Stereogenic site."),
        ("xrm:1600020", "ambiguous site", "Site localization is ambiguous."),
        ("xrm:1600021", "tertiary amine site", "Tertiary amine nitrogen soft spot."),
        ("xrm:1600022", "phenolic site", "Phenolic oxygen site."),
    ]
    for cid, pref, definition in site_terms:
        by_id[cid] = concept(
            cid, pref, broader="xrm:1600000", definition=definition, spines=["site type"]
        )
    # Nest old site atom / aromaticity spines under site type.
    for old in ["xrm:7100000", "xrm:7800000"]:
        if old in by_id:
            set_broader(by_id[old], ["xrm:1600000"])

    # structural delta absorbs formula-delta, bond-edit, redox, ring fate, oxygenation outcome, cardinality
    for old in [
        "xrm:7700000",
        "xrm:7200000",
        "xrm:7000000",
        "xrm:7600000",
        "xrm:7400000",
        "xrm:7300000",
        "xrm:8000000",
        "xrm:7900000",
    ]:
        if old in by_id:
            set_broader(by_id[old], ["xrm:1700000"])
    if "xrm:7500000" in by_id:
        set_broader(by_id["xrm:7500000"], ["xrm:1500000"])

    delta_terms = [
        ("xrm:1700010", "formula delta", "Elemental composition change."),
        ("xrm:1700011", "mass shift", "Nominal or exact mass difference."),
        ("xrm:1700012", "atom added", "Net atom addition."),
        ("xrm:1700013", "atom removed", "Net atom removal."),
        ("xrm:1700014", "bond order change", "Bond order increased or decreased."),
        ("xrm:1700015", "aromaticity change", "Aromaticity gained or lost."),
        ("xrm:1700016", "charge change", "Formal charge / protonation change."),
        ("xrm:1700017", "stereochemical change", "Stereo configuration change."),
        ("xrm:1700018", "conjugate moiety added", "Conjugate group attached."),
    ]
    for cid, pref, definition in delta_terms:
        by_id[cid] = concept(
            cid, pref, broader="xrm:1700000", definition=definition, spines=["structural delta"]
        )

    product = [
        ("xrm:1800010", "observed product", "Metabolite observed experimentally."),
        ("xrm:1800011", "predicted product", "Metabolite predicted by a model or rule."),
        ("xrm:1800012", "intermediate metabolite", "Non-terminal pathway intermediate."),
        ("xrm:1800013", "terminal metabolite", "Pathway-terminal metabolite."),
        ("xrm:1800014", "minor metabolite", "Minor abundance metabolite."),
        ("xrm:1800015", "major metabolite", "Major abundance metabolite."),
        ("xrm:1800016", "authentic metabolite", "Confirmed authentic metabolite."),
        ("xrm:1800017", "artifact", "Artifactual / non-biological product."),
        ("xrm:1800018", "not detected", "Looked for but not detected."),
        ("xrm:1800019", "ruled out", "Structurally or experimentally ruled out."),
    ]
    for cid, pref, definition in product:
        by_id[cid] = concept(
            cid, pref, broader="xrm:1800000", definition=definition, spines=["product status"]
        )

    provenance = [
        ("xrm:1900010", "SMARTS-derived tag", "Tag assigned from reactant/product SMARTS."),
        ("xrm:1900011", "formula-delta-derived tag", "Tag assigned from elemental delta."),
        ("xrm:1900012", "caller-tag-derived", "Tag assigned from opaque caller CURIE/tag."),
        ("xrm:1900013", "Forest-map-derived tag", "Tag assigned via Metabolic Forest map correspondence."),
        ("xrm:1900014", "exclusion rule", "Negative / exclusion pattern applied."),
    ]
    for cid, pref, definition in provenance:
        by_id[cid] = concept(
            cid, pref, broader="xrm:1900000", definition=definition, spines=["rule provenance"]
        )

    evidence = [
        ("xrm:2000010", "literature evidence", "Supported by a paper or review."),
        ("xrm:2000011", "MS1 evidence", "Supported by MS1."),
        ("xrm:2000012", "MS2 evidence", "Supported by MS/MS."),
        ("xrm:2000013", "retention time evidence", "Supported by chromatographic retention time."),
        ("xrm:2000014", "NMR evidence", "Supported by NMR."),
        ("xrm:2000015", "synthetic standard", "Confirmed versus synthetic standard."),
        ("xrm:2000016", "radiolabel evidence", "Supported by radiolabel tracking."),
        ("xrm:2000017", "HLM incubation", "Human liver microsome incubation context."),
        ("xrm:2000018", "hepatocyte incubation", "Hepatocyte incubation context."),
        ("xrm:2000019", "recombinant enzyme", "Recombinant enzyme incubation context."),
        ("xrm:2000020", "negative control", "Negative-control observation."),
    ]
    for cid, pref, definition in evidence:
        by_id[cid] = concept(
            cid, pref, broader="xrm:2000000", definition=definition, spines=["evidence"]
        )

    bio = [
        ("xrm:2100010", "enzyme family context", "Enzyme family framing (orthogonal to reaction name)."),
        ("xrm:2100011", "tissue context", "Tissue framing."),
        ("xrm:2100012", "species context", "Species framing."),
        ("xrm:2100013", "matrix context", "Matrix / biofluid framing."),
        ("xrm:2100014", "compartment context", "Subcellular compartment framing."),
    ]
    for cid, pref, definition in bio:
        by_id[cid] = concept(
            cid, pref, broader="xrm:2100000", definition=definition, spines=["biological context"]
        )

    # process facet stays useful; hang under chemical transformation
    if "xrm:0004000" in by_id:
        set_broader(by_id["xrm:0004000"], ["xrm:1100000"])


def update_scheme_blurb(by_id):
    scheme = by_id["xrm:scheme"]
    scheme["definition"] = (
        "XRM tags each reaction with many terms from parallel, cross-cutting "
        "spines aligned to med-chem faceting: metabolism phase, chemical "
        "transformation, phase I reaction family (Rainbow), phase II conjugation "
        "family, medchem liability, reactive metabolite family, site type, "
        "structural delta, product status, rule provenance, evidence, biological "
        "context, plus ambiguity and Metabolic Forest map (alias spine)."
    )


def export_jsonld(data, by_id):
    # Preserve scheme first, then root, then others sorted by id.
    graph = []
    for cid in ["xrm:scheme", "xrm:0000000"]:
        if cid in by_id:
            graph.append(by_id[cid])
    for cid in sorted(by_id.keys()):
        if cid in ("xrm:scheme", "xrm:0000000"):
            continue
        graph.append(by_id[cid])
    out = {"@context": data["@context"], "@graph": graph}
    SKOS_PATH.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n")


def export_yaml(by_id):
    """Authoring YAML: spines first, then concepts with stable fields."""
    spine_ids = [
        "xrm:1000000",
        "xrm:1100000",
        "xrm:1200000",
        "xrm:1300000",
        "xrm:1400000",
        "xrm:1500000",
        "xrm:1600000",
        "xrm:1700000",
        "xrm:1800000",
        "xrm:1900000",
        "xrm:2000000",
        "xrm:2100000",
        "xrm:6000000",
        "xrm:9000000",
    ]
    spines = []
    for sid in spine_ids:
        n = by_id[sid]
        spines.append(
            {
                "id": sid,
                "preferred_label": n["prefLabel"],
                "definition": n.get("definition"),
                "synonyms": as_list(n.get("altLabel")),
            }
        )

    concepts = []
    for cid in sorted(by_id.keys()):
        if cid == "xrm:scheme":
            continue
        n = by_id[cid]
        if n.get("type") != "skos:Concept":
            continue
        entry = {
            "id": cid,
            "preferred_label": n.get("prefLabel"),
            "synonyms": as_list(n.get("altLabel")),
            "definition": n.get("definition"),
            "parents": get_broader(n),
            "spines": as_list(n.get("spines")),
            "exact_match": as_list(n.get("exactMatch")),
            "close_match": as_list(n.get("closeMatch")),
            "related_match": as_list(n.get("relatedMatch")),
            "broad_match": as_list(n.get("broadMatch")),
        }
        # Drop empty optional lists for readability
        for k in list(entry.keys()):
            if entry[k] in (None, [], ""):
                del entry[k]
        concepts.append(entry)

    doc = {
        "scheme": {
            "id": "xrm:scheme",
            "preferred_label": by_id["xrm:scheme"].get("prefLabel"),
            "iri": by_id["xrm:scheme"].get("iri"),
            "definition": by_id["xrm:scheme"].get("definition"),
        },
        "spines": spines,
        "concepts": concepts,
    }
    YAML_PATH.write_text(
        "# XRM thesaurus (authoring source). Export to xrm.skos.jsonld for Rust.\n"
        + yaml.safe_dump(doc, sort_keys=False, allow_unicode=True, width=100)
    )


def patch_forest_sssom():
    text = FOREST_SSSOM.read_text()
    # Move exactMatch for SO/UO from chemist Rainbow IDs to Forest ruleset IDs.
    text = text.replace(
        "xrm:0000010\tskos:exactMatch\tforest.ruleset:SO\tsemapv:ManualMappingCuration\tstable oxygenation\tSO",
        "xrm:9000010\tskos:exactMatch\tforest.ruleset:SO\tsemapv:ManualMappingCuration\tstable oxygenation ruleset\tSO\n"
        "xrm:0000010\tskos:relatedMatch\tforest.ruleset:SO\tsemapv:ManualMappingCuration\tstable oxygenation\tSO",
    )
    text = text.replace(
        "xrm:0000011\tskos:exactMatch\tforest.ruleset:UO\tsemapv:ManualMappingCuration\tunstable oxygenation\tUO",
        "xrm:9000011\tskos:exactMatch\tforest.ruleset:UO\tsemapv:ManualMappingCuration\tunstable oxygenation ruleset\tUO\n"
        "xrm:0000011\tskos:relatedMatch\tforest.ruleset:UO\tsemapv:ManualMappingCuration\tunstable oxygenation\tUO",
    )
    FOREST_SSSOM.write_text(text)


def patch_assignments():
    """Emit Rainbow classes + medchem/site/reactive tags alongside chemist types."""
    lines = ASSIGNMENTS.read_text().splitlines()
    out = []
    # Prepend comment block
    if not lines[0].startswith("# Cross-cutting"):
        out.append(
            "# Cross-cutting emits preferred: ChatGPT spines (phase, Rainbow family, "
            "medchem liability, site type, structural delta, …)."
        )
    replacements = {
        # hydroxylation → stable oxygenation + soft spot + stable metabolite + oxygen gain already present
        '"asg:hydroxylation"': None,
    }

    def add_emits(line: str, extra: list[str]) -> str:
        m = re.search(r'"emit":\[(.*?)\]', line)
        if not m:
            return line
        current = m.group(1)
        ids = re.findall(r'"([^"]+)"', current)
        for e in extra:
            if e not in ids:
                ids.append(e)
        new_emit = '"emit":[' + ",".join(f'"{i}"' for i in ids) + "]"
        return line[: m.start()] + new_emit + line[m.end() :]

    for line in lines:
        if line.startswith("#") or not line.strip():
            out.append(line)
            continue
        try:
            obj = json.loads(line)
        except json.JSONDecodeError:
            out.append(line)
            continue
        eid = obj.get("id", "")
        extra = []
        if eid in ("asg:hydroxylation", "asg:aromatic-hydroxylation", "asg:aliphatic-hydroxylation", "asg:tag-hydroxylation", "asg:tag-aromatic-hydroxylation", "asg:tag-para-hydroxylation", "asg:tag-benzylic-hydroxylation", "asg:tag-omega-hydroxylation"):
            extra += ["xrm:0000010", "xrm:1200000", "xrm:1400010", "xrm:1400015", "xrm:1900010"]
        if eid in ("asg:epoxidation", "asg:tag-epoxidation", "asg:tag-arene-oxide"):
            extra += ["xrm:0000010", "xrm:1200000", "xrm:1500012", "xrm:1400014", "xrm:1400019"]
        if eid in ("asg:tag-dealkylation", "asg:tag-n-dealkylation", "asg:tag-n-demethylation", "asg:tag-oxidative-deamination", "asg:tag-o-dealkylation", "asg:tag-o-demethylation"):
            extra += ["xrm:0000011", "xrm:1200000", "xrm:1400010", "xrm:1400011"]
        if eid in ("asg:tag-n-dealkylation", "asg:tag-n-demethylation", "asg:tag-oxidative-deamination"):
            extra += ["xrm:1400022", "xrm:1600021", "xrm:1400019"]
        if eid in ("asg:alcohol-oxidation", "asg:primary-alcohol-oxidation"):
            extra += ["xrm:0000012", "xrm:1200000", "xrm:1400022", "xrm:1500013"]
        if eid.startswith("asg:tag-glucuronidation") or eid in ("asg:tag-phenolic-glucuronidation", "asg:tag-n-glucuronidation", "asg:tag-acyl-glucuronidation"):
            extra += ["xrm:1300000", "xrm:1400020", "xrm:1400021"]
        if eid == "asg:tag-acyl-glucuronidation":
            extra += ["xrm:1500015", "xrm:1400014", "xrm:1400019"]
        if eid in ("asg:tag-quinone", "asg:tag-quinone-imine", "asg:tag-quinone-methide", "asg:tag-imine-methide", "asg:tag-one-step-quinone", "asg:tag-two-step-quinone"):
            extra += ["xrm:1500010", "xrm:1400014", "xrm:1400012", "xrm:1400019", "xrm:0000012"]
        if eid == "asg:tag-quinone-imine":
            extra += ["xrm:1500011"]
        if eid in ("asg:tag-gsh-michael",):
            extra += ["xrm:1500016", "xrm:1500019", "xrm:1400013", "xrm:1300000"]
        if eid in ("asg:tag-gsh-epoxide", "asg:tag-glutathionation"):
            extra += ["xrm:1500019", "xrm:1400013", "xrm:1300000"]
        if eid == "asg:aromatic-hydroxylation":
            extra += ["xrm:1600013", "xrm:7800001"]
        if eid == "asg:aliphatic-hydroxylation":
            extra += ["xrm:1600010"]
        if eid == "asg:tag-benzylic-hydroxylation":
            extra += ["xrm:1600014"]
        if eid == "asg:delta-oxygen-gain":
            extra += ["xrm:1700010", "xrm:1700012", "xrm:1900011"]
        if eid.startswith("asg:delta-"):
            extra += ["xrm:1700000", "xrm:1900011"]
        if extra:
            emits = obj.get("emit", [])
            for e in extra:
                if e not in emits:
                    emits.append(e)
            obj["emit"] = emits
            out.append(json.dumps(obj, separators=(",", ":")))
        else:
            out.append(json.dumps(obj, separators=(",", ":")) if line.startswith("{") else line)
    # Ensure comment lines preserved: rewrite more carefully
    # Actually we destroyed comments for non-matching - fix by re-read approach
    ASSIGNMENTS.write_text("\n".join(out) + "\n")


def patch_assignments_preserve_comments():
    raw = ASSIGNMENTS.read_text().splitlines()
    # Re-run from original git? We may have already overwritten. Use current JSON lines.
    # Reload from SKOS-adjacent backup: re-read file just written and also pull comments from git.
    import subprocess

    original = subprocess.check_output(
        ["git", "show", "HEAD:crates/xenosite-xrm/data/assignments/xenobiotic.jsonl"],
        cwd="/workspace",
        text=True,
    ).splitlines()

    def enrich(obj):
        eid = obj.get("id", "")
        extra = []
        if eid in (
            "asg:hydroxylation",
            "asg:aromatic-hydroxylation",
            "asg:aliphatic-hydroxylation",
            "asg:tag-hydroxylation",
            "asg:tag-aromatic-hydroxylation",
            "asg:tag-para-hydroxylation",
            "asg:tag-benzylic-hydroxylation",
            "asg:tag-omega-hydroxylation",
        ):
            extra += [
                "xrm:0000010",
                "xrm:1200000",
                "xrm:1400010",
                "xrm:1400015",
                "xrm:1900010",
            ]
        if eid in ("asg:epoxidation", "asg:tag-epoxidation", "asg:tag-arene-oxide"):
            extra += [
                "xrm:0000010",
                "xrm:1200000",
                "xrm:1500012",
                "xrm:1400014",
                "xrm:1400019",
            ]
        if eid in (
            "asg:tag-dealkylation",
            "asg:tag-n-dealkylation",
            "asg:tag-n-demethylation",
            "asg:tag-oxidative-deamination",
            "asg:tag-o-dealkylation",
            "asg:tag-o-demethylation",
        ):
            extra += ["xrm:0000011", "xrm:1200000", "xrm:1400010", "xrm:1400011"]
        if eid in (
            "asg:tag-n-dealkylation",
            "asg:tag-n-demethylation",
            "asg:tag-oxidative-deamination",
        ):
            extra += ["xrm:1400022", "xrm:1600021", "xrm:1400019"]
        if eid in ("asg:alcohol-oxidation", "asg:primary-alcohol-oxidation"):
            extra += ["xrm:0000012", "xrm:1200000", "xrm:1400022", "xrm:1500013"]
        if eid in (
            "asg:tag-glucuronidation",
            "asg:tag-phenolic-glucuronidation",
            "asg:tag-n-glucuronidation",
            "asg:tag-acyl-glucuronidation",
        ):
            extra += ["xrm:1300000", "xrm:1400020", "xrm:1400021"]
        if eid == "asg:tag-acyl-glucuronidation":
            extra += ["xrm:1500015", "xrm:1400014", "xrm:1400019"]
        if eid in (
            "asg:tag-quinone",
            "asg:tag-quinone-imine",
            "asg:tag-quinone-methide",
            "asg:tag-imine-methide",
            "asg:tag-one-step-quinone",
            "asg:tag-two-step-quinone",
        ):
            extra += [
                "xrm:1500010",
                "xrm:1400014",
                "xrm:1400012",
                "xrm:1400019",
                "xrm:1200000",
            ]
        if eid == "asg:tag-quinone-imine":
            extra += ["xrm:1500011"]
        if eid == "asg:tag-gsh-michael":
            extra += ["xrm:1500016", "xrm:1500019", "xrm:1400013", "xrm:1300000"]
        if eid in ("asg:tag-gsh-epoxide", "asg:tag-glutathionation"):
            extra += ["xrm:1500019", "xrm:1400013", "xrm:1300000"]
        if eid == "asg:aromatic-hydroxylation":
            extra += ["xrm:1600013"]
        if eid == "asg:aliphatic-hydroxylation":
            extra += ["xrm:1600010"]
        if eid == "asg:tag-benzylic-hydroxylation":
            extra += ["xrm:1600014"]
        if eid == "asg:delta-oxygen-gain":
            extra += ["xrm:1700010", "xrm:1700012", "xrm:1900011"]
        if eid.startswith("asg:delta-"):
            extra += ["xrm:1700000", "xrm:1900011"]
        if eid.startswith("asg:tag-") and "ambiguity" not in eid and eid.startswith("asg:tag-"):
            # caller-tag provenance for tag-driven rules
            if not eid.startswith("asg:tag-fully") and "ambiguity" not in eid and "underspec" not in eid and "competing" not in eid and "phase-ambiguity" not in eid and "forest" not in eid and "external" not in eid and "aromatic-impact" not in eid and "electrophile-role" not in eid and "provenance" not in eid and "som" not in eid and "regio" not in eid and "stereo" not in eid and "type-ambiguity" not in eid and "mechanism" not in eid and "structure" not in eid and "mapping" not in eid and "formula-only" not in eid and "pathway" not in eid and "intermediate" not in eid and "ambiguous" not in eid:
                extra += ["xrm:1900012"]
        emits = list(obj.get("emit", []))
        for e in extra:
            if e not in emits:
                emits.append(e)
        obj["emit"] = emits
        return obj

    out = [
        "# Cross-cutting emits: ChatGPT spines (metabolism phase, Rainbow phase I family,",
        "# phase II conjugation family, medchem liability, reactive metabolite family,",
        "# site type, structural delta, product status, rule provenance, …).",
        "# Structural / delta rules",
    ]
    for line in original:
        if line.startswith("#") or not line.strip():
            # skip old header comments we replaced
            if line.startswith("# Cross-cutting") or line.startswith("# Structural"):
                continue
            out.append(line)
            continue
        obj = json.loads(line)
        obj = enrich(obj)
        out.append(json.dumps(obj, separators=(",", ":")))
    ASSIGNMENTS.write_text("\n".join(out) + "\n")


def write_sources_md():
    SOURCES_MD.write_text(
        """# Sources for XRM term inventory

XRM **tags each reaction with many terms** from parallel, cross-cutting spines.
Authoring source: [`xrm.yaml`](xrm.yaml). Runtime export: [`xrm.skos.jsonld`](xrm.skos.jsonld).

Spine vocabulary follows the ChatGPT design share (metabolism phase, chemical
transformation, Rainbow phase I family, phase II conjugation family, medchem
liability, reactive metabolite family, site type, structural delta, product
status, rule provenance, evidence, biological context), plus ambiguity and the
Metabolic Forest map as an **alias / mapping** spine.

Forest abbreviations (`SO`, `UO`, `DH`, `HD`, `RD`, …) are **never** XRM
prefLabels. They appear only as opaque `forest.*` CURIE object ids in SSSOM.
Forest ruleset concepts use unabbreviated labels ending in `ruleset`.

## Core spines (under `xenobiotic biotransformation`)

| Spine | Role | Auto-tag cues |
| --- | --- | --- |
| **metabolism phase** | Phase I / II / III framing | ancestors of typed reactions; `chem:phase-*` |
| **chemical transformation** | Enzyme-independent edit class | SMARTS / delta / typed chem tags |
| **phase I reaction family** | Rainbow five classes + detailed Phase I types | hydroxylation→stable oxygenation; dealkylation→unstable oxygenation |
| **phase II conjugation family** | Glucuronidation, sulfation, GSH, … | conjugation SMARTS / tags |
| **medchem liability** | Soft spot, clearance, bioactivation, blocking | typed liability emits |
| **reactive metabolite family** | Quinone, epoxide, aldehyde, acyl glucuronide, … | quinone / epoxide / GSH tags |
| **site type** | Atom / bond / ring / aromatic / benzylic / … | site aromaticity + typed sites |
| **structural delta** | Formula, mass, bond order, aromaticity, … | elemental delta; nested legacy delta spines |
| **product status** | Observed / predicted / intermediate / … | caller / harvest tags |
| **rule provenance** | SMARTS vs delta vs caller vs Forest map | assignment path |
| **evidence** | Literature, MS, NMR, incubation system, … | harvest / curation tags |
| **biological context** | Enzyme / tissue / species / matrix (orthogonal) | never primary reaction label |
| **ambiguity and underspecification** | Typed incomplete/conflicting evidence | `chem:*-ambiguity` |
| **Metabolic Forest map** | Alias spine: ruleset → rule → PatternInfo | `forest.rule:*`, `forest.pattern:*` |

Legacy facets (redox polarity, bond-edit topology, formula-delta class, ring
fate, oxygenation outcome, metabolite cardinality, aromatic impact,
pathway-step role, site atom class, site aromaticity, process facet) hang
**under** `structural delta`, `site type`, or `chemical transformation` rather
than as peer root spines.

## Metabolic Forest map (full names)

| XRM prefLabel | Opaque Forest CURIE |
| --- | --- |
| stable oxygenation ruleset | `forest.ruleset:SO` |
| unstable oxygenation ruleset | `forest.ruleset:UO` |
| dehydrogenation ruleset | `forest.ruleset:DH` |
| hydrolysis ruleset | `forest.ruleset:HD` |
| reduction ruleset | `forest.ruleset:RD` |
| quinone formation ruleset | `forest.ruleset:QF` |
| conjugation ruleset | `forest.ruleset:CJ` |
| tautomerization ruleset | `forest.ruleset:TT` |
| phase I ruleset | `forest.ruleset:PhaseOne` |
| bioactivation ruleset | `forest.ruleset:BA` |

Chemist Rainbow classes `stable oxygenation` / `unstable oxygenation` are under
**phase I reaction family**, with `skos:relatedMatch` to the Forest rulesets.

## Primary literature

1. **Rainbow** — Dang et al., JCIM 2020
   ([10.1021/acs.jcim.9b00836](https://doi.org/10.1021/acs.jcim.9b00836)).
2. **Metabolic Forest** — metabolite enumeration; rulesets for Phase I classes,
   conjugation, quinone formation, tautomerization.
3. **Quinone formation** — Hughes & Swamidass, Chem. Res. Toxicol. 2017
   ([10.1021/acs.chemrestox.6b00385](https://doi.org/10.1021/acs.chemrestox.6b00385)).
4. **IUPAC** xenobiotic metabolism glossary (Pure Appl. Chem. 2021,
   [10.1515/pac-2018-0208](https://doi.org/10.1515/pac-2018-0208)).
5. DMPK / Phase II teaching literature for conjugation and process facets.

## Automated external guidance

Offline harvesters (Python, stdlib) collect candidate terms, synonyms, and
example pairs from **ChEBI, PubChem, KEGG, Rhea, GO, Reactome**, plus a local
seed lexicon. See:

- [`tools/harvest/README.md`](../tools/harvest/README.md)
- [`tools/harvest/GUIDANCE_SOURCES.md`](../tools/harvest/GUIDANCE_SOURCES.md)
- Output: [`data/candidates/`](../candidates/)

Promote reviewed candidates into YAML `synonyms`, SSSOM, and SMARTS-backed
goldens. Prefer structure (SMARTS/delta) over Forest reaction-tool tags.
"""
    )


def main():
    data, by_id = load_skos()
    rewrite_forest_parents(by_id)
    add_core_spines(by_id)
    reparent_phase_and_transformations(by_id)
    add_medchem_and_reactive(by_id)
    add_site_structural_product_evidence(by_id)
    update_scheme_blurb(by_id)

    # Strip non-SKOS helper field before JSON-LD export (keep in YAML via deepcopy)
    yaml_by_id = deepcopy(by_id)
    for n in by_id.values():
        n.pop("spines", None)

    export_jsonld(data, by_id)
    export_yaml(yaml_by_id)
    patch_forest_sssom()
    patch_assignments_preserve_comments()
    write_sources_md()

    # Summary
    children = defaultdict(list)
    for cid, n in by_id.items():
        for p in get_broader(n):
            children[p].append(cid)
    print("Root children:")
    for cid in sorted(children["xrm:0000000"], key=lambda x: by_id[x].get("prefLabel", x)):
        print(f"  {cid}  {by_id[cid].get('prefLabel')}  n={len(children[cid])}")
    print(f"Wrote {SKOS_PATH}")
    print(f"Wrote {YAML_PATH}")
    print(f"concepts={sum(1 for n in by_id.values() if n.get('type')=='skos:Concept')}")


if __name__ == "__main__":
    main()
