#!/usr/bin/env python3
"""Expand underfilled ChatGPT spines toward ~425 canonical chemist concepts.

Does not add Forest-map rules/patterns. Does not enumerate site-localized names.
"""

from __future__ import annotations

import json
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
YAML_PATH = ROOT / "data/ontology/xrm.yaml"
SKOS_PATH = ROOT / "data/ontology/xrm.skos.jsonld"

# Next free id ranges (avoid collisions with existing).
NEXT = {
    "site": 1600100,
    "reactive": 1500100,
    "medchem": 1400100,
    "evidence": 2000100,
    "bio": 2100100,
    "p2": 1300100,
    "p1": 1200100,
    "product": 1800100,
    "delta": 1700100,
    "phase": 1000100,
    "chem": 1100100,
    "prov": 1900100,
}


def next_id(bucket: str) -> str:
    n = NEXT[bucket]
    NEXT[bucket] = n + 1
    return f"xrm:{n:07d}"


def main():
    doc = yaml.safe_load(YAML_PATH.read_text())
    existing_labels = {c["preferred_label"].lower() for c in doc["concepts"]}
    existing_ids = {c["id"] for c in doc["concepts"]}

    def add(pref, parents, definition, synonyms=None, bucket="chem", related=None):
        if pref.lower() in existing_labels:
            return
        cid = next_id(bucket)
        while cid in existing_ids:
            cid = next_id(bucket)
        entry = {
            "id": cid,
            "preferred_label": pref,
            "definition": definition,
            "parents": parents if isinstance(parents, list) else [parents],
        }
        if synonyms:
            entry["synonyms"] = synonyms
        if related:
            entry["related_match"] = related if isinstance(related, list) else [related]
        doc["concepts"].append(entry)
        existing_labels.add(pref.lower())
        existing_ids.add(cid)

    # --- metabolism phase / pathway role ---
    add(
        "sequential metabolism",
        "xrm:1000000",
        "Multi-step metabolic sequence rather than a single transformation.",
        ["sequential biotransformation"],
        "phase",
    )
    add(
        "intermediate metabolite formation",
        "xrm:1000000",
        "Formation of a non-terminal metabolic intermediate.",
        bucket="phase",
    )
    add(
        "first-pass metabolism",
        "xrm:1000000",
        "Pathway role associated with first-pass metabolic clearance.",
        bucket="phase",
    )

    # --- phase I detailed (Rainbow-compatible finer types) ---
    for pref, definition, syn in [
        ("aromatic epoxidation", "Epoxidation on an aromatic system.", ["arene epoxidation"]),
        ("aliphatic epoxidation", "Epoxidation of an aliphatic alkene.", []),
        ("N-oxidation", "Oxidation at nitrogen.", ["nitrogen oxidation class"]),
        ("S-oxidation", "Oxidation at sulfur.", ["sulfur oxidation class"]),
        ("C-dealkylation", "Oxidative cleavage of a carbon–carbon alkyl attachment.", []),
        ("oxidative deamination class", "Rainbow UO class member: oxidative deamination.", []),
        ("alcohol dehydrogenation", "Dehydrogenation of an alcohol to carbonyl.", []),
        ("amine dehydrogenation", "Dehydrogenation at an amine (iminium-forming path).", []),
        ("ester hydrolysis class", "Hydrolytic cleavage of an ester.", []),
        ("amide hydrolysis class", "Hydrolytic cleavage of an amide.", []),
        ("carbonyl reduction class", "Reduction of a carbonyl to alcohol.", []),
        ("nitro reduction class", "Reduction of a nitro group.", []),
        ("reductive dehalogenation class", "Reductive removal of halogen.", []),
        ("hydrogenation class", "Addition of hydrogen across unsaturation.", []),
        ("aldehyde formation class", "Phase I path that forms an aldehyde.", ["aldehyde-forming class"]),
    ]:
        add(pref, "xrm:1200000", definition, syn or None, "p1")

    # --- phase II conjugation detail ---
    for pref, definition, syn in [
        ("O-glucuronidation class", "Glucuronidation at oxygen.", ["phenolic or alcoholic O-glucuronidation"]),
        ("N-glucuronidation class", "Glucuronidation at nitrogen.", []),
        ("acyl glucuronidation class", "Glucuronidation of a carboxylic acid.", []),
        ("O-sulfation class", "Sulfation at oxygen.", []),
        ("N-sulfation class", "Sulfation at nitrogen.", []),
        ("GSH conjugation class", "Glutathione conjugation family.", ["glutathione conjugation class"]),
        ("cysteine conjugation", "Conjugation with cysteine.", []),
        ("N-acetylcysteine conjugation", "Formation of an N-acetylcysteine (mercapturate precursor) conjugate.", []),
        ("O-methylation class", "Methylation at oxygen.", []),
        ("N-methylation class", "Methylation at nitrogen.", []),
        ("S-methylation class", "Methylation at sulfur.", []),
        ("N-acetylation class", "Acetylation at nitrogen.", []),
        ("O-acetylation class", "Acetylation at oxygen.", []),
        ("glycine conjugation class", "Amino-acid conjugation with glycine.", []),
        ("taurine conjugation class", "Amino-acid conjugation with taurine.", []),
        ("glutamine conjugation class", "Amino-acid conjugation with glutamine.", []),
        ("phosphorylation conjugation", "Phase II-like phosphorylation of a xenobiotic.", ["xenobiotic phosphorylation"]),
        ("carbamoyl glucuronidation class", "Carbamoyl glucuronide formation.", []),
        ("quaternary N-glucuronidation class", "Quaternary ammonium N-glucuronide formation.", []),
        ("epoxide GSH conjugation", "GSH attack on an epoxide.", []),
        ("Michael GSH conjugation", "GSH attack on a Michael acceptor.", []),
        ("halide-displacement GSH conjugation", "GSH displacement of a halide.", []),
    ]:
        parents = ["xrm:1300000", "xrm:0000002"]
        add(pref, parents, definition, syn or None, "p2")

    # --- reactive / product metabolite classes ---
    for pref, definition, syn in [
        ("quinone-like reactive metabolite", "Broad class of quinoid reactive metabolites.", ["quinoid metabolite"]),
        ("ortho-quinone", "Ortho-quinone reactive metabolite.", []),
        ("para-quinone", "Para-quinone reactive metabolite.", []),
        ("quinone methide", "Quinone-methide reactive metabolite.", []),
        ("imine methide", "Imine-methide reactive metabolite.", []),
        ("arene oxide", "Arene oxide reactive intermediate.", []),
        ("thiophene epoxide", "Thiophene S-oxide/epoxide-like reactive intermediate.", []),
        ("acyl glucuronide reactive", "Reactive acyl glucuronide.", []),
        ("isocyanate", "Isocyanate reactive metabolite.", []),
        ("ketene", "Ketene reactive metabolite.", []),
        ("nitroso metabolite", "Nitroso metabolite class.", []),
        ("hydroxylamine metabolite", "Hydroxylamine metabolite class.", []),
        ("nitrene", "Nitrene reactive species.", []),
        ("free radical metabolite", "Radical metabolite species.", []),
        ("GSH adduct", "Glutathione adduct product class.", ["glutathione adduct"]),
        ("mercapturic acid", "Mercapturic acid product class.", []),
        ("glucuronide conjugate", "Glucuronide conjugate product class.", []),
        ("sulfate conjugate", "Sulfate conjugate product class.", []),
        ("phenol metabolite", "Phenolic metabolite product class.", []),
        ("catechol metabolite", "Catechol metabolite product class.", []),
        ("alcohol metabolite", "Alcohol metabolite product class.", []),
        ("carboxylic acid metabolite", "Carboxylic acid metabolite product class.", []),
        ("N-oxide metabolite", "N-oxide metabolite product class.", []),
        ("sulfoxide metabolite", "Sulfoxide metabolite product class.", []),
        ("sulfone metabolite", "Sulfone metabolite product class.", []),
        ("imine metabolite", "Imine metabolite product class.", []),
        ("carbinolamine intermediate", "Carbinolamine intermediate class.", []),
        ("hemiacetal intermediate", "Hemiacetal intermediate class.", []),
        ("epoxide metabolite", "Epoxide metabolite product class.", ["epoxide product"]),
        ("aldehyde metabolite", "Aldehyde metabolite product class.", []),
        ("iminium metabolite", "Iminium metabolite product class.", []),
        ("Michael acceptor metabolite", "Michael-acceptor metabolite class.", []),
        ("aziridine metabolite", "Aziridine reactive metabolite.", []),
        ("sulfenic acid metabolite", "Sulfenic acid intermediate/product.", []),
        ("nitro anion radical", "Nitro anion radical intermediate.", []),
    ]:
        add(pref, "xrm:1500000", definition, syn or None, "reactive")

    # --- site environment ---
    for pref, definition, syn in [
        ("benzylic carbon", "Benzylic carbon environment.", ["benzylic site environment"]),
        ("allylic carbon", "Allylic carbon environment.", []),
        ("phenol", "Phenolic site environment.", ["phenolic oxygen"]),
        ("aniline", "Aniline / aromatic amine site.", ["aromatic amine site"]),
        ("tertiary amine", "Tertiary amine nitrogen environment.", []),
        ("secondary amine", "Secondary amine nitrogen environment.", []),
        ("primary amine", "Primary amine nitrogen environment.", []),
        ("electron-rich arene", "Electron-rich aromatic ring environment.", []),
        ("electron-poor arene", "Electron-poor aromatic ring environment.", []),
        ("heteroarene", "Heteroaromatic ring environment.", ["heteroaryl site"]),
        ("furan", "Furan ring site environment.", []),
        ("thiophene", "Thiophene ring site environment.", []),
        ("imidazole", "Imidazole ring site environment.", []),
        ("indole", "Indole ring site environment.", []),
        ("pyridine", "Pyridine ring site environment.", []),
        ("pyrimidine", "Pyrimidine ring site environment.", []),
        ("alkyl ether", "Alkyl ether oxygen environment.", []),
        ("aryl ether", "Aryl ether oxygen environment.", []),
        ("thioether", "Thioether sulfur environment.", []),
        ("thiol", "Thiol sulfur environment.", []),
        ("carboxylic acid", "Carboxylic acid site environment.", []),
        ("ester", "Ester site environment.", []),
        ("amide", "Amide site environment.", []),
        ("carbamate", "Carbamate site environment.", []),
        ("urea", "Urea site environment.", []),
        ("nitrile", "Nitrile site environment.", []),
        ("nitroarene", "Nitroaromatic site environment.", []),
        ("alkyl halide", "Alkyl halide site environment.", []),
        ("aryl halide", "Aryl halide site environment.", []),
        ("terminal methyl", "Terminal methyl carbon environment.", []),
        ("omega carbon", "ω-carbon of an alkyl chain.", ["ω-carbon"]),
        ("omega-1 carbon", "ω-1 carbon of an alkyl chain.", ["ω-1 carbon"]),
        ("para aromatic carbon", "Para position on an aromatic ring.", []),
        ("ortho aromatic carbon", "Ortho position on an aromatic ring.", []),
        ("meta aromatic carbon", "Meta position on an aromatic ring.", []),
        ("attachment atom O", "Conjugation attachment at oxygen.", ["O-attachment"]),
        ("attachment atom N", "Conjugation attachment at nitrogen.", ["N-attachment"]),
        ("attachment atom S", "Conjugation attachment at sulfur.", ["S-attachment"]),
        ("attachment atom C", "Conjugation attachment at carbon.", ["C-attachment"]),
        ("acyl attachment", "Conjugation via acyl attachment.", []),
        ("alpha to nitrogen", "Carbon alpha to nitrogen.", []),
        ("alpha to oxygen", "Carbon alpha to oxygen.", []),
        ("alpha to sulfur", "Carbon alpha to sulfur.", []),
        # Bare hybridization / vinyl / alkyne carbons are out of scope
        # (general chem); see data/ontology/SCOPE.md.
    ]:
        add(pref, "xrm:1600000", definition, syn or None, "site")

    # --- medchem interpretation ---
    for pref, definition, syn in [
        ("aromatic soft spot", "Aromatic site that is a metabolic soft spot.", []),
        ("benzylic soft spot", "Benzylic metabolic soft spot.", []),
        ("amine soft spot", "Amine metabolic soft spot.", []),
        ("ether soft spot", "Ether metabolic soft spot (O-dealkylation risk).", []),
        ("thioether soft spot", "Thioether metabolic soft spot.", []),
        ("ester soft spot", "Ester metabolic soft spot (hydrolysis risk).", []),
        ("amide soft spot", "Amide metabolic soft spot.", []),
        ("clearance hotspot", "Site dominating metabolic clearance.", []),
        ("bioactivation hotspot", "Site associated with bioactivation risk.", []),
        ("detoxification pathway", "Pathway framed as detoxifying.", ["detoxication pathway"]),
        ("GSH trapping opportunity", "Site/path amenable to GSH trapping assays.", []),
        ("metabolic blocking vector", "Suggested vector for metabolic blocking.", []),
        ("fluorine blocking", "Fluorine substitution as a blocking strategy.", []),
        ("methyl blocking", "Methyl / alkyl blocking strategy.", []),
        ("scaffold hop opportunity", "Liability suggesting scaffold change.", []),
        ("polarity decreasing", "Transformation that decreases polarity.", []),
        ("lipophilicity decreasing", "Transformation that decreases lipophilicity.", []),
        ("reactive intermediate risk", "Risk of forming a reactive intermediate.", []),
        ("idiosyncratic toxicity risk", "Liability framed toward idiosyncratic toxicity.", []),
        ("time-dependent inhibition risk", "Liability framed toward TDI.", ["TDI risk"]),
        ("quinone risk", "Risk of quinone / quinone-imine formation.", []),
        ("acyl glucuronide risk", "Risk of reactive acyl glucuronide.", []),
        ("aldehyde risk", "Risk of aldehyde formation.", []),
        ("iminium risk", "Risk of iminium formation.", []),
        ("epoxide risk", "Risk of epoxide formation.", []),
        ("stable circulating metabolite", "Stable metabolite likely to circulate.", []),
        ("active metabolite opportunity", "Path that may yield an active metabolite.", []),
        ("prodrug activation", "Metabolic activation of a prodrug.", []),
        ("metabolic switching risk", "Blocking one site may switch metabolism elsewhere.", []),
        ("species-difference risk", "Liability sensitive to species differences.", []),
    ]:
        add(pref, "xrm:1400000", definition, syn or None, "medchem")

    # --- evidence / assertion ---
    for pref, definition, syn in [
        ("predicted assertion", "Assertion from a predictive model.", ["model-predicted"]),
        ("observed assertion", "Assertion from experimental observation.", []),
        ("curated assertion", "Manually curated assertion.", []),
        ("conflicting evidence", "Evidence sources conflict.", []),
        ("literature figure evidence", "Supported by a literature figure.", []),
        ("literature table evidence", "Supported by a literature table.", []),
        ("text-span evidence", "Supported by a cited text span.", []),
        ("exact-mass evidence", "Supported by accurate mass.", []),
        ("isotope pattern evidence", "Supported by isotopic pattern.", []),
        ("fragmentation evidence", "Supported by MS fragmentation.", []),
        ("chromatographic evidence", "Supported by chromatography beyond RT alone.", []),
        ("in vitro evidence", "Supported by in vitro experiment.", []),
        ("in vivo evidence", "Supported by in vivo experiment.", []),
        ("clinical evidence", "Supported by clinical observation.", []),
        ("radiolabel recovery", "Supported by radiolabel recovery/ADME.", []),
        ("trapping assay evidence", "Supported by trapping assay (e.g. GSH).", []),
        ("synthetic standard match", "Matches a synthetic authentic standard.", []),
        ("recombinant isoform evidence", "Supported by recombinant enzyme isoform data.", []),
        ("S9 incubation", "S9 fraction incubation context.", []),
        ("microsome incubation", "Microsomal incubation (non-HLM-specific).", []),
        ("hepatocyte evidence", "Hepatocyte-based evidence.", []),
        ("plasma observation", "Observed in plasma.", []),
        ("urine observation", "Observed in urine.", []),
        ("bile observation", "Observed in bile.", []),
        ("low confidence assertion", "Low-confidence assertion.", []),
        ("high confidence assertion", "High-confidence assertion.", []),
    ]:
        add(pref, "xrm:2000000", definition, syn or None, "evidence")

    # --- biological context (orthogonal) ---
    for pref, definition, syn in [
        ("CYP enzyme family context", "Cytochrome P450 family context (orthogonal).", ["P450 family context"]),
        ("UGT enzyme family context", "UGT family context (orthogonal).", []),
        ("SULT enzyme family context", "Sulfotransferase family context (orthogonal).", []),
        ("GST enzyme family context", "Glutathione S-transferase context (orthogonal).", []),
        ("NAT enzyme family context", "N-acetyltransferase context (orthogonal).", []),
        ("FMO enzyme family context", "Flavin monooxygenase context (orthogonal).", []),
        ("AO enzyme family context", "Aldehyde oxidase context (orthogonal).", []),
        ("MAO enzyme family context", "Monoamine oxidase context (orthogonal).", []),
        ("CES enzyme family context", "Carboxylesterase context (orthogonal).", []),
        ("liver tissue context", "Liver tissue framing.", []),
        ("intestine tissue context", "Intestine tissue framing.", []),
        ("kidney tissue context", "Kidney tissue framing.", []),
        ("lung tissue context", "Lung tissue framing.", []),
        ("brain tissue context", "Brain tissue framing.", []),
        ("human species context", "Human species framing.", []),
        ("rat species context", "Rat species framing.", []),
        ("mouse species context", "Mouse species framing.", []),
        ("dog species context", "Dog species framing.", []),
        ("monkey species context", "Nonhuman primate framing.", ["NHP species context"]),
        ("plasma matrix context", "Plasma matrix framing.", []),
        ("urine matrix context", "Urine matrix framing.", []),
        ("bile matrix context", "Bile matrix framing.", []),
        ("feces matrix context", "Feces matrix framing.", []),
        ("microsome matrix context", "Microsome incubation matrix.", []),
        ("cytosol matrix context", "Cytosol incubation matrix.", []),
        ("hepatocyte matrix context", "Hepatocyte incubation matrix.", []),
        ("endoplasmic reticulum compartment", "ER compartment framing.", []),
        ("cytosol compartment", "Cytosol compartment framing.", []),
        ("mitochondria compartment", "Mitochondrial compartment framing.", []),
        ("blood compartment", "Blood / systemic compartment framing.", []),
    ]:
        add(pref, "xrm:2100000", definition, syn or None, "bio")

    # --- structural delta extras ---
    for pref, definition, syn in [
        ("monooxygenation", "Net addition of one oxygen atom.", ["+O"]),
        ("dioxygenation", "Net addition of two oxygen atoms.", ["+2O"]),
        ("dehydrogenation delta", "Net loss of H2.", ["-2H"]),
        ("hydrogenation delta", "Net gain of H2.", ["+2H"]),
        ("demethylation delta", "Net loss of CH2 (demethylation mass shift).", ["-CH2"]),
        ("deethylation delta", "Net loss of C2H4.", ["-C2H4"]),
        ("glucuronide mass shift", "Mass shift for glucuronide conjugation.", ["+C6H8O6"]),
        ("sulfate mass shift", "Mass shift for sulfate conjugation.", ["+SO3"]),
        ("GSH mass shift", "Mass shift for glutathione conjugation.", ["+GSH"]),
        ("aromaticity loss delta", "Loss of aromaticity.", []),
        ("aromaticity gain delta", "Gain of aromaticity.", []),
        ("bond order increase", "Increase in bond order.", []),
        ("bond order decrease", "Decrease in bond order.", []),
        ("ring expansion delta", "Ring size increase.", []),
        ("ring contraction delta", "Ring size decrease.", []),
    ]:
        add(pref, "xrm:1700000", definition, syn or None, "delta")

    # --- product status extras ---
    for pref, definition, syn in [
        ("circulating metabolite", "Metabolite observed in circulation.", []),
        ("excreted metabolite", "Metabolite observed in excreta.", []),
        ("primary metabolite", "First-generation metabolite.", []),
        ("secondary metabolite", "Second-generation metabolite.", []),
        ("trace metabolite", "Trace-level metabolite.", []),
        ("putative metabolite", "Putative / unconfirmed metabolite.", []),
        ("structure elucidated", "Structure fully elucidated.", []),
        ("structure partially elucidated", "Structure only partially elucidated.", []),
    ]:
        add(pref, "xrm:1800000", definition, syn or None, "product")

    # --- rule provenance extras ---
    for pref, definition, syn in [
        ("Rainbow class derived", "Tag derived from Rainbow Phase I class mapping.", []),
        ("XenoNet path derived", "Tag derived from a multi-step XenoNet/Forest path.", []),
        ("legacy model output", "Tag mapped from a legacy model output label.", []),
        ("harvest-derived synonym", "Label/synonym promoted from offline harvest.", []),
        ("manual curation", "Tag asserted by manual curation.", []),
        ("priority override", "Higher-priority rule overrode a broader tag.", []),
        ("site template expansion", "Localized label generated from a site template (not a concept).", []),
    ]:
        add(pref, "xrm:1900000", definition, syn or None, "prov")

    # Localization / display-name templates are NOT SKOS concepts.
    # They live in data/ontology/site_templates.yml (see SCOPE.md).

    # Cross-spine related links (SKOS relatedMatch) for a few key pairs
    def find_id(label):
        for c in doc["concepts"]:
            if c["preferred_label"] == label:
                return c["id"]
        return None

    pairs = [
        ("N-dealkylation", "aldehyde forming"),
        ("N-dealkylation", "aldehyde metabolite"),
        ("aromatic hydroxylation", "phenol metabolite"),
        ("quinone formation", "quinone"),
        ("quinone-imine formation", "quinone imine"),
        ("epoxidation", "epoxide"),
        ("acyl glucuronidation", "acyl glucuronide"),
        ("glucuronidation", "glucuronide conjugate"),
        ("sulfation", "sulfate conjugate"),
        ("glutathionation", "GSH adduct"),
    ]
    by_label = {c["preferred_label"]: c for c in doc["concepts"]}
    for a, b in pairs:
        ca, cb = by_label.get(a), by_label.get(b)
        if not ca or not cb:
            continue
        for src, tgt in ((ca, cb["id"]), (cb, ca["id"])):
            rel = list(src.get("related_match") or [])
            if tgt not in rel:
                rel.append(tgt)
            src["related_match"] = rel

    doc["concepts"].sort(key=lambda c: c["id"])
    YAML_PATH.write_text(
        "# XRM thesaurus (authoring source). Export to xrm.skos.jsonld for Rust.\n"
        "# Site-localized names (e.g. C4 hydroxylation) are NOT concepts; generate from templates.\n"
        + yaml.safe_dump(doc, sort_keys=False, allow_unicode=True, width=100)
    )
    non_forest = sum(
        1
        for c in doc["concepts"]
        if not c["id"].startswith("xrm:9")
        and c["id"] != "xrm:0000000"
    )
    alts = sum(len(c.get("synonyms") or []) for c in doc["concepts"])
    print(f"concepts={len(doc['concepts'])} non_forest≈{non_forest} synonyms={alts}")


if __name__ == "__main__":
    main()
