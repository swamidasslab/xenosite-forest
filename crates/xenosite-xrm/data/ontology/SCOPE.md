# Scope: metabolism-specific only

XRM is **not** a full chemical ontology. Do not import or mirror ChEBI/RXNO/MOP
hierarchies wholesale. Keep terms that a medicinal chemist uses when reading
xenobiotic metabolism:

- reaction / transformation names
- Rainbow / Phase I–II families
- site environments and leaving groups (metabolism-facing motifs)
- reactive metabolite / conjugate product classes
- med-chem liability and pharmacological role (active / inactive / prodrug)
- whether a tag is about the **parent**, the **product**, or the **reaction**
- evidence, provenance, biological context (orthogonal)

ChEBI, Rhea, KEGG, RXNO, MOP, GO, MeSH, ECO, and CHMO are **mapping targets**
via SSSOM, not parents to extend.

## Out of scope (do not mint as SKOS)

- Bare hybridization / bonding atoms (`sp2 carbon`, `sp3 carbon`, `vinyl carbon`,
  `alkyne carbon`) without a metabolism handle
- General organic-chemistry hierarchies (functional groups, named reactions,
  spectroscopy, reagents) that are not used when reading xenobiotic metabolism
- Localization / display-name **templates** — those live in
  [`site_templates.yml`](site_templates.yml) and feed
  `AnnotationBundle::site_label`, not the concept inventory

If a candidate is general organic chemistry without a metabolism handle, leave
it out or keep it only as an external `closeMatch` / `relatedMatch`.
Site environments that chemists actually use for Soft Spots / bioactivation
(benzylic, allylic, aniline N, phenol O, ortho/meta/para aromatic, ω / ω−1)
stay in scope.
