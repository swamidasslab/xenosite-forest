# Scope: metabolism-specific only

XRM is **not** a full chemical ontology. Do not import or mirror ChEBI/RXNO/MOP
hierarchies wholesale. Keep terms that a medicinal chemist uses when reading
xenobiotic metabolism:

- reaction / transformation names
- Rainbow / Phase I–II families
- site environments and leaving groups
- reactive metabolite / conjugate product classes
- med-chem liability and pharmacological role (active / inactive / prodrug)
- whether a tag is about the **parent**, the **product**, or the **reaction**
- evidence, provenance, biological context (orthogonal)

ChEBI, Rhea, KEGG, RXNO, MOP, GO, MeSH, ECO, and CHMO are **mapping targets**
via SSSOM, not parents to extend.

If a candidate is general organic chemistry without a metabolism handle, leave
it out or keep it only as an external `closeMatch` / `relatedMatch`.
