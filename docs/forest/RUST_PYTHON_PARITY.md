# Rust ↔ Python leaf-rule product parity

Tracking sheet for chematic (Rust) vs RDKit (Python) leaf parity. Product
identity is **RDKit CSMI** after round-trip. Site identity is **topological**
(engine-local ranks / pair orbits) — not raw atom indexes across engines.

Do not put exception name lists in test bodies. Unpaired or intentionally
divergent leaves use **rule attributes**:

- Python: ``ReactionRule.rust_parity_exception`` (non-empty reason string)
- Rust: ``RuleSet::parity_exception`` / ``with_parity_exception``

Corpus completeness (every PatternInfo possibility / every ``when``):
``tests/forest/test_rule_parity_corpus.py`` + ``PARITY_FUZZ_MOLS``.

Pairing completeness: ``tests/forest/test_rule_parity_pairs.py``.

Product/site parity: ``tests/forest/test_rule_parity_fuzz.py`` (parametric over
``PARITY_FUZZ_MOLS`` × paired leaves).

Native doors: ``xenosite_forest.RuleSet`` (``leaf`` / ``phase_one`` /
``default_ruleset`` / ``all_rules`` / ``catalog_names``), ``find_path``
(``ruleset=`` override), ``bfs`` / ``dfs`` / ``enumerate``. Python wrappers:
``xenosite.forest.find_path_rust``.

---

## Work order (tackle in this sequence)

| # | Issue | Status | Notes |
|---|--------|--------|-------|
| 0 | Corpus meta-test covers every pattern / every ``when`` | **done** | ``test_rule_parity_corpus``; grow ``PARITY_FUZZ_MOLS`` only |
| 1 | Discover all Rust + Python leaves; pair or attribute-exception | **done** | ``rule_parity_pairs``; EpoxideHydration / TautomerRule / ConjugationRule excused on attributes |
| 2 | Expose rulesets + ``find_path`` + enumerate to Python | **in progress** | ``leaf`` / catalogs / ``find_path(ruleset=)`` / ``bfs``/``dfs``/``enumerate`` wired; finish wrappers + rebuild |
| 3 | Parametric parity (bounded corpus × paired leaves), not Hypothesis | **done** | ``test_leaf_rule_product_parity`` |
| 4 | **Python pair unique-edit / topology dedup is wrong for pairs** | **open** | ResonancePair ``pair_site_signature`` can disagree with Rust ``atom_pair_orbit_id``. Fix Python (calling a Rust orbit/dedup helper is OK) **or** parity iteration without unique-edit dedup until fixed. Blocks fair site-count compares for DH/QF/Hydrogenation/…. |
| 5 | Acetylation: chematic SMIRKS apply miss | **open** | Match + unique-edit fire; ``[*:1][#6](=[#8])[#6]`` yields nothing; ``[*:1]C(=O)C`` applies. Align product SMIRKS with chematic dialect (least ad hoc: same spelling that applies). |
| 6 | Dephosphorylation: Rust emits nothing on ``COP(=O)(O)O`` | **open** | Likely recursive SMARTS ``$([#8][#6])`` / P valence on chematic. |
| 7 | AzoSplitting / ThiopheneSulfurOxidation: Rust empty on aromatic examples | **open** | ResonanceRule ``=,:`` path; Kekulé parent selection / aromatic match gap. |
| 8 | NitrogenReduction: Rust empty on ``CCNO`` (hydroxylamine) | **open** | Resonance / pattern arm not firing. |
| 9 | BenzodioxoleReduction: wrong products | **open** | Python catechol ``Oc1ccccc1O`` (+ C); Rust ring-opened junk. Adjudicate chemically; fix apply / aromatic handling. |
| 10 | SulfurOxidation: form / set mismatch on ``CCS`` | **open** | Python ``CCSO`` + ``CC[S+][O-]``; Rust ``CCS[O]`` (and missing zwitterion). Valence / canonical form. |
| 11 | Dehydration: site-count mismatch with matching products | **open** | Python both arms share atom site → one topo key; Rust alcohol vs beta ``site_atoms`` differ. May clear once #4 / site_map reporting aligned. |
| 12 | Drive parametric suite green; no silent skips | **blocked on 4–11** | Exceptions only via rule attributes when chemically justified. |

Depth-1 PhaseOne product diffs (separate probe, not leaf-only):
``tests/forest/probe_d1_diff.py`` / ``artifacts/d1_*`` — Dealkylation-heavy;
EpoxideHydration Rust-only (attribute exception above).

---

## Adjudication rule

When sets disagree: prefer **chemical correctness**, then **least ad hoc**
(shared SMARTS/effect data, not a one-off branch). Record lasting archive↔live
disagreements in [DIVERGENCES](DIVERGENCES.md). Record lasting Rust↔Python
gaps here until fixed or excused on the rule attribute.
