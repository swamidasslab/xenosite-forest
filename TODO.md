# TODO

- Forest WAE: Epoxidation (and peers) still emit `CsmiDedupWarning` after aromatic `incident_orders` on ResonanceRule — chem suite is green with CSMI ignored; unique-edit not fully collapsing symmetric epoxides
- **Gate (PR #15 / pre-v0.7.0):** do not merge or tag v0.7.0 until `bond_atom` is gone, DH uses atom-pair Sites, `site_kind` + `_example_substrates` + `test_site_kind` are green, and PAIR_ORBITS/HEURISTICS/LOG docs match the atom-pair taxonomy
- Expand each rule's `_example_substrates` so examples cover all patterns/whens (meta-test currently only checks site_kind shape)

- Later (optional): hardcode `canonical_emitted_sites=True` on flaky forest tests (topo / guided gold / path hard / unique-edit) if lex site identity helps — fuzz already draws the flip via Hypothesis
- PatternInfo coverage gaps: product-side SMARTS parse; `pin` maps present in reactant; `skip_same_rings` / `edit` checked beyond presence; `breaks_ring` and `partner`/`partner_h` filled by `resolve_effect` (not only declared survival)
- `find_path_to_MS1`: condensed StepPlan-like form of reactions + possible sites to reach an MS1 m/z (z/charge) within tolerance; isotope-aware; all in-tolerance resolutions; alongside PhaseOne + conjugation, rules for common MS changes and MS-relevant input adjustments (e.g. C12→C14, add H+)
- After MS1 is solid/tested: `find_path_to_MS2` same idea plus MS2 fragmentation rules to develop
- Optimize `alternating_path`: quadratic because each BFS step copies the path; odd rings can leave the search (parity to the other side is fixed by whether the path skips parity), and those decomposed subproblems are then 2-colorable so all-pairs odd-distance BFS per anchor is not required
- xf: hang remaining helpers still free (`resonance_bond_maps`, MCS pairwise, `carry_forest` / split, `reordered_forest_labels`); drop thin shims `stamp_forest_labels` / `install_forest` / `ensure_forest` / `forest_trace` once tests use `mol.xf` only
- Deferred: tautomer SMARTS matching / TautomerQuery (preferred direction; Status: not decided) — orthogonal to landed `TautomerRule` stub (`NotImplementedError`; see docs/forest/HEURISTICS)
- Later: SMARTS least-common atoms first
- forest: phase I names still come from the search, not from `atom_trace["additions"][id]["phase1"]` (the field is reserved; the id, rules, effect, and formula delta are stored)
- Guided `and_cleave_plan` SMARTS: optional per-SMARTS flags naming atoms made dirty on match, so match-before/after deps are less leaky (fewer false precedes; less rely on apply-replay)
- Rank ``PathOutcome`` hits by likelihood via ``xenosite-predict`` using cleavage-side bags (without enumerating side reactions) — not this pass
- Guided-search heuristics from fuzz cases that hit ``max_expansions`` / miss under budget
- Use ``sanitize_dropped`` to find remaining rules/SMARTS that still clean-fail after include_sites / QF plan-only enum; further pre-clean rejects if needed
- Short-circuit mol edits that would yield invalid structures before sanitize/clean (edits are costly; reject early)
- Speed up `clean` / `_sanitize_kekulize` (dominant Full wall cost: sanitize + kekulize + SMILES round-trip per product)
