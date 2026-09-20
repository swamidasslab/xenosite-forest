# TODO

- Defer epoxidation and N-dealkylation phase-I shapes; those rules are not in the default set
- Later: SMARTS least-common atoms first, and the isotope atom first when a tag is used
- Tautomer matching (TautomerQuery versus a Tautomerization rule) is deferred; it is out of this phase because it is not needed for Phase I + quinone formation + conjugation
- refactor_poc PhaseOne omits Tautomerization, Sulfation, Glucuronidation, and Glutathionation
- refactor_poc: phase I names still come from the search, not from `atom_trace["additions"][id]["phase1"]` (the field is reserved; the id, rules, effect, and formula delta are stored)
- refactor_poc search: consult `src/xenosite/refactor_poc/HEURISTICS.md` before changing find_path or explaining a miss
- Guided `and_cleave_plan` SMARTS: optional per-SMARTS flags naming atoms made dirty on match, so match-before/after deps are less leaky (fewer false precedes; less rely on apply-replay)
- Rank ``PathOutcome`` hits by likelihood via ``xenosite-predict`` using cleavage-side bags (without enumerating side reactions) — not this pass
- Guided-search heuristics from fuzz cases that hit ``max_expansions`` / miss under budget
- Use ``sanitize_dropped`` to find remaining rules/SMARTS that still clean-fail after include_sites / QF plan-only enum; further pre-clean rejects if needed
- Short-circuit mol edits that would yield invalid structures before sanitize/clean (edits are costly; reject early)
- Speed up `clean` / `_sanitize_kekulize` (dominant Full wall cost: sanitize + kekulize + SMILES round-trip per product)
- Speed up metabolite enumeration by skipping topologically equivalent SMARTS matches in the reaction loop (before RunReactants / tagging), not only at emit time
