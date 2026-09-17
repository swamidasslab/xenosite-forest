# TODO

- Short-circuit mol edits that would yield invalid structures before sanitize/clean (edits are costly; reject early)
- Speed up `clean` / `_sanitize_kekulize` (dominant Full wall cost: sanitize + kekulize + SMILES round-trip per product)
- Speed up `AtomTracker.tags` (avoid `ast.literal_eval` / deepcopy on stringified atom-trace props)
- Consider recording `mol._forest["atom_refs"]` inside rules (optional) so BFS/`metabolize` share creation tracking with `Step.apply` — defer to atom-ref / tags refactor
- Speed up metabolite enumeration by skipping topologically equivalent SMARTS matches in the reaction loop (before RunReactants / tagging), not only at emit time
