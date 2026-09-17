# TODO

- Speed up `clean` / `_sanitize_kekulize` (dominant Full wall cost: sanitize + kekulize + SMILES round-trip per product)
- Speed up `AtomTracker.tags` (avoid `ast.literal_eval` / deepcopy on stringified atom-trace props)
- Speed up metabolite enumeration by skipping topologically equivalent SMARTS matches in the reaction loop (before RunReactants / tagging), not only at emit time
