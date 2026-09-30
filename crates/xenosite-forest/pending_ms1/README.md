# Pending MS1 port (not compiled)

Parked from `cursor/find-path-ms1-7f58` until `ApplyN` lands on HEAD
`canonical_plan`, then move into `src/` / `tests/` / `examples/`.

| File | Role |
|------|------|
| `mass.rs` | **Already on branch** (`src/mass.rs`) — 10 unit tests green |
| `find_path_ms1.rs` | ~30 unit tests + internal fuzz; needs `ApplyN` + Deps::with_apply_n |
| `ms1_apply_fuzz.rs` | Integration fuzz (apply→recover by m/z); DeferredSite APIs adapted; needs `failure_persistence` regressions file |
| `apply_n_bench.rs` | Bench for ApplyN emit |

## Gate properties (do not weaken)

- One-hop: every materialized leaf product CSMI appears in `find_path_ms1` hits at that product's [M+H]⁺.
- Hits: last step / hit smiles within tol; path mz-error monotonic; plan linearizations replay.
- No redundant plans among isobars (`same_linearizations` / same-product skeleton twins).
- ApplyN count=N recovers N-hop chains; refuse when count exceeds reachable at mass.
- Mixed two-hop corpus (`MIXED_TWO_HOP_CASES`) deterministic + fuzz.

