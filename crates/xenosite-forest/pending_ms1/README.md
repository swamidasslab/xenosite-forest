# Pending MS1 port (not compiled)

Parked from `cursor/find-path-ms1-7f58` until wired onto HEAD.

| File | Role |
|------|------|
| `mass.rs` | **On branch** (`src/mass.rs`) |
| `ApplyN` | **On HEAD** (`canonical_plan`) |
| `find_path_ms1.rs` | **Moved to** `src/find_path_ms1.rs` |
| `ms1_apply_fuzz.rs` | **Moved to** `tests/ms1_apply_fuzz.rs` (+ `proptest-regressions/ms1_apply_fuzz.txt`) — green |
| `apply_n_bench.rs` | Bench for ApplyN emit (still parked) |

## Gate properties (do not weaken)

- One-hop: every materialized leaf product CSMI appears in `find_path_ms1` hits at that product's [M+H]⁺.
- Hits: last step / hit smiles within tol; path mz-error monotonic; plan linearizations replay.
- No redundant plans for the **same** product CSMI (`same_linearizations` / same-product skeleton twins). Distinct isobar CSMIs may share a rule/site plan (smile-gated yield).
- ApplyN count=N recovers N-hop chains; refuse when count exceeds reachable at mass.
- Mixed two-hop corpus (`MIXED_TWO_HOP_CASES`) deterministic + fuzz (`max_paths` ≥ 128 for dimethoxy Dealk×2). Omit BzdRed+OH (non-monotonic mz) and anisole Dealk+DH (no DH sites after dealk).
