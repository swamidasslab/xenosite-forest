# Path-search performance (archive vs promoted forest)

Head-to-head `find_path` after the POC→`xenosite.forest` swap.

| Side | Package | Ruleset | API |
| --- | --- | --- | --- |
| **archive** | `xenosite._archive_forest` | PhaseOneQF (guided default) | rule-attached `guided_path.find_path` |
| **live** | `xenosite.forest` (ex-POC) | PhaseOne | `find_path.find_path` |

Tree identity for this table: post-swap working tree on `feature/rule-refactor` (includes `516d66f` ResonancePair path fix). Raw log: `artifacts/bench_find_path_h2h_post_swap.out`.

## Summary

- **10/10** live hits; **9/10** both hit. Archive exhausts on MeOPhOH→hydroxyquinone; live still finds it.
- Wall (all cases): archive **4.961s**, live **1.229s**.
- Wall (both_ok only, n=9): archive **2.684s**, live **0.207s** (~**13×**).

Billed units differ by design (archive: linearizations+site_applies; live: mol_edits+nodes) — compare wall and hit/miss first.

## Re-run

```bash
uv run python tests/forest/bench_find_path_h2h.py
```

Writes `artifacts/bench_find_path_h2h_post_swap.{out,live.log}`.
