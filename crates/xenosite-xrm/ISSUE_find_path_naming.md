# Proposal: Optional find_path adapter for reaction naming tags

**Status:** open for discussion (not implemented in `xenosite-xrm`)  
**GitHub issue:** https://github.com/swamidasslab/xenosite-forest/issues/31  
**Related PR:** https://github.com/swamidasslab/xenosite-forest/pull/30

## Summary

Keep `xenosite-xrm` free of forest imports. Optionally add a **separate** thin adapter (new crate or `xenosite-forest` helper) that calls Rust `find_path` on a reactant→product pair and feeds opaque tags into `Namer::name`.

## Why this is desirable

Structural SMARTS/delta assignment covers many Phase I transforms, but:

- Multi-step paths (e.g. quinone formation as OH + DH plan) are hard to encode as a single reactant/product SMARTS pair.
- Phase II conjugations often use star adducts / large groups; tag-based naming is more reliable when the emitter is known.
- `find_path` already returns the identity the thesaurus maps to.

## Evidence from Rust `find_path` API

`PathOutcome` / `PathStep` (in `xenosite-forest`) already expose:

| Field | Use as opaque tag |
| --- | --- |
| `PathStep.leaf_rule()` | `forest.rule:{name}` |
| `PathStep.pattern_name` | `forest.pattern:{leaf}/{pattern}` |
| `PathStep.namespace()` / `rule_path` | `forest.ruleset:{name}` when present |
| `PathOutcome.plan` elementary `Step.rule` | same as leaf rule strings |

Python binding shape today: `{"smiles", "steps": [{"rule", "site"}]}` from chematic `find_path`.

Suggested adapter contract (outside `xenosite-xrm`):

```text
find_path(reactant, product) → tags: ["forest.rule:Hydroxylation", "forest.pattern:Hydroxylation/h2", ...]
Namer::name(MappedReaction { reactant, product, tags, ... })
```

## Boundary rule

- **Do not** add `xenosite-forest` as a dependency of `xenosite-xrm`.
- Adapter may depend on both; core namer remains config-only.
- SSSOM file `data/mappings/xrm-forest.sssom.tsv` already declares the CURIE vocabulary.

## Open questions

1. Single best path vs all paths within `max_paths` when tagging?
2. Should plan elementary steps (QF → Hydroxylation + Dehydrogenation) emit multiple chemical terms or the composite `quinone formation` only?
3. Conjugation / `Full` catalog availability on the Rust default ruleset vs PhaseOne-only `find_path_default`.

## Acceptance for a follow-up PR

- Adapter crate or forest helper with no changes to xrm assignment engine semantics.
- Golden tests: ethane→ethanol tags include `forest.rule:Hydroxylation`; namer still returns Phase I → SO → hydroxylation → aliphatic….
- Boundary test in xrm continues to fail if forest is added to `Cargo.toml`.
