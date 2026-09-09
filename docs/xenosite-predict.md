# Notes from xenosite-predict

Downstream package: sibling [`xenosite-predict`](https://github.com/swamidasslab/xenosite-predict) (`../xenosite-predict`). It enumerates Phase I and conjugate structures with this library. Two adapters live there until this package grows the same behavior.

## RDKit 2026 valence caches

Python-2 XenoNet and older RDKit compute implicit/explicit hydrogens lazily inside `RunReactants` and `MolToSmiles`. RDKit 2026 asserts those caches already exist (`getNumImplicitHs` / `getValence` without `calcImplicitValence`). Resonance copies from `join_fragments` look like valid mols but have an empty cache, so `Dehydrogenation` dies on the first aromatic drug.

This checkout already calls `UpdatePropertyCache(strict=False)` in `SmartsReactionRule.metabolites` and skips the rule on `RuntimeError`. The **installed wheel** xenosite-predict was using did not. Predict therefore shims the call sites it hits:

[`../xenosite-predict/src/xenosite/predict/forest_rdkit.py`](../../xenosite-predict/src/xenosite/predict/forest_rdkit.py)

That file replaces `SmartsReactionRule.metabolites`, `_kekulize`, `clean`, `can_smi`, and `RuleSet.metabolites` (unique-key SMILES). Please fold the same refreshes into this tree (including `clean` fragments and unique SMILES), and **do not `return` the whole rule** when one resonance copy fails — skip that copy and continue.

### Examples that crashed on the unpatched wheel, then enumerated

These are from the 327-molecule descriptor suite vs Python-2 XenoNet. On RDKit 2026.03.5, `PhaseOneRS.metabolites(..., unique=True)` raised a valence precondition until the shim. After the shim they produce metabolites (diphenhydramine and ibuprofen are in predict’s `tests/test_forest_rdkit.py`).

| Name / note | SMILES |
| --- | --- |
| Diphenhydramine | `CN(C)CCOC(c1ccccc1)c1ccccc1` |
| Ibuprofen | `CC(C)Cc1ccc(C(C)C(=O)O)cc1` |
| Atenolol | `CC(C)NCC(O)COc1ccc(CC(N)=O)cc1` |
| Cinnarizine | `C(=Cc1ccccc1)CN1CCN(C(c2ccccc2)c2ccccc2)CC1` |
| Warfarin | `CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O` |
| Trimethoprim | `COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC` |
| Omeprazole-like | `COc1cc2nc(SCc3ccccc3C)[nH]c2cc1OC` |
| Fluconazole | `OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F` |
| Sulfisoxazole | `Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1` |
| Phenobarbital | `CCC1(c2ccccc2)C(=O)NC(=O)NC1=O` |
| Chloramphenicol | `O=C(NCC(O)c1ccc([N+](=O)[O-])cc1)C(Cl)Cl` |
| Atropine | `CN1C2CCC1CC(OC(=O)C(CO)c1ccccc1)C2` |
| Penicillin G core | `CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N21` |
| Sudoxicam-like | `CN1C(C(=O)Nc2nccs2)=C(O)c2ccccc2S1(=O)=O` |

On the 327-molecule suite the shim took **171 crashes → 0**. Shared Phase I edge *weights* already matched Python-2; the crash was structure enumeration, not ONNX.

## Conjugation rules (DNA / cyanide vs GSH)

Forest ships `Glucuronidation` and `Glutathionation` under `CJ`. Predict maps:

| Predict head | Forest ruleset |
| --- | --- |
| `ugt` | `CJ.Glucuronidation` |
| `reactivity.gsh` | `CJ.Glutathionation` |
| `reactivity.protein` | `CJ.Glutathionation` (star label `Protein`) |
| `reactivity.dna` | glutathionation **without** the thiol-disulfide SMARTS |
| `reactivity.cyanide` | same as DNA |

`Glutathionation.smarts` includes `[#16h1:1]>>…` (substrate thiol → mixed disulfide). That is a GSH reaction, not DNA or cyanide. Predict subclasses:

```python
class GlutathionationNoThiol(Glutathionation):
    smarts = [s for s in Glutathionation.smarts if "[#16h1" not in s]
```

Adapter: [`../xenosite-predict/src/xenosite/predict/conjugates.py`](../../xenosite-predict/src/xenosite/predict/conjugates.py).

Please add a built-in variant (or a flag) that keeps epoxide, C–Cl, and terminal alkene and drops `[#16h1`. Dummy `*` atoms in predict are labeled with CXSMILES (`GlcA` / `GSH` / `Protein` / `DNA` / `CN`); forest can keep emitting a bare `*`.

### Examples where GSH vs DNA/CN output differs

Thiol substrates match `[#16h1` and must **not** appear as DNA/CN conjugates:

| Note | SMILES |
| --- | --- |
| Aliphatic thiol | `CCS` |
| Thiophenol | `Sc1ccccc1` |
| Cysteine-like thiol | `SC[C@H](N)C(=O)O` |

Electrophiles that **should** still conjugate for GSH, DNA, and CN (no thiol SMARTS):

| Note | SMILES |
| --- | --- |
| Styrene oxide | `c1ccccc1C1OC1` |
| Benzyl chloride | `ClCc1ccccc1` |
| Terminal alkene | `C=CC` |
