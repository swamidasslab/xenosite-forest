# Divergences

Canonical-SMILES disagreements between this proof of concept and the old
library. A listing is the record. It is not permission to drop a real
metabolite the old code emits, and it is not permission to hide a product
only the new code emits. Sites may differ when `topol_equiv` puts them in
the same atom class. Order does not matter.

## Hydroxylation of butylbenzene

The old second pattern writes a carbonyl. Both patterns add OH and remove
one H. The `>>[*:1]=O` products are not metabolites.

- Reactant: `c1ccc(CCCC)cc1`
- Only old: `CC(=O)CCc1ccccc1`, `CCC(=O)Cc1ccccc1`, `CCCC(=O)c1ccccc1`
- Only new: none
- Shared: `CC(O)CCc1ccccc1`, `CCC(O)Cc1ccccc1`, `CCCC(O)c1ccccc1`, `CCCCc1ccc(O)cc1`, `CCCCc1cccc(O)c1`, `CCCCc1ccccc1O`, `OCCCCc1ccccc1`

Ethane (`CC`) does not hit that carbonyl. Both libraries yield `CCO` only.

More correct. Not a problem in the proof of concept. The three old-only strings are the `>>[*:1]=O` products, and those are not metabolites. Both patterns here add OH and remove one H, which is the shared set. Ethane does not hit that carbonyl, and both libraries already yield `CCO` only.

## Azo splitting of pyridazine

The old library kekulizes once. That form of pyridazine has an N-N single bond, so the split does not run. Every Kekulé form is searched here, and the N=N writing matches. RDKit's `[#7:1]=[#7:2]>>[*:1].[*:2]` then drops one nitrogen. `C=CC=CN` is not a metabolite: pyridazine is C4H4N2.

- Reactant: `c1ccnnc1`
- Only old: none
- Only new: `C=CC=CN`
- Shared: none

The same SMARTS on `c1ccc(N=Nc2ccnnc2)cc1` is what both libraries emit, including the ring-opened fragments. Those sets match.

Problem. The proof of concept is less correct on pyridazine. `C=CC=CN` is not a metabolite: the SMARTS drops one nitrogen, and the rule says both fragments stay. The old library emits nothing only because its one Kekulé form has an N-N single bond. The phenylazo-pyridazine sets already match, so the bad product is this ring writing.

## Nitro reduction of 1,4-dinitrobenzene

Kekulé copies rewrite one nitro as `[N+](=[O-])O` before the reaction. The oxygen then has valence 2. Nitrobenzene has no second nitro, and its fragments match.

- Reactant: `O=[N+]([O-])c1ccc([N+](=O)[O-])cc1`
- Only old: `O=Nc1ccc([N+](=O)[O-])cc1`
- Only new: `O=NC1=CC=C([N+](=[O-])O)C=C1` (does not parse)
- Shared: `O`

Problem. The proof of concept is less correct. The rule cleaves one N-O and leaves the nitroso, and the old-only SMILES is that nitroso with the second nitro still `[N+](=O)[O-]`. The new-only SMILES does not parse, because a Kekulé copy rewrote that nitro as `[N+](=[O-])O` and the oxygen has valence 2. Nitrobenzene has no second nitro, and its fragments match.

## Thiophene S-oxidation of benzothiophene

The old library kekulizes once. That form is `C1=CC=C2SC=CC2=C1`, so the thiophene ring is not `C=C-C=C-S` and the oxidation does not run. The other Kekulé form is that pattern, and the same SMARTS writes the S-oxide. Thiophene itself has no second writing to miss: both libraries emit `[O-][s+]1cccc1`. Dibenzothiophene is the same miss.

- Reactant: `c1ccc2sccc2c1`
- Only old: none
- Only new: `[O-][s+]1ccc2ccccc21`
- Shared: none

Dibenzothiophene (`c1ccc2c(c1)sc1ccccc12`): only new `[O-][s+]1c2ccccc2c2ccccc21`.

More correct. Not a problem in the proof of concept. The other Kekulé form is the `C=C-C=C-S` pattern, and the same SMARTS writes the S-oxide that only the proof of concept emits. The old library stops on `C1=CC=C2SC=CC2=C1`, which is not that pattern. Thiophene has no second writing, and both libraries emit `[O-][s+]1cccc1`. Dibenzothiophene is the same miss.
