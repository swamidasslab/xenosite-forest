# Divergences

Canonical-SMILES disagreements between this proof of concept and the old
library. A listing is the record. It is not permission to drop a real
metabolite the old code emits, and it is not permission to hide a product
only the new code emits. Sites may differ when `topol_equiv` puts them in
the same atom class. Order does not matter.

Also give an assessment of whether the divergence is a problem or not, whether it is more or less correct and why, or if it is neutral.

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

The old library kekulizes once. That form of pyridazine has an N-N single bond, so the split does not run. Every Kekulé form is searched here, and the N=N writing matches. RDKit's `[#7:1]=[#7:2]>>[*:1].[*:2]` then drops one nitrogen. `C=CC=CN` is not a metabolite: pyridazine is C4H4N2. Both Kekulé forms parse, and neither atom is charged, so the dropped nitrogen is the reaction, not a charge left on the wrong atom.

- Reactant: `c1ccnnc1`
- Only old: none
- Only new: `C=CC=CN`
- Shared: none

The same SMARTS on `c1ccc(N=Nc2ccnnc2)cc1` is what both libraries emit, including the ring-opened fragments. Those sets match.

Problem. The proof of concept is less correct on pyridazine. `C=CC=CN` is not a metabolite: the SMARTS drops one nitrogen, and the rule says both fragments stay. The old library emits nothing only because its one Kekulé form has an N-N single bond. The phenylazo-pyridazine sets already match, so the bad product is this ring writing. Moving formal charge with the bond flip does not change it.

## Nitro reduction of 1,4-dinitrobenzene

A Kekulé copy used to rewrite one nitro as `[N+](=[O-])O` and leave the charge on the oxygen that had become double-bonded. That oxygen has valence 2, so the SMILES did not parse. The charge now follows the bond order. Both libraries yield the nitroso, and the leaving oxygen.

- Reactant: `O=[N+]([O-])c1ccc([N+](=O)[O-])cc1`
- Only old: none
- Only new: none
- Shared: `O=Nc1ccc([N+](=O)[O-])cc1`, `O`

Not a problem. The nitroso is the metabolite the rule describes, and it parses. Nitrobenzene was already that pair.

## Thiophene S-oxidation of benzothiophene

The old library kekulizes once. That form is `C1=CC=C2SC=CC2=C1`, so the thiophene ring is not `C=C-C=C-S` and the oxidation does not run. The other Kekulé form is that pattern, and the same SMARTS writes the S-oxide. Thiophene itself has no second writing to miss: both libraries emit `[O-][s+]1cccc1`. Dibenzothiophene is the same miss.

- Reactant: `c1ccc2sccc2c1`
- Only old: none
- Only new: `[O-][s+]1ccc2ccccc21`
- Shared: none

Dibenzothiophene (`c1ccc2c(c1)sc1ccccc12`): only new `[O-][s+]1c2ccccc2c2ccccc21`.

More correct. Not a problem in the proof of concept. The other Kekulé form is the `C=C-C=C-S` pattern, and the same SMARTS writes the S-oxide that only the proof of concept emits. The old library stops on `C1=CC=C2SC=CC2=C1`, which is not that pattern. Thiophene has no second writing, and both libraries emit `[O-][s+]1cccc1`. Dibenzothiophene is the same miss.

## Arene-oxide methyl sulfone

The old second SMARTS deletes the epoxide oxygen and writes `-S(C)(=O)(=O)`. The methyl carbon is not in the reactant. The product is not a sulfate conjugate.

- Reactant: `C1=CC2OC2C=C1`
- Only old: `CS(=O)(=O)c1ccccc1`
- Only new: none
- Shared: none

The same writing on `Clc1ccc(c(c1)Cl)C1=CC2(Cl)OC2C=C1Cl` is only old `CS(=O)(=O)c1cc(Cl)c(-c2ccc(Cl)cc2Cl)cc1Cl`.

`OC1=CC2OC2C=C1` still shares the sulfate `O=S(=O)(O)OC1=CC2OC2C=C1`. Only old there: `CS(=O)(=O)c1ccc(O)cc1`, `CS(=O)(=O)c1cccc(O)c1`.

More correct. Not a problem in the proof of concept. Those strings are methyl sulfones. A sulfate conjugate keeps the oxygen and adds `S(=O)(=O)O`. The epoxide oxygen has no hydrogen, so it is not a sulfate site.

## Carbonyl heteroatom rewritten as oxygen

The old first SMARTS is `[#8:1][#6:2](=[O,N,P,S:3])[#6:4]` and the product writes `=[#8:3]`. Nitrogen, phosphorus, and sulfur become oxygen. That product is a carboxylic glucuronide, not the glucuronide of the matched group.

- Reactant: `N=C(O)c1ccccc1`
- Only old: `O=C(OC1OC(C(=O)O)C(O)C(O)C1O)c1ccccc1`
- Only new: none
- Shared: `N=C(OC1OC(C(=O)O)C(O)C(O)C1O)c1ccccc1`

The same writing on `S=C(O)c1ccccc1` is only old `O=C(OC1OC(C(=O)O)C(O)C(O)C1O)c1ccccc1`. Shared: `O=C(O)C1OC(OC(=S)c2ccccc2)C(O)C(O)C1O`.

On `N=C([O-])C` the alcohol pattern does not match. Only old: `CC(=O)OC1OC(C(=O)O)C(O)C(O)C1O`. Shared: none. `S=C([O-])C` and `P=C([O-])C` are the same only-old carboxylate.

More correct. Not a problem in the proof of concept. Those strings replace the heteroatom with oxygen. A glucuronide of `N=C(O)-` or `S=C(O)-` keeps that atom. The pattern here matches `=[#8]` on an OH or an anion, so benzoate still gives the acyl glucuronide. An ester oxygen is not that site: the old product does not sanitize, and this pattern does not match it.
