"""RULE MODIFICATION EXAMPLE

>>> from rdkit import Chem
>>> mol = Chem.MolFromSmiles("CN")

>>> OriginalRule = rules.Dealkylation()
>>> ModifiedRule = rules.Dealkylation(rxns=rules.Dealkylation.smarts[1:])

>>> len(list(ModifiedRule.metabolites(mol)))
2

>>> len(list(OriginalRule.metabolites(mol)))
3


To see all the rule names of a rule set constructed from other rule sets.
>>> print('\\n'.join(sorted([rule.name for ruleset in  PhaseOneRS for rule in ruleset])))
Dealkylation
Dehydration
Dehydrogenation
Dephosphorylation
Epoxidation
EpoxideOpening
Hydrogenation
Hydrolysis
Hydroxylation
NitrogenOxidation
NitrogenReduction
OxidativeDehalogenation
OxygenReduction
ReductiveDehalogenation
SulfurOxidation
SulfurReduction

Current phase1 site labeling paradigms

>>> print('\\n'.join(sorted([rule.name+'\t'+rule.phase1_sites_on for ruleset in  PhaseOneRS for rule in ruleset])))
Dealkylation bonds
Dehydration bonds
Dehydrogenation atom_hydrogen
Dephosphorylation bonds
Epoxidation bonds
EpoxideOpening bonds
Hydrogenation atoms
Hydrolysis bonds
Hydroxylation atom_hydrogen
NitrogenOxidation atoms
NitrogenReduction bonds
OxidativeDehalogenation bonds
OxygenReduction bonds
ReductiveDehalogenation bonds
SulfurOxidation atoms
SulfurReduction bonds

>>> from xenosite.metabolite.bfs import bfs
>>> next(bfs(['CCN','CCO'],phase1=True,outmols=False))
"(['CCN', 'CCO'], [('Dealkylation', frozenset({'2.3'}))])"

>>> next(bfs(['CCO','CC'],phase1=True,outmols=False))
"(['CCO', 'CC'], [('Dehydration', frozenset({'2.3'}))])"

# >>> next(bfs(['CCO','C=CO'],phase1=True,outmols=False))
# "(['CCO', 'C=CO'], [('Dehydrogenation', frozenset({'2.h', '1.h'}))])"

>>> next(bfs(['COP(=O)(O)O','CO'],phase1=True,outmols=False))
"(['COP(=O)(O)O', 'CO'], [('Dephosphorylation', frozenset({'2.3'}))])"

>>> next(bfs(['C=C','C1OC1'],phase1=True,outmols=False))
"(['C=C', 'C1CO1'], [('Epoxidation', frozenset({'1.2'}))])"

>>> next(bfs(['C1OC1','CCO'],phase1=True,outmols=False))
"(['C1CO1', 'CCO'], [('EpoxideOpening', frozenset({'1.2'}))])"

# >>> next(bfs(['C=CO','CCO'],phase1=True,outmols=False))
# "(['C=CO', 'CCO'], [('Hydrogenation', frozenset({'2.2', '1.1'}))])"

>>> next(bfs(['CC','CCO'],phase1=True,outmols=False))
"(['CC', 'CCO'], [('Hydroxylation', frozenset({'1.h'}))])"

>>> next(bfs(['O=C(O)C','CC=O'],phase1=True,outmols=False))
"(['CC(=O)O', 'CC=O'], [('Hydrolysis', frozenset({'2.3'}))])"

>>> next(bfs(['CCN','CCNO'],phase1=True,outmols=False))
"(['CCN', 'CCNO'], [('NitrogenOxidation', frozenset({'3.3'}))])"

>>> next(bfs(['CCNO','CCN'],phase1=True,outmols=False))
"(['CCNO', 'CCN'], [('Dehydration', frozenset({'3.4'}))])"
>>> next(bfs(['[O-]-[N+](C1=CC=C(O1)C=NN2C(=O)NC(=O)C2)=O','N(C1=CC=C(O1)C=NN2C(=O)NC(=O)C2)=O'],phase1=True,outmols=False))
"(['O=C1CN(N=CC2=CC=C([N+](=O)[O-])O2)C(=O)N1', 'O=Nc1ccc(C=NN2CC(=O)NC2=O)o1'], [('NitrogenReduction', frozenset({'1.2'}))])"

>>> next(bfs(['CC=O','CCO'],phase1=True,outmols=False)) # doctest: +SKIP
"(['CC=O', 'CCO'], [('Hydrogenation', frozenset({'2.2', '3.3'}))])"

>>> next(bfs(['CCCl','CCO'],phase1=True,outmols=False))
"(['CCCl', 'CCO'], [('OxidativeDehalogenation', frozenset({'2.3'}))])"

>>> next(bfs(['CCCl','CC'],phase1=True,outmols=False))
"(['CCCl', 'CC'], [('ReductiveDehalogenation', frozenset({'2.3'}))])"

>>> next(bfs(['CCS','CCSO'],phase1=True,outmols=False))
"(['CCS', 'CCSO'], [('SulfurOxidation', frozenset({'3.3'}))])"

>>> next(bfs(['CCSO','CCS'],phase1=True,outmols=False))
"(['CCSO', 'CCS'], [('SulfurReduction', frozenset({'3.4'}))])"

"""

# First Party
from .rulesets import RuleSet
from . import rules as rules

StableOxygenation_PhaseOne = RuleSet(name='SO', longname='StableOxygenation')
Dehydrogenation_PhaseOne = RuleSet(name='DH', longname='Dehydrogenation')
Hydrolysis_PhaseOne = RuleSet(name='HD', longname='Hydrolysis')
Reduction_PhaseOne = RuleSet(name='RD', longname='Reduction')
UnstableOxygenation_PhaseOne = RuleSet(
    name='UO', longname='Unstable Oxygenation')

Dehydrogenation_PhaseOne.add_rule(rules.Dehydrogenation())

Hydrolysis_PhaseOne.add_rule(rules.Dephosphorylation())
Hydrolysis_PhaseOne.add_rule(rules.EpoxideOpening())
Hydrolysis_PhaseOne.add_rule(rules.Hydrolysis())

Reduction_PhaseOne.add_rule(rules.Dehydration())
Reduction_PhaseOne.add_rule(rules.Hydrogenation())
Reduction_PhaseOne.add_rule(rules.NitrogenReduction())
Reduction_PhaseOne.add_rule(rules.OxygenReduction())
Reduction_PhaseOne.add_rule(rules.ReductiveDehalogenation())
Reduction_PhaseOne.add_rule(rules.SulfurReduction())

StableOxygenation_PhaseOne.add_rule(rules.Hydroxylation())
StableOxygenation_PhaseOne.add_rule(rules.Epoxidation())
StableOxygenation_PhaseOne.add_rule(rules.SulfurOxidation())
StableOxygenation_PhaseOne.add_rule(rules.NitrogenOxidation())

UnstableOxygenation_PhaseOne.add_rule(rules.Dealkylation())
UnstableOxygenation_PhaseOne.add_rule(rules.OxidativeDehalogenation())

PhaseOneRS = RuleSet(rules=[
    Dehydrogenation_PhaseOne, Hydrolysis_PhaseOne, Reduction_PhaseOne,
    StableOxygenation_PhaseOne, UnstableOxygenation_PhaseOne
])

PhaseOneRS.rulenames = sorted([x.name for x in PhaseOneRS])


if __name__ == '__main__':
    import doctest
    doctest.testmod()
