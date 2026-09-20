"""bfs and dfs enumerate one ruleset. They do not search toward a target."""

from rdkit import Chem

from xenosite.forest.find_path import bfs, dfs
from xenosite.forest.rules import Hydroxylation
from xenosite.forest.rulesets import RuleSet


def _oxygens(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return sum(atom.GetAtomicNum() == 8 for atom in mol.GetAtoms())


def test_depth_one_hydroxylation_of_ethane_is_ethanol():
    ruleset = RuleSet((Hydroxylation,), name="Poc")
    smiles = {info["csmi"] for _mol, info in bfs("CC", ruleset, depth=1)}
    assert smiles == {"CCO"}


def test_dfs_reaches_depth_two_before_the_frontier_is_done():
    ruleset = RuleSet((Hydroxylation,), name="Poc")
    reactant = "c1ccc(CCCC)cc1"
    dfs_oxygen = []
    for _mol, info in dfs(reactant, ruleset, depth=2):
        dfs_oxygen.append(_oxygens(info["csmi"]))
        if len(dfs_oxygen) == 2:
            break
    assert dfs_oxygen == [1, 2]

    bfs_oxygen = [
        _oxygens(info["csmi"]) for _mol, info in bfs(reactant, ruleset, depth=2)
    ]
    first_diol = bfs_oxygen.index(2)
    assert first_diol > 1
    assert set(bfs_oxygen[:first_diol]) == {1}


def _canon(smiles):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles), canonical=True)


def test_filter_sites_keeps_the_chain_alcohol_and_drops_the_ring():
    reactant = "c1ccc(CCCC)cc1"
    mol = Chem.MolFromSmiles(reactant)
    ring = {atom.GetIdx() for atom in mol.GetAtoms() if atom.IsInRing()}

    def filter_sites(mol, site, info):
        atoms = site if isinstance(site, frozenset) else frozenset((site,))
        return atoms.isdisjoint(ring)

    ruleset = RuleSet((Hydroxylation,), name="Poc")
    smiles = {
        _canon(info["csmi"])
        for _mol, info in bfs(
            reactant, ruleset, filter_sites=filter_sites, depth=1
        )
    }
    assert _canon("OCCCCc1ccccc1") in smiles
    assert _canon("CCCCc1ccc(O)cc1") not in smiles
    assert _canon("CCCCc1ccccc1O") not in smiles


def test_filter_rules_can_refuse_the_only_child():
    ruleset = RuleSet((Hydroxylation,), name="Poc")
    products = list(
        bfs("CC", ruleset, filter_rules=lambda mol, rule, info: False, depth=1)
    )
    assert products == []
