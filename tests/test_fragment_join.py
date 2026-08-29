import pytest
from rdkit.Chem.rdmolfiles import MolFromSmiles, MolToSmiles
from rdkit.Chem.rdmolops import RemoveStereochemistry
from xenosite.metabolite.base import AromaticSystems, Resonate

from test_quinone import examples


def can_mol(m):
    
    m = MolFromSmiles(MolToSmiles(m))

    if m is None:
        return ""
    RemoveStereochemistry(m)

    return MolToSmiles(m)  



RES = Resonate()

@pytest.mark.parametrize(
    "reactant",
    argvalues=[x[1] for x in examples],
    ids=[x[0] for x in examples],
)
def test_fragmentation_and_join(reactant):
    mol = MolFromSmiles(reactant)
    AS = AromaticSystems()

    for fragment_to_resonante, the_other_fragments, system, bond_types in AS.fragments(
            mol):

        joined = RES.join_fragments(
            [fragment_to_resonante] + the_other_fragments, bond_types)
        
        assert can_mol(mol) == can_mol(joined)