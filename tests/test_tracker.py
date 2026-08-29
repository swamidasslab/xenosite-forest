from rdkit.Chem.rdchem import Mol
from xenosite.forest.base import AtomTracker
from xenosite.forest.bfs import bfs
from collections import defaultdict




def test_tracker():
    ruleset = "Full"
    reactant = "CC"
    product = "OCCO"
    depth = 2
    correct_phase1_sites = (frozenset({"1.h"}), frozenset({"3.h"}))


    result = list(
        bfs(
            [reactant, product],
            ruleset,
            depth=depth,
            all_paths=False,
            outmols=True,
            phase1=True,
            ismi=True,
        )
    )

    path_mols: list[Mol] = result[0][2] #type: ignore
    sites = [r[1] for r in result[0][1]]
    record: dict = AtomTracker.tags(path_mols[-1]) #type: ignore
    depths = list(range(len(path_mols)))
    # original_atoms = {k: v for k, v in list(record.items()) if "Depth0" in v}
    depth = len(result[0][1])

    # rulenames = [r[0] for r in result[0][1]]

    assert depth > 0

    # all depths represented
    assert depths == AtomTracker.depths(path_mols[-1]) 

    # identifies correct sites of metabolism
    if correct_phase1_sites:
        assert correct_phase1_sites == tuple(sites)

    # all atoms in molecules tagged
    for mol, d in zip(path_mols, depths):
        tags: dict = AtomTracker.tags(mol, depth=d) #type: ignore
        # print(d, tags)
        assert mol.GetNumHeavyAtoms() == len(tags)

    # atom traces are all same atom type
    # this loop shows key look up logic required to trace atoms
    tagged_atom_types = defaultdict(set)
    for mol, d in zip(path_mols, depths):
        r: dict = AtomTracker.tags(record, depth=d) # type: ignore

        for tag, info in list(r.items()):
            idx = info['idx'][info['depth'].index(d)]

            if depth in r:
                tagged_atom_types[tag].add(
                    mol.GetAtomWithIdx(idx).GetAtomicNum())

    for k,v in tagged_atom_types.items():
        assert len(v) == 1


    max_record_length = max([
        len(x['depth'])
        for x in list(AtomTracker.tags(record, depth=0).values()) #type: ignore
    ])

    assert max_record_length == len(depths)
    assert tuple(depths) == tuple(AtomTracker.depths(path_mols[-1]))
        


