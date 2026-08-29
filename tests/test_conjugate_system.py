import pytest
from xenosite.forest.base import AromaticSystems
from rdkit.Chem.rdmolfiles import MolFromSmiles


@pytest.mark.parametrize(
    "smi, expected",
    [
        (
            "c1ccc2c(c1)c3cccc4ccc5cccc2c5c34",
            [{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}],
        ),
        ("C(C=CC=CCC(C=CC=CC1=CC=CC=C1)=CC)C", [{11, 12, 13, 14, 15, 16}]),
        (
            "C1=C5C=CC(=C6C=CC4=CC=CC(=C1CC2=CC=CC3=C2[N](C=C3)[H])C4=C56)CC7=CC=CC8=C7C=CC=C8",
            [
                {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 24, 25},
                {15, 16, 17, 18, 19, 20, 21, 22, 23},
                {27, 28, 29, 30, 31, 32, 33, 34, 35, 36},
            ],
        ),
        (
            "c1ccc2c(Cc3ccc4cc(Cc5cccc6cc[nH]c56)c5cccc6ccc3c4c65)cccc2c1",
            [
                {0, 1, 2, 3, 4, 32, 33, 34, 35, 36},
                {6, 7, 8, 9, 10, 11, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31},
                {13, 14, 15, 16, 17, 18, 19, 20, 21},
            ],
        ),
        (
            "C1=C3C=CC(=C4C=CC2=CC=CC(=C1)C2=C34)CC5=CC=CC6=C5C=CC=C6",
            [
                {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15},
                {17, 18, 19, 20, 21, 22, 23, 24, 25, 26},
            ],
        ),
        ("CC(=O)Nc1ccc(O)cc1", [{4, 5, 6, 7, 9, 10}]),
        (
            "O=C(NC(C(=O)NCC(=O)O)CSC1=C(N(C(C(=N1)NC2=C(C=C(C=C2Cl)OC(F)F)Cl)=O)C(C3CC3)C)O)CCC(C(=O)O)N",
            [{13, 14, 15, 16, 17, 18}, {20, 21, 22, 23, 24, 25}],
        ),
        (
            "COc1cc(ccc1O)c1oc2cc(O)cc(c2c(=O)c1O)O",
            [{2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 15, 16, 17, 18, 20}],
        ),
        (
            "Nc1nc(N)ncc1Cc(cc2OC)cc(OC)c2OC",
            [{1, 2, 3, 5, 6, 7}, {9, 10, 11, 14, 15, 18}],
        ),
        (
            "NC1=C(C=NC(=N1)N)CC2=CC(=C(C(=C2NC)OC)OC)OC",
            [{1, 2, 3, 4, 5, 6}, {9, 10, 11, 12, 13, 14}],
        ),
    ],
)
def test_perception(smi, expected):
    mol = MolFromSmiles(smi)
    aromatic_systems = AromaticSystems().systems(mol)
    assert aromatic_systems == expected
