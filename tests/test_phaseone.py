from rdkit.Chem.rdmolfiles import MolFromSmiles, MolToSmiles
from xenosite.forest import rules, bfs
import pytest


def test_rule_modification():
    mol = MolFromSmiles("CN")
    OriginalRule = rules.Dealkylation()
    ModifiedRule = rules.Dealkylation(rxns=rules.Dealkylation.smarts[1:])

    assert len(list(ModifiedRule.metabolites(mol))) == 2
    assert len(list(OriginalRule.metabolites(mol))) == 3

    # OriginalRule = rules.Dealkylation()

    # ModifiedRule = rules.Dealkylation()

    # s = list(ModifiedRule.smarts)
    # s.remove(
    #     "[#6H3:1][#7,#8H0,#16:2]>>([*:2].[*:1](=O)O)"
    # )  # removes a carboxylic acid reaction
    # ModifiedRule.smarts = s

    # assert len(list(OriginalRule.metabolites(mol))) == 3



def cmp_pathway(pathway, target):
    smi, sites, path = pathway

    assert smi == target[0]
    assert sites == target[1]
    assert [MolToSmiles(m) for m in path] == [MolToSmiles(MolFromSmiles(s)) for s in target[2]]  # type: ignore


def test_dehydrogenation():
    result = next(bfs(["CCO", "C=CO"], phase1=True))

    cmp_pathway(
        result,
        ("C=CO", [("Dehydrogenation", {"1.h", "2.h"})], ["CCO", "C=CO"]),
    )


def test_dealkylation1():
    result = next(bfs(["CCN", "CCO"], phase1=True))

    cmp_pathway(
        result,
        ("CCO", [("Dealkylation", {"2.3"})], ["CCN", "CCO"]),
    )


def test_dehydration1():
    result = next(bfs(["CCO", "CC"], phase1=True))

    cmp_pathway(
        result,
        ("CC", [("Dehydration", {"2.3"})], ["CCO", "CC"]),
    )


def test_dehydrogenation1():
    result = next(bfs(["CCO", "C=CO"], phase1=True))

    cmp_pathway(
        result,
        ("C=CO", [("Dehydrogenation", {"2.h", "1.h"})], ["CCO", "C=CO"]),
    )


def test_dephosphorylation1():
    result = next(bfs(["COP(=O)(O)O", "CO"], phase1=True))

    cmp_pathway(
        result,
        ("CO", [("Dephosphorylation", {"2.3"})], ["COP(=O)(O)O", "CO"]),
    )


def test_epoxidation1():
    result = next(bfs(["C=C", "C1OC1"], phase1=True))

    cmp_pathway(
        result,
        ("C1CO1", [("Epoxidation", {"1.2"})], ["C=C", "C1CO1"]),
    )


def test_epoxidation_opening1():
    result = next(bfs(["C1OC1", "CCO"], phase1=True))

    cmp_pathway(
        result,
        ("CCO", [("EpoxideOpening", {"1.2"})], ["C1CO1", "CCO"]),
    )


def test_hydrogenation1():
    result = next(bfs(["CC=O", "CCO"], phase1=True))

    cmp_pathway(
        result,
        ("CCO", [("Hydrogenation", {"2.2", "3.3"})], ["CC=O", "CCO"]),
    )


def test_hydroxylation1():
    result = next(bfs(["CC", "CCO"], phase1=True))

    cmp_pathway(
        result,
        ("CCO", [("Hydroxylation", {"1.h"})], ["CC", "CCO"]),
    )


def test_hydrolysis1():
    result = next(bfs(["O=C(O)C", "CC=O"], phase1=True))

    cmp_pathway(
        result,
        ("CC=O", [("Hydrolysis", {"2.3"})], ["O=C(O)C", "CC=O"]),
    )


def test_nitrogen_oxidation1():
    result = next(bfs(["CCN", "CCNO"], phase1=True))

    cmp_pathway(
        result,
        ("CCNO", [("NitrogenOxidation", {"3.3"})], ["CCN", "CCNO"]),
    )


def test_dehydration2():
    result = next(bfs(["CCNO", "CCN"], phase1=True))

    cmp_pathway(
        result,
        ("CCN", [("Dehydration", {"3.4"})], ["CCNO", "CCN"]),
    )


@pytest.mark.xfail(reason="canonical SMILES mismatch for this nitro reduction path")
def test_nitrogen_reduction1():
    result = next(
        bfs(
            [
                "[O-]-[N+](C1=CC=C(O1)C=NN2C(=O)NC(=O)C2)=O",
                "N(C1=CC=C(O1)C=NN2C(=O)NC(=O)C2)=O",
            ],
            phase1=True,
        )
    )

    cmp_pathway(
        result,
        (
            "O=Nc1ccc(C=NN2CC(=O)NC2=O)o1",
            [("NitrogenReduction", {"1.2"})],
            [
                "O=C1CN(N=CC2=CC=C([N+](=O)[O-])O2)C(=O)N1",
                "O=Nc1ccc(C=NN2CC(=O)NC2=O)o1",
            ],
        ),
    )


def test_hydrogenation2():
    result = next(bfs(["CC=O", "CCO"], phase1=True))

    cmp_pathway(
        result,
        ("CCO", [("Hydrogenation", {"2.2", "3.3"})], ["CC=O", "CCO"]),
    )


def test_oxidative_dehalogenation():
    result = next(bfs(["CCCl", "CCO"], phase1=True))

    cmp_pathway(
        result,
        ("CCO", [("OxidativeDehalogenation", {"2.3"})], ["CCCl", "CCO"]),
    )


def test_reductive_dehalogenation():
    result = next(bfs(["CCCl", "CC"], phase1=True))

    cmp_pathway(
        result,
        ("CC", [("ReductiveDehalogenation", {"2.3"})], ["CCCl", "CC"]),
    )


def test_sulfur_oxidation():
    result = next(bfs(["CCS", "CCSO"], phase1=True))

    cmp_pathway(
        result,
        ("CCSO", [("SulfurOxidation", {"3.3"})], ["CCS", "CCSO"]),
    )


def test_sulfer_reduction():
    result = next(bfs(["CCSO", "CCS"], phase1=True))

    cmp_pathway(
        result,
        ("CCS", [("SulfurReduction", {"3.4"})], ["CCSO", "CCS"]),
    )
