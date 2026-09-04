from rdkit.Chem.rdmolfiles import MolFromSmiles, MolToSmiles
from xenosite.forest import rules, bfs, PhaseOneRS
from xenosite.forest.utils import canon_smi, unmapped_smiles, is_rdkit_valid, clean


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

    assert canon_smi(smi) == canon_smi(target[0])
    assert sites == target[1]
    assert canon_smi(path) == canon_smi(target[2])


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


HISTIDINE = "NC(Cc1cnc[nH]1)C(=O)O"


def _phase1_histidine_products():
    mol = MolFromSmiles(HISTIDINE)
    out = []
    for (rxn, site), mets in PhaseOneRS.metabolites(mol):
        for m in mets:
            out.append((rxn, site, unmapped_smiles(m), m))
    return out


def test_histidine_phase1_does_not_emit_invalid_metabolites():
    products = _phase1_histidine_products()
    assert products
    invalid = [
        (rxn, site, smi)
        for rxn, site, smi, mol in products
        if MolFromSmiles(smi) is None or not is_rdkit_valid(mol)
    ]
    assert invalid == []


def test_histidine_nitrogen_oxidation_is_chemically_valid():
    mol = MolFromSmiles(HISTIDINE)
    smiles = []
    for _, mets in rules.NitrogenOxidation().metabolize(mol):
        smiles.extend(canon_smi(m) for m in mets)

    # hydroxylamine and nitroso on the amino group; N-OH on pyrrole N;
    # N-oxide on the pyridine-like imidazole N. No pentavalent N.
    assert canon_smi("O=C(O)C(CC1=CN=CN1)NO") in smiles
    assert canon_smi("O=NC(CC1=CN=CN1)C(=O)O") in smiles
    assert canon_smi("NC(CC1=CN=CN1O)C(=O)O") in smiles
    assert canon_smi("NC(CC1=C[N+]([O-])=CN1)C(=O)O") in smiles
    assert canon_smi("NC(CC1=CN(O)=CN1)C(=O)O") not in smiles
    assert canon_smi("NC(CC1=CN=CN1=O)C(=O)O") not in smiles
    assert all(MolFromSmiles(s) is not None for s in smiles)


def test_clean_preserves_imidazole_nh():
    mol = MolFromSmiles(HISTIDINE)
    cleaned = clean(mol)
    assert len(cleaned) == 1
    assert MolFromSmiles(MolToSmiles(cleaned[0])) is not None

