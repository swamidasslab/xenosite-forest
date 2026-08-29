import pytest
from xenosite.metabolite import rules
from rdkit.Chem.rdmolfiles import MolFromSmiles, MolToSmiles
from rdkit.Chem.rdmolops import RemoveStereochemistry


examples = {
    "Glucuronidation": [
        (
            "AcetaminophenGlucuronide_SRT",
            "CC(=O)Nc1ccc(O)cc1",
            "CC(=O)NC1=CC=C(C=C1)OC2OC(C(C(O)C2O)O)C(=O)O",
        )
    ],
    "Epoxidation": [
        (
            "Amitriptyline_Epoxidation_Res_SRT",
            "c3cc2c(/C(c1c(cccc1)CC2)=C\\CCN(C)C)cc3",
            "C1=CC=C2C(=C1)CCC3=C(C2=CCCN(C)C)C=CC4C3O4",
        ),
        (
            "Amitriptyline_Epoxidation_SR",
            "c3cc2c(/C(c1c(cccc1)CC2)=C\\CCN(C)C)cc3",
            "C1=CC=C3C(=C1)CCC2=CC4C(C=C2C3=CCCN(C)C)O4",
        ),
        (
            "Amitriptyline_Epoxidation_SRT",
            "c3cc2c(/C(c1c(cccc1)CC2)=C\\CCN(C)C)cc3",
            "C1=CC=C2C(=C1)CCC3=C(C2=CCCN(C)C)C=CC4C3O4",
        ),
        ("EpoxideFormation_SRT", "C=Cc1ccccc1", "C1OC1c1ccccc1"),
        ("Example1_ResonanceStructure1_SRT", "c1c(cc(cc1)CC)C", "C1=C(C=C(C2C1O2)CC)C"),
        ("Example1_ResonanceStructure2_SRT", "c1c(cc(cc1)CC)C", "C12C(C=C(C=C1C)CC)O2"),
        (
            "Reaction3401_29363_Step1",
            "c1ccc2c(c1)ccc1c2cc2c(n1)ccc1c2cccc1",
            "C1=CC2=C(C=C1)C3=C(C=C2)N=C4C(=C3)C5=C(C=C4)C6C(C=C5)O6",
        ),
        (
            "Reaction3401_29363_Step1_aromatic_product",
            "c1ccc2c(c1)ccc1c2cc2c(n1)ccc1c2cccc1",
            "C1=CC2OC2c2ccc3nc4ccc5ccccc5c4cc3c21",
        ),
        (
            "Reaction3401_29363_Step1_kek_product",
            "c1ccc2c(c1)ccc1c2cc2c(n1)ccc1c2cccc1",
            "C1=CC=C2C(=C1)C=CC3=NC4=C(C=C23)C5=C(C=C4)C6C(C=C5)O6",
        ),
        ("Reaction49071_Epoxidation_SRT", "CC=C", "CC1CO1"),
        (
            "Reaction52978_SRT",
            "Clc1cc(Cl)c(cc1c1cc(Cl)c(cc1Cl)Cl)Cl",
            "ClC1=C(C=C(C2(C1O2)Cl)Cl)c1cc(Cl)c(cc1Cl)Cl",
        ),
        (
            "Reaction52978_SRT_site1",
            "Clc1cc(Cl)c(cc1c1cc(Cl)c(cc1Cl)Cl)Cl",
            "ClC1=C(Cl)C2OC2(Cl)C(c2cc(Cl)c(Cl)cc2Cl)=C1",
        ),
        (
            "Reaction52978_SRT_site2",
            "Clc1cc(Cl)c(cc1c1cc(Cl)c(cc1Cl)Cl)Cl",
            "ClC1=CC(=C(C2C1(Cl)O2)Cl)C3=C(C=C(C(=C3)Cl)Cl)Cl",
        ),
    ],
    "Glutathionation": [
        (
            "Amitriptyline_Epoxide_GSH_Conjugation",
            "C1=CC=C2C(=C1)CCC3=C(C2=CCCN(C)C)C=CC4C3O4",
            "C1=CC=C2C(=C1)CCC3=C(C2=CCCN(C)C)C=CC(C3O)SC[C@H](NC(CC[C@@H](C(=O)O)N)=O)C(=O)NCC(=O)O",
        ),
        (
            "ButadieneReactiveMetabolite",
            "C=CC(=O)CO",
            "C(CC(=O)N[C@@H](CSCCC(=O)CO)C(=O)NCC(=O)O)[C@@H](C(=O)O)N",
        ),
        # (
        #     "NotButadieneReactiveMetaboliteSoWillNotConjugateToGSH",
        #     "CC=CC(=O)CO",
        #     "C(CC(=O)N[C@@H](CSCCC(=O)CO)C(=O)NCC(=O)O)[C@@H](C(=O)O)N",
        # ),
    ],
    "Hydrogenation": [
        (
            "BigMolHydrogenation_SRT",
            "C1=CC=C6C(=C1)C2C(C(=NC(=C2)C(=O)N3CCCC3)C4=CC=C5C(=N4)C(C(=CC5=O)N)=O)N6",
            "C1=CC=C6C(=C1)C2C(C(=NC(=C2)C(=O)N3CCCC3)C4=CC=C5C(=N4)C(=C(C=C5O)N)O)N6",
        ),
        ("NAPQI_Reduction_SR", "CC(=O)N=C1C=CC(=O)C=C1", "CC(=O)Nc1ccc(O)cc1"),
        ("NAPQI_Reduction_SRT", "CC(=O)N=C1C=CC(=O)C=C1", "CC(=O)Nc1ccc(O)cc1"),
        (
            "Reaction1224_SRT",
            "Fc1ccc(c(c1)Br)C1N=C(NC2=C1C(=O)NC2)c1nccs1",
            "Fc1ccc(c(c1)Br)C1NC(NC2=C1C(=O)NC2)c1nccs1",
        ),
        (
            "Reaction52711_SRT",
            "O=C1CC[C@]2(C(=C1)C(O)C[C@@H]1[C@]2(F)[C@@H](O)CC2=C1CCC2(C)C)C",
            "O=C1CCC2(C(C1)C(O)CC1C2(F)C(O)CC2=C1CCC2(C)C)C",
        ),
    ],
    "QuinoneFormation": [
        (
            "BigMolQuinoneFormation_SRT",
            "C1=CC=C6C(=C1)C2C(C(=NC(=C2)C(=O)N3CCCC3)C4=CC=C5C(=N4)C(=C(C=C5O)N)O)N6",
            "C1=CC=C6C(=C1)C2C(C(=NC(=C2)C(=O)N3CCCC3)C4=CC=C5C(=N4)C(C(=CC5=O)N)=O)N6",
        ),
        (
            "Delavirdine_QuinoneFormation_SR",
            "CC(C)Nc1cccnc1N2CCN(CC2)C(=O)c3[nH]c4ccc(N[S](C)(=O)=O)cc4c3",
            "CC(C)N=C1C=CC(N=C1N2CCN(CC2)C(=O)C3=CC4=C([N]3[H])C=CC(=C4)N[S](C)(=O)=O)=O",
        ),
        (
            "Delavirdine_SR",
            "CC(C)Nc1cccnc1N2CCN(CC2)C(=O)c3[nH]c4ccc(N[S](C)(=O)=O)cc4c3",
            "CC(C)N=C1C=CC(N=C1N2CCN(CC2)C(=O)C3=CC4=C([N]3[H])C=CC(=C4)N[S](C)(=O)=O)=O",
        ),
        ("Epoxidation_SR_Num", "CC(=O)Nc1ccc(O)cc1", "CC(=O)N=C1C=CC(=O)C=C1"),
        (
            "Hydroxydelavirdine_QuinoneFormation",
            "CC(C)Nc1ccc(O)nc1N1CCN(C(=O)c2cc3cc(NS(C)(=O)=O)ccc3[nH]2)CC1",
            "CC(C)N=C1C=CC(N=C1N2CCN(CC2)C(=O)C3=CC4=C([N]3[H])C=CC(=C4)N[S](C)(=O)=O)=O",
        ),
        (
            "LongRangeQuinone_SRT",
            "C12=C(C=CC(=C1)O)N([H])C(=[N]2)C",
            "O=C2C=CC1=NC(=NC1=C2)C",
        ),
        ("NAPQI_Formation_SRT", "CC(=O)Nc1ccc(O)cc1", "O=C\\1/C=C\\C(=N/C(=O)C)/C=C/1"),
        (
            "NAPQI_QuinoneFormation_Site_Debug_SR",
            "CC(=O)Nc1ccc(O)cc1",
            "CC(=O)N=C1C=CC(=O)C=C1",
        ),
        (
            "Raloxifene_SRT",
            "O=C(c1c3ccc(O)cc3sc1c2ccc(O)cc2)c5ccc(OCCN4CCCCC4)cc5",
            "O=C(C2=C1C=CC(=O)C=C1SC2=C3C=CC(=O)C=C3)C4=CC=C(C=C4)OCCN5CCCCC5",
        ),
        (
            "Reaction26565_SRT",
            "CCc1nn(c(=O)n1CCOc1ccccc1)CCCN1CCN(CC1)c1ccc(c(c1)Cl)O",
            "CCc1nn(c(=O)n1CCOc1ccccc1)CCCN1CC[N+](=C2C=CC(=O)C(=C2)Cl)CC1",
        ),
        (
            "Reaction545_SRT",
            "NC(=O)N1c2ccccc2C=Cc2c1cccc2",
            "O=c1ccc2c(c1)ccc1c(n2)cccc1",
        ),
        (
            "Reaction61536_SRT",
            "O=c1[nH]c2c(C)ccnc2n(c2c1cccn2)C1CC1",
            "C=c1ccnc2c1nc(=O)c1c(n2C2CC2)nccc1",
        ),
        (
            "Reaction74945_SRT",
            "O=C1Nc2ccc(c(c2[C@@](N1)(C#CC1CC1)C(F)(F)F)F)F",
            "O=C1N=C2C=CC(=O)C(=C2[C@@](N1)(C#CC1CC1)C(F)(F)F)F",
        ),
        ("Reaction7751_SRT", "CC(=O)Nc1ccc(cc1)O", "N=C1C=CC(=O)C=C1"),
        (
            "Reaction9370_SRT",
            "CN(CCOc1cc(C)c(cc1C(C)C)OC(=O)C)C",
            "CC(C1=CC(=O)C(=CC1=O)C)C",
        ),
        (
            "Reaction96622_SRT",
            "Cc1cc(c(c(c1Cc1ncc[nH]1)C)O)C(C)(C)C",
            "C=C1C=C(C(=O)C(=C1Cc1ncc[nH]1)C)C(C)(C)C",
        ),
        ("Reaction9676_SRT", "COc1ccc(c(c1)C(C)(C)C)O", "O=C1C=CC(=O)C(=C1)C(C)(C)C"),
        ("Toluene_SR1", "c1ccc(C)cc1", "C=C1C(C=CC=C1)=O"),
    ],
    "Dealkylation": [
        ("DealkyNextToSulfur_SR", "C1=C(C=CN=C1CC)C(=S)N", "C1=C(C=CN=C1CC)C(=S)O"),
        ("DealkylationEx1_SRT", "c1ccncc1", "N=CC=CC=C=O"),
        (
            "Dealkylation_Ex1_SRT",
            "OC(=O)CCc1ccc(cc1C(=O)N[C@@H](c1cc(C)cc(c1)C)CC(C)C)C(Oc1cc(C)ccc1C)O",
            "CC1=CC(C)=CC(=C1)C(CC(C)C)NC(=O)C1=CC(C=O)=CC=C1CCC(=O)O",
        ),
        (
            "Reaction3401_29363_Step2_aromatic_reactant",
            "C1=CC2OC2c2ccc3nc4ccc5ccccc5c4cc3c21",
            "[C@H]1(O)c2ccc3c(c2C=C[C@@H]1O)cc1c(n3)ccc2c1cccc2",
        ),
        (
            "Reaction6127_23201_Metabolite_Dealkylation",
            "Clc1ccc(cc1)C(c1ccccc1)C",
            "Clc1ccc(cc1)C(=O)c1ccccc1",
        ),
        ("Reaction8560_37327_SRT", "CCOP(=O)(C(N1C(=O)CNC1=O)C)OCC", "CC(=O)O"),
        (
            "Reaction94457_SRT",
            "OCc1c(F)cccc1N(C(=O)OC([n+]1cnn(c1)C[C@@]([C@H](c1scc(n1)c1ccc(cc1)C#N)C)(c1cc(F)ccc1F)O)C)C",
            "N#Cc1ccc(cc1)c1csc(n1)C(C(c1cc(F)ccc1F)(Cn1cncn1)O)C",
        ),
        ("RingOpening2_SRT", "C1CCNC1", "O=CCCCN"),
    ],
    "ReductiveDehalogenation": [
        ("DebugRedutiveDehalo", "BrCCCCCBr", "CCCCCBr"),
        (
            "Reaction9310_110536_SRT",
            "Br/C/1=C/CCC(Br)C(Br)CCC(C(CC1)Br)Br",
            "Br/C/1=C/CC/C=C(/Br)\\CCC(C(CC1)Br)Br",
        ),
    ],
    "Dehydration": [
        ("DehydrationEx1_SRT", "Oc1ccccc1", "c1ccccc1"),
        (
            "Reaction119294_SRT",
            "ONC(=O)/C=C/c1cccc(c1)S(=O)(=O)Nc1ccccc1",
            "NC(=O)C=Cc1cccc(c1)S(=O)(=O)Nc1ccccc1",
        ),
        # (
        #     "Reaction6127_23201_Dehydration",
        #     "Clc1ccc(cc1)[C@@](c1ccccc1)(O)C",
        #     "Clc1ccc(cc1)[C@@](c1ccccc1)C",
        # ),
        (
            "Reaction7389_11214_Dehydration_SRT",
            "O-CCN(c1cc2c(cc1F)c(=O)c(cn2C1CC1)C(=O)O)CCNCC#N",
            "CCN(c1cc2c(cc1F)c(=O)c(cn2C1CC1)C(=O)O)CCNCC#N",
        ),
        (
            "Reaction7389_11214_Step2",
            "OCCN(c1cc2c(cc1F)c(=O)c(cn2C1CC1)C(=O)O)CCNCC#N",
            "C=CN(c1cc2c(cc1F)c(=O)c(cn2C1CC1)C(=O)O)CCNCC#N",
        ),
    ],
    "Dehydrogenation": [
        ("DehydrogenationGeneralEx1_SRT", "OC=CC=CO", "O=C-C=C-C=O"),
        ("DehydrogenationGeneralEx2_SRT", "OC=CO", "O=CC=O"),
        ("Dehydrogenation_Ex1", "OC=CC=CC=CC=CN", "N=CC=CC=CC=CC=O"),
        (
            "Reaction3183_3757_Already_Hydroxylated1_SRT",
            "c1ccc2c(c1)Sc1c(N2CCC2CCCCN2C)cc(cc1)S(O)C",
            "c1ccc2c(c1)Sc1c(N2CCC2CCCCN2C)cc(cc1)S(=O)C",
        ),
        (
            "Reaction3183_3757_Already_Hydroxylated2_SRT",
            "c1ccc2c(c1)Sc1c(N2CCC2CCCCN2C)cc(cc1)S(O)(C)C",
            "c1ccc2c(c1)Sc1c(N2CCC2CCCCN2C)cc(cc1)S(=O)(C)C",
        ),
        ("Reaction5963_62119_SR", "ONCCc1ccc(c(c1)O)O", "O/N=C\\Cc1ccc(c(c1)O)O"),
        (
            "Reaction59840_SRT",
            "CC(CCC[C@H]([C@H]1CC[C@@H]2[C@]1(C)CC[C@H]1C2=CC[C@@H]2[C@]1(C)CC[C@@H](C2)O)C)C",
            "CC(CCCC(C1CCC2C1(C)CCC1C2=CC=C2C1(C)CCC(C2)O)C)C",
        ),
        (
            "Reaction7389_11214_Dehydrogenation",
            "CCN(c1cc2c(cc1F)c(=O)c(cn2C1CC1)C(=O)O)CCNCC#N",
            "C=CN(c1cc2n(cc(c(=O)c2cc1F)C(=O)O)C1CC1)CCNCC#N",
        ),
    ],
    "Hydroxylation": [
        (
            "Delavirdine_Hydroxylation",
            "CC(C)Nc1cccnc1N2CCN(CC2)C(=O)c3[nH]c4ccc(N[S](C)(=O)=O)cc4c3",
            "CC(C)Nc1ccc(O)nc1N1CCN(C(=O)c2cc3cc(NS(C)(=O)=O)ccc3[nH]2)CC1",
        ),
        (
            "Reaction54714_H_SRT",
            "O=C([C@@H]1CCCCN1)Nc1c(C)cccc1C",
            "O=C(C1CCCCN1)Nc1c(C)ccc(c1C)O",
        ),
        (
            "Reaction54714_SRT",
            "O=C(C1CCCCN1)Nc1c(C)cccc1C",
            "O=C(C1CCCCN1)Nc1c(C)ccc(c1C)O",
        ),
        (
            "Reaction8492_99558_SRT",
            "OCC(c1nnc2n1cc(cc2)c1ocnc1c1cc(F)ccc1F)C",
            "OCC(c1nnc2n1cc(cc2)c1ocnc1c1cc(F)ccc1F)(O)C",
        ),
        ("Reaction932_37343_SRT", "[CH2+]CCl", "OCCCl"),
        ("Reaction932_37343_SRT_2", "Cl[CH+]C", "CC(O)Cl"),
        (
            "Reaction_4051_70821_SRT",
            "O=N(=O)c1ccc(c(c1)Oc1ccccc1)NS(=O)(=O)C",
            "Oc1ccc(cc1)Oc1cc(ccc1NS(=O)(=O)C)N(=O)=O",
        ),
    ],
    "Dephosphorylation": [("Dephosphorylation_SRT1", "O=P(O)(O)OC", "CO")],
    "EpoxideOpening": [
        ("EpoxideOpeningToDiol", "C1=CC=CC2OC12", "C1=CC=CC(O)C1O"),
        (
            "Reaction77559_SRT",
            "Oc1ccc(cc1)C1OC1c1cc(O)cc(c1)O",
            "OC1=CC=C(C=C1)CC(O)C2=CC(=CC(=C2)O)O",
        ),
    ],
    "Hydrolysis": [("HydrolyisEx1_SRT", "C1=CC=CC=C1C(=O)OC(C)(C)C", "O=C(O)c1ccccc1")],
    "Tautomerization": [
        ("LongRangeTautomerizationTest1", "ClCC=CC=CC=CC=CO", "ClCCC=CC=CC=CC=O"),
        ("LongRangeTautomerizationTest2", "ClCCC=CC=CC=CC=O", "ClCC=CC=CC=CC=CO"),
        (
            "LongRangeTautomerizationTest2_Rule_Set",
            "ClCCC=CC=CC=CC=O",
            "ClCC=CC=CC=CC=CO",
        ),
        ("TautomerizationTest1", "O=C1CCCCC1", "OC1=CCCCC1"),
        ("TautomerizationTest1_Rule_Set", "O=C1CCCCC1", "OC1=CCCCC1"),
        ("TautomerizationTest2", "OC1=CCCCC1", "O=C1CCCCC1"),
        ("TautomerizationTest3", "O=C1C=CCCC1", "OC1=CC=CCC1"),
        ("TautomerizationTest4", "OC1=CC=CCC1", "O=C1C=CCCC1"),
        # ("TautomerizationTest5", "Oc1ccccc1", "O=C1CC=CC=C1"),
    ],
    "NitrogenOxidation": [
        ("NitrogenOxidation_Demo1_SRT", "CCCN", "CCCNO"),
        (
            "Reaction3183_3757_SulfurOxidation5_SRT",
            "c1ccc2c(c1)Sc1c(N2CCC2CCCCN2)cc(cc1)S(C)(C)C",
            "c1ccc2c(c1)Sc1c(N2CCC2CCCCN2O)cc(cc1)S(C)(C)C",
        ),
        (
            "Reaction3183_3757_SulfurOxidation6_SRT",
            "c1ccc2c(c1)Sc1c(N2CCC2CCCCN2)cc(cc1)S(C)(C)C",
            "CS(C)(C)C1C=CC2SC3=CC=CC=C3N(O)(CCC3CCCCN3)C=2C=1",
        ),
    ],
    "NitrogenReduction": [
        (
            "NitrogenReductionEx1_SRT",
            "CC[C@](c1cc2c3nc4cccc(c4cc3Cn2c(=O)c1CO)N(=O)=O)(C(=O)O)O",
            "CCC(O)(C(=O)O)c1cc2n(c(=O)c1CO)Cc1cc3c(N)cccc3nc1-2",
        ),
        (
            "Reaction1323_61994_SRT",
            "O=Nc1c(C)c(N(=O)=O)c(c(c1C)N(=O)=O)C(C)(C)C",
            "[O-][N+](=O)c1c(C)c(N)c(c(c1C(C)(C)C)[N+](=O)[O-])C",
        ),
        (
            "Reaction375_38982_Step1",
            "Fc1ccc2c(c1)onc2C1CCN(CC1)CCc1c(C)nc2n(c1=O)CCCC2",
            "FC1=CC(=C(C=C1)C(=N)C2CCN(CC2)CCC3=C(C)N=C4N(C3=O)CCCC4)O",
        ),
        (
            "Reaction85420_SRT",
            "CC[C@](c1cc2c3nc4cccc(c4cc3Cn2c(=O)c1CO)N(=O)=O)(C(=O)O)O",
            "CCC(c1cc2c3nc4cccc(c4cc3Cn2c(=O)c1CO)N)(C(=O)O)O",
        ),
        ("RingOpeningEx1", "c1nocc1", "OC=CC=N"),
    ],
    "BenzodioxoleReduction": [
        (
            "Reaction112612_SRT",
            "COc1cc(cc(c1OC)OC)[C@@H]1c2cc3OCOc3cc2C[C@@H]2[C@@H]1C(=O)OC2",
            "COc1cc(cc(c1OC)OC)C1C2C(=O)OCC2Cc2c1cc(O)c(c2)O",
        )
    ],
    "OxidativeDehalogenation": [
        (
            "Reaction1176_11623_SRT",
            "ClC(C(c1ccccc1Cl)c1ccc(cc1)Cl)Cl",
            "OC(=O)C(c1ccccc1Cl)c1ccc(cc1)Cl",
        ),
        (
            "Reaction1176_11623_SRT_DoubleBondOxygen",
            "ClC(C(c1ccccc1Cl)c1ccc(cc1)Cl)Cl",
            "ClC(C(c1ccccc1Cl)c1ccc(cc1)Cl)=O",
        ),
        (
            "Reaction1176_11623_SRT_SingleBondOxygen",
            "ClC(C(c1ccccc1Cl)c1ccc(cc1)Cl)Cl",
            "ClC(C(c1ccccc1Cl)c1ccc(cc1)Cl)O",
        ),
        ("Reaction2790_28049_SRT", "FCC(F)(F)F", "OC(=O)C(F)(F)F"),
        (
            "Reaction443_69616_SRT",
            "Brc1ccc(c(c1)Br)Oc1ccc(cc1Br)Br",
            "Brc1ccc(c(c1)Br)Oc1ccc(c(c1Br)Br)O",
        ),
    ],
    "SulfurOxidation": [
        (
            "Reaction3183_3757_1Step_SRT",
            "CSc1ccc2c(c1)N(CCC1CCCCN1C)c1c(S2)cccc1",
            "c1ccc2c(c1)Sc1c(N2CCC2CCCCN2C)cc(cc1)S(=O)C",
        ),
        (
            "Reaction3183_3757_SulfurOxidation2_SRT",
            "c1ccc2c(c1)Sc1c(N2CCC2CCCCN2C)cc(cc1)SC",
            "c1ccc2c(c1)Sc1c(N2CCC2CCCCN2C)cc(cc1)S(O)C",
        ),
        (
            "Reaction3183_3757_SulfurOxidation3_SRT",
            "c1ccc2c(c1)Sc1c(N2CCC2CCCCN2C)cc(cc1)S(C)C",
            "c1ccc2c(c1)Sc1c(N2CCC2CCCCN2C)cc(cc1)S(O)(C)C",
        ),
        (
            "Reaction3183_3757_SulfurOxidation4_SRT",
            "c1ccc2c(c1)Sc1c(N2CCC2CCCCN2C)cc(cc1)S(C)(C)C",
            "c1ccc2c(c1)Sc1c(N2CCC2CCCCN2C)cc(cc1)S(O)(C)(C)C",
        ),
        (
            "Reaction3183_3757_SulfurOxidation_SRT",
            "c1ccc2c(c1)Sc1c(N2CCC2CCCCN2C)cc(cc1)S",
            "c1ccc2c(c1)Sc1c(N2CCC2CCCCN2C)cc(cc1)SO",
        ),
    ],
    "OxygenReduction": [
        (
            "Reaction7389_11214_Reduction_SRT",
            "O=CCN(c1cc2c(cc1F)c(=O)c(cn2C1CC1)C(=O)O)CCNCC#N",
            "O-CCN(c1cc2c(cc1F)c(=O)c(cn2C1CC1)C(=O)O)CCNCC#N",
        )
    ],
    "AzoSplitting": [
        (
            "Reaction7673_72241_SR",
            "OC(=O)c1cc(/N=N/c2ccc(c(c2)C(=O)O)O)ccc1O",
            "Nc1ccc(c(c1)C(=O)O)O",
        )
    ],
    "Sulfation": [
        (
            "Reaction_162_35034_SRT",
            "Clc1ccc(c(c1)Cl)C1=CC2(Cl)OC2C=C1Cl",
            "Clc1ccc(c(c1)Cl)c1cc(Cl)c(cc1Cl)S(=O)(=O)C",
        )
    ],
}


def can_smi(smi):
    m = MolFromSmiles(smi)
    if m is None:
        return ""
    RemoveStereochemistry(m)

    return MolToSmiles(m)  # type: ignore


def can_mol(m):
    if m is None:
        return ""
    RemoveStereochemistry(m)

    return MolToSmiles(m)  # type: ignore


def pytest_collection_modifyitems(items):
    for item in items:
        print(item)
        # check that we are altering a test named `test_xxx`
        # and it accepts the `value` arg
        if item.originalname == "test_xxx" and "value" in item.fixturenames:
            item._nodeid = item.nodeid.replace("]", "").replace("xxx[", "")


@pytest.mark.parametrize("rule", examples.keys())
def test_rule_exists(rule):
    rules.__dict__[rule]


@pytest.mark.parametrize(
    "rule, reactant, product",
    [(rule, e[1], e[2]) for rule, exs in examples.items() for e in exs],
    ids=[f"{rule.lower()}-{e[0]}" for rule, exs in examples.items() for e in exs],
)
def test_rule(rule, reactant, product):
    R = rules.__dict__[rule]()

    canonical_product = can_smi(product)
    r = MolFromSmiles(reactant)
    RemoveStereochemistry(r)

    for _, mols in R.metabolites(r):
        for p in mols:
            predicted = MolToSmiles(p)
            if predicted and canonical_product == can_smi(predicted):
                return
            
    assert False, f"Failed to find {product} in {reactant}"
