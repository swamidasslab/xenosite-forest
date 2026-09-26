"""Hand-curated seed lexicon: reaction types + classic example pairs."""

from __future__ import annotations

from .base import candidate

# Classic structure-backed examples (PubChem CIDs noted where known).
SEED = [
    {
        "pref_label": "aromatic hydroxylation",
        "synonyms": ["aryl hydroxylation", "benzene hydroxylation"],
        "spines": ["chemist reaction type", "site aromaticity", "aromatic and conjugated-system impact"],
        "examples": [
            {
                "reactant_smiles": "c1ccccc1",
                "product_smiles": "Oc1ccccc1",
                "label": "benzene → phenol",
                "refs": [
                    {"system": "pubchem", "id": "CID:241"},
                    {"system": "pubchem", "id": "CID:996"},
                ],
            }
        ],
    },
    {
        "pref_label": "aliphatic hydroxylation",
        "synonyms": ["alkyl hydroxylation"],
        "spines": ["chemist reaction type", "site aromaticity"],
        "examples": [
            {
                "reactant_smiles": "CC",
                "product_smiles": "CCO",
                "label": "ethane → ethanol",
                "refs": [
                    {"system": "pubchem", "id": "CID:6324"},
                    {"system": "pubchem", "id": "CID:702"},
                ],
            },
            {
                "reactant_smiles": "C",
                "product_smiles": "CO",
                "label": "methane → methanol (KEGG R01142 motif)",
                "refs": [{"system": "kegg", "id": "R01142"}],
            },
        ],
    },
    {
        "pref_label": "benzylic hydroxylation",
        "synonyms": ["toluene side-chain hydroxylation"],
        "spines": ["chemist reaction type"],
        "examples": [
            {
                "reactant_smiles": "Cc1ccccc1",
                "product_smiles": "OCc1ccccc1",
                "label": "toluene → benzyl alcohol",
                "refs": [
                    {"system": "pubchem", "id": "CID:1140"},
                    {"system": "pubchem", "id": "CID:244"},
                ],
            }
        ],
    },
    {
        "pref_label": "epoxidation",
        "synonyms": ["alkene epoxidation", "oxirane formation"],
        "spines": ["chemist reaction type", "ring fate", "electrophile role"],
        "examples": [
            {
                "reactant_smiles": "C=C",
                "product_smiles": "C1CO1",
                "label": "ethene → oxirane",
                "refs": [
                    {"system": "pubchem", "id": "CID:6325"},
                    {"system": "pubchem", "id": "CID:6354"},
                ],
            }
        ],
    },
    {
        "pref_label": "O-demethylation",
        "synonyms": ["oxidative O-demethylation", "aryl methyl ether cleavage"],
        "spines": ["chemist reaction type", "oxygenation outcome", "metabolite cardinality"],
        "examples": [
            {
                "reactant_smiles": "COc1ccccc1",
                "product_smiles": "Oc1ccccc1",
                "label": "anisole → phenol",
                "refs": [
                    {"system": "pubchem", "id": "CID:7519"},
                    {"system": "pubchem", "id": "CID:996"},
                ],
            }
        ],
    },
    {
        "pref_label": "N-demethylation",
        "synonyms": ["oxidative N-demethylation"],
        "spines": ["chemist reaction type", "oxygenation outcome"],
        "examples": [
            {
                "reactant_smiles": "CN(C)C",
                "product_smiles": "CNC",
                "label": "trimethylamine motif → dimethylamine motif",
                "refs": [],
            }
        ],
    },
    {
        "pref_label": "alcohol oxidation",
        "synonyms": ["alcohol dehydrogenation", "alcohol to aldehyde"],
        "spines": ["chemist reaction type", "redox polarity"],
        "examples": [
            {
                "reactant_smiles": "CCO",
                "product_smiles": "CC=O",
                "label": "ethanol → acetaldehyde",
                "refs": [
                    {"system": "pubchem", "id": "CID:702"},
                    {"system": "pubchem", "id": "CID:177"},
                ],
            }
        ],
    },
    {
        "pref_label": "phenol ortho-hydroxylation",
        "synonyms": ["phenol 2-hydroxylation", "catechol formation from phenol"],
        "spines": ["chemist reaction type", "aromatic and conjugated-system impact"],
        "examples": [
            {
                "reactant_smiles": "Oc1ccccc1",
                "product_smiles": "Oc1ccccc1O",
                "label": "phenol → catechol (KEGG R00815)",
                "refs": [{"system": "kegg", "id": "R00815"}],
            }
        ],
    },
    {
        "pref_label": "glucuronidation",
        "synonyms": ["glucuronosylation", "UDP-glucuronosyl transfer"],
        "spines": ["chemist reaction type", "electrophile role"],
        "examples": [
            {
                "reactant_smiles": "C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc43)[C@@H]1CCC2=O",
                "product_smiles": None,
                "label": "estrone → estrone 3-O-glucuronide (ChEBI parent/child)",
                "refs": [
                    {"system": "chebi", "id": "CHEBI:17263"},
                    {"system": "chebi", "id": "CHEBI:28919"},
                    {"system": "reactome", "id": "R-HSA-156588"},
                ],
            }
        ],
    },
    {
        "pref_label": "glycine conjugation",
        "synonyms": ["hippurate formation", "benzoyl glycine conjugation"],
        "spines": ["chemist reaction type"],
        "examples": [
            {
                "reactant_smiles": "O=C(O)c1ccccc1",
                "product_smiles": "O=C(O)CNC(=O)c1ccccc1",
                "label": "benzoic acid → hippuric acid",
                "refs": [
                    {"system": "pubchem", "id": "CID:243"},
                    {"system": "pubchem", "id": "CID:464"},
                ],
            }
        ],
    },
    {
        "pref_label": "quinone-imine formation",
        "synonyms": ["quinone imine formation", "NAPQI-type formation"],
        "spines": ["chemist reaction type", "electrophile role", "aromatic and conjugated-system impact"],
        "examples": [
            {
                "reactant_smiles": "CC(=O)Nc1ccc(O)cc1",
                "product_smiles": "CC(=O)N=C1C=CC(=O)C=C1",
                "label": "acetaminophen → NAPQI (quinone-imine)",
                "refs": [{"system": "pubchem", "id": "CID:1983"}],
            }
        ],
    },
    {
        "pref_label": "sarcosine demethylation",
        "synonyms": ["N-methylglycine demethylation"],
        "spines": ["chemist reaction type"],
        "examples": [
            {
                "reactant_smiles": "CNCC(=O)O",
                "product_smiles": "NCC(=O)O",
                "label": "sarcosine → glycine (KEGG R00610)",
                "refs": [{"system": "kegg", "id": "R00610"}],
            }
        ],
    },
]


def harvest(_cfg=None):
    out = []
    for i, row in enumerate(SEED):
        out.append(
            candidate(
                candidate_id=f"cand:seed:{i}:{row['pref_label'].replace(' ', '_')}",
                pref_label=row["pref_label"],
                synonyms=row.get("synonyms"),
                examples=row.get("examples"),
                suggested_spines=row.get("spines"),
                sources=[{"system": "seed", "id": f"seed:{i}", "url": ""}],
                notes="Hand-curated seed example for XRM review",
            )
        )
    return out
