"""Pluggable harvest sources for XRM candidate terms."""

from .seed import harvest as harvest_seed
from .chebi import harvest as harvest_chebi
from .pubchem import harvest as harvest_pubchem
from .kegg import harvest as harvest_kegg
from .rhea import harvest as harvest_rhea
from .go import harvest as harvest_go
from .reactome import harvest as harvest_reactome

REGISTRY = {
    "seed": harvest_seed,
    "chebi": harvest_chebi,
    "pubchem": harvest_pubchem,
    "kegg": harvest_kegg,
    "rhea": harvest_rhea,
    "go": harvest_go,
    "reactome": harvest_reactome,
}

__all__ = ["REGISTRY"]
