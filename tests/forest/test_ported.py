"""Every target placement stays available to the search. See HEURISTICS.md."""

from xenosite.forest.find_path import atom_diff
from xenosite.forest.rdkitutil import as_mol, mcs_matches


def _rings(mol, embeddings):
    rings = set()
    for embedding in embeddings:
        rings.add(
            frozenset(
                index
                for index in embedding
                if mol.GetAtomWithIdx(index).GetIsAromatic()
            )
        )
    rings.discard(frozenset())
    return rings


def test_multi_mcs_union_would_prune_both_benzyls():
    """Two PhCHO placements. One core must not be the union of both rings."""

    reactant = as_mol("c1ccccc1CNCc1ccccc1")
    target = as_mol("O=Cc1ccccc1")
    embeddings = mcs_matches(reactant, target).embeddings
    rings = _rings(reactant, embeddings)
    assert len(rings) >= 2
    assert len({len(embedding) for embedding in embeddings}) == 1
    union = set().union(*rings)
    assert all(set(ring) != union for ring in rings)

    diff = atom_diff(reactant, target)
    mappings = getattr(diff, "mappings", (diff.mapping,))
    assert len(tuple(mappings)) >= 2
    for ring in rings:
        assert any(
            ring <= set(mapping) and not union <= set(mapping) for mapping in mappings
        )


def test_smaller_secondary_embedding_alone_misses_naphthaldehyde():
    """The phenyl remainder must not be the only embedding kept."""

    reactant = as_mol("c1ccc2ccccc2c1CNCc1ccccc1")
    target = as_mol("O=Cc1cccc2ccccc12")
    embeddings = mcs_matches(reactant, target).embeddings
    sizes = sorted({len(embedding) for embedding in embeddings})
    assert len(sizes) >= 2 and sizes[0] < sizes[-1], sizes
    assert any(len(embedding) == sizes[-1] for embedding in embeddings)
