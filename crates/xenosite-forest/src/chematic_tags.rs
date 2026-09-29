//! Sync forest [`Tag`] with chematic caller tags (`1..=u16::MAX`).
//!
//! Chematic: `set_tag(..., None)` or `Some(0)` clears; `atom_tag` returns only
//! non-zero labels. Our [`Tag`] is always non-zero, so clear is `Option::None`.

use chematic::core::{AtomIdx, Molecule};

use crate::labels::Tag;

/// Push a sidecar label onto chematic (`None` clears; never writes `0`).
pub fn set_label(mol: &mut Molecule, idx: AtomIdx, tag: Option<Tag>) {
    mol.set_tag(idx, tag.map(Tag::get));
}

/// Read a chematic tag (`None` if missing or cleared).
pub fn get_label(mol: &Molecule, idx: AtomIdx) -> Option<Tag> {
    mol.atom_tag(idx).map(Tag::from_nonzero)
}
