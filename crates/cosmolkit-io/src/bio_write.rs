//! IO-owned options and typed failures for detached BIO coordinate mmCIF output.
//!
//! These options drive a **coordinate-only** serializer: the emitted document
//! contains the block name, `_entry.id`, `_atom_site` and (when atoms carry
//! anisotropic tensors) `_atom_site_anisotrop` — nothing else. Crystal,
//! symmetry/space-group, NCS, assembly, connection, cis-peptide, refinement
//! and all other source categories are preserved on input objects but are
//! never serialized here; no lossless roundtrip is claimed.

pub mod atoms;
pub mod tags;
pub mod value;

use std::fmt;

use cosmolkit_bio::BioStructureError;

use crate::cif::CifReadError;

/// The exact seven writer controls for BIO coordinate mmCIF output.
///
/// Source profile: pinned Gemmi `make_mmcif_document` invoked with
/// `MmcifOutputGroups(false)` then `atoms=block_name=entry=true`; only the
/// `group_pdb` and `auth_all` group bits remain user-selectable, plus the
/// five `cif::WriteOptions` layout controls. No other group switch exists in
/// this scope and no symmetry parameter is exposed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BioMmcifWriteParams {
    /// Emit the `_atom_site.group_PDB` column (ATOM/HETATM records).
    pub group_pdb: bool,
    /// Emit `_atom_site.auth_atom_id`/`auth_comp_id` in addition to label ids.
    pub auth_all: bool,
    /// Write single-row loops as tag/value pairs.
    pub prefer_pairs: bool,
    /// Omit blank lines between categories (blocks still separated).
    pub compact: bool,
    /// Put `#` empty comments before/after categories.
    pub misuse_hash: bool,
    /// Column at which pair values start (0 = single space).
    pub align_pairs: u16,
    /// Max loop column width used for value alignment (0 = no alignment).
    pub align_loops: u16,
}

impl Default for BioMmcifWriteParams {
    fn default() -> Self {
        // Gemmi✔️✔️: MmcifOutputGroups groups = MmcifOutputGroups(false);
        // Gemmi✔️✔️: groups.atoms = groups.block_name = groups.entry = true;
        // Gemmi✔️✔️: // group_pdb and auth_all stay caller-selected bits
        // Gemmi✔️✔️: struct WriteOptions {
        // Gemmi✔️✔️:   bool prefer_pairs = false;
        // Gemmi✔️✔️:   bool compact = false;
        // Gemmi✔️✔️:   bool misuse_hash = false;
        // Gemmi✔️✔️:   std::uint16_t align_pairs = 0;
        // Gemmi✔️✔️:   std::uint16_t align_loops = 0;
        // Gemmi✔️✔️: };
        // Behavior: fixed coordinate profile defaults; every group other than
        // atoms/block_name/entry is permanently off in this scope.
        // Complexity: constant-size construction, no allocation.
        Self {
            group_pdb: true,
            auth_all: false,
            prefer_pairs: false,
            compact: false,
            misuse_hash: false,
            align_pairs: 0,
            align_loops: 0,
        }
    }
}

/// A structural invariant, CIF mutation, invalid text or file IO failure.
#[derive(Debug)]
pub enum BioMmcifWriteError {
    Structure(BioStructureError),
    Cif(CifReadError),
    InvalidText { field: &'static str },
    Io(std::io::Error),
}

impl fmt::Display for BioMmcifWriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Structure(error) => write!(f, "BIO structure invalid for mmCIF writing: {error}"),
            Self::Cif(error) => write!(f, "mmCIF document error: {error}"),
            Self::InvalidText { field } => write!(f, "invalid UTF-8 in {field}"),
            Self::Io(error) => write!(f, "mmCIF IO error: {error}"),
        }
    }
}

impl std::error::Error for BioMmcifWriteError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Structure(error) => Some(error),
            Self::Cif(error) => Some(error),
            Self::Io(error) => Some(error),
            Self::InvalidText { .. } => None,
        }
    }
}

impl From<BioStructureError> for BioMmcifWriteError {
    fn from(error: BioStructureError) -> Self {
        Self::Structure(error)
    }
}
impl From<CifReadError> for BioMmcifWriteError {
    fn from(error: CifReadError) -> Self {
        Self::Cif(error)
    }
}
impl From<std::io::Error> for BioMmcifWriteError {
    fn from(error: std::io::Error) -> Self {
        Self::Io(error)
    }
}

#[cfg(test)]
mod scope_tests {
    use super::{BioMmcifWriteError as Error, BioMmcifWriteParams as Params};
    use crate::cif::{CifCheckLevel, read_cif_document};
    use cosmolkit_bio::BioStructureError;
    use std::error::Error as _;

    #[test]
    fn bio_pdbscope_scope_exact_seven_defaults_and_independence() {
        let defaults = Params::default();
        assert!(defaults.group_pdb);
        assert!(!defaults.auth_all);
        assert!(!defaults.prefer_pairs);
        assert!(!defaults.compact);
        assert!(!defaults.misuse_hash);
        assert_eq!(defaults.align_pairs, 0);
        assert_eq!(defaults.align_loops, 0);
        // Exactly seven controls exist; construction names every field.
        let all_set = Params {
            group_pdb: false,
            auth_all: true,
            prefer_pairs: true,
            compact: true,
            misuse_hash: true,
            align_pairs: 33,
            align_loops: 30,
        };
        assert_ne!(all_set, defaults);
        // Each control flips independently of the others.
        let mut one = defaults;
        one.group_pdb = false;
        assert_eq!(
            (one.group_pdb, one.auth_all, one.prefer_pairs),
            (false, false, false)
        );
        let mut two = defaults;
        two.auth_all = true;
        assert_eq!(
            (two.group_pdb, two.auth_all, two.align_pairs),
            (true, true, 0)
        );
        let mut three = defaults;
        three.align_loops = 30;
        assert_eq!(
            (three.compact, three.misuse_hash, three.align_loops),
            (false, false, 30)
        );
    }

    #[test]
    fn bio_pdbscope_scope_error_chains_retained() {
        let io = Error::from(std::io::Error::new(
            std::io::ErrorKind::PermissionDenied,
            "source path",
        ));
        assert_eq!(io.source().unwrap().to_string(), "source path");
        assert!(matches!(io, Error::Io(_)));
        let cif = Error::from(
            read_cif_document("invalid", "original.cif", CifCheckLevel::Syntax).unwrap_err(),
        );
        assert!(cif.source().unwrap().to_string().contains("original.cif"));
        assert!(matches!(cif, Error::Cif(_)));
        let structure = Error::from(BioStructureError::RowIndexTooLarge { value: usize::MAX });
        assert!(structure.source().is_some());
        assert!(matches!(structure, Error::Structure(_)));
        assert!(Error::InvalidText { field: "metadata" }.source().is_none());
    }
}
