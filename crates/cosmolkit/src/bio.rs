//! Public structural objects and their thin algorithm adapters.
//! See dev/bio_architecture.md for ownership and lightweight operation semantics.
mod runtime;
use crate::{
    BioMmcifReadError, BioPdbReadError, BioReadError, BioStructureError, ProteinProjectionError,
};
pub use runtime::{BioStructure, Protein};

/// Reading and explicit amino-acid projection failures retain their source.
#[derive(Debug)]
pub enum ProteinReadError {
    Structure(BioReadError),
    Pdb(BioPdbReadError),
    Mmcif(BioMmcifReadError),
    Projection(ProteinProjectionError),
}
impl std::fmt::Display for ProteinReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Structure(e) => write!(f, "{e}"),
            Self::Pdb(e) => write!(f, "{e}"),
            Self::Mmcif(e) => write!(f, "{e}"),
            Self::Projection(e) => write!(f, "{e}"),
        }
    }
}
impl std::error::Error for ProteinReadError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(match self {
            Self::Structure(e) => e,
            Self::Pdb(e) => e,
            Self::Mmcif(e) => e,
            Self::Projection(e) => e,
        })
    }
}

/// Final data validation failure; a failed operation never replaces its input.
#[derive(Debug)]
pub enum BioOperationError {
    Structure(BioStructureError),
    Protein(ProteinProjectionError),
}
impl std::fmt::Display for BioOperationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Structure(e) => write!(f, "{e}"),
            Self::Protein(e) => write!(f, "{e}"),
        }
    }
}
impl std::error::Error for BioOperationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(match self {
            Self::Structure(e) => e,
            Self::Protein(e) => e,
        })
    }
}

fn translate_coordinates_impl(
    access: runtime::TranslateCoordinatesAccess<'_>,
    offset: [f64; 3],
) -> Result<(), BioOperationError> {
    cosmolkit_bio::translate_coordinates(access.coordinates, offset);
    Ok(())
}

// Compile probes live beside real operation bodies, outside private storage.
#[cfg(cosmolkit_bio_privacy_probe)]
#[allow(dead_code, unused_variables)]
fn bio_privacy_probe(
    access: runtime::TranslateCoordinatesAccess<'_>,
    structure: BioStructure,
    protein: Protein,
) {
    let _ = access.coordinates.positions_mut();
    #[cfg(cosmolkit_bio_privacy_case = "undeclared")]
    let _ = access.atoms;
    #[cfg(cosmolkit_bio_privacy_case = "storage")]
    {
        let _ = structure.data;
        let _ = protein.structure;
    }
}
