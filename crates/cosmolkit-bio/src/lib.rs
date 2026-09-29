//! Detached structural-biology value and algorithm boundaries.

mod hierarchy;
mod metadata;
mod protein;
mod relationships;
mod residue;
mod secondary_structure;
mod source_ids;
mod structure_metadata;

pub use hierarchy::{
    AltLocRequest, BioAltLocGroupId, BioAssembly, BioAssemblyGenerator, BioAssemblyId,
    BioAssemblyOperator, BioAssemblySpecialKind, BioAtomId, BioAtomRow, BioCalcFlag, BioChainId,
    BioChainRow, BioCoordinateBlock, BioCoordinateFormat, BioCrystalCell, BioCrystalInfo,
    BioEntityDbRef, BioEntityId, BioEntityRow, BioModelId, BioModelRow, BioNcsOperator,
    BioNearestImage, BioResidueId, BioResidueRow, BioRowSpan, BioSiftsUnpResidue, BioStructureData,
    BioStructureError, BioStructureParts, BioTransform, ChainKind, EntityKind, PolymerKind,
    ResidueKind, altloc_matches, find_nearest_image, is_same_conformer, set_crystal_cell,
    set_crystal_fractional_transform, set_crystal_space_group_hm, set_crystal_z_pdb_if_nonempty,
    setup_cell_images,
};
pub use metadata::{
    BioBasicRefinementInfo, BioDiffractionInfo, BioExperimentInfo, BioExperimentalCrystalInfo,
    BioMetadata, BioRefinementInfo, BioRefinementRestraint, BioReflectionsInfo,
    BioSoftwareClassification, BioSoftwareItem, BioTlsGroup, BioTlsSelection,
};
pub use protein::{
    ProteinAtomIter, ProteinAtomRef, ProteinChainIter, ProteinChainRef, ProteinData,
    ProteinProjectionError, ProteinResidueIter, ProteinResidueRef, ProteinSelectionSummary,
    protein_atoms, protein_chain, protein_chains, protein_residues, protein_selection_summary,
    validate_protein_structure,
};
pub use relationships::{
    AtomAddress, BioAsu, BioCisPep, BioConnection, BioConnectionKind, BioModRes, ResidueAddress,
};
pub use residue::{
    ResidueCode, ResidueInfo, ResidueInfoKind, ResidueSequenceError,
    UNKNOWN_TABULATED_RESIDUE_INDEX, expand_one_letter, expand_one_letter_sequence,
    find_residue_info, find_residue_info_index, residue_code, residue_info, residue_info_checked,
};
pub use secondary_structure::{BioHelix, BioHelixClass, BioSheet, BioStrand};
pub use source_ids::{
    AltLocLabel, AtomName, AtomSourceIds, ChainSourceIds, EntitySourceIds, PdbAtomSerial,
    PdbChainId, PdbSeqId, ResidueName, ResidueSourceIds,
};
pub use structure_metadata::BioStructureSourceState;

use cosmolkit_model::TopologyBlock;

/// Translate coordinate rows only; metadata and anisotropic tensors are unchanged.
pub fn translate_coordinates(coordinates: &mut BioCoordinateBlock, offset: [f64; 3]) {
    // Project-defined coordinate-only translation. No full Gemmi structure
    // transformation or metadata rewrite is claimed. One pass, O(n), no allocation.
    for position in coordinates.positions_mut() {
        for axis in 0..3 {
            position[axis] += offset[axis];
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BioError {
    Unsupported,
}

pub fn select_residues(topology: &TopologyBlock, query: &str) -> Result<TopologyBlock, BioError> {
    let _ = (topology, query);
    Err(BioError::Unsupported)
}
