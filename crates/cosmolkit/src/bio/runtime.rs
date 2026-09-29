//! Private storage and generated transactions; sibling operation bodies see only access structs.
use crate::*;
use cosmolkit_bio::BioStructureData;

#[derive(Debug, Clone, PartialEq)]
pub struct BioStructure {
    data: BioStructureData,
}
/// An amino-acid-only projection with one owned, already filtered [`BioStructure`].
///
/// `as_bio_structure` borrows that structure; `into_bio_structure` consumes
/// the wrapper. Neither operation restores non-protein rows removed when the
/// projection was created.
///
/// ```
/// use cosmolkit::{BioStructure, Protein};
/// let pdb = "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \n";
/// let protein = Protein::from_pdb(pdb).unwrap();
/// let view: &BioStructure = protein.as_bio_structure();
/// assert_eq!(view.num_atoms(), protein.num_atoms());
/// assert_eq!(protein.selection_summary().atoms, protein.num_atoms());
/// let owned: BioStructure = protein.into_bio_structure();
/// assert_eq!(owned.num_atoms(), 1);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Protein {
    structure: BioStructure,
}

impl BioStructure {
    /// Validate and construct a structure from detached BIO parts.
    pub fn from_parts(parts: BioStructureParts) -> Result<Self, BioStructureError> {
        Ok(Self {
            data: BioStructureData::from_parts(parts)?,
        })
    }
    pub fn validate_parts(parts: &BioStructureParts) -> Result<(), BioStructureError> {
        BioStructureData::validate_parts(parts)
    }
    pub fn validate(&self) -> Result<(), BioStructureError> {
        self.data.validate()
    }
    pub fn into_parts(self) -> BioStructureParts {
        self.data.into_parts()
    }
    pub fn input_format(&self) -> BioCoordinateFormat {
        self.data.input_format()
    }
    /// Read structural text with the IO owner's explicit dispatch and source context.
    pub fn from_text_with_params(text: &str, params: &BioReadParams) -> Result<Self, BioReadError> {
        Ok(Self {
            data: cosmolkit_io::read_bio_structure(text, params)?,
        })
    }
    /// Read structural text using the default content-detection parameters.
    pub fn from_text(text: &str) -> Result<Self, BioReadError> {
        Self::from_text_with_params(text, &BioReadParams::default())
    }
    /// Read a regular structural file using the selected coordinate format.
    pub fn read_with_format(
        path: &std::path::Path,
        format: BioCoordinateFormat,
    ) -> Result<Self, BioReadError> {
        Ok(Self {
            data: cosmolkit_io::read_bio_structure_file(path, format)?,
        })
    }
    /// Read a regular structural file, selecting its format from its extension.
    pub fn read(path: &std::path::Path) -> Result<Self, BioReadError> {
        Self::read_with_format(path, BioCoordinateFormat::Unknown)
    }
    /// Experimental Gemmi-aligned PDB text reader; source scope and limits are inherited.
    pub fn from_pdb(text: &str) -> Result<Self, BioPdbReadError> {
        Self::from_pdb_with_params(text, &BioPdbReadParams::default())
    }
    pub fn from_pdb_with_params(
        text: &str,
        params: &BioPdbReadParams,
    ) -> Result<Self, BioPdbReadError> {
        Ok(Self {
            data: cosmolkit_io::read_pdb_bio_structure(text, "<string>", params)?,
        })
    }
    /// Experimental coordinate mmCIF reader; no chemical Molecule conversion.
    pub fn from_mmcif(text: &str) -> Result<Self, BioMmcifReadError> {
        Ok(Self {
            data: cosmolkit_io::read_mmcif_bio_structure(text, "<string>")?,
        })
    }
    pub fn protein(&self) -> Result<Protein, ProteinProjectionError> {
        Ok(Protein {
            structure: BioStructure {
                data: self.data.protein()?.into_structure(),
            },
        })
    }
    fn operation_data_mut(&mut self) -> &mut BioStructureData {
        &mut self.data
    }
    fn validate_operation(&self) -> Result<(), BioOperationError> {
        self.data.validate().map_err(BioOperationError::Structure)
    }
    pub fn models(&self) -> &[cosmolkit_bio::BioModelRow] {
        self.data.models()
    }
    pub fn num_models(&self) -> usize {
        self.data.models().len()
    }
    pub fn chains(&self) -> &[cosmolkit_bio::BioChainRow] {
        self.data.chains()
    }
    pub fn num_chains(&self) -> usize {
        self.data.chains().len()
    }
    pub fn residues(&self) -> &[cosmolkit_bio::BioResidueRow] {
        self.data.residues()
    }
    pub fn num_residues(&self) -> usize {
        self.data.residues().len()
    }
    pub fn atoms(&self) -> &[cosmolkit_bio::BioAtomRow] {
        self.data.atoms()
    }
    pub fn atom_position(&self, atom: BioAtomId) -> Option<[f64; 3]> {
        self.data.atom_position(atom)
    }
    pub fn residue_atoms(&self, residue: BioResidueId) -> Option<&[BioAtomRow]> {
        self.data.residue_atoms(residue)
    }
    pub fn num_atoms(&self) -> usize {
        self.data.atoms().len()
    }
    pub fn entities(&self) -> &[cosmolkit_bio::BioEntityRow] {
        self.data.entities()
    }
    pub fn num_entities(&self) -> usize {
        self.data.entities().len()
    }
    pub fn connections(&self) -> &[cosmolkit_bio::BioConnection] {
        self.data.connections()
    }
    pub fn cispeps(&self) -> &[cosmolkit_bio::BioCisPep] {
        self.data.cispeps()
    }
    pub fn mod_residues(&self) -> &[cosmolkit_bio::BioModRes] {
        self.data.mod_residues()
    }
    pub fn helices(&self) -> &[cosmolkit_bio::BioHelix] {
        self.data.helices()
    }
    pub fn sheets(&self) -> &[cosmolkit_bio::BioSheet] {
        self.data.sheets()
    }
    pub fn metadata(&self) -> &cosmolkit_bio::BioMetadata {
        self.data.metadata()
    }
    pub fn source_state(&self) -> &cosmolkit_bio::BioStructureSourceState {
        self.data.source_state()
    }
    pub fn name(&self) -> &str {
        &self.data.source_state().name
    }
    pub fn has_origx(&self) -> bool {
        self.data.source_state().has_origx
    }
    pub fn origx(&self) -> &BioTransform {
        &self.data.source_state().origx
    }
    /// Borrow the last identity NCS operation ID retained in source metadata.
    pub fn ncs_oper_identity_id(&self) -> Option<&str> {
        self.data.ncs_oper_identity_id()
    }
    pub fn resolution(&self) -> f64 {
        self.data.source_state().resolution
    }
    pub fn ter_status(&self) -> u8 {
        self.data.source_state().ter_status
    }
    pub fn coordinates(&self) -> &cosmolkit_bio::BioCoordinateBlock {
        self.data.coordinates()
    }
    pub fn crystal(&self) -> Option<&cosmolkit_bio::BioCrystalInfo> {
        self.data.crystal()
    }
    pub fn ncs_operators(&self) -> &[cosmolkit_bio::BioNcsOperator] {
        self.data.ncs_operators()
    }
    pub fn assemblies(&self) -> &[cosmolkit_bio::BioAssembly] {
        self.data.assemblies()
    }
    pub fn find_entity(&self, source_id: &str) -> Option<(BioEntityId, &BioEntityRow)> {
        self.data.find_entity(source_id)
    }
    pub fn find_entity_of_subchain(&self, subchain: &str) -> Option<(BioEntityId, &BioEntityRow)> {
        self.data.find_entity_of_subchain(subchain)
    }
    pub fn find_atom(
        &self,
        residue_id: BioResidueId,
        name: AtomName,
        request: AltLocRequest,
        element: Option<Element>,
    ) -> Option<(BioAtomId, &BioAtomRow)> {
        self.data.find_atom(residue_id, name, request, element)
    }
    pub fn atom_by_altloc(
        &self,
        residue_id: BioResidueId,
        name: AtomName,
        altloc: Option<AltLocLabel>,
    ) -> Result<(BioAtomId, &BioAtomRow), BioStructureError> {
        self.data.atom_by_altloc(residue_id, name, altloc)
    }
}

impl Protein {
    /// Return the original input format retained by the underlying structure.
    ///
    /// This reads source metadata; it does not infer a format from the coordinates.
    pub fn input_format(&self) -> BioCoordinateFormat {
        self.structure.input_format()
    }

    /// Borrow the already filtered structural hierarchy without cloning it.
    ///
    /// The view is immutable; it cannot become a mutable structure reference.
    ///
    /// ```compile_fail
    /// use cosmolkit::{BioStructure, Protein};
    /// fn invalid(protein: &mut Protein) {
    ///     let _: &mut BioStructure = protein.as_bio_structure();
    /// }
    /// ```
    pub fn as_bio_structure(&self) -> &BioStructure {
        &self.structure
    }

    /// Consume this protein and return its already filtered structural hierarchy.
    pub fn into_bio_structure(self) -> BioStructure {
        self.structure
    }

    /// Read a structure and project its amino-acid hierarchy.
    pub fn from_text_with_params(
        text: &str,
        params: &BioReadParams,
    ) -> Result<Self, ProteinReadError> {
        BioStructure::from_text_with_params(text, params)
            .map_err(ProteinReadError::Structure)?
            .protein()
            .map_err(ProteinReadError::Projection)
    }
    /// Detect structural content and retain only its amino-acid hierarchy.
    pub fn from_text(text: &str) -> Result<Self, ProteinReadError> {
        BioStructure::from_text(text)
            .map_err(ProteinReadError::Structure)?
            .protein()
            .map_err(ProteinReadError::Projection)
    }
    /// Read a structural file with explicit format and project its amino acids.
    pub fn read_with_format(
        path: &std::path::Path,
        format: BioCoordinateFormat,
    ) -> Result<Self, ProteinReadError> {
        BioStructure::read_with_format(path, format)
            .map_err(ProteinReadError::Structure)?
            .protein()
            .map_err(ProteinReadError::Projection)
    }
    /// Read a structural file by extension and project its amino acids.
    pub fn read(path: &std::path::Path) -> Result<Self, ProteinReadError> {
        BioStructure::read(path)
            .map_err(ProteinReadError::Structure)?
            .protein()
            .map_err(ProteinReadError::Projection)
    }
    pub fn from_pdb(text: &str) -> Result<Self, ProteinReadError> {
        Self::from_pdb_with_params(text, &BioPdbReadParams::default())
    }
    pub fn from_pdb_with_params(
        text: &str,
        params: &BioPdbReadParams,
    ) -> Result<Self, ProteinReadError> {
        BioStructure::from_pdb_with_params(text, params)
            .map_err(ProteinReadError::Pdb)?
            .protein()
            .map_err(ProteinReadError::Projection)
    }
    pub fn from_mmcif(text: &str) -> Result<Self, ProteinReadError> {
        BioStructure::from_mmcif(text)
            .map_err(ProteinReadError::Mmcif)?
            .protein()
            .map_err(ProteinReadError::Projection)
    }
    fn operation_data_mut(&mut self) -> &mut BioStructureData {
        &mut self.structure.data
    }
    fn validate_operation(&self) -> Result<(), BioOperationError> {
        cosmolkit_bio::validate_protein_structure(&self.structure.data)
            .map_err(BioOperationError::Protein)
    }
    pub fn num_models(&self) -> usize {
        self.structure.models().len()
    }
    pub fn num_chains(&self) -> usize {
        self.structure.chains().len()
    }
    pub fn num_residues(&self) -> usize {
        self.structure.residues().len()
    }
    pub fn num_atoms(&self) -> usize {
        self.structure.atoms().len()
    }
    pub fn selection_summary(&self) -> ProteinSelectionSummary {
        cosmolkit_bio::protein_selection_summary(&self.structure.data)
    }
    pub fn chains(&self) -> ProteinChainIter<'_> {
        cosmolkit_bio::protein_chains(&self.structure.data)
    }
    pub fn residues(&self) -> ProteinResidueIter<'_> {
        cosmolkit_bio::protein_residues(&self.structure.data)
    }
    pub fn atoms(&self) -> ProteinAtomIter<'_> {
        cosmolkit_bio::protein_atoms(&self.structure.data)
    }
    pub fn chain(&self, index: usize) -> Option<ProteinChainRef<'_>> {
        cosmolkit_bio::protein_chain(&self.structure.data, index)
    }
}

#[allow(dead_code)]
#[derive(Debug)]
pub(super) struct BioOperationSpec {
    pub id: &'static str,
    pub target: &'static str,
    pub value_method: &'static str,
    pub inplace_method: &'static str,
    pub read: &'static [&'static str],
    pub write: &'static [&'static str],
}
cosmolkit_macros::bio_structure_ops! {
    translate_coordinates(offset: [f64; 3]) {
        targets: [BioStructure, Protein],
        value: with_translated_coordinates,
        inplace: translate_,
        body: super::translate_coordinates_impl,
        access: TranslateCoordinatesAccess,
        read: [],
        write: [coordinates: cosmolkit_bio::BioCoordinateBlock],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    const ATOM: &str =
        "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \n";

    mod probes {
        use super::*;
        cosmolkit_macros::bio_structure_ops! {
            failure(mode: u8) {
                targets: [BioStructure, Protein],
                value: probe_failure,
                inplace: probe_failure_,
                body: super::fail_body,
                access: FailureAccess,
                read: [atoms: Vec<cosmolkit_bio::BioAtomRow>],
                write: [coordinates: cosmolkit_bio::BioCoordinateBlock],
            }
            protein_guard(rows: Vec<cosmolkit_bio::BioResidueRow>) {
                targets: [Protein],
                value: probe_residues,
                inplace: probe_residues_,
                body: super::replace_residues,
                access: ResidueAccess,
                read: [],
                write: [residues: Vec<cosmolkit_bio::BioResidueRow>],
            }
        }
    }
    fn fail_body(access: probes::FailureAccess<'_>, mode: u8) -> Result<(), BioOperationError> {
        assert_eq!(access.atoms.len(), access.coordinates.len());
        if mode == 0 {
            return Ok(());
        } // Declared writable does not mean mandatory writing.
        access.coordinates.positions_mut()[0][0] = 99.0;
        match mode {
            1 => Err(BioOperationError::Structure(
                BioStructureError::AtomNotFound,
            )),
            2 => {
                *access.coordinates = BioCoordinateBlock::default();
                Ok(())
            }
            _ => panic!("intentional rollback probe"),
        }
    }
    fn replace_residues(
        access: probes::ResidueAccess<'_>,
        rows: Vec<cosmolkit_bio::BioResidueRow>,
    ) -> Result<(), BioOperationError> {
        *access.residues = rows;
        Ok(())
    }

    #[test]
    fn bio_read_support_defaults_root_types_and_older_error_paths() {
        use std::error::Error;
        let params: crate::BioReadParams = Default::default();
        assert_eq!(params.format, BioCoordinateFormat::Unknown);
        assert_eq!(params.source_name, "<string>");
        let _: Option<crate::BioReadError> = None;
        let structure_error = crate::BioReadError::WrongFormat {
            source_name: "named.cif".into(),
            format: BioCoordinateFormat::ChemComp,
        };
        let projected = ProteinReadError::Structure(structure_error);
        assert!(projected.to_string().contains("named.cif"));
        assert!(
            projected
                .source()
                .unwrap()
                .downcast_ref::<crate::BioReadError>()
                .is_some()
        );
        assert!(matches!(
            Protein::from_pdb("ATOM\n"),
            Err(ProteinReadError::Pdb(_))
        ));
        assert!(matches!(
            Protein::from_mmcif("invalid cif"),
            Err(ProteinReadError::Mmcif(_))
        ));
    }

    #[test]
    fn bio_ops_both_objects_share_impl_and_only_copy_writable_blocks() {
        let source = BioStructure::from_pdb(ATOM).unwrap();
        let snapshot = source.clone();
        assert!(Arc::ptr_eq(
            &source.data.coordinates,
            &snapshot.data.coordinates
        ));
        let transformed = source.with_translated_coordinates([2.0, 3.0, 4.0]).unwrap();
        let mut inplace = source.clone();
        inplace.translate_([2.0, 3.0, 4.0]).unwrap();
        assert_eq!(transformed, inplace);
        assert_eq!(source, snapshot);
        assert_eq!(source.coordinates().positions(), &[[1.0, 2.0, 3.0]]);
        assert_eq!(transformed.coordinates().positions(), &[[3.0, 5.0, 7.0]]);
        assert!(Arc::ptr_eq(&source.data.atoms, &transformed.data.atoms));
        assert!(Arc::ptr_eq(
            &source.data.residues,
            &transformed.data.residues
        ));
        assert!(Arc::ptr_eq(
            &source.data.metadata,
            &transformed.data.metadata
        ));
        assert!(!Arc::ptr_eq(
            &source.data.coordinates,
            &transformed.data.coordinates
        ));
        let source = source.protein().unwrap();
        let snapshot = source.clone();
        let transformed = source.with_translated_coordinates([2.0, 3.0, 4.0]).unwrap();
        let mut inplace = source.clone();
        inplace.translate_([2.0, 3.0, 4.0]).unwrap();
        assert_eq!(inplace, transformed);
        assert_eq!(source, snapshot);
        assert_eq!(
            transformed.atoms().next().unwrap().position(),
            [3.0, 5.0, 7.0]
        );
        assert!(Arc::ptr_eq(
            &source.structure.data.atoms,
            &transformed.structure.data.atoms
        ));
        assert!(!Arc::ptr_eq(
            &source.structure.data.coordinates,
            &transformed.structure.data.coordinates
        ));
        assert_eq!(BIO_STRUCTURE_OPS.len(), 2);
        for spec in BIO_STRUCTURE_OPS {
            assert_eq!(spec.write, ["coordinates"]);
            assert!(spec.read.is_empty());
            for method in [spec.value_method, spec.inplace_method] {
                let id = format!("{}.{}", spec.target, method);
                let row = BINDING_CONTRACT
                    .iter()
                    .find(|row| row.semantic_id == id)
                    .unwrap();
                assert_eq!(row.callable.unwrap().operation_semantic_id, Some(method));
                assert_eq!(
                    row.callable.unwrap().state_model,
                    if method == spec.value_method {
                        StateModel::ValueReturning
                    } else {
                        StateModel::InPlace
                    }
                );
            }
        }
    }

    #[test]
    fn bio_ops_failure_validation_and_unwind_are_atomic_on_both_targets() {
        let source = BioStructure::from_pdb(ATOM).unwrap();
        let protein = source.protein().unwrap();
        assert_eq!(source.probe_failure(0).unwrap(), source);
        assert_eq!(protein.probe_failure(0).unwrap(), protein);
        for mode in [1, 2] {
            let mut value = source.clone();
            assert!(value.probe_failure_(mode).is_err());
            assert_eq!(value, source);
            assert!(Arc::ptr_eq(
                &value.data.coordinates,
                &source.data.coordinates
            ));
            let mut value = protein.clone();
            assert!(value.probe_failure_(mode).is_err());
            assert_eq!(value, protein);
        }
        let mut value = source.clone();
        assert!(
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| value.probe_failure_(3)))
                .is_err()
        );
        assert_eq!(value, source);
        let mut value = protein.clone();
        assert!(
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| value.probe_failure_(3)))
                .is_err()
        );
        assert_eq!(value, protein);
    }

    #[test]
    fn bio_ops_protein_cannot_commit_non_amino_acid_rows() {
        let source = BioStructure::from_pdb(ATOM).unwrap().protein().unwrap();
        let water = BioStructure::from_pdb(&ATOM.replace("ALA", "HOH")).unwrap();
        let mut value = source.clone();
        let error = value
            .probe_residues_(water.residues().to_vec())
            .unwrap_err();
        assert!(matches!(
            error,
            BioOperationError::Protein(ProteinProjectionError::NonAminoAcidResidue { index: 0 })
        ));
        assert_eq!(value, source);
    }

    #[test]
    fn protein_recovery_storage_single_wrapper_preserves_projection_and_cow() {
        let mixed = format!(
            "{ATOM}HETATM    2  O   HOH A   2       4.000   5.000   6.000  1.00 20.00           O  \n"
        );
        let source = BioStructure::from_pdb(&mixed).unwrap();
        let protein = source.protein().unwrap();
        assert_eq!(source.atoms().len(), 2);
        assert_eq!(protein.atoms().len(), 1);
        assert_eq!(protein.structure.atoms().len(), 1);
        let snapshot = protein.clone();
        assert!(Arc::ptr_eq(
            &protein.structure.data.atoms,
            &snapshot.structure.data.atoms
        ));
        let translated = protein
            .with_translated_coordinates([1.0, 0.0, 0.0])
            .unwrap();
        let mut inplace = protein.clone();
        inplace.translate_([1.0, 0.0, 0.0]).unwrap();
        assert_eq!(translated, inplace);
        assert_eq!(protein, snapshot);
        assert_eq!(
            translated.atoms().next().unwrap().position(),
            [2.0, 2.0, 3.0]
        );
        assert!(Arc::ptr_eq(
            &protein.structure.data.atoms,
            &translated.structure.data.atoms
        ));
        assert!(!Arc::ptr_eq(
            &protein.structure.data.coordinates,
            &translated.structure.data.coordinates
        ));
    }

    #[test]
    fn protein_recovery_storage_failed_finalization_and_unwind_preserve_receiver() {
        let source = BioStructure::from_pdb(ATOM).unwrap().protein().unwrap();
        let snapshot = source.clone();
        let water = BioStructure::from_pdb(&ATOM.replace("ALA", "HOH")).unwrap();
        let mut candidate = source.clone();
        assert!(matches!(
            candidate.probe_residues_(water.residues().to_vec()),
            Err(BioOperationError::Protein(
                ProteinProjectionError::NonAminoAcidResidue { index: 0 }
            ))
        ));
        assert_eq!(candidate, snapshot);
        assert!(Arc::ptr_eq(
            &candidate.structure.data.atoms,
            &snapshot.structure.data.atoms
        ));
        let mut candidate = source.clone();
        assert!(
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                candidate.probe_failure_(3)
            }))
            .is_err()
        );
        assert_eq!(candidate, snapshot);
    }
}
