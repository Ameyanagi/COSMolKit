use cosmolkit_core::{
    AtomSpec, BondOrder, BondSpec, Conformer3D, Element, Hybridization, MmffMolProperties,
    MmffPublicApiError, MmffVariant, MolecularForceField, MolecularForceFieldError, Molecule,
    MoleculeBuilder, UffPublicApiError, mmff_get_molecule_force_field,
    uff_get_molecule_force_field,
};

fn ethanol() -> Molecule {
    let mol = Molecule::from_smiles("CCO")
        .unwrap()
        .with_hydrogens()
        .unwrap();
    assert_eq!(mol.num_atoms(), 9);
    let mut builder = mol.to_builder();
    builder
        .add_conformer(Conformer3D::new(
            7,
            vec![
                [0., 0., 0.],
                [1.7, 0.1, 0.1],
                [2.2, 1.5, 0.3],
                [-0.5, 0.9, 0.2],
                [-0.4, -0.8, 0.7],
                [-0.4, -0.3, -0.9],
                [2.1, -0.8, 0.8],
                [2.1, -0.3, -0.9],
                [3.1, 1.4, 0.5],
            ],
            true,
        ))
        .unwrap();
    builder.build().unwrap()
}

fn fields(mol: &Molecule) -> Vec<MolecularForceField> {
    vec![
        mmff_get_molecule_force_field(mol, MmffVariant::Mmff94, f64::INFINITY, 7, true)
            .unwrap()
            .unwrap(),
        mmff_get_molecule_force_field(mol, MmffVariant::Mmff94s, f64::INFINITY, 7, true)
            .unwrap()
            .unwrap(),
        uff_get_molecule_force_field(mol, f64::INFINITY, 7, true)
            .unwrap()
            .unwrap(),
    ]
}

#[test]
fn owned_fields_move_safely_and_physical_gradients_match_energy_differences() {
    let mol = ethanol();
    let original = mol.clone();
    // Each opaque handle moves from its factory through Vec and Iterator before
    // evaluation; no builder-local or former container address may be retained.
    for mut field in fields(&mol) {
        let positions = field.positions();
        assert_eq!(positions, mol.conformers_3d()[0].coordinates());
        assert!(field.energy().is_finite());
        field.set_fixed_points(&[0, 8]).unwrap();
        let gradient = field.gradient();
        assert_eq!(gradient.len(), positions.len());
        assert!(gradient.iter().flatten().all(|value| value.is_finite()));
        assert!(
            gradient[0].iter().any(|value| value.abs() > 1e-5),
            "physical gradient must not be masked at a fixed atom"
        );
        let step = 1e-5;
        for atom in 0..positions.len() {
            for axis in 0..3 {
                let mut displaced = positions.clone();
                displaced[atom][axis] += step;
                field.set_positions(&displaced).unwrap();
                let plus = field.energy();
                displaced[atom][axis] -= 2. * step;
                field.set_positions(&displaced).unwrap();
                let minus = field.energy();
                let numeric = (plus - minus) / (2. * step);
                let analytical = gradient[atom][axis];
                assert!(
                    (numeric - analytical).abs() <= 1e-4 * (1. + analytical.abs()),
                    "gradient {atom}/{axis}: analytical={analytical}, finite difference={numeric}"
                );
            }
        }
        field.set_positions(&positions).unwrap();
        let before = field.energy();
        let status = field.minimize(60, 1e-4, 1e-6).unwrap();
        assert!(matches!(status, 0 | 1));
        let after = field.energy();
        assert!(
            after.is_finite() && after < before,
            "minimization must reduce energy"
        );
        let final_positions = field.positions();
        for atom in [0, 8] {
            for axis in 0..3 {
                assert_eq!(
                    positions[atom][axis].to_bits(),
                    final_positions[atom][axis].to_bits()
                );
            }
        }
        assert_eq!(field.fixed_points(), [0, 8]);
        assert_eq!(
            mol, original,
            "factories and relaxation preserve the source molecule"
        );
    }
}

#[test]
fn factories_require_complete_parameters_and_never_return_partial_fields() {
    let mut builder = MoleculeBuilder::new();
    builder.add_atom(AtomSpec::new(Element::H));
    builder.add_3d_conformer(vec![[0., 0., 0.]]).unwrap();
    let hydrogen = builder.build().unwrap();
    for variant in [MmffVariant::Mmff94, MmffVariant::Mmff94s] {
        assert!(
            mmff_get_molecule_force_field(&hydrogen, variant, 100., -1, true)
                .unwrap()
                .is_none()
        );
    }

    let mut builder = MoleculeBuilder::new();
    let c = builder.add_atom(AtomSpec::new(Element::C).with_hybridization(Hybridization::S));
    let f = builder.add_atom(AtomSpec::new(Element::F).with_hybridization(Hybridization::Sp3));
    builder
        .add_bond(BondSpec::new(c, f, BondOrder::Single))
        .unwrap();
    builder
        .add_3d_conformer(vec![[0., 0., 0.], [1.4, 0., 0.]])
        .unwrap();
    let missing_type = builder.build().unwrap();
    assert!(
        uff_get_molecule_force_field(&missing_type, 10., -1, true)
            .unwrap()
            .is_none()
    );
}

#[test]
fn repeated_mmff_preparation_preserves_aromatic_types_and_same_xyz_force_field() {
    for (smiles, aromatic) in [("c1ccccc1", true), ("c1ccc2ccccc2c1", true), ("CCO", false)] {
        let mol = Molecule::from_smiles(smiles)
            .unwrap()
            .with_hydrogens()
            .unwrap();
        // Use fixed finite, noncoincident XYZ so the two factory paths evaluate
        // identical geometry without another conformer-generation dependency.
        let coordinates: Vec<_> = (0..mol.num_atoms())
            .map(|i| {
                [
                    (i % 5) as f64 * 1.2,
                    (i / 5) as f64 * 1.3,
                    ((i * i) % 7) as f64 * 0.23,
                ]
            })
            .collect();
        let mut builder = mol.to_builder();
        builder.add_3d_conformer(coordinates.clone()).unwrap();
        let mol = builder.build().unwrap();
        let original = mol.clone();

        for variant in [MmffVariant::Mmff94, MmffVariant::Mmff94s] {
            let first = MmffMolProperties::new(&mol, variant.as_rdkit_str(), 0).unwrap();
            assert!(first.is_valid(), "{smiles}/{variant:?}");
            assert!(first.molecule.prop("_MMFFSanitized").is_some());
            let first_snapshot = first.molecule.clone();
            let repeated =
                MmffMolProperties::new(&first.molecule, variant.as_rdkit_str(), 0).unwrap();
            assert!(repeated.is_valid(), "repeated {smiles}/{variant:?}");
            assert_eq!(first.molecule, first_snapshot);
            assert_eq!(first.atom_properties, repeated.atom_properties);
            assert_eq!(first.aromaticity, repeated.aromaticity);
            assert_eq!(first.molecule.conformers_3d()[0].coordinates(), coordinates);
            assert_eq!(
                repeated.molecule.conformers_3d()[0].coordinates(),
                coordinates
            );

            for (idx, atom) in first.molecule.atoms().iter().enumerate() {
                let expected_aromatic = aromatic && atom.atomic_number() == 6;
                assert_eq!(atom.is_aromatic(), expected_aromatic);
                assert_eq!(
                    repeated.molecule.atoms()[idx].is_aromatic(),
                    expected_aromatic
                );
                assert_eq!(
                    first.get_mmff_partial_charge(idx).unwrap(),
                    repeated.get_mmff_partial_charge(idx).unwrap()
                );
                if aromatic {
                    let expected_type = if atom.atomic_number() == 6 { 37 } else { 5 };
                    assert_eq!(first.get_mmff_atom_type(idx).unwrap(), expected_type);
                    assert_eq!(repeated.get_mmff_atom_type(idx).unwrap(), expected_type);
                }
            }
            for (idx, bond) in first.molecule.bonds().iter().enumerate() {
                let expected_aromatic = aromatic
                    && first.molecule.atoms()[bond.begin().index()].atomic_number() == 6
                    && first.molecule.atoms()[bond.end().index()].atomic_number() == 6;
                assert_eq!(bond.is_aromatic(), expected_aromatic);
                assert_eq!(
                    repeated.molecule.bonds()[idx].is_aromatic(),
                    expected_aromatic
                );
                assert_eq!(bond.order(), repeated.molecule.bonds()[idx].order());
            }
            if aromatic {
                let (bond_type, parameters) =
                    first.get_mmff_bond_stretch_params(0, 1).unwrap().unwrap();
                assert_eq!(bond_type, 0);
                assert_eq!((parameters.kb, parameters.r0), (5.573, 1.374));
                assert_eq!(
                    repeated.get_mmff_bond_stretch_params(0, 1).unwrap(),
                    first.get_mmff_bond_stretch_params(0, 1).unwrap()
                );
            }

            let mut direct = mmff_get_molecule_force_field(&mol, variant, f64::INFINITY, -1, true)
                .unwrap()
                .unwrap();
            assert_eq!(direct.positions(), coordinates);
            let direct_energy = direct.energy();
            assert!(direct_energy.is_finite());
            let direct_gradient = direct.gradient();
            assert!(
                direct_gradient
                    .iter()
                    .flatten()
                    .all(|value| value.is_finite())
            );
            // These include the adapter's preparation-then-factory path. Each
            // factory prepares an already _MMFFSanitized graph one more time.
            for prepared_molecule in [&first.molecule, &repeated.molecule] {
                let mut prepared = mmff_get_molecule_force_field(
                    prepared_molecule,
                    variant,
                    f64::INFINITY,
                    -1,
                    true,
                )
                .unwrap()
                .unwrap();
                assert_eq!(prepared.positions(), coordinates);
                let prepared_energy = prepared.energy();
                assert!(prepared_energy.is_finite());
                assert_eq!(
                    direct_energy, prepared_energy,
                    "energy {smiles}/{variant:?}"
                );
                let prepared_gradient = prepared.gradient();
                assert!(
                    prepared_gradient
                        .iter()
                        .flatten()
                        .all(|value| value.is_finite())
                );
                assert_eq!(
                    direct_gradient, prepared_gradient,
                    "gradient {smiles}/{variant:?}"
                );
            }
            assert_eq!(mol, original);
        }
    }
}

#[test]
fn malformed_inputs_return_structured_errors_and_setters_are_atomic() {
    let mol = ethanol();
    for threshold in [f64::NAN, -1., f64::NEG_INFINITY] {
        assert!(matches!(
            mmff_get_molecule_force_field(&mol, MmffVariant::Mmff94, threshold, 7, true),
            Err(MmffPublicApiError::ForceField(
                MolecularForceFieldError::InvalidThreshold { .. }
            ))
        ));
        assert!(matches!(
            uff_get_molecule_force_field(&mol, threshold, 7, true),
            Err(UffPublicApiError::ForceField(
                MolecularForceFieldError::InvalidThreshold { .. }
            ))
        ));
    }
    assert!(matches!(
        mmff_get_molecule_force_field(&mol, MmffVariant::Mmff94, 100., -2, true),
        Err(MmffPublicApiError::ForceField(
            MolecularForceFieldError::InvalidConformerId { conf_id: -2 }
        ))
    ));
    assert!(matches!(
        uff_get_molecule_force_field(&mol, 10., -2, true),
        Err(UffPublicApiError::ForceField(
            MolecularForceFieldError::InvalidConformerId { conf_id: -2 }
        ))
    ));
    assert!(matches!(
        uff_get_molecule_force_field(&mol, 10., 99, true),
        Err(UffPublicApiError::Missing3dConformer { conf_id: 99 })
    ));
    let no_coordinates = Molecule::from_smiles("CCO").unwrap();
    assert!(matches!(
        mmff_get_molecule_force_field(&no_coordinates, MmffVariant::Mmff94, 100., -1, true),
        Err(MmffPublicApiError::Builder(_))
    ));
    assert!(matches!(
        uff_get_molecule_force_field(&no_coordinates, 10., -1, true),
        Err(UffPublicApiError::Missing3dConformer { conf_id: -1 })
    ));
    assert!(matches!(
        mmff_get_molecule_force_field(&Molecule::new(), MmffVariant::Mmff94, 100., -1, true),
        Err(MmffPublicApiError::Builder(_))
    ));
    assert!(matches!(
        uff_get_molecule_force_field(&Molecule::new(), 10., -1, true),
        Err(UffPublicApiError::EmptyMolecule)
    ));
    for mut field in fields(&mol) {
        let positions = field.positions();
        assert!(matches!(
            field.set_positions(&positions[..2]),
            Err(MolecularForceFieldError::CoordinateCount { .. })
        ));
        let mut invalid = positions.clone();
        invalid[1][2] = f64::NAN;
        assert_eq!(
            field.set_positions(&invalid),
            Err(MolecularForceFieldError::NonfiniteCoordinate { atom: 1, axis: 2 })
        );
        assert_eq!(field.positions(), positions);
        field.set_fixed_points(&[0]).unwrap();
        assert_eq!(
            field.set_fixed_points(&[0, 9]),
            Err(MolecularForceFieldError::InvalidFixedPoint { atom: 9, atoms: 9 })
        );
        assert_eq!(
            field.set_fixed_points(&[0, 0]),
            Err(MolecularForceFieldError::DuplicateFixedPoint { atom: 0 })
        );
        assert_eq!(field.fixed_points(), [0]);
        assert_eq!(
            field.minimize(10, 0., 1e-6),
            Err(MolecularForceFieldError::InvalidTolerances)
        );
        assert_eq!(
            field.minimize(10, 1e-4, f64::NAN),
            Err(MolecularForceFieldError::InvalidTolerances)
        );
        assert_eq!(field.positions(), positions);
    }
}
