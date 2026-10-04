use cosmolkit_core::{EmbedParameters, Molecule, embed_multiple_confs};

fn parameters() -> EmbedParameters {
    let mut params = EmbedParameters::etkdg_v3();
    params.random_seed = 7;
    params.num_threads = 1;
    // This selects the actual fourth-dimension stage. The bonded carbonyl also
    // supplies nonempty distance and improper-planarity contributions.
    params.use_random_coords = true;
    params.timeout = 0;
    params
}

#[test]
fn public_etkdg_random_coordinates_execute_owned_3d_and_4d_fields_reproducibly() {
    let molecule = Molecule::from_smiles("CC=O")
        .unwrap()
        .with_hydrogens()
        .unwrap();
    let original = molecule.clone();
    let (embedded, ids) = embed_multiple_confs(&molecule, 1, &mut parameters()).unwrap();
    assert_eq!(ids, [0]);
    assert_eq!(embedded.conformers_3d().len(), 1);
    let conformer = &embedded.conformers_3d()[0];
    assert!(conformer.is_3d());
    assert_eq!(conformer.coordinates().len(), molecule.num_atoms());
    assert!(
        conformer
            .coordinates()
            .iter()
            .flatten()
            .all(|value| value.is_finite())
    );
    assert_eq!(embedded.atoms(), molecule.atoms());
    assert_eq!(embedded.bonds(), molecule.bonds());
    assert_eq!(molecule, original);

    let (repeated, repeated_ids) = embed_multiple_confs(&molecule, 1, &mut parameters()).unwrap();
    assert_eq!(repeated_ids, ids);
    assert_eq!(
        repeated.conformers_3d()[0].coordinates(),
        conformer.coordinates()
    );
    assert_eq!(molecule, original);
}
