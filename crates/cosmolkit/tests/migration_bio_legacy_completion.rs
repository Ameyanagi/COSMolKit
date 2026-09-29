#![cfg(feature = "bio")]
use std::alloc::GlobalAlloc;

use cosmolkit::{
    BINDING_CONTRACT, BindingKind, BindingSupport, BioStructure, Protein, ProteinAtomIter,
    ProteinChainIter, ProteinChainRef, ProteinResidueIter, ProteinResidueRef,
};

const CIF: &str = include_str!("../../../testdata/bio/fixtures/gemmi_full_feature_sample.cif");

thread_local! {
    static COUNTING: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    static ALLOCATIONS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

struct CurrentThreadAllocator;

// SAFETY: All allocation/deallocation is delegated unchanged to System;
// the TLS probe only counts allocations on the calling thread and never allocates.
unsafe impl std::alloc::GlobalAlloc for CurrentThreadAllocator {
    unsafe fn alloc(&self, layout: std::alloc::Layout) -> *mut u8 {
        COUNTING.with(|enabled| {
            if enabled.get() {
                ALLOCATIONS.with(|count| count.set(count.get() + 1));
            }
        });
        // SAFETY: same layout and allocation contract as the caller.
        unsafe { std::alloc::System.alloc(layout) }
    }

    unsafe fn alloc_zeroed(&self, layout: std::alloc::Layout) -> *mut u8 {
        COUNTING.with(|enabled| {
            if enabled.get() {
                ALLOCATIONS.with(|count| count.set(count.get() + 1));
            }
        });
        // SAFETY: same layout and allocation contract as the caller.
        unsafe { std::alloc::System.alloc_zeroed(layout) }
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: std::alloc::Layout, size: usize) -> *mut u8 {
        COUNTING.with(|enabled| {
            if enabled.get() {
                ALLOCATIONS.with(|count| count.set(count.get() + 1));
            }
        });
        // SAFETY: same pointer, layout and size contract as the caller.
        unsafe { std::alloc::System.realloc(ptr, layout, size) }
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: std::alloc::Layout) {
        // SAFETY: same pointer and layout contract as the caller.
        unsafe { std::alloc::System.dealloc(ptr, layout) }
    }
}

#[global_allocator]
static ALLOCATOR: CurrentThreadAllocator = CurrentThreadAllocator;

fn measure_allocations<T>(work: impl FnOnce() -> T) -> (T, usize) {
    // Initialize both TLS slots before entering the measured region.
    COUNTING.with(|_| {});
    ALLOCATIONS.with(|count| count.set(0));
    COUNTING.with(|enabled| enabled.set(true));
    struct Reset;
    impl Drop for Reset {
        fn drop(&mut self) {
            COUNTING.with(|enabled| enabled.set(false));
        }
    }
    let reset = Reset;
    let result = work();
    drop(reset);
    (result, ALLOCATIONS.with(std::cell::Cell::get))
}

#[test]
fn bio_legacy_i06_six_traversals_create_step_and_exhaust_without_allocations() {
    let structure = BioStructure::from_mmcif(CIF).unwrap();
    let protein = structure.protein().unwrap();
    let chain = protein.chain(0).unwrap();
    let residue = chain.residues().next().unwrap();
    let empty = BioStructure::from_mmcif("data_empty\n")
        .unwrap()
        .protein()
        .unwrap();

    // Fixed source has one chain, two CYS residues and two atoms (IDs 0, 1).
    // Results are held in stack arrays; no Vec is constructed in a measured closure.
    macro_rules! check_ids {
        ($iterator:expr, $expected:expr) => {{
            let ((ids, count), allocations) = measure_allocations(|| {
                let mut iter = $iterator;
                let first = iter.next().map(|item| item.id().value());
                let mut ids = [u32::MAX; 8];
                let mut count = 0;
                if let Some(id) = first {
                    ids[count] = id;
                    count += 1;
                }
                for item in iter.by_ref() {
                    ids[count] = item.id().value();
                    count += 1;
                }
                assert!(iter.next().is_none());
                (ids, count)
            });
            assert_eq!(allocations, 0);
            let expected: &[u32] = $expected;
            assert_eq!(&ids[..count], expected);
        }};
    }
    check_ids!(protein.chains(), &[0]);
    check_ids!(protein.residues(), &[0, 1]);
    check_ids!(protein.atoms(), &[0, 1]);
    check_ids!(chain.residues(), &[0, 1]);
    check_ids!(chain.atoms(), &[0, 1]);
    check_ids!(residue.atoms(), &[0]);
    check_ids!(empty.chains(), &[]);
    check_ids!(empty.residues(), &[]);
    check_ids!(empty.atoms(), &[]);
    let (early, count) = measure_allocations(|| {
        let mut iter = protein.atoms();
        iter.next().unwrap().id().value()
    });
    assert_eq!((early, count), (0, 0));
}

#[test]
fn bio_legacy_i05_canonical_borrowed_iterators_and_registry() {
    let _: for<'a> fn(&'a Protein) -> ProteinChainIter<'a> = Protein::chains;
    let _: for<'a> fn(&'a Protein) -> ProteinResidueIter<'a> = Protein::residues;
    let _: for<'a> fn(&'a Protein) -> ProteinAtomIter<'a> = Protein::atoms;
    fn borrowed<'a>(chain: ProteinChainRef<'a>, residue: ProteinResidueRef<'a>) {
        let _: ProteinResidueIter<'a> = chain.residues();
        let _: ProteinAtomIter<'a> = chain.atoms();
        let _: ProteinAtomIter<'a> = residue.atoms();
    }
    let _: for<'a> fn(ProteinChainRef<'a>, ProteinResidueRef<'a>) = borrowed;
    for (id, rust_type, python, javascript) in [
        (
            "types.ProteinChainIter",
            "ProteinChainIter",
            "ProteinChainIter",
            "ProteinChainIter",
        ),
        (
            "types.ProteinResidueIter",
            "ProteinResidueIter",
            "ProteinResidueIter",
            "ProteinResidueIter",
        ),
        (
            "types.ProteinAtomIter",
            "ProteinAtomIter",
            "ProteinAtomIter",
            "ProteinAtomIter",
        ),
    ] {
        let matches = BINDING_CONTRACT
            .iter()
            .filter(|row| row.semantic_id == id)
            .collect::<Vec<_>>();
        assert_eq!(matches.len(), 1, "{id}");
        assert!(matches[0].rust_path.contains(rust_type));
        assert_eq!(matches[0].python_name, python);
        assert_eq!(matches[0].javascript_name, javascript);
        assert_eq!(matches[0].feature, "bio");
        assert_eq!(matches[0].support, BindingSupport::Experimental);
    }
    for id in [
        "Protein.chains",
        "Protein.residues",
        "Protein.atoms",
        "ProteinChainRef.residues",
        "ProteinChainRef.atoms",
        "ProteinResidueRef.atoms",
    ] {
        let matches = BINDING_CONTRACT
            .iter()
            .filter(|row| row.semantic_id == id)
            .collect::<Vec<_>>();
        assert_eq!(matches.len(), 1, "{id}");
        assert_eq!(matches[0].feature, "bio");
        assert_eq!(matches[0].support, BindingSupport::Experimental);
    }
    let structure = BioStructure::from_mmcif(CIF).unwrap();
    let protein = structure.protein().unwrap();
    let mut chains = protein.chains();
    let first = chains.next().unwrap();
    assert_eq!(first.id().value(), 0);
    let first_residue = first.residues().next().unwrap();
    assert_eq!(first_residue.id(), protein.residues().next().unwrap().id());
    assert_eq!(
        first_residue.atoms().next().unwrap().id(),
        protein.atoms().next().unwrap().id()
    );
    assert_eq!(
        first.atoms().next().unwrap().id(),
        protein.atoms().next().unwrap().id()
    );
    assert_eq!(protein.chains().count(), protein.num_chains());
    assert_eq!(protein.residues().count(), protein.num_residues());
    assert_eq!(protein.atoms().count(), protein.num_atoms());
    assert_eq!(
        first_residue.atoms().next().unwrap().position(),
        protein.atoms().next().unwrap().position()
    );
}

#[test]
fn bio_legacy_n03_public_borrow_registry_and_protein_projection() {
    let _: for<'a> fn(&'a BioStructure) -> Option<&'a str> = BioStructure::ncs_oper_identity_id;
    let row = BINDING_CONTRACT
        .iter()
        .find(|entry| entry.semantic_id == "BioStructure.ncs_oper_identity_id")
        .unwrap();
    assert_eq!(row.feature, "bio");
    assert_eq!(row.support, BindingSupport::Experimental);
    assert_eq!(row.python_name, "ncs_oper_identity_id");
    assert_eq!(row.javascript_name, "ncsOperIdentityId");
    assert_eq!(row.callable.unwrap().kind, BindingKind::Instance);
    assert!(row.callable.unwrap().parameters.is_empty());
    let absent = BioStructure::from_mmcif("data_demo\n_entry.id DEMO\n").unwrap();
    assert_eq!(absent.ncs_oper_identity_id(), None);
    let structure = BioStructure::from_mmcif(CIF).unwrap();
    assert_eq!(structure.ncs_oper_identity_id(), Some("I"));
    assert_eq!(
        structure.ncs_oper_identity_id().unwrap().as_ptr(),
        structure.source_state().info["_struct_ncs_oper.id"].as_ptr()
    );
    let protein = structure.protein().unwrap();
    assert_eq!(protein.as_bio_structure().ncs_oper_identity_id(), Some("I"));
}
