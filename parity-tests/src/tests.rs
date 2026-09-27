use super::*;
use registry::Pair;
use std::cell::Cell;

fn fingerprint_tasks() -> Result<Vec<&'static Task>> {
    Ok(registry::TASKS
        .iter()
        .filter(|t| !matches!(t.operation, registry::Operation::Molecular(_)))
        .collect())
}
fn fingerprint_corpus() -> Corpus {
    Corpus {
        fingerprints: registry::builtin(),
        molecules: vec![],
    }
}
fn fingerprint(input: &Input) -> &registry::FingerprintInput {
    let Input::Fingerprint(input) = input else {
        panic!("fingerprint expected")
    };
    input
}

#[test]
fn public_input_adapter_preserves_extreme_counts_and_stored_zero() {
    let entries = vec![(1, i32::MIN), (2, i32::MAX), (3, 0)];
    let cases = Corpus {
        fingerprints: vec![Pair {
            id: "adapter_boundaries".into(),
            length: 4,
            left: entries.clone(),
            right: entries.clone(),
        }],
        molecules: vec![],
    };
    let tasks = fingerprint_tasks().unwrap();
    registry::validate(&cases, &tasks).unwrap();
    for task in tasks {
        for input in registry::expand(&cases, task) {
            assert_eq!(
                execute::run(&input).unwrap().output,
                registry::Value::Fingerprint(registry::FingerprintValue {
                    length: 4,
                    entries: entries.clone()
                })
            );
        }
    }
}

// Framework tests intentionally use synthetic reference values. They test
// validation/comparison machinery, NOT RDKit parity; the CLI uses the oracle.
fn fixture(data: &Path, task: &Task, cases: &Corpus) -> PathBuf {
    let inputs = registry::expand(cases, task);
    let records: Vec<_> = inputs
        .iter()
        .map(|input| Record {
            input: input.clone(),
            output: registry::Value::Fingerprint(registry::FingerprintValue {
                length: fingerprint(input).case.length,
                entries: vec![],
            }),
        })
        .collect();
    let input = encode(&inputs).unwrap();
    let reference = encode(&records).unwrap();
    let directory = generation(data, task, &input);
    fs::create_dir_all(&directory).unwrap();
    fs::write(directory.join("input.json"), &input).unwrap();
    fs::write(directory.join("reference.json"), &reference).unwrap();
    fs::write(
        directory.join("manifest.json"),
        encode(&identity(task, &input, &reference, inputs.len())).unwrap(),
    )
    .unwrap();
    directory
}

#[test]
fn missing_last_task_stops_before_first_operation() {
    let temp = tempfile::tempdir().unwrap();
    let tasks = fingerprint_tasks().unwrap();
    let cases = fingerprint_corpus();
    fixture(temp.path(), tasks[0], &cases);
    let calls = Cell::new(0);
    let result = run_with(
        &tasks,
        &cases,
        temp.path(),
        |_| Err("oracle unavailable".into()),
        |i| {
            calls.set(calls.get() + 1);
            execute::run(i)
        },
    );
    assert!(result.unwrap_err().contains("fuzzy_or"));
    assert_eq!(calls.get(), 0);
}

fn synthetic(inputs: &[Input]) -> Result<Vec<Record>> {
    Ok(inputs
        .iter()
        .map(|input| Record {
            input: input.clone(),
            output: registry::Value::Fingerprint(registry::FingerprintValue {
                length: fingerprint(input).case.length,
                entries: vec![],
            }),
        })
        .collect())
}

#[test]
fn fresh_run_prepares_every_task_before_execution_then_reuses_cache() {
    let temp = tempfile::tempdir().unwrap();
    let tasks = fingerprint_tasks().unwrap();
    let cases = fingerprint_corpus();
    let generated = Cell::new(0);
    let executed = Cell::new(0);
    let report = run_with(
        &tasks,
        &cases,
        temp.path(),
        |inputs| {
            generated.set(generated.get() + 1);
            synthetic(inputs)
        },
        |input| {
            assert_eq!(generated.get(), tasks.len());
            executed.set(executed.get() + 1);
            Ok(synthetic(std::slice::from_ref(input))?.remove(0))
        },
    )
    .unwrap();
    assert_eq!(executed.get(), cases.fingerprints.len() * 4);
    assert!(report.iter().all(|row| row.matches));
    let (_, preparation) = ensure_ready_with(&tasks, &cases, temp.path(), |_| {
        panic!("valid cache must not invoke oracle")
    })
    .unwrap();
    assert_eq!(
        preparation,
        Preparation {
            reused_tasks: 2,
            generated_tasks: 0,
            rows: report.len()
        }
    );
}

#[test]
fn repairs_only_corrupt_task_and_preserves_evidence() {
    let temp = tempfile::tempdir().unwrap();
    let tasks = fingerprint_tasks().unwrap();
    let cases = fingerprint_corpus();
    fixture(temp.path(), tasks[0], &cases);
    let broken = fixture(temp.path(), tasks[1], &cases);
    fs::write(broken.join("reference.json"), b"[]").unwrap();
    let (_, preparation) = ensure_ready_with(&tasks, &cases, temp.path(), |inputs| {
        assert_eq!(
            fingerprint(&inputs[0]).operation,
            registry::Operation::FuzzyOr
        );
        synthetic(inputs)
    })
    .unwrap();
    assert_eq!(preparation.generated_tasks, 1);
    assert_eq!(preparation.reused_tasks, 1);
    let backup = fs::read_dir(temp.path())
        .unwrap()
        .map(|e| e.unwrap().path())
        .find(|p| {
            p.file_name()
                .unwrap()
                .to_string_lossy()
                .starts_with(".invalid-")
        })
        .unwrap();
    assert_eq!(read(&backup.join("reference.json")).unwrap(), b"[]");
    assert_eq!(
        preflight(&tasks, &cases, temp.path()).unwrap().len(),
        cases.fingerprints.len() * 4
    );
}

#[test]
fn malformed_oracle_never_publishes_or_executes() {
    let temp = tempfile::tempdir().unwrap();
    let tasks = fingerprint_tasks().unwrap();
    let cases = fingerprint_corpus();
    let result = run_with(
        &tasks,
        &cases,
        temp.path(),
        |_| Ok(vec![]),
        |_| panic!("must not execute"),
    );
    assert!(result.unwrap_err().contains("reference row count"));
    assert!(
        !generation(
            temp.path(),
            tasks[0],
            &encode(&registry::expand(&cases, tasks[0])).unwrap()
        )
        .exists()
    );
}

#[test]
fn failed_final_generation_does_not_execute_prepared_first_task() {
    let temp = tempfile::tempdir().unwrap();
    let tasks = fingerprint_tasks().unwrap();
    let cases = fingerprint_corpus();
    let result = run_with(
        &tasks,
        &cases,
        temp.path(),
        |inputs| {
            if fingerprint(&inputs[0]).operation == registry::Operation::FuzzyOr {
                Err("oracle failure".into())
            } else {
                synthetic(inputs)
            }
        },
        |_| panic!("global barrier must prevent execution"),
    );
    assert!(result.unwrap_err().contains("0 Rust operation calls"));
    assert!(preflight(&tasks[..1], &cases, temp.path()).is_ok());
}

#[test]
fn stale_manifest_is_repaired_and_unselected_tasks_are_not_generated() {
    let temp = tempfile::tempdir().unwrap();
    let tasks = registry::select(Some("fuzzy_and")).unwrap();
    let cases = fingerprint_corpus();
    let directory = fixture(temp.path(), tasks[0], &cases);
    let mut manifest: Manifest =
        serde_json::from_slice(&read(&directory.join("manifest.json")).unwrap()).unwrap();
    manifest.rdkit_version = "stale".into();
    fs::write(directory.join("manifest.json"), encode(&manifest).unwrap()).unwrap();
    let (_, preparation) = ensure_ready_with(&tasks, &cases, temp.path(), |inputs| {
        assert_eq!(
            fingerprint(&inputs[0]).operation,
            registry::Operation::FuzzyAnd
        );
        synthetic(inputs)
    })
    .unwrap();
    assert_eq!(preparation.generated_tasks, 1);
    assert_eq!(preparation.rows, cases.fingerprints.len() * 2);
    assert!(preflight(&fingerprint_tasks().unwrap(), &cases, temp.path()).is_err());
}

#[test]
fn comparison_failures_never_rewrite_reference_data() {
    let temp = tempfile::tempdir().unwrap();
    let tasks = registry::select(Some("fuzzy_and")).unwrap();
    let cases = fingerprint_corpus();
    let directory = fixture(temp.path(), tasks[0], &cases);
    let before = read(&directory.join("reference.json")).unwrap();
    let report = run_with(
        &tasks,
        &cases,
        temp.path(),
        |_| panic!("no oracle for valid cache"),
        |_| Err("CK failed".into()),
    )
    .unwrap();
    assert!(report.iter().all(|row| !row.matches));
    assert_eq!(read(&directory.join("reference.json")).unwrap(), before);
}

#[test]
fn corrupted_reference_stops_before_execution() {
    let temp = tempfile::tempdir().unwrap();
    let tasks = fingerprint_tasks().unwrap();
    let cases = fingerprint_corpus();
    fixture(temp.path(), tasks[0], &cases);
    let last = fixture(temp.path(), tasks[1], &cases);
    fs::write(last.join("reference.json"), b"[]").unwrap();
    assert!(
        preflight(&tasks, &cases, temp.path())
            .err()
            .unwrap()
            .contains("corrupted")
    );
}

#[test]
fn selected_task_does_not_require_unselected_task() {
    let temp = tempfile::tempdir().unwrap();
    let tasks = registry::select(Some("fuzzy_and")).unwrap();
    let cases = fingerprint_corpus();
    fixture(temp.path(), tasks[0], &cases);
    assert_eq!(
        preflight(&tasks, &cases, temp.path()).unwrap().len(),
        cases.fingerprints.len() * 2
    );
}

#[test]
fn wrong_case_identity_is_rejected_even_with_updated_checksum() {
    let temp = tempfile::tempdir().unwrap();
    let tasks = registry::select(Some("fuzzy_or")).unwrap();
    let cases = fingerprint_corpus();
    let directory = fixture(temp.path(), tasks[0], &cases);
    let mut records: Vec<Record> =
        serde_json::from_slice(&read(&directory.join("reference.json")).unwrap()).unwrap();
    if let Input::Fingerprint(input) = &mut records[0].input {
        input.width = registry::Width::U64;
    }
    let reference = encode(&records).unwrap();
    let input = read(&directory.join("input.json")).unwrap();
    fs::write(directory.join("reference.json"), &reference).unwrap();
    fs::write(
        directory.join("manifest.json"),
        encode(&identity(tasks[0], &input, &reference, records.len())).unwrap(),
    )
    .unwrap();
    assert!(
        preflight(&tasks, &cases, temp.path())
            .err()
            .unwrap()
            .contains("case/operation/width")
    );
}

#[test]
fn comparison_checks_values_and_does_not_drop_failures() {
    let temp = tempfile::tempdir().unwrap();
    let tasks = fingerprint_tasks().unwrap();
    let cases = fingerprint_corpus();
    for task in &tasks {
        fixture(temp.path(), task, &cases);
    }
    let ready = preflight(&tasks, &cases, temp.path()).unwrap();
    let expected = ready.len();
    let report = compare(ready, |input| {
        Ok(Record {
            input: input.clone(),
            output: registry::Value::Fingerprint(registry::FingerprintValue {
                length: fingerprint(input).case.length,
                entries: vec![(1, 99)],
            }),
        })
    });
    assert_eq!(report.len(), expected);
    assert!(report.iter().all(|row| !row.matches));
}

#[test]
fn registry_rejects_unknown_task_and_invalid_corpus() {
    assert!(registry::select(Some("fuzzy_annd")).is_err());
    let tasks = fingerprint_tasks().unwrap();
    assert!(registry::validate(&Corpus::default(), &tasks).is_err());
    let mut cases = fingerprint_corpus();
    cases.fingerprints.push(cases.fingerprints[0].clone());
    assert!(registry::validate(&cases, &tasks).is_err());
    let mut cases = fingerprint_corpus();
    cases.fingerprints[0].left = vec![(16, 3)];
    assert!(registry::validate(&cases, &tasks).is_err());
}

#[test]
fn molecular_missing_final_reference_blocks_all_rust_calls() {
    let temp = tempfile::tempdir().unwrap();
    let tasks = registry::select(None).unwrap();
    let mut cases = fingerprint_corpus();
    cases.molecules.push(registry::SmilesCase {
        id: "ethanol".into(),
        smiles: "CCO".into(),
    });
    let result = run_with(
        &tasks,
        &cases,
        temp.path(),
        |inputs| {
            if let Input::Molecular {
                profile: registry::molecule_plan::Profile::RemoveHydrogens { .. },
                ..
            } = &inputs[0]
            {
                return Err("last reference missing".into());
            }
            inputs
                .iter()
                .map(|input| {
                    let output = match input {
                        Input::Fingerprint(input) => {
                            registry::Value::Fingerprint(registry::FingerprintValue {
                                length: input.case.length,
                                entries: vec![],
                            })
                        }
                        Input::Molecular { .. } => {
                            registry::Value::Molecular(molecular::Outcome::Error {
                                stage: molecular::Stage::Parse,
                                detail: "synthetic diagnostic".into(),
                            })
                        }
                    };
                    Ok(Record {
                        input: input.clone(),
                        output,
                    })
                })
                .collect()
        },
        |_| panic!("global molecular barrier bypassed"),
    );
    assert!(result.unwrap_err().contains("last reference missing"));
}

#[test]
fn molecular_schema_and_comparison_do_not_accept_wrong_types_or_equal_errors() {
    use registry::molecule_plan::Profile;
    assert!(
        molecular::validate_output(
            &Profile::MolecularWeight { only_heavy: false },
            &molecular::Outcome::Text("18".into())
        )
        .is_err()
    );
    let input = Input::Molecular {
        case: registry::SmilesCase {
            id: "invalid".into(),
            smiles: "C(".into(),
        },
        profile: Profile::SanitizeAll,
    };
    let record = Record {
        input,
        output: registry::Value::Molecular(molecular::Outcome::Error {
            stage: molecular::Stage::Parse,
            detail: "same text".into(),
        }),
    };
    let ready = Ready {
        records: vec![record.clone()],
    };
    assert!(!compare(ready, |_| Ok(record.clone()))[0].matches);
}

#[test]
fn molecular_corpus_preserves_blank_records_cx_whitespace_and_duplicates() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("input.smi");
    fs::write(&path, "CCO\n\nC |$label$| name\nCCO\n").unwrap();
    let cases = molecular::read_corpus(&path).unwrap();
    assert_eq!(cases.len(), 4);
    assert_eq!(cases[1].smiles, "");
    assert_eq!(cases[2].smiles, "C |$label$| name");
    assert_ne!(cases[0].id, cases[3].id);
}

#[test]
fn selected_input_family_does_not_load_an_unselected_corpus() {
    let tasks = registry::select(Some("fuzzy_and")).unwrap();
    let cases = corpus(None, &tasks).unwrap();
    assert!(cases.molecules.is_empty());
    assert_eq!(cases.fingerprints.len(), 6);
    assert!(
        corpus(Some(Path::new("missing.smi")), &tasks)
            .unwrap_err()
            .contains("input family")
    );
}

#[test]
fn coordinate_comparison_checks_shape_nonfinite_topology_and_declared_tolerance() {
    let sample = |x: f64, count| molecular::Outcome::Coordinates2d {
        topology: molecular::Topology {
            atoms: vec![],
            bonds: vec![],
        },
        xy_bits: vec![[x.to_bits(), 0]; count],
    };
    assert!(molecular::matches(&sample(0.0, 1), &sample(0.5e-8, 1)));
    assert!(!molecular::matches(&sample(0.0, 1), &sample(2e-8, 1)));
    assert!(!molecular::matches(&sample(0.0, 1), &sample(0.0, 0)));
    assert!(!molecular::matches(
        &sample(f64::NAN, 1),
        &sample(f64::NAN, 1)
    ));
    assert!(!molecular::matches(
        &sample(f64::INFINITY, 1),
        &sample(f64::INFINITY, 1)
    ));
}

#[test]
fn distance_matrix_profiles_and_schema_preserve_all_entries_and_float_bits() {
    use registry::molecule_plan::{Profile, TaskId};
    let profiles = TaskId::DistanceMatrix.profiles();
    assert_eq!(profiles.len(), 4);
    let profile = Profile::DistanceMatrix {
        use_bond_order: true,
        use_atom_weights: true,
    };
    assert!(profiles.contains(&profile));
    let matrix = molecular::Outcome::Matrix {
        dimension: 1,
        values_bits: vec![f64::INFINITY.to_bits()],
    };
    // Atomic number zero has an infinite weighted diagonal in the source.
    assert!(molecular::validate_output(&profile, &matrix).is_ok());
    assert!(molecular::matches(&matrix, &matrix));
    let wrong_shape = molecular::Outcome::Matrix {
        dimension: 2,
        values_bits: vec![0],
    };
    assert!(molecular::validate_output(&profile, &wrong_shape).is_err());
    let positive_zero = molecular::Outcome::Matrix {
        dimension: 1,
        values_bits: vec![0],
    };
    let negative_zero = molecular::Outcome::Matrix {
        dimension: 1,
        values_bits: vec![(-0.0f64).to_bits()],
    };
    assert!(!molecular::matches(&positive_zero, &negative_zero));
}
