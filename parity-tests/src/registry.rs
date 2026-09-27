//! Executable tasks, typed input families and parameter matrices.
//! Future catalog rows are not executable registration or parity claims.
pub mod molecule_plan;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Operation {
    FuzzyAnd,
    FuzzyOr,
    Molecular(molecule_plan::TaskId),
}

impl Operation {
    pub fn name(self) -> &'static str {
        match self {
            Self::FuzzyAnd => "fuzzy_and",
            Self::FuzzyOr => "fuzzy_or",
            Self::Molecular(id) => id.name(),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Width {
    U32,
    U64,
}

pub struct Task {
    pub operation: Operation,
    pub widths: &'static [Width],
}

pub const TASKS: &[Task] = &[
    Task {
        operation: Operation::FuzzyAnd,
        widths: &[Width::U32, Width::U64],
    },
    Task {
        operation: Operation::FuzzyOr,
        widths: &[Width::U32, Width::U64],
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::SmilesRead),
        widths: &[],
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::Sanitize),
        widths: &[],
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::Kekulize),
        widths: &[],
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::MolecularWeight),
        widths: &[],
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::ExactMolecularWeight),
        widths: &[],
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::MolecularFormula),
        widths: &[],
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::AddHydrogens),
        widths: &[],
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::RemoveHydrogens),
        widths: &[],
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::Coordinates2d),
        widths: &[],
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::DistanceMatrix),
        widths: &[],
    },
];

pub const RDKIT_VERSION: &str = "2026.03.1";

// Only typed integer indices/counts cross the oracle boundary.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Pair {
    pub id: String,
    pub length: u64,
    pub left: Vec<(u64, i32)>,
    pub right: Vec<(u64, i32)>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FingerprintInput {
    pub case: Pair,
    pub operation: Operation,
    pub width: Width,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FingerprintValue {
    pub length: u64,
    pub entries: Vec<(u64, i32)>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Record {
    pub input: Input,
    pub output: Value,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SmilesCase {
    pub id: String,
    pub smiles: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Input {
    Fingerprint(FingerprintInput),
    Molecular {
        case: SmilesCase,
        profile: molecule_plan::Profile,
    },
}

impl Input {
    pub fn task_name(&self) -> &'static str {
        match self {
            Self::Fingerprint(input) => input.operation.name(),
            Self::Molecular { profile, .. } => {
                use molecule_plan::Profile::*;
                match profile {
                    SmilesRead { .. } => "smiles_read",
                    SanitizeAll => "sanitize",
                    Kekulize { .. } => "kekulize",
                    MolecularWeight { .. } => "molecular_weight",
                    ExactMolecularWeight { .. } => "exact_molecular_weight",
                    MolecularFormula { .. } => "molecular_formula",
                    AddHydrogens { .. } => "add_hydrogens",
                    RemoveHydrogens { .. } => "remove_hydrogens",
                    Coordinates2dDefault => "coordinates_2d",
                    CipLabels { .. } => "cip_labels",
                    PotentialStereo { .. } => "potential_stereo",
                    Valence { .. } => "valence",
                    DistanceMatrix { .. } => "distance_matrix",
                }
            }
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Value {
    Fingerprint(FingerprintValue),
    Molecular(crate::molecular::Outcome),
}

#[derive(Clone, Debug, Default)]
pub struct Corpus {
    pub fingerprints: Vec<Pair>,
    pub molecules: Vec<SmilesCase>,
}

impl Task {
    pub fn count(&self, cases: &Corpus) -> usize {
        match self.operation {
            Operation::Molecular(id) => cases.molecules.len() * id.profiles().len(),
            _ => cases.fingerprints.len() * self.widths.len(),
        }
    }
}

pub fn select(name: Option<&str>) -> Result<Vec<&'static Task>, String> {
    let selected: Vec<_> = TASKS
        .iter()
        .filter(|t| name.is_none_or(|n| n == t.operation.name()))
        .collect();
    if selected.is_empty() {
        return Err(format!("unknown task: {name:?}"));
    }
    Ok(selected)
}

pub fn validate(corpus: &Corpus, tasks: &[&Task]) -> Result<(), String> {
    use std::collections::BTreeSet;
    if tasks.is_empty() {
        return Err("empty task selection".into());
    }
    for task in tasks {
        if task.count(corpus) == 0 {
            return Err(format!(
                "{}: corpus/profile selection is empty",
                task.operation.name()
            ));
        }
    }
    let mut molecule_ids = BTreeSet::new();
    for case in &corpus.molecules {
        if case.id.is_empty() || !molecule_ids.insert(&case.id) {
            return Err("empty/duplicate molecular case ID".into());
        }
    }
    let mut ids = BTreeSet::new();
    for c in &corpus.fingerprints {
        if c.id.is_empty() || !ids.insert(&c.id) {
            return Err("empty/duplicate case ID".into());
        }
        for width in tasks.iter().flat_map(|t| t.widths) {
            let max = match width {
                Width::U32 => u32::MAX as u64,
                Width::U64 => u64::MAX,
            };
            if c.length > max {
                return Err(format!("{}: length does not fit {width:?}", c.id));
            }
            for values in [&c.left, &c.right] {
                let mut keys = BTreeSet::new();
                for &(key, _) in values {
                    if !keys.insert(key) || key >= c.length {
                        return Err(format!("{}: invalid/duplicate index", c.id));
                    }
                }
            }
        }
    }
    Ok(())
}

pub fn expand(cases: &Corpus, task: &Task) -> Vec<Input> {
    if let Operation::Molecular(id) = task.operation {
        return cases
            .molecules
            .iter()
            .flat_map(|case| {
                id.profiles()
                    .into_iter()
                    .map(move |profile| Input::Molecular {
                        case: case.clone(),
                        profile,
                    })
            })
            .collect();
    }
    cases
        .fingerprints
        .iter()
        .flat_map(|case| {
            task.widths.iter().map(move |&width| {
                Input::Fingerprint(FingerprintInput {
                    case: case.clone(),
                    operation: task.operation,
                    width,
                })
            })
        })
        .collect()
}

pub fn builtin() -> Vec<Pair> {
    // Each row names a source branch, not an arbitrary sampling size.
    [
        ("empty_both", vec![], vec![]),
        ("left_empty_right_tail", vec![], vec![(2, 3)]),
        ("right_empty_remove_left", vec![(2, -3)], vec![]),
        (
            "shared_signed_min_max",
            vec![(1, 5), (3, -2), (8, 4)],
            vec![(1, 3), (3, -4), (9, 7)],
        ),
        (
            "disjoint_interleaved",
            vec![(1, 2), (3, -1)],
            vec![(0, 4), (2, -3), (4, 5)],
        ),
        (
            "explicit_zero_shared_and_exclusive",
            vec![(1, 0), (3, 0)],
            vec![(1, -2), (4, 0)],
        ),
    ]
    .into_iter()
    .map(|(id, left, right)| Pair {
        id: id.into(),
        length: 16,
        left,
        right,
    })
    .collect()
}
