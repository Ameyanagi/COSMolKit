//! Owned molecular force fields with checked coordinates and exact fixed atoms.
//!
//! This Rust ownership adapter reuses the existing source-backed MMFF/UFF
//! contributions and BFGS driver. It does not expose the movable, self-referential
//! primitive `ForceField` or its contributions.

use std::collections::HashSet;

use super::{
    OwnedForceField,
    core::{ForceField, ForceFieldVec3},
};

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum MolecularForceFieldError {
    #[error("force-field threshold must be nonnegative and not NaN, got {threshold}")]
    InvalidThreshold { threshold: f64 },
    #[error("force-field conformer id must be >= -1, got {conf_id}")]
    InvalidConformerId { conf_id: isize },
    #[error("force-field coordinates require {expected} atoms, got {actual}")]
    CoordinateCount { expected: usize, actual: usize },
    #[error("nonfinite force-field coordinate at atom {atom}, axis {axis}")]
    NonfiniteCoordinate { atom: usize, axis: usize },
    #[error("fixed atom {atom} is out of range for {atoms} atoms")]
    InvalidFixedPoint { atom: usize, atoms: usize },
    #[error("fixed atom {atom} appears more than once")]
    DuplicateFixedPoint { atom: usize },
    #[error("force-field tolerances must be finite and positive")]
    InvalidTolerances,
}

/// A movable, opaque force field owning its XYZ coordinates in Angstroms.
///
/// Coordinates follow the input molecule's atom order. Energies are in kcal/mol
/// and gradients in kcal/(mol Angstrom). Fixed atoms constrain minimization;
/// [`Self::gradient`] always returns the physical, unconstrained derivative.
/// Setters validate the whole input before changing the field. Neither building
/// nor minimizing this value changes the source [`crate::Molecule`].
pub struct MolecularForceField {
    field: OwnedForceField,
}

impl std::fmt::Debug for MolecularForceField {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("MolecularForceField")
            .field("positions", &self.positions())
            .field("fixed_points", &self.fixed_points())
            .finish()
    }
}

impl MolecularForceField {
    pub(super) fn from_force_field(owner: ForceField) -> Self {
        Self {
            field: OwnedForceField::from_force_field(owner),
        }
    }

    #[must_use]
    pub fn positions(&self) -> Vec<[f64; 3]> {
        self.field
            .positions()
            .iter()
            .map(|p| [p.x, p.y, p.z])
            .collect()
    }

    pub fn set_positions(
        &mut self,
        positions: &[[f64; 3]],
    ) -> Result<(), MolecularForceFieldError> {
        validate_positions(positions, self.field.num_points())?;
        let points: Vec<_> = positions
            .iter()
            .map(|&[x, y, z]| ForceFieldVec3::new(x, y, z))
            .collect();
        self.field.set_positions(&points);
        Ok(())
    }

    #[must_use]
    pub fn fixed_points(&self) -> Vec<usize> {
        self.field.fixed_points().to_vec()
    }

    pub fn set_fixed_points(
        &mut self,
        fixed_points: &[usize],
    ) -> Result<(), MolecularForceFieldError> {
        let atoms = self.field.num_points();
        let mut seen = HashSet::with_capacity(fixed_points.len());
        for &atom in fixed_points {
            if atom >= atoms {
                return Err(MolecularForceFieldError::InvalidFixedPoint { atom, atoms });
            }
            if !seen.insert(atom) {
                return Err(MolecularForceFieldError::DuplicateFixedPoint { atom });
            }
        }
        self.field.set_fixed_points(fixed_points);
        Ok(())
    }

    #[must_use]
    pub fn energy(&mut self) -> f64 {
        self.field.energy_current(None)
    }

    #[must_use]
    pub fn gradient(&mut self) -> Vec<[f64; 3]> {
        self.field
            .gradient_current()
            .chunks_exact(3)
            .map(|p| [p[0], p[1], p[2]])
            .collect()
    }

    /// Minimize with the existing source-backed BFGS implementation. Returns
    /// zero on convergence and one when more iterations are needed.
    pub fn minimize(
        &mut self,
        max_its: usize,
        force_tol: f64,
        energy_tol: f64,
    ) -> Result<i32, MolecularForceFieldError> {
        if !force_tol.is_finite() || force_tol <= 0. || !energy_tol.is_finite() || energy_tol <= 0.
        {
            return Err(MolecularForceFieldError::InvalidTolerances);
        }
        let original = self.field.positions().to_vec();
        let status = self.field.minimize(max_its, force_tol, energy_tol);
        if let Err(error) = validate_positions(&self.positions(), self.field.num_points()) {
            self.field.set_positions(&original);
            return Err(error);
        }
        Ok(status)
    }
}

pub(super) fn validate_options(
    threshold: f64,
    conf_id: isize,
) -> Result<(), MolecularForceFieldError> {
    // Positive infinity explicitly requests all nonbonded pairs.
    if threshold.is_nan() || threshold < 0. {
        return Err(MolecularForceFieldError::InvalidThreshold { threshold });
    }
    if conf_id < -1 {
        return Err(MolecularForceFieldError::InvalidConformerId { conf_id });
    }
    Ok(())
}

pub(super) fn validate_positions(
    positions: &[[f64; 3]],
    expected: usize,
) -> Result<(), MolecularForceFieldError> {
    if positions.len() != expected {
        return Err(MolecularForceFieldError::CoordinateCount {
            expected,
            actual: positions.len(),
        });
    }
    for (atom, point) in positions.iter().enumerate() {
        for (axis, coordinate) in point.iter().enumerate() {
            if !coordinate.is_finite() {
                return Err(MolecularForceFieldError::NonfiniteCoordinate { atom, axis });
            }
        }
    }
    Ok(())
}
