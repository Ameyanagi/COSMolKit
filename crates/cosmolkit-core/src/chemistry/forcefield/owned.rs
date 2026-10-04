//! Dimension-generic evaluator for source-built force fields.
//!
//! Context and contributions stay separate. The primitive field is used only
//! as a coordinate/distance context; it never evaluates its own contributions.

use std::cell::RefCell;

use super::core::{ForceField, ForceFieldContrib, ForceFieldVec3, bfgs_minimize, scale_gradient};

pub(crate) struct OwnedForceField {
    // Terms are dropped before their context, and neither can escape.
    contribs: Vec<Box<dyn ForceFieldContrib>>,
    owner: ForceField,
    fixed_points: Vec<usize>,
}

impl OwnedForceField {
    pub(crate) fn from_force_field(mut owner: ForceField) -> Self {
        let mut contribs = std::mem::take(owner.contribs_mut());
        for contrib in &mut contribs {
            contrib.set_force_field(std::ptr::null());
        }
        let fixed_points = std::mem::take(owner.fixed_points_mut());
        owner.initialize();
        Self {
            contribs,
            owner,
            fixed_points,
        }
    }

    pub(crate) fn dimension(&self) -> usize {
        self.owner.dimension()
    }

    pub(crate) fn num_points(&self) -> usize {
        self.owner.num_points()
    }

    pub(crate) fn positions(&self) -> &[ForceFieldVec3] {
        self.owner.positions()
    }

    pub(crate) fn set_positions(&mut self, positions: &[ForceFieldVec3]) {
        assert_eq!(positions.len(), self.num_points());
        self.owner.positions_mut().clone_from_slice(positions);
        self.owner.init_distance_matrix();
    }

    pub(crate) fn fixed_points(&self) -> &[usize] {
        &self.fixed_points
    }

    pub(crate) fn set_fixed_points(&mut self, points: &[usize]) {
        assert!(points.iter().all(|&point| point < self.num_points()));
        self.fixed_points = points.to_vec();
    }

    pub(crate) fn energy_current(&mut self, energies: Option<&mut Vec<f64>>) -> f64 {
        // Match the primitive's early return, including the caller's existing
        // per-term buffer when a field has no contributions.
        if self.contribs.is_empty() {
            return 0.;
        }
        let pos = self.flat_positions();
        self.evaluate(&pos, None, energies)
    }

    pub(crate) fn gradient_current(&mut self) -> Vec<f64> {
        let pos = self.flat_positions();
        self.gradient_at(&pos)
    }

    pub(crate) fn gradient_at(&mut self, pos: &[f64]) -> Vec<f64> {
        assert_eq!(pos.len(), self.dimension() * self.num_points());
        let mut grad = vec![0.; pos.len()];
        self.evaluate(pos, Some(&mut grad), None);
        grad
    }

    pub(crate) fn minimize(&mut self, max_its: usize, force_tol: f64, energy_tol: f64) -> i32 {
        if self.contribs.is_empty() {
            return 0;
        }
        let dimension = self.dimension();
        let original = self.flat_positions();
        let pins = self.fixed_points.clone();
        let (status, mut positions) = {
            let evaluator = RefCell::new(&mut *self);
            bfgs_minimize(
                original.clone(),
                force_tol,
                |pos| evaluator.borrow_mut().evaluate(pos, None, None),
                |pos, grad| {
                    grad.fill(0.);
                    evaluator.borrow_mut().evaluate(pos, Some(&mut *grad), None);
                    for &atom in &pins {
                        grad[dimension * atom..dimension * (atom + 1)].fill(0.);
                    }
                    scale_gradient(grad)
                },
                0,
                None,
                energy_tol,
                max_its,
            )
        };
        for &atom in &pins {
            let range = dimension * atom..dimension * (atom + 1);
            positions[range.clone()].copy_from_slice(&original[range]);
        }
        self.owner.gather(&positions);
        self.owner.init_distance_matrix();
        status
    }

    fn flat_positions(&self) -> Vec<f64> {
        let mut pos = vec![0.; self.num_points() * self.dimension()];
        self.owner.scatter(&mut pos);
        pos
    }

    fn evaluate(
        &mut self,
        pos: &[f64],
        mut grad: Option<&mut [f64]>,
        mut energies: Option<&mut Vec<f64>>,
    ) -> f64 {
        let Self {
            owner, contribs, ..
        } = self;
        owner.gather(pos);
        owner.init_distance_matrix();
        // Soundness boundary: the owner contains no terms. Rebind each private
        // term from a fresh shared owner borrow after all mutable context work.
        // Evaluation only shares the context (its cache uses RefCell). No mutable
        // owner borrow occurs until every term call has ended. Clear the pointers
        // afterward; moves or coordinate edits cannot leave a usable stale pointer.
        // Contributions and owner references never escape this evaluator.
        let shared: &ForceField = owner;
        let pointer = shared as *const ForceField;
        for contrib in contribs.iter_mut() {
            contrib.set_force_field(pointer);
        }
        if let Some(energies) = &mut energies {
            energies.clear();
            energies.reserve(contribs.len());
        }
        let mut energy = 0.;
        for contrib in contribs.iter() {
            if let Some(grad) = grad.as_deref_mut() {
                contrib.get_grad(pos, grad);
            } else {
                let term = contrib.get_energy(pos);
                energy += term;
                if let Some(energies) = &mut energies {
                    energies.push(term);
                }
            }
        }
        for contrib in contribs.iter_mut() {
            contrib.set_force_field(std::ptr::null());
        }
        energy
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::distgeom::{DistViolationContribs, FourthDimContribs};
    use crate::chemistry::forcefield::uff::{
        bond_stretch::BondStretchContrib, params::AtomicParams,
    };

    fn four_dimensional_field() -> OwnedForceField {
        let mut primitive = ForceField::new(4);
        primitive.positions_mut().extend([
            ForceFieldVec3::new4(0., -0., 0., -0.),
            ForceFieldVec3::new4(2.5, 0.4, -0.3, 1.),
        ]);
        let mut distance = DistViolationContribs::new(&primitive);
        distance.add_contrib(0, 1, 1.5, 1.5, 1.);
        primitive.add_contrib(Box::new(distance));
        let mut fourth = FourthDimContribs::new(&primitive);
        fourth.add_contrib(0, 1.);
        fourth.add_contrib(1, 1.);
        primitive.add_contrib(Box::new(fourth));
        primitive.fixed_points_mut().push(0);
        OwnedForceField::from_force_field(primitive)
    }

    #[test]
    fn owned_evaluator_preserves_four_dimensional_terms_and_fixed_point_stride() {
        let mut field = vec![four_dimensional_field()].pop().unwrap();
        assert_eq!(field.dimension(), 4);
        assert_eq!(field.fixed_points(), [0]);
        let original = field.positions().to_vec();
        let mut energies = vec![999.];
        let before = field.energy_current(Some(&mut energies));
        assert_eq!(energies.len(), 2);
        assert_eq!(before, energies.iter().sum::<f64>());
        let gradient = field.gradient_current();
        assert_eq!(gradient.len(), 8);
        assert!(gradient[7].abs() > 1.);
        assert!(gradient.iter().all(|v| v.is_finite()));
        assert!(matches!(field.minimize(40, 1e-4, 1e-6), 0 | 1));
        assert!(field.energy_current(None) < before);
        for axis in 0..4 {
            assert_eq!(
                field.positions()[0][axis].to_bits(),
                original[0][axis].to_bits()
            );
        }
        assert!(field.positions()[1].w.abs() < original[1].w.abs());
    }

    #[test]
    fn owned_evaluator_preserves_source_fourth_dimension_gradient() {
        let mut field = four_dimensional_field();
        field.energy_current(None);
        let pos = vec![0., 0., 0., 0., 1.8, 0.2, 0., 0.4];
        let gradient = field.gradient_at(&pos);
        // The distance term is strictly inside its smooth upper-bound branch:
        // d² = 3.44 > 1.5². RDKit FourthDimContribs.h returns energy weight*w²
        // but source getGrad returns weight*w, preserving a half derivative.
        // Account for that source convention explicitly instead of changing it.
        assert!((gradient[7] - 0.7760987654320989).abs() < 1e-12);
        let h = 1e-5;
        for axis in 0..pos.len() {
            let mut plus = pos.clone();
            plus[axis] += h;
            let plus_energy = field.evaluate(&plus, None, None);
            let mut minus = pos.clone();
            minus[axis] -= h;
            let minus_energy = field.evaluate(&minus, None, None);
            let numeric = (plus_energy - minus_energy) / (2. * h);
            let source_correction = if axis % 4 == 3 { pos[axis] } else { 0. };
            let expected = gradient[axis] + source_correction;
            assert!(
                (numeric - expected).abs() < 1e-5,
                "4D axis {axis}: finite difference={numeric}, source gradient={}, source fourth correction={source_correction}",
                gradient[axis]
            );
        }
    }

    #[test]
    fn owned_evaluator_refreshes_trial_coordinates_before_gradient() {
        let mut primitive = ForceField::new(3);
        primitive.positions_mut().extend([
            ForceFieldVec3::new(0., 0., 0.),
            ForceFieldVec3::new(2.5, 0.4, -0.3),
        ]);
        // Only the source bond constructor's used parameters are needed here.
        let params = AtomicParams {
            r1: 0.757,
            theta0: 0.,
            x1: 0.,
            d1: 0.,
            zeta: 0.,
            z1: 1.912,
            v1: 0.,
            u1: 0.,
            gmp_xi: 5.343,
            gmp_hardness: 0.,
            gmp_radius: 0.,
        };
        let bond = BondStretchContrib::new(&primitive, 0, 1, 1., &params, &params).unwrap();
        primitive.add_contrib(Box::new(bond));
        let mut field = OwnedForceField::from_force_field(primitive);
        // UFF BondStretchContrib calls cached ForceField::distance. Populate the
        // cache at A, then request gradient at distinct, nonzero distance B.
        field.energy_current(None);
        let pos = vec![0., 0., 0., 1.8, 0.2, 0.4];
        let gradient = field.gradient_at(&pos);
        let h = 1e-5;
        for axis in 0..pos.len() {
            let mut plus = pos.clone();
            plus[axis] += h;
            let plus_energy = field.evaluate(&plus, None, None);
            let mut minus = pos.clone();
            minus[axis] -= h;
            let minus_energy = field.evaluate(&minus, None, None);
            let numeric = (plus_energy - minus_energy) / (2. * h);
            assert!(
                (numeric - gradient[axis]).abs() < 1e-5,
                "cached UFF axis {axis}: finite difference={numeric}, gradient={}",
                gradient[axis]
            );
        }
    }
}
