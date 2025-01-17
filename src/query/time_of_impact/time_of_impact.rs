use crate::math::{Isometry, Real, UnitVector, Vector};
use crate::query::{DefaultQueryDispatcher, QueryDispatcher, Unsupported};
use crate::shape::Shape;

/// The status of the time-of-impact computation algorithm.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum TOIStatus {
    /// The TOI algorithm ran out of iterations before achieving convergence.
    ///
    /// The content of the `TOI` will still be a conservative approximation of the actual result so
    /// it is often fine to interpret this case as a success.
    OutOfIterations,
    /// The TOI algorithm converged successfully.
    Converged,
    /// Something went wrong during the TOI computation, likely due to numerical instabilities.
    ///
    /// The content of the `TOI` will still be a conservative approximation of the actual result so
    /// it is often fine to interpret this case as a success.
    Failed,
    /// The two shape already overlap at the time 0.
    ///
    /// The witness points and normals provided by the `TOI` will have undefined values.
    Penetrating,
}

/// The result of a time-of-impact (TOI) computation.
#[derive(Copy, Clone, Debug)]
pub struct TOI {
    /// The time at which the objects touch.
    pub toi: Real,
    /// The local-space closest point on the first shape at the time of impact.
    ///
    /// Undefined if `status` is `Penetrating`.
    pub witness1: Vector,
    /// The local-space closest point on the second shape at the time of impact.
    ///
    /// Undefined if `status` is `Penetrating`.
    pub witness2: Vector,
    /// The local-space outward normal on the first shape at the time of impact.
    ///
    /// Undefined if `status` is `Penetrating`.
    pub normal1: UnitVector,
    /// The local-space outward normal on the second shape at the time of impact.
    ///
    /// Undefined if `status` is `Penetrating`.
    pub normal2: UnitVector,
    /// The way the time-of-impact computation algorithm terminated.
    pub status: TOIStatus,
}

impl TOI {
    /// Swaps every data of this TOI result such that the role of both shapes are inverted.
    ///
    /// In practice, this makes it so that `self.witness1` and `self.normal1` become `self.witness2` and `self.normal2` and vice-versa.
    pub fn swapped(self) -> Self {
        Self {
            toi: self.toi,
            witness1: self.witness2,
            witness2: self.witness1,
            normal1: self.normal2,
            normal2: self.normal1,
            status: self.status,
        }
    }

    /// Transform `self.witness1` and `self.normal1` by `pos`.
    pub fn transform1_by(&self, pos: Isometry) -> Self {
        Self {
            toi: self.toi,
            witness1: pos * self.witness1,
            witness2: self.witness2,
            normal1: pos * self.normal1,
            normal2: self.normal2,
            status: self.status,
        }
    }
}

/// Computes the smallest time when two shapes under translational movement are separated by a
/// distance smaller or equal to `distance`.
///
/// Returns `0.0` if the objects are touching or penetrating.
pub fn time_of_impact(
    pos1: Isometry,
    vel1: Vector,
    g1: &dyn Shape,
    pos2: Isometry,
    vel2: Vector,
    g2: &dyn Shape,
    max_toi: Real,
    stop_at_penetration: bool,
) -> Result<Option<TOI>, Unsupported> {
    let pos12 = pos1.inv_mul(pos2);
    let vel12 = pos1.rotation.inverse() * (vel2 - vel1);
    DefaultQueryDispatcher.time_of_impact(pos12, vel12, g1, g2, max_toi, stop_at_penetration)
}
