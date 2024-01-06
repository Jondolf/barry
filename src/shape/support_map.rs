//! Traits for support mapping based shapes.

use crate::math::{Isometry, UnitVector, Vector};

/// Traits of convex shapes representable by a support mapping function.
///
/// # Parameters:
///   * V - type of the support mapping direction argument and of the returned point.
pub trait SupportMap {
    // Evaluates the support function of this shape.
    //
    // A support function is a function associating a vector to the shape point which maximizes
    // their dot product.
    fn local_support_point(&self, dir: Vector) -> Vector;

    /// Same as `self.local_support_point` except that `dir` is normalized.
    fn local_support_point_toward(&self, dir: UnitVector) -> Vector {
        self.local_support_point(*dir)
    }

    // Evaluates the support function of this shape transformed by `transform`.
    //
    // A support function is a function associating a vector to the shape point which maximizes
    // their dot product.
    fn support_point(&self, transform: Isometry, dir: Vector) -> Vector {
        let local_dir = transform.rotation.inverse() * dir;
        transform.transform_point(self.local_support_point(local_dir))
    }

    /// Same as `self.support_point` except that `dir` is normalized.
    fn support_point_toward(&self, transform: Isometry, dir: UnitVector) -> Vector {
        let local_dir = transform.rotation.inverse() * dir;
        transform.transform_point(self.local_support_point_toward(local_dir))
    }
}
