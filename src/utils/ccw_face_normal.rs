use crate::math::{UnitVector, Vector};
use bevy_math::primitives::InvalidDirectionError;

/// Computes the direction pointing toward the right-hand-side of an oriented segment.
///
/// Returns `None` if the segment is degenerate.
#[inline]
#[cfg(feature = "dim2")]
pub fn ccw_face_normal(pts: [Vector; 2]) -> Result<UnitVector, InvalidDirectionError> {
    let ab = pts[1] - pts[0];
    let res = Vector::new(ab[1], -ab[0]);
    UnitVector::new(res)
}

/// Computes the normal of a counter-clock-wise triangle.
///
/// Returns `None` if the triangle is degenerate.
#[inline]
#[cfg(feature = "dim3")]
pub fn ccw_face_normal(pts: [Vector; 3]) -> Result<UnitVector, InvalidDirectionError> {
    let ab = pts[1] - pts[0];
    let ac = pts[2] - pts[0];
    let res = ab.cross(ac);
    UnitVector::new(res)
}
