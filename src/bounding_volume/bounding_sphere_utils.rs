use crate::math::{Real, Vector};
use crate::utils;

/// Computes the bounding sphere of a set of point, given its center.
// FIXME: return a bounding sphere?
#[inline]
pub fn point_cloud_bounding_sphere_with_center(pts: &[Vector], center: Vector) -> (Vector, Real) {
    let mut sqradius = 0.0;

    for pt in pts.iter() {
        let distance_squared = pt.distance_squared(center);

        if distance_squared > sqradius {
            sqradius = distance_squared
        }
    }

    (center, sqradius.sqrt())
}

/// Computes a bounding sphere of the specified set of point.
// FIXME: return a bounding sphere?
#[inline]
pub fn point_cloud_bounding_sphere(pts: &[Vector]) -> (Vector, Real) {
    point_cloud_bounding_sphere_with_center(pts, utils::center(pts))
}
