use std::iter::IntoIterator;

use crate::bounding_volume::Aabb;
use crate::math::{Isometry, Vector, DIM};
use crate::shape::SupportMap;

/// Computes the [`Aabb`] of an [support mapped shape](SupportMap).
#[cfg(feature = "dim3")]
pub fn support_map_aabb<G>(m: Isometry, i: &G) -> Aabb
where
    G: SupportMap,
{
    let mut min = Vector::ZERO;
    let mut max = Vector::ZERO;
    let mut basis = Vector::ZERO;

    for d in 0..DIM {
        // FIXME: this could be further improved iterating on `m`'s columns, and passing
        // Id as the transformation matrix.
        basis[d] = 1.0;
        max[d] = i.support_point(m, basis)[d];

        basis[d] = -1.0;
        min[d] = i.support_point(m, basis)[d];

        basis[d] = 0.0;
    }

    Aabb::new(min, max)
}

/// Computes the [`Aabb`] of an [support mapped shape](SupportMap).
pub fn local_support_map_aabb<G>(i: &G) -> Aabb
where
    G: SupportMap,
{
    let mut min = Vector::ZERO;
    let mut max = Vector::ZERO;
    let mut basis = Vector::ZERO;

    for d in 0..DIM {
        // FIXME: this could be further improved iterating on `m`'s columns, and passing
        // Id as the transformation matrix.
        basis[d] = 1.0;
        max[d] = i.local_support_point(basis)[d];

        basis[d] = -1.0;
        min[d] = i.local_support_point(basis)[d];

        basis[d] = 0.0;
    }

    Aabb::new(min, max)
}

/// Computes the [`Aabb`] of a set of points transformed by `m`.
pub fn point_cloud_aabb<'a, I>(m: Isometry, pts: I) -> Aabb
where
    I: IntoIterator<Item = &'a Vector>,
{
    let mut it = pts.into_iter().copied();

    let p0 = it.next().expect(
        "Point cloud Aabb construction: the input iterator should yield at least one point.",
    );
    let wp0 = m.transform_point(p0);
    let mut min = wp0;
    let mut max = wp0;

    for pt in it {
        let wpt = m * pt;
        min = min.min(wpt);
        max = max.max(wpt);
    }

    Aabb::new(min, max)
}

/// Computes the [`Aabb`] of a set of points.
pub fn local_point_cloud_aabb<'a, I>(pts: I) -> Aabb
where
    I: IntoIterator<Item = &'a Vector>,
{
    let mut it = pts.into_iter().copied();

    let p0 = it.next().expect(
        "Point cloud Aabb construction: the input iterator should yield at least one point.",
    );
    let mut min = p0;
    let mut max = p0;

    for pt in it {
        min = min.min(pt);
        max = max.max(pt);
    }

    Aabb::new(min, max)
}
