use crate::math::{AnyReal, AnyVector};
#[cfg(feature = "dim3")]
use crate::{
    bounding_volume,
    math::{Real, Vector},
};

/// Returns the index of the support point of a list of points.
pub fn support_point_id<V: AnyVector>(direction: V, points: &[V]) -> Option<usize> {
    let mut argmax = None;
    let mut max = -V::Real::MAX;

    for (id, pt) in points.iter().copied().enumerate() {
        let dot = direction.dot(pt);

        if dot > max {
            argmax = Some(id);
            max = dot;
        }
    }

    argmax
}

/// Returns the index of the support point of an indexed list of points.
pub fn indexed_support_point_id<I, V: AnyVector>(
    direction: V,
    points: &[V],
    idx: I,
) -> Option<usize>
where
    I: Iterator<Item = usize>,
{
    let mut argmax = None;
    let mut max = -V::Real::MAX;

    for i in idx.into_iter() {
        let dot = direction.dot(points[i]);

        if dot > max {
            argmax = Some(i);
            max = dot;
        }
    }

    argmax
}

/// Returns the number `n` such that `points[idx.nth(n)]` is the support point.
#[cfg(feature = "dim3")] // We only use this in 3D right now.
pub fn indexed_support_point_nth<I, V: AnyVector>(
    direction: V,
    points: &[V],
    idx: I,
) -> Option<usize>
where
    I: Iterator<Item = usize>,
{
    let mut argmax = None;
    let mut max = -V::Real::MAX;

    for (k, i) in idx.into_iter().enumerate() {
        let dot = direction.dot(points[i]);

        if dot > max {
            argmax = Some(k);
            max = dot;
        }
    }

    argmax
}

/// Scale and center the given set of point depending on their Aabb.
#[cfg(feature = "dim3")]
pub fn normalize(coords: &mut [Vector]) -> (Vector, Real) {
    let aabb = bounding_volume::details::local_point_cloud_aabb(&*coords);
    let diag = aabb.mins.distance(aabb.maxs);
    let center = aabb.center();

    for c in coords.iter_mut() {
        *c = (*c - center) / diag;
    }

    (center, diag)
}
