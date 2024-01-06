use crate::bounding_volume::SimdAabb;
use crate::math::{SimdVector, Vector, SIMD_WIDTH};
use crate::partitioning::{SimdVisitStatus, SimdVisitor};
use simba::simd::SimdBool as _;
use std::marker::PhantomData;

// FIXME: add a point cost fn.

/// Spatial partitioning structure visitor collecting nodes that may contain a given point.
pub struct PointIntersectionsVisitor<'a, T, F> {
    simd_point: SimdVector,
    /// Callback executed for each leaf which Aabb contains `self.point`.
    callback: &'a mut F,
    _phantom: PhantomData<T>,
}

impl<'a, T, F> PointIntersectionsVisitor<'a, T, F>
where
    F: FnMut(&T) -> bool,
{
    /// Creates a new `PointIntersectionsVisitor`.
    #[inline]
    pub fn new(point: Vector, callback: &'a mut F) -> PointIntersectionsVisitor<'a, T, F> {
        PointIntersectionsVisitor {
            simd_point: SimdVector::splat(point),
            callback,
            _phantom: PhantomData,
        }
    }
}

impl<'a, T, F> SimdVisitor<T, SimdAabb> for PointIntersectionsVisitor<'a, T, F>
where
    F: FnMut(&T) -> bool,
{
    #[inline]
    fn visit(&mut self, bv: &SimdAabb, b: Option<[Option<&T>; SIMD_WIDTH]>) -> SimdVisitStatus {
        let mask = bv.contains_local_point(self.simd_point);

        if let Some(data) = b {
            let bitmask = mask.bitmask();
            for ii in 0..SIMD_WIDTH {
                if (bitmask & (1 << ii)) != 0
                    && data[ii].is_some()
                    && !(self.callback)(data[ii].unwrap())
                {
                    return SimdVisitStatus::ExitEarly;
                }
            }
        }

        SimdVisitStatus::MaybeContinue(mask)
    }
}
