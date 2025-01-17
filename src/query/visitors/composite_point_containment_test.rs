use crate::bounding_volume::SimdAabb;
use crate::math::{SimdVector, Vector, SIMD_WIDTH};
use crate::partitioning::{SimdVisitStatus, SimdVisitor};
use crate::query::point::point_query::PointQuery;
use crate::shape::TypedSimdCompositeShape;
use crate::utils::IsometryOpt;
use simba::simd::SimdBool as _;

/// Visitor for checking if a composite shape contains a specific point.
pub struct CompositePointContainmentTest<'a, S: 'a> {
    /// The composite shape on which the point containment test should be performed.
    pub shape: &'a S,
    /// The point to be tested.
    pub point: Vector,
    /// A traversal will set this to `true` if the point is inside of `self.shape`.
    pub found: bool,
}

impl<'a, S> CompositePointContainmentTest<'a, S> {
    /// Creates a new visitor for the testing containment of the given `point`
    /// into the given `shape`.
    pub fn new(shape: &'a S, point: Vector) -> Self {
        Self {
            shape,
            point,
            found: false,
        }
    }
}

impl<'a, S: TypedSimdCompositeShape> SimdVisitor<S::PartId, SimdAabb>
    for CompositePointContainmentTest<'a, S>
{
    #[inline]
    fn visit(
        &mut self,
        bv: &SimdAabb,
        b: Option<[Option<&S::PartId>; SIMD_WIDTH]>,
    ) -> SimdVisitStatus {
        let point = SimdVector::splat(self.point);
        let mask = bv.contains_local_point(point);

        if let Some(data) = b {
            let bitmask = mask.bitmask();

            for ii in 0..SIMD_WIDTH {
                if (bitmask & (1 << ii)) != 0 && data[ii].is_some() {
                    self.shape
                        .map_typed_part_at(*data[ii].unwrap(), |part_pos, obj| {
                            if obj
                                .contains_local_point(part_pos.inverse_transform_point(self.point))
                            {
                                self.found = true;
                            }
                        });

                    if self.found {
                        return SimdVisitStatus::ExitEarly;
                    }
                }
            }
        }

        SimdVisitStatus::MaybeContinue(mask)
    }
}
