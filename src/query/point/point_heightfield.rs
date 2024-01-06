use crate::bounding_volume::Aabb;
use crate::math::{Real, Vector};
use crate::query::{PointProjection, PointQuery, PointQueryWithLocation};
use crate::shape::{FeatureId, GenericHeightField, HeightFieldStorage, TrianglePointLocation};

impl<Storage: HeightFieldStorage> PointQuery for GenericHeightField<Storage> {
    fn project_local_point_with_max_dist(
        &self,
        point: Vector,
        solid: bool,
        max_dist: Real,
    ) -> Option<PointProjection> {
        let aabb = Aabb::new(
            point - Vector::splat(max_dist),
            point + Vector::splat(max_dist),
        );
        let mut sq_smallest_dist = Real::MAX;
        let mut best_proj = None;

        self.map_elements_in_local_aabb(&aabb, &mut |_, triangle| {
            let proj = triangle.project_local_point(point, solid);
            let sq_dist = point.distance_squared(proj.point);

            if sq_dist < sq_smallest_dist {
                sq_smallest_dist = sq_dist;

                if sq_dist.sqrt() <= max_dist {
                    best_proj = Some(proj);
                }
            }
        });

        best_proj
    }

    #[inline]
    fn project_local_point(&self, point: Vector, _: bool) -> PointProjection {
        let mut smallest_dist = Real::MAX;
        let mut best_proj = PointProjection::new(false, point);

        #[cfg(feature = "dim2")]
        let iter = self.segments();
        #[cfg(feature = "dim3")]
        let iter = self.triangles();
        for elt in iter {
            let proj = elt.project_local_point(point, false);
            let dist = point.distance_squared(proj.point);

            if dist < smallest_dist {
                smallest_dist = dist;
                best_proj = proj;
            }
        }

        best_proj
    }

    #[inline]
    fn project_local_point_and_get_feature(&self, point: Vector) -> (PointProjection, FeatureId) {
        // FIXME: compute the feature properly.
        (self.project_local_point(point, false), FeatureId::Unknown)
    }

    // FIXME: implement distance_to_point too?

    #[inline]
    fn contains_local_point(&self, _point: Vector) -> bool {
        false
    }
}

impl<Storage: HeightFieldStorage> PointQueryWithLocation for GenericHeightField<Storage> {
    type Location = (usize, TrianglePointLocation);

    #[inline]
    fn project_local_point_and_get_location(
        &self,
        _point: Vector,
        _: bool,
    ) -> (PointProjection, Self::Location) {
        unimplemented!()
    }
}
