use crate::math::{self, Real, Vector2};
use crate::shape::Capsule;
use crate::transformation::utils;

impl Capsule {
    /// Discretize the boundary of this capsule as a polygonal line.
    pub fn to_polyline(&self, nsubdiv: u32) -> Vec<Vector2> {
        let pi = math::real_consts::PI;
        let dtheta = pi / (nsubdiv as Real);

        let mut points: Vec<Vector2> = Vec::with_capacity(nsubdiv as usize);

        utils::push_xy_arc(self.radius, nsubdiv, dtheta, &mut points);

        let npoints = points.len();

        for i in 0..npoints {
            let new_point = points[i] + Vector2::new(0.0, self.half_height());

            points.push(-new_point);
            points[i] = new_point;
        }

        utils::transformed(points, self.canonical_transform())
    }
}
