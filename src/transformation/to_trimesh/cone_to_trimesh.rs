use crate::math::{self, Real, Vector3};
use crate::shape::Cone;
use crate::transformation::utils;

impl Cone {
    /// Discretize the boundary of this cone as a triangle-mesh.
    pub fn to_trimesh(&self, nsubdiv: u32) -> (Vec<Vector3>, Vec<[u32; 3]>) {
        let diameter = self.radius * 2.0;
        let height = self.half_height * 2.0;
        let scale = Vector3::new(diameter, height, diameter);
        let (vtx, idx) = unit_cone(nsubdiv);
        (utils::scaled(vtx, scale), idx)
    }
}

/// Generates a cone with unit height and diameter.
fn unit_cone(nsubdiv: u32) -> (Vec<Vector3>, Vec<[u32; 3]>) {
    let dtheta = math::real_consts::TAU / (nsubdiv as Real);
    let mut coords = Vec::new();
    let mut indices = Vec::new();

    utils::push_circle(0.5, nsubdiv, dtheta, -0.5, &mut coords);

    coords.push(Vector3::new(0.0, 0.5, 0.0));

    utils::push_degenerate_top_ring_indices(0, coords.len() as u32 - 1, nsubdiv, &mut indices);
    utils::push_filled_circle_indices(0, nsubdiv, &mut indices);

    (coords, indices)
}
