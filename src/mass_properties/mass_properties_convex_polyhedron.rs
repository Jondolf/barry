use crate::mass_properties::MassProperties;
use crate::math::{Real, Vector, DIM};

impl MassProperties {
    /// Computes the mass properties of a convex polyhedron.
    pub fn from_convex_polyhedron(
        density: Real,
        vertices: &[Vector],
        indices: &[[u32; DIM]],
    ) -> MassProperties {
        Self::from_trimesh(density, vertices, indices)
    }
}
