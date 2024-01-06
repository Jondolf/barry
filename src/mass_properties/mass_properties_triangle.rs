use crate::mass_properties::MassProperties;
use crate::math::{Real, Vector};
use crate::shape::Triangle;

impl MassProperties {
    /// Computes the mass properties of a triangle.
    pub fn from_triangle(density: Real, a: Vector, b: Vector, c: Vector) -> MassProperties {
        let triangle = Triangle::new(a, b, c);
        let area = triangle.area();
        let com = triangle.center();

        if area == 0.0 {
            return MassProperties::new(com, 0.0, 0.0);
        }

        let ipart = triangle.unit_angular_inertia();

        Self::new(com, area * density, ipart * area * density)
    }
}
