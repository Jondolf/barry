#![allow(dead_code)] // TODO: remove this

use crate::mass_properties::MassProperties;
use crate::math::{Real, Vector};
use crate::shape::Triangle;

impl MassProperties {
    /// Computes the mass properties of a convex polygon.
    pub fn from_convex_polygon(density: Real, vertices: &[Vector]) -> MassProperties {
        let (area, com) = convex_polygon_area_and_center_of_mass(vertices);

        if area == 0.0 {
            return MassProperties::new(com, 0.0, 0.0);
        }

        let mut itot = 0.0;

        let mut iterpeek = vertices.iter().peekable();
        let first_element = *iterpeek.peek().unwrap(); // store first element to close the cycle in the end with unwrap_or
        while let Some(elem) = iterpeek.next() {
            let triangle = Triangle::new(com, *elem, **iterpeek.peek().unwrap_or(&first_element));
            let area = triangle.area();
            let ipart = triangle.unit_angular_inertia();
            itot += ipart * area;
        }

        Self::new(com, area * density, itot * density)
    }
}

/// Computes the area and center-of-mass of a convex polygon.
pub fn convex_polygon_area_and_center_of_mass(convex_polygon: &[Vector]) -> (Real, Vector) {
    let geometric_center =
        convex_polygon.iter().fold(Vector::ZERO, |e1, e2| e1 + *e2) / convex_polygon.len() as Real;
    let mut res = Vector::ZERO;
    let mut areasum = 0.0;

    let mut iterpeek = convex_polygon.iter().copied().peekable();
    let first_element = *iterpeek.peek().unwrap(); // Stores first element to close the cycle in the end with unwrap_or.
    while let Some(elem) = iterpeek.next() {
        let (a, b, c) = (
            elem,
            iterpeek.peek().copied().unwrap_or(first_element),
            geometric_center,
        );
        let area = Triangle::new(a, b, c).area();
        let center = (a + b + c) / 3.0;

        res += center * area;
        areasum += area;
    }

    if areasum == 0.0 {
        (areasum, geometric_center)
    } else {
        (areasum, res / areasum)
    }
}
