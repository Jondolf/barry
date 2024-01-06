use crate::math::UnitVector;
use crate::math::{Isometry, Real};
use crate::query::Contact;
use crate::shape::Ball;

/// Contact between balls.
#[inline]
pub fn contact_ball_ball(
    pos12: Isometry,
    b1: &Ball,
    b2: &Ball,
    prediction: Real,
) -> Option<Contact> {
    let r1 = b1.radius;
    let r2 = b2.radius;
    let center2_1 = pos12.translation;
    let distance_squared = center2_1.length();
    let sum_radius = r1 + r2;
    let sum_radius_with_error = sum_radius + prediction;

    if distance_squared < sum_radius_with_error * sum_radius_with_error {
        let normal1 = if distance_squared != 0.0 {
            UnitVector::new(center2_1).unwrap()
        } else {
            UnitVector::X
        };
        let normal2 = -(pos12.rotation.inverse() * normal1);
        let point1 = *normal1 * r1;
        let point2 = *normal2 * r2;

        Some(Contact::new(
            point1,
            point2,
            normal1,
            normal2,
            distance_squared.sqrt() - sum_radius,
        ))
    } else {
        None
    }
}
