use crate::math::{Isometry, Real, UnitVector, Vector};
use crate::query::{self, Ray, TOIStatus, TOI};
use crate::shape::Ball;
use num::Zero;

/// Time Of Impact of two balls under translational movement.
#[inline]
pub fn time_of_impact_ball_ball(
    pos12: Isometry,
    vel12: Vector,
    b1: &Ball,
    b2: &Ball,
    max_toi: Real,
) -> Option<TOI> {
    let rsum = b1.radius + b2.radius;
    let radius = rsum;
    let center = Vector::from(-pos12.translation);
    let ray = Ray::new(Vector::ZERO, vel12);

    if let (inside, Some(toi)) = query::details::ray_toi_with_ball(center, radius, &ray, true) {
        if toi > max_toi {
            return None;
        }

        let dpt = ray.point_at(toi) - center;
        let normal1;
        let normal2;
        let witness1;
        let witness2;

        if radius.is_zero() {
            normal1 = UnitVector::X;
            normal2 = pos12.rotation.inverse() * -UnitVector::X;
            witness1 = Vector::ZERO;
            witness2 = Vector::ZERO;
        } else {
            normal1 = UnitVector::new_unchecked(dpt / radius);
            normal2 = pos12.rotation.inverse() * -normal1;
            witness1 = Vector::from(*normal1 * b1.radius);
            witness2 = Vector::from(*normal2 * b2.radius);
        }

        let status = if inside && center.length_squared() < rsum * rsum {
            TOIStatus::Penetrating
        } else {
            TOIStatus::Converged
        };

        Some(TOI {
            toi,
            normal1,
            normal2,
            witness1,
            witness2,
            status,
        })
    } else {
        None
    }
}
