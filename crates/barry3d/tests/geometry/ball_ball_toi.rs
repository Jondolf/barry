// Issue #35

use na::{self, Isometry3, Vector3};
use barry3d::math::Real;
use barry3d::query;
use barry3d::shape::Ball;

#[test]
fn test_ball_ball_toi() {
    let b = Ball::new(0.5);
    let m1 = Isometry3::identity();
    let m2 = Isometry3::translation(0.0, 10.0, 0.0);
    let vel1 = Vector3::new(0.0, 10.0, 0.0);
    let vel2 = Vector3::zeros();

    let cast = query::time_of_impact(&m1, &vel1, &b, &m2, &vel2, &b, Real::MAX, true).unwrap();

    assert_eq!(cast.unwrap().toi, 0.9);
}
