// Issue #35

use barry2d::math::Real;
use barry2d::math::{Isometry2, Vector2};
use barry2d::query;
use barry2d::shape::Ball;

#[test]
fn test_ball_ball_toi() {
    let b = Ball::new(0.5);
    let m1 = Isometry2::IDENTITY;
    let m2 = Isometry2::from_xy(0.0, 10.0);
    let v1 = Vector2::new(0.0, 10.0);
    let v2 = Vector2::ZERO;

    let cast = query::time_of_impact(m1, v1, &b, m2, v2, &b, Real::MAX, true).unwrap();

    assert_eq!(cast.unwrap().toi, 0.9);
}
