// Issue #123

use barry3d::math::{Isometry3, Vector3};
use barry3d::query;
use barry3d::shape::{Ball, Triangle};

#[test]
fn ball_triangle_toi_infinite_loop_issue() {
    let b = Ball::new(0.375f32);
    let t = Triangle::new(
        Vector3::new(0.5, -0.5, 0.0),
        Vector3::new(-0.5, -0.5, 0.0),
        Vector3::new(-0.5, 0.5, 0.0),
    );

    let m1 = Isometry3::from_xyz(0.0, 0.0, 0.0);
    let m2 = Isometry3::from_xyz(11.5, 5.5, 0.0);
    let vel1 = Vector3::new(0.0, 0.000000000000000000000000000000000000000006925, 0.0);
    let vel2 = Vector3::ZERO;

    let cast = query::time_of_impact(m1, vel1, &b, m2, vel2, &t, std::f32::MAX, true).unwrap();

    println!("TOI: {:?}", cast);
    assert!(cast.is_none()); // The provided velocity is too small.
}
