use barry3d::math::{Isometry3, Vector3};
use barry3d::query;
use barry3d::shape::{Ball, Cuboid};

fn main() {
    let cuboid = Cuboid::new(Vector3::new(1.0, 1.0, 1.0));
    let ball = Ball::new(1.0);
    let cuboid_pos = Isometry3::IDENTITY;
    let ball_pos_intersecting = Isometry3::from_xyz(1.0, 1.0, 1.0);
    let ball_pos_disjoint = Isometry3::from_xyz(3.0, 3.0, 3.0);

    let intersecting =
        query::intersection_test(ball_pos_intersecting, &ball, cuboid_pos, &cuboid).unwrap();
    let not_intersecting =
        !query::intersection_test(ball_pos_disjoint, &ball, cuboid_pos, &cuboid).unwrap();

    assert!(intersecting);
    assert!(not_intersecting);
}
