use barry2d::math::{Isometry2, Vector2};
use barry2d::query;
use barry2d::shape::{Ball, Cuboid};

fn main() {
    let cuboid = Cuboid::new(Vector2::new(1.0, 1.0));
    let ball = Ball::new(1.0);

    let cuboid_pos = Isometry2::IDENTITY;
    let ball_pos_intersecting = Isometry2::from_xy(1.0, 1.0);
    let ball_pos_disjoint = Isometry2::from_xy(3.0, 3.0);

    assert!(query::intersection_test(ball_pos_intersecting, &ball, cuboid_pos, &cuboid).unwrap());
    assert!(!query::intersection_test(ball_pos_disjoint, &ball, cuboid_pos, &cuboid).unwrap());
}
