#[macro_use]
extern crate approx; // for relative_eq!

use barry2d::math::{Isometry2, Vector2};
use barry2d::query;
use barry2d::shape::{Ball, Cuboid};

fn main() {
    let cuboid = Cuboid::new(Vector2::new(1.0, 1.0));
    let ball = Ball::new(1.0);

    let cuboid_pos = Isometry2::IDENTITY;
    let ball_pos_intersecting = Isometry2::from_xy(0.0, 1.0);
    let ball_pos_disjoint = Isometry2::from_xy(0.0, 3.0);

    let dist_intersecting =
        query::distance(ball_pos_intersecting, &ball, cuboid_pos, &cuboid).unwrap();
    let dist_disjoint = query::distance(ball_pos_disjoint, &ball, cuboid_pos, &cuboid).unwrap();

    assert_eq!(dist_intersecting, 0.0);
    assert!(relative_eq!(dist_disjoint, 1.0, epsilon = 1.0e-7));
}
