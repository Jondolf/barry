use barry3d::math::{Isometry3, Vector3};
use barry3d::query;
use barry3d::shape::{Ball, Cuboid};

fn main() {
    let cuboid = Cuboid::new(Vector3::new(1.0, 1.0, 1.0));
    let ball = Ball::new(1.0);
    let prediction = 1.0;

    let cuboid_pos = Isometry3::IDENTITY;
    let ball_pos_penetrating = Isometry3::from_xyz(1.0, 1.0, 1.0);
    let ball_pos_in_prediction = Isometry3::from_xyz(2.0, 2.0, 2.0);
    let ball_pos_too_far = Isometry3::from_xyz(3.0, 3.0, 3.0);

    let ctct_penetrating =
        query::contact(ball_pos_penetrating, &ball, cuboid_pos, &cuboid, prediction).unwrap();
    let ctct_in_prediction = query::contact(
        ball_pos_in_prediction,
        &ball,
        cuboid_pos,
        &cuboid,
        prediction,
    )
    .unwrap();
    let ctct_too_far =
        query::contact(ball_pos_too_far, &ball, cuboid_pos, &cuboid, prediction).unwrap();

    assert!(ctct_penetrating.unwrap().dist <= 0.0);
    assert!(ctct_in_prediction.unwrap().dist >= 0.0);
    assert_eq!(ctct_too_far, None);
}
