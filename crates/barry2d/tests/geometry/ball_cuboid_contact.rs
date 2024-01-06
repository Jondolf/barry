use barry2d::math::{Isometry2, Vector2};
use barry2d::query;
use barry2d::shape::{Ball, Cuboid};
#[cfg(feature = "improved_fixed_point_support")]
use simba::scalar::FixedI40F24;

#[test]
fn test_ball_cuboid_query_contact() {
    let cuboid = Cuboid::new(Vector2::new(0.5, 0.5));
    let cuboid_pos = Isometry2::from_xy(0.0, 4.0);
    let ball = Ball::new(0.5);
    let ball_pos = Isometry2::from_xy(0.0517938, 3.05178815);
    let ct = query::contact(cuboid_pos, &cuboid, ball_pos, &cuboid, 0.0).unwrap();
    assert!(ct.is_none());
    let ct = query::contact(ball_pos, &ball, cuboid_pos, &cuboid, 0.0).unwrap();
    assert!(ct.is_some());
}

#[test]
fn test_false_negative() {
    let contact = query::contact(
        Isometry2::from_xy(1.0, 1.0),
        &Ball::new(1.0),
        Isometry2::IDENTITY,
        &Cuboid::new(Vector2::new(1.0, 1.0)),
        1.0,
    )
    .unwrap()
    .unwrap();

    assert!(contact.dist < 0.0);
}
