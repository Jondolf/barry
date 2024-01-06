use barry2d::math::{Isometry2, Vector2};
use barry2d::query::PointQuery;
use barry2d::shape::Cuboid;

fn main() {
    let cuboid = Cuboid::new(Vector2::new(1.0, 2.0));
    let pt_inside = Vector2::ZERO;
    let pt_outside = Vector2::new(2.0, 2.0);

    // Solid projection.
    assert_eq!(
        cuboid.distance_to_point(Isometry2::IDENTITY, pt_inside, true),
        0.0
    );

    // Non-solid projection.
    assert_eq!(
        cuboid.distance_to_point(Isometry2::IDENTITY, pt_inside, false),
        -1.0
    );

    // The other point is outside of the cuboid so the `solid` flag has no effect.
    assert_eq!(
        cuboid.distance_to_point(Isometry2::IDENTITY, pt_outside, false),
        1.0
    );
    assert_eq!(
        cuboid.distance_to_point(Isometry2::IDENTITY, pt_outside, true),
        1.0
    );
}
