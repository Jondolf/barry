use barry3d::math::{Isometry3, Vector3};
use barry3d::query::PointQuery;
use barry3d::shape::Cuboid;

fn main() {
    let cuboid = Cuboid::new(Vector3::new(1.0, 2.0, 2.0));
    let pt_inside = Vector3::ZERO;
    let pt_outside = Vector3::new(2.0, 2.0, 2.0);

    // Solid projection.
    assert_eq!(
        cuboid.distance_to_point(Isometry3::IDENTITY, pt_inside, true),
        0.0
    );

    // Non-solid projection.
    assert_eq!(
        cuboid.distance_to_point(Isometry3::IDENTITY, pt_inside, false),
        -1.0
    );

    // The other point is outside of the cuboid so the `solid` flag has no effect.
    assert_eq!(
        cuboid.distance_to_point(Isometry3::IDENTITY, pt_outside, false),
        1.0
    );
    assert_eq!(
        cuboid.distance_to_point(Isometry3::IDENTITY, pt_outside, true),
        1.0
    );
}