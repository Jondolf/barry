use barry3d::math::{Isometry3, Vector3};
use barry3d::query::{Ray, RayCast};
use barry3d::shape::Cuboid;

fn main() {
    let cuboid = Cuboid::new(Vector3::new(1.0, 2.0, 1.0));
    let ray_inside = Ray::new(Vector3::ZERO, Vector3::Y);
    let ray_miss = Ray::new(Vector3::new(2.0, 2.0, 2.0), Vector3::new(1.0, 1.0, 1.0));

    // Solid cast.
    assert_eq!(
        cuboid
            .cast_ray(Isometry3::IDENTITY, &ray_inside, std::f32::MAX, true)
            .unwrap(),
        0.0
    );

    // Non-solid cast.
    assert_eq!(
        cuboid
            .cast_ray(Isometry3::IDENTITY, &ray_inside, std::f32::MAX, false)
            .unwrap(),
        2.0
    );

    // The other ray does not intersect this shape.
    assert!(cuboid
        .cast_ray(Isometry3::IDENTITY, &ray_miss, std::f32::MAX, false)
        .is_none());
    assert!(cuboid
        .cast_ray(Isometry3::IDENTITY, &ray_miss, std::f32::MAX, true)
        .is_none());
}
