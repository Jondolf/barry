use barry2d::math::Vector2;
use barry2d::shape::Cuboid;

fn main() {
    let cuboid = Cuboid::new(Vector2::new(2.0f32, 1.0));

    assert!(cuboid.half_extents.x == 2.0);
    assert!(cuboid.half_extents.y == 1.0);
}
