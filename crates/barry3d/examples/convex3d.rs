use barry3d::math::Vector3;
use barry3d::shape::ConvexPolyhedron;

fn main() {
    let points = [
        Vector3::new(0.0f32, 0.0, 1.0),
        Vector3::new(0.0, 0.0, -1.0),
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(0.0, -1.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(-1.0, 0.0, 0.0),
        Vector3::new(0.0, 0.0, 0.0),
    ];

    let convex = ConvexPolyhedron::from_convex_hull(&points).expect("Invalid convex shape.");
    convex.check_geometry();
}
