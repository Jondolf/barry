use barry2d::math::Vector2;
use barry2d::shape::ConvexPolygon;

fn main() {
    let points = vec![
        Vector2::new(-1.0f32, 1.0),
        Vector2::new(-0.5, -0.5),
        Vector2::new(0.5, -0.5),
        Vector2::new(1.0, 1.0),
    ];

    let convex = ConvexPolygon::from_convex_polyline(points).expect("Invalid convex polygon.");
    assert!(convex.points().len() == 4);
}
