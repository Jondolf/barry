extern crate num_traits as num;

use barry2d::{math::Vector2, shape::ConvexPolygon};

fn main() {
    let points = [
        Vector2::new(-1.0f32, 1.0),
        Vector2::new(-0.5, -0.5),
        Vector2::new(0.0, 0.5),
        Vector2::new(0.5, -0.5),
        Vector2::new(1.0, 1.0),
    ];

    let convex = ConvexPolygon::from_convex_hull(&points).expect("Invalid convex polygon.");
    assert!(convex.points().len() == 4);
}
