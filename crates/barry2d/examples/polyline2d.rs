use barry2d::math::Vector2;
use barry2d::shape::Polyline;

fn main() {
    let points = vec![
        Vector2::new(0.0, 1.0),
        Vector2::new(-1.0, -1.0),
        Vector2::new(0.0, -0.5),
        Vector2::new(1.0, -1.0),
    ];

    let indices = vec![
        [0, 1],
        [1, 2],
        [2, 3],
        [3, 0], // This forms a loop.
    ];

    // Build the polyline.
    let polyline = Polyline::new(points, Some(indices));

    assert_eq!(polyline.vertices().len(), 4);
}
