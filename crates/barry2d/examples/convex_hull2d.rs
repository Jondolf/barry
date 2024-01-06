use barry2d::math::Vector2;
use barry2d::transformation;

fn main() {
    let points = [
        Vector2::new(0.77705324, 0.05374551),
        Vector2::new(0.35096353, 0.9873069),
        Vector2::new(0.09537989, 0.44411153),
        Vector2::new(0.108208835, 0.72445065),
        Vector2::new(0.7661844, 0.86163324),
        Vector2::new(0.5185994, 0.66594696),
        Vector2::new(0.768981, 0.23657233),
        Vector2::new(0.058607936, 0.09037298),
        Vector2::new(0.8818559, 0.3804205),
        Vector2::new(0.9571466, 0.17664945),
    ];

    let _ = transformation::convex_hull(&points[..]);
}
