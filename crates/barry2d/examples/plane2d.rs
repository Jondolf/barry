use barry2d::math::UnitVector2;
use barry2d::shape::HalfSpace;

fn main() {
    let halfspace = HalfSpace::new(UnitVector2::Y);

    assert!(halfspace.normal.x == 0.0);
    assert!(halfspace.normal.y == 1.0);
}
