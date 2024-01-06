use barry3d::math::UnitVector3;
use barry3d::shape::HalfSpace;

fn main() {
    let halfspace = HalfSpace::new(UnitVector3::Y);

    assert!(halfspace.normal.x == 0.0);
    assert!(halfspace.normal.y == 1.0);
    assert!(halfspace.normal.z == 0.0);
}
