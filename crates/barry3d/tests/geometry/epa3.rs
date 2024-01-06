use barry3d::math::{Isometry3, UnitVector3, Vector3};
use barry3d::query;
use barry3d::shape::Cuboid;

#[test]
#[allow(non_snake_case)]
fn cuboid_cuboid_EPA() {
    let c = Cuboid::new(Vector3::new(2.0, 1.0, 1.0));
    let m1 = Isometry3::from_xyz(3.5, 0.0, 0.0);
    let m2 = Isometry3::IDENTITY;

    let res = query::details::contact_support_map_support_map(m1.inv_mul(m2), &c, &c, 10.0)
        .expect("Penetration not found.");
    assert_eq!(res.dist, -0.5);
    assert_eq!(res.normal1, -UnitVector3::X);

    let m1 = Isometry3::from_xyz(0.0, 0.2, 0.0);
    let res = query::details::contact_support_map_support_map(m1.inv_mul(m2), &c, &c, 10.0)
        .expect("Penetration not found.");
    assert_eq!(res.dist, -1.8);
    assert_eq!(res.normal1, -UnitVector3::Y);
}
