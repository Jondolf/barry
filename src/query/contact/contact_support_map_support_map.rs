use crate::math::{Isometry, Real, UnitVector};
use crate::query::epa::EPA;
use crate::query::gjk::{self, CSOPoint, GJKResult, VoronoiSimplex};
use crate::query::Contact;
use crate::shape::SupportMap;

/// Contact between support-mapped shapes (`Cuboid`, `ConvexHull`, etc.)
pub fn contact_support_map_support_map<G1: ?Sized, G2: ?Sized>(
    pos12: Isometry,
    g1: &G1,
    g2: &G2,
    prediction: Real,
) -> Option<Contact>
where
    G1: SupportMap,
    G2: SupportMap,
{
    let simplex = &mut VoronoiSimplex::new();
    match contact_support_map_support_map_with_params(pos12, g1, g2, prediction, simplex, None) {
        GJKResult::ClosestPoints(point1, point2_1, normal1) => {
            let dist = (point2_1 - point1).dot(*normal1);
            let point2 = pos12.inverse_transform_point(point2_1);
            let normal2 = pos12.rotation.inverse() * -normal1;
            Some(Contact::new(point1, point2, normal1, normal2, dist))
        }
        GJKResult::NoIntersection(_) => None,
        GJKResult::Intersection => unreachable!(),
        GJKResult::Proximity(_) => unreachable!(),
    }
}

/// Contact between support-mapped shapes (`Cuboid`, `ConvexHull`, etc.)
///
/// This allows a more fine grained control other the underlying GJK algorigtm.
/// The vector-typed result is the vector that should be passed as `init` for
/// subsequent executions of the algorithm. It is also the contact
/// normal (that points toward the outside of the first solid).
pub fn contact_support_map_support_map_with_params<G1: ?Sized, G2: ?Sized>(
    pos12: Isometry,
    g1: &G1,
    g2: &G2,
    prediction: Real,
    simplex: &mut VoronoiSimplex,
    init_dir: Option<UnitVector>,
) -> GJKResult
where
    G1: SupportMap,
    G2: SupportMap,
{
    let dir = if let Some(init_dir) = init_dir {
        init_dir
    } else if let Ok(init_dir) =
        UnitVector::new_with_min(pos12.translation, crate::math::DEFAULT_EPSILON)
    {
        init_dir
    } else {
        UnitVector::X
    };

    //println!("{:?}", simplex);
    simplex.reset(CSOPoint::from_shapes(pos12, g1, g2, dir));

    let cpts = gjk::closest_points(pos12, g1, g2, prediction, true, simplex);
    if cpts != GJKResult::Intersection {
        return cpts;
    }

    // The point is inside of the CSO: use the fallback algorithm
    let mut epa = EPA::new();
    if let Some((p1, p2, n)) = epa.closest_points(pos12, g1, g2, simplex) {
        return GJKResult::ClosestPoints(p1, p2, n);
    }

    // Everything failed
    GJKResult::NoIntersection(UnitVector::X)
}
