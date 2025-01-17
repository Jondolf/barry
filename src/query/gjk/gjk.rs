//! The Gilbert–Johnson–Keerthi distance algorithm.

use crate::math::{Isometry, Real, UnitVector, Vector, DIM};
use crate::query::gjk::{CSOPoint, ConstantOrigin, VoronoiSimplex};
use crate::query::{self, Ray};
use crate::shape::SupportMap;

use num::{Bounded, Zero};

/// Results of the GJK algorithm.
#[derive(Clone, Debug, PartialEq)]
pub enum GJKResult {
    /// Result of the GJK algorithm when the origin is inside of the polytope.
    Intersection,
    /// Result of the GJK algorithm when a projection of the origin on the polytope is found.
    ///
    /// Both points and vector are expressed in the local-space of the first geometry involved
    /// in the GJK execution.
    ClosestPoints(Vector, Vector, UnitVector),
    /// Result of the GJK algorithm when the origin is too close to the polytope but not inside of it.
    ///
    /// The returned vector is expressed in the local-space of the first geometry involved in the
    /// GJK execution.
    Proximity(UnitVector),
    /// Result of the GJK algorithm when the origin is too far away from the polytope.
    ///
    /// The returned vector is expressed in the local-space of the first geomety involved in the
    /// GJK execution.
    NoIntersection(UnitVector),
}

/// The absolute tolerence used by the GJK algorithm.
pub const EPS_TOLERANCE: Real = 10.0 * crate::math::DEFAULT_EPSILON;

/// Projects the origin on the boundary of the given shape.
///
/// The origin is assumed to be outside of the shape. If it is inside,
/// use the EPA algorithm instead.
/// Return `None` if the origin is not inside of the shape or if
/// the EPA algorithm failed to compute the projection.
///
/// Return the projected point in the local-space of `g`.
pub fn project_origin<G: ?Sized>(m: Isometry, g: &G, simplex: &mut VoronoiSimplex) -> Option<Vector>
where
    G: SupportMap,
{
    match closest_points(
        m.inverse(),
        g,
        &ConstantOrigin,
        Real::max_value(),
        true,
        simplex,
    ) {
        GJKResult::Intersection => None,
        GJKResult::ClosestPoints(p, _, _) => Some(p),
        _ => unreachable!(),
    }
}

/*
 * Separating Axis GJK
 */
/// Projects the origin on a shape using the Separating Axis GJK algorithm.
/// The algorithm will stop as soon as the polytope can be proven to be at least `max_dist` away
/// from the origin.
///
/// # Arguments:
/// * simplex - the simplex to be used by the GJK algorithm. It must be already initialized
///             with at least one point on the shape boundary.
/// * exact_dist - if `false`, the gjk will stop as soon as it can prove that the origin is at
/// a distance smaller than `max_dist` but not inside of `shape`. In that case, it returns a
/// `GJKResult::Proximity(sep_axis)` where `sep_axis` is a separating axis. If `false` the gjk will
/// compute the exact distance and return `GJKResult::Projection(point)` if the origin is closer
/// than `max_dist` but not inside `shape`.
pub fn closest_points<G1: ?Sized, G2: ?Sized>(
    pos12: Isometry,
    g1: &G1,
    g2: &G2,
    max_dist: Real,
    exact_dist: bool,
    simplex: &mut VoronoiSimplex,
) -> GJKResult
where
    G1: SupportMap,
    G2: SupportMap,
{
    let _eps = crate::math::DEFAULT_EPSILON;
    let _eps_tol: Real = EPS_TOLERANCE;
    let _eps_rel: Real = _eps_tol.sqrt();

    // FIXME: reset the simplex if it is empty?
    let mut proj = simplex.project_origin_and_reduce();

    let mut old_dir;

    if let Ok(proj_dir) = UnitVector::new(proj) {
        old_dir = -proj_dir;
    } else {
        return GJKResult::Intersection;
    }

    let mut max_bound = Real::max_value();
    let mut dir;
    let mut niter = 0;

    loop {
        let old_max_bound = max_bound;

        if let Ok((new_dir, dist)) = UnitVector::new_and_length(-proj) {
            dir = new_dir;
            max_bound = dist;
        } else {
            // The origin is on the simplex.
            return GJKResult::Intersection;
        }

        if max_bound >= old_max_bound {
            if exact_dist {
                let (p1, p2) = result(simplex, true);
                return GJKResult::ClosestPoints(p1, p2, old_dir); // upper bounds inconsistencies
            } else {
                return GJKResult::Proximity(old_dir);
            }
        }

        let cso_point = CSOPoint::from_shapes(pos12, g1, g2, dir);
        let min_bound = -dir.dot(cso_point.point);

        assert!(min_bound.is_finite());

        if min_bound > max_dist {
            return GJKResult::NoIntersection(dir);
        } else if !exact_dist && min_bound > 0.0 && max_bound <= max_dist {
            return GJKResult::Proximity(old_dir);
        } else if max_bound - min_bound <= _eps_rel * max_bound {
            if exact_dist {
                let (p1, p2) = result(simplex, false);
                return GJKResult::ClosestPoints(p1, p2, dir); // the distance found has a good enough precision
            } else {
                return GJKResult::Proximity(dir);
            }
        }

        if !simplex.add_point(cso_point) {
            if exact_dist {
                let (p1, p2) = result(simplex, false);
                return GJKResult::ClosestPoints(p1, p2, dir);
            } else {
                return GJKResult::Proximity(dir);
            }
        }

        old_dir = dir;
        proj = simplex.project_origin_and_reduce();

        if simplex.dimension() == DIM {
            if min_bound >= _eps_tol {
                if exact_dist {
                    let (p1, p2) = result(simplex, true);
                    return GJKResult::ClosestPoints(p1, p2, old_dir);
                } else {
                    // NOTE: previous implementation used old_proj here.
                    return GJKResult::Proximity(old_dir);
                }
            } else {
                return GJKResult::Intersection; // Point inside of the cso.
            }
        }
        niter += 1;
        if niter == 10000 {
            return GJKResult::NoIntersection(UnitVector::X);
        }
    }
}

/// Casts a ray on a support map using the GJK algorithm.
pub fn cast_local_ray<G: ?Sized>(
    shape: &G,
    simplex: &mut VoronoiSimplex,
    ray: &Ray,
    max_toi: Real,
) -> Option<(Real, Vector)>
where
    G: SupportMap,
{
    let g2 = ConstantOrigin;
    minkowski_ray_cast(Isometry::IDENTITY, shape, &g2, ray, max_toi, simplex)
}

/// Compute the normal and the distance that can travel `g1` along the direction
/// `dir` so that `g1` and `g2` just touch.
///
/// The `dir` vector must be expressed in the local-space of the first shape.
pub fn directional_distance<G1: ?Sized, G2: ?Sized>(
    pos12: Isometry,
    g1: &G1,
    g2: &G2,
    dir: Vector,
    simplex: &mut VoronoiSimplex,
) -> Option<(Real, Vector, Vector, Vector)>
where
    G1: SupportMap,
    G2: SupportMap,
{
    let ray = Ray::new(Vector::ZERO, dir);
    minkowski_ray_cast(pos12, g1, g2, &ray, Real::max_value(), simplex).map(|(toi, normal)| {
        let witnesses = if !toi.is_zero() {
            result(simplex, simplex.dimension() == DIM)
        } else {
            // If there is penetration, the witness points
            // are undefined.
            (Vector::ZERO, Vector::ZERO)
        };

        (toi, normal, witnesses.0, witnesses.1)
    })
}

// Ray-cast on the Minkowski Difference `g1 - pos12 * g2`.
fn minkowski_ray_cast<G1: ?Sized, G2: ?Sized>(
    pos12: Isometry,
    g1: &G1,
    g2: &G2,
    ray: &Ray,
    max_toi: Real,
    simplex: &mut VoronoiSimplex,
) -> Option<(Real, Vector)>
where
    G1: SupportMap,
    G2: SupportMap,
{
    let _eps = crate::math::DEFAULT_EPSILON;
    let _eps_tol: Real = EPS_TOLERANCE;
    let _eps_rel: Real = _eps_tol.sqrt();

    let ray_length = ray.dir.length();

    if relative_eq!(ray_length, 0.0) {
        return None;
    }

    let mut ltoi = 0.0;
    let mut curr_ray = Ray::new(ray.origin, ray.dir / ray_length);
    let dir = -curr_ray.dir;
    let mut ldir = dir;

    // Initialize the simplex.
    let support_point = CSOPoint::from_shapes(pos12, g1, g2, UnitVector::new_unchecked(dir));
    simplex.reset(support_point.translate(-curr_ray.origin));

    // FIXME: reset the simplex if it is empty?
    let mut proj = simplex.project_origin_and_reduce();
    let mut max_bound = Real::MAX;
    let mut dir;
    let mut niter = 0;
    let mut last_chance = false;

    loop {
        let old_max_bound = max_bound;

        if let Ok((new_dir, dist)) = UnitVector::new_and_length(-proj) {
            dir = new_dir;
            max_bound = dist;
        } else {
            return Some((ltoi / ray_length, ldir));
        }

        let support_point = if max_bound >= old_max_bound {
            // Upper bounds inconsistencies. Consider the projection as a valid support point.
            last_chance = true;
            CSOPoint::single_point(proj + curr_ray.origin)
        } else {
            CSOPoint::from_shapes(pos12, g1, g2, dir)
        };

        if last_chance && ltoi > 0.0 {
            // last_chance && ltoi > 0.0 && (support_point.point - curr_ray.origin).dot(ldir) >= 0.0 {
            return Some((ltoi / ray_length, ldir));
        }

        // Clip the ray on the support halfspace (None <=> t < 0)
        // The configurations are:
        //   dir.dot(curr_ray.dir)  |   t   |               Action
        // −−−−−−−−−−−−−−−−−−−−-----+−−−−−−−+−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
        //          < 0             |  < 0  | Continue.
        //          < 0             |  > 0  | New lower bound, move the origin.
        //          > 0             |  < 0  | Miss. No intersection.
        //          > 0             |  > 0  | New higher bound.
        match query::details::ray_toi_with_halfspace(support_point.point, dir, &curr_ray) {
            Some(t) => {
                if dir.dot(curr_ray.dir) < 0.0 && t > 0.0 {
                    // new lower bound
                    ldir = *dir;
                    ltoi += t;

                    // NOTE: we divide by ray_length instead of doing max_toi * ray_length
                    // because the multiplication may cause an overflow if max_toi is set
                    // to Real::max_value() by users that want to have an infinite ray.
                    if ltoi / ray_length > max_toi {
                        return None;
                    }

                    let shift = curr_ray.dir * t;
                    curr_ray.origin += shift;
                    max_bound = Real::max_value();
                    simplex.modify_pnts(&|pt| pt.translate_mut(-shift));
                    last_chance = false;
                }
            }
            None => {
                if dir.dot(curr_ray.dir) > _eps_tol {
                    // miss
                    return None;
                }
            }
        }

        if last_chance {
            return None;
        }

        let min_bound = -dir.dot(support_point.point - curr_ray.origin);

        assert!(min_bound.is_finite());

        if max_bound - min_bound <= _eps_rel * max_bound {
            // This is needed when using fixed-points to avoid missing
            // some castes.
            // FIXME: I feel like we should always return `Some` in
            // this case, even with floating-point numbers. Though it
            // has not been sufficiently tested with floats yet to be sure.
            if cfg!(feature = "improved_fixed_point_support") {
                return Some((ltoi / ray_length, ldir));
            } else {
                return None;
            }
        }

        let _ = simplex.add_point(support_point.translate(-curr_ray.origin));
        proj = simplex.project_origin_and_reduce();

        if simplex.dimension() == DIM {
            if min_bound >= _eps_tol {
                return None;
            } else {
                return Some((ltoi / ray_length, ldir)); // Point inside of the cso.
            }
        }

        niter += 1;
        if niter == 10000 {
            return None;
        }
    }
}

fn result(simplex: &VoronoiSimplex, prev: bool) -> (Vector, Vector) {
    let mut res = (Vector::ZERO, Vector::ZERO);
    if prev {
        for i in 0..simplex.prev_dimension() + 1 {
            let coord = simplex.prev_proj_coord(i);
            let point = simplex.prev_point(i);
            res.0 += point.orig1 * coord;
            res.1 += point.orig2 * coord;
        }

        res
    } else {
        for i in 0..simplex.dimension() + 1 {
            let coord = simplex.proj_coord(i);
            let point = simplex.point(i);
            res.0 += point.orig1 * coord;
            res.1 += point.orig2 * coord;
        }

        res
    }
}
