use crate::math::{AnyReal, AnyVector, Isometry, Real};
use crate::query::ClosestPoints;
use crate::shape::{Segment, SegmentPointLocation};

/// Closest points between segments.
#[inline]
pub fn closest_points_segment_segment(
    pos12: Isometry,
    seg1: &Segment,
    seg2: &Segment,
    margin: Real,
) -> ClosestPoints {
    let (loc1, loc2) = closest_points_segment_segment_with_locations(pos12, seg1, seg2);
    let p1 = seg1.point_at(&loc1);
    let p2 = seg2.point_at(&loc2);

    if p1.distance_squared(pos12 * p2) <= margin * margin {
        ClosestPoints::WithinMargin(p1, p2)
    } else {
        ClosestPoints::Disjoint
    }
}

// FIXME: use this specialized procedure for distance/interference/contact determination as well.
/// Closest points between two segments.
#[inline]
pub fn closest_points_segment_segment_with_locations(
    pos12: Isometry,
    seg1: &Segment,
    seg2: &Segment,
) -> (SegmentPointLocation, SegmentPointLocation) {
    let seg2_1 = seg2.transformed(pos12);
    closest_points_segment_segment_with_locations_nD((seg1.a, seg1.b), (seg2_1.a, seg2_1.b))
}

/// Segment-segment closest points computation in an arbitrary dimension.
#[allow(non_snake_case)]
#[inline]
pub fn closest_points_segment_segment_with_locations_nD<V: AnyVector>(
    seg1: (V, V),
    seg2: (V, V),
) -> (SegmentPointLocation, SegmentPointLocation) {
    let zero = V::Real::ZERO;
    let one = V::Real::ONE;

    // Inspired by RealField-time collision detection by Christer Ericson.
    let d1 = seg1.1 - seg1.0;
    let d2 = seg2.1 - seg2.0;
    let r = seg1.0 - seg2.0;

    let a = d1.length_squared();
    let e = d2.length_squared();
    let f = d2.dot(r);

    let mut s;
    let mut t;

    let _eps = V::Real::from_real(crate::math::DEFAULT_EPSILON);
    if a <= _eps && e <= _eps {
        s = zero;
        t = zero;
    } else if a <= _eps {
        s = zero;
        t = (f / e).clamp(zero, one);
    } else {
        let c = d1.dot(r);
        if e <= _eps {
            t = zero;
            s = (-c / a).clamp(zero, one);
        } else {
            let b = d1.dot(d2);
            let ae = a * e;
            let bb = b * b;
            let denom = ae - bb;

            // Use absolute and ulps error to test collinearity.
            if denom > _eps && !ulps_eq!(ae.to_real(), bb.to_real()) {
                s = ((b * f - c * e) / denom).clamp(zero, one);
            } else {
                s = zero;
            }

            t = (b * s + f) / e;

            if t < zero {
                t = zero;
                s = (-c / a).clamp(zero, one);
            } else if t > one {
                t = one;
                s = ((b - c) / a).clamp(zero, one);
            }
        }
    }

    let loc1 = if s == zero {
        SegmentPointLocation::OnVertex(0)
    } else if s == one {
        SegmentPointLocation::OnVertex(1)
    } else {
        SegmentPointLocation::OnEdge([(one - s).to_real(), s.to_real()])
    };

    let loc2 = if t == zero {
        SegmentPointLocation::OnVertex(0)
    } else if t == one {
        SegmentPointLocation::OnVertex(1)
    } else {
        SegmentPointLocation::OnEdge([(one - t).to_real(), t.to_real()])
    };

    (loc1, loc2)
}
