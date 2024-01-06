use crate::math::{Real, Vector};

/// Computes the center of a set of point.
#[inline]
pub fn center(pts: &[Vector]) -> Vector {
    assert!(
        pts.len() >= 1,
        "Cannot compute the center of less than 1 point."
    );

    let denom: Real = 1.0 / (pts.len() as Real);

    let mut piter = pts.iter();
    let mut res = *piter.next().unwrap() * denom;

    for pt in piter.copied() {
        res += pt * denom;
    }

    res
}
