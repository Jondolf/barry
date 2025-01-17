//! Utilities useful for various generations tasks.

use crate::math::{Isometry, Real, Rotation, Vector};
#[cfg(feature = "dim3")]
use {crate::math::DIM, num::Zero};

/// Applies in-place a transformation to an array of points.
pub fn transform(points: &mut [Vector], m: Isometry) {
    points.iter_mut().for_each(|p| *p = m * *p);
}

/// Returns the transformed version of a vector of points.
pub fn transformed(mut points: Vec<Vector>, m: Isometry) -> Vec<Vector> {
    transform(&mut points, m);
    points
}

/// Returns the transformed version of a vector of points.
pub fn scaled(mut points: Vec<Vector>, scale: Vector) -> Vec<Vector> {
    points.iter_mut().for_each(|p| *p *= scale);
    points
}

// FIXME: remove that in favor of `push_xy_circle` ?
/// Pushes a discretized counterclockwise circle to a buffer.
#[cfg(feature = "dim3")]
#[inline]
pub fn push_circle(radius: Real, nsubdiv: u32, dtheta: Real, y: Real, out: &mut Vec<Vector>) {
    let mut curr_theta = Real::zero();

    for _ in 0..nsubdiv {
        out.push(Vector::new(
            curr_theta.cos() * radius,
            y.clone(),
            curr_theta.sin() * radius,
        ));
        curr_theta = curr_theta + dtheta;
    }
}

/// Pushes a discretized counterclockwise circle to a buffer.
/// The circle is contained on the plane spanned by the `x` and `y` axis.
#[inline]
#[cfg(feature = "dim2")]
pub fn push_xy_arc(radius: Real, nsubdiv: u32, dtheta: Real, out: &mut Vec<Vector>) {
    let mut curr_theta: Real = 0.0;

    for _ in 0..nsubdiv {
        let mut pt_coords = Vector::ZERO;

        pt_coords[0] = curr_theta.cos() * radius;
        pt_coords[1] = curr_theta.sin() * radius;
        out.push(Vector::from(pt_coords));

        curr_theta = curr_theta + dtheta;
    }
}

/// Creates the faces from two circles with the same discretization.
#[cfg(feature = "dim3")]
#[inline]
pub fn push_ring_indices(
    base_lower_circle: u32,
    base_upper_circle: u32,
    nsubdiv: u32,
    out: &mut Vec<[u32; DIM]>,
) {
    push_open_ring_indices(base_lower_circle, base_upper_circle, nsubdiv, out);

    // adjust the last two triangles
    push_rectangle_indices(
        base_upper_circle,
        base_upper_circle + nsubdiv - 1,
        base_lower_circle,
        base_lower_circle + nsubdiv - 1,
        out,
    );
}

/// Creates the faces from two circles with the same discretization.
#[cfg(feature = "dim3")]
#[inline]
pub fn push_open_ring_indices(
    base_lower_circle: u32,
    base_upper_circle: u32,
    nsubdiv: u32,
    out: &mut Vec<[u32; DIM]>,
) {
    assert!(nsubdiv > 0);

    for i in 0..nsubdiv - 1 {
        let bli = base_lower_circle + i;
        let bui = base_upper_circle + i;
        push_rectangle_indices(bui + 1, bui, bli + 1, bli, out);
    }
}

/// Creates the faces from a circle and a point that is shared by all triangle.
#[cfg(feature = "dim3")]
#[inline]
pub fn push_degenerate_top_ring_indices(
    base_circle: u32,
    point: u32,
    nsubdiv: u32,
    out: &mut Vec<[u32; DIM]>,
) {
    push_degenerate_open_top_ring_indices(base_circle, point, nsubdiv, out);

    out.push([base_circle + nsubdiv - 1, point, base_circle]);
}

/// Creates the faces from a circle and a point that is shared by all triangle.
#[cfg(feature = "dim3")]
#[inline]
pub fn push_degenerate_open_top_ring_indices(
    base_circle: u32,
    point: u32,
    nsubdiv: u32,
    out: &mut Vec<[u32; DIM]>,
) {
    assert!(nsubdiv > 0);

    for i in 0..nsubdiv - 1 {
        out.push([base_circle + i, point, base_circle + i + 1]);
    }
}

/// Pushes indices so that a circle is filled with triangles. Each triangle will have the
/// `base_circle` point in common.
/// Pushes `nsubdiv - 2` elements to `out`.
#[cfg(feature = "dim3")]
#[inline]
pub fn push_filled_circle_indices(base_circle: u32, nsubdiv: u32, out: &mut Vec<[u32; DIM]>) {
    for i in base_circle + 1..base_circle + nsubdiv - 1 {
        out.push([base_circle, i, i + 1]);
    }
}

/// Given four corner points, pushes to two counterclockwise triangles to `out`.
///
/// # Arguments:
/// * `ul` - the up-left point.
/// * `dl` - the down-left point.
/// * `dr` - the down-right point.
/// * `ur` - the up-right point.
#[cfg(feature = "dim3")]
#[inline]
pub fn push_rectangle_indices(ul: u32, ur: u32, dl: u32, dr: u32, out: &mut Vec<[u32; DIM]>) {
    out.push([ul.clone(), dl, dr.clone()]);
    out.push([dr, ur, ul]);
}

/// Reverses the clockwising of a set of faces.
#[cfg(feature = "dim3")]
#[inline]
pub fn reverse_clockwising(indices: &mut [[u32; DIM]]) {
    indices.iter_mut().for_each(|idx| idx.swap(0, 1));
}

/// Pushes the index buffer of a closed loop.
#[cfg(feature = "dim3")]
#[inline]
pub fn push_circle_outline_indices(indices: &mut Vec<[u32; 2]>, range: std::ops::Range<u32>) {
    indices.extend((range.start..range.end - 1).map(|i| [i, i + 1]));
    indices.push([range.end - 1, range.start]);
}

/// Pushes the index buffer of an open chain.
#[cfg(feature = "dim3")]
#[inline]
pub fn push_open_circle_outline_indices(indices: &mut Vec<[u32; 2]>, range: std::ops::Range<u32>) {
    indices.extend((range.start..range.end - 1).map(|i| [i, i + 1]));
}

/// Pushes to `out_vtx` a set of points forming an arc starting at `start`, ending at `end` with
/// revolution center at `center`. The curve is approximated by pushing `nsubdivs` points.
/// The `start` and `end` point are not pushed to `out_vtx`.
///
/// ALso pushes to `out_idx` the appropriate index buffer to form the arc (including attaches to
/// the `start` and `end` points).
#[cfg(feature = "dim3")]
pub fn push_arc_and_idx(
    center: Vector,
    start: u32,
    end: u32,
    nsubdivs: u32,
    out_vtx: &mut Vec<Vector>,
    out_idx: &mut Vec<[u32; 2]>,
) {
    let base = out_vtx.len() as u32;
    push_arc(
        center,
        out_vtx[start as usize],
        out_vtx[end as usize],
        nsubdivs,
        out_vtx,
    );
    push_arc_idx(start, base..base + nsubdivs - 1, end, out_idx);
}

/// Pushes to `out` a set of points forming an arc starting at `start`, ending at `end` with
/// revolution center at `center`. The curve is approximated by pushing `nsubdivs` points.
/// The `start` and `end` point are not pushed to `out`.
pub fn push_arc(center: Vector, start: Vector, end: Vector, nsubdivs: u32, out: &mut Vec<Vector>) {
    assert!(nsubdivs > 0);

    let start_len = (start - center).length();
    let end_len = (end - center).length();

    let start_dir = (start - center) / start_len;
    let end_dir = (end - center) / end_len;

    if !start_dir.is_finite() || !end_dir.is_finite() {
        return;
    }

    let len_inc = (end_len - start_len) / nsubdivs as Real;

    #[cfg(feature = "dim2")]
    let rot = Some(Rotation::from_scaled_rotation_arc_colinear(
        start_dir,
        end_dir,
        1.0 / nsubdivs as Real,
    ));

    #[cfg(feature = "dim3")]
    let rot =
        Rotation::from_scaled_rotation_arc_colinear(start_dir, end_dir, 1.0 / nsubdivs as Real);

    if let Some(rot) = rot {
        let mut curr_dir = start_dir;
        let mut curr_len = start_len;

        for _ in 0..nsubdivs - 1 {
            curr_dir = rot * curr_dir;
            curr_len += len_inc;

            out.push(center + curr_dir * curr_len);
        }
    }
}

/// Pushes the index buffer for an arc between `start` and `end` and intermediate points in the
/// range `arc`.
#[cfg(feature = "dim3")]
pub fn push_arc_idx(start: u32, arc: std::ops::Range<u32>, end: u32, out: &mut Vec<[u32; 2]>) {
    if arc.is_empty() {
        out.push([start, end]);
    } else {
        out.push([start, arc.start]);
        for i in arc.start..arc.end - 1 {
            out.push([i, i + 1])
        }
        out.push([arc.end - 1, end])
    }
}

/// Applies a revolution, using the Y symmetry axis passing through the origin.
#[cfg(feature = "dim3")]
pub fn apply_revolution(
    collapse_bottom: bool,
    collapse_top: bool,
    circle_ranges: &[std::ops::Range<u32>],
    nsubdivs: u32,
    out_vtx: &mut Vec<Vector>, // Must be set to the half-profile.
    out_idx: &mut Vec<[u32; 2]>,
) {
    use crate::math;
    use bevy_math::Vec3Swizzles;

    let ang_increment = math::real_consts::TAU / (nsubdivs as Real);
    let angs = [
        ang_increment * (nsubdivs / 4) as Real,
        ang_increment * (nsubdivs / 2) as Real,
        ang_increment * ((3 * nsubdivs) / 4) as Real,
    ];

    let half_profile_len = out_vtx.len();

    for k in 0..half_profile_len as u32 - 1 {
        out_idx.push([k, k + 1]);
    }

    let mut range = 0..half_profile_len;

    if collapse_bottom {
        range.start += 1;
    }
    if collapse_top {
        range.end -= 1;
    }

    // Push rotated profiles.
    for i in 0..3 {
        let base = out_vtx.len() as u32;
        let rot = Rotation::from_scaled_axis(Vector::Y * angs[i]);

        if collapse_bottom {
            out_idx.push([0, base]);
        }

        for k in range.clone() {
            out_vtx.push(rot * out_vtx[k]);
        }

        for k in 0..range.len() as u32 - 1 {
            out_idx.push([base + k, base + k + 1]);
        }

        if collapse_top {
            out_idx.push([base + range.len() as u32 - 1, half_profile_len as u32 - 1]);
        }
    }

    // Push circles.
    // TODO: right now, this duplicates some points, to simplify the index
    //       buffer construction.
    for circle_range in circle_ranges {
        for i in circle_range.clone() {
            let pt = out_vtx[i as usize];
            let base = out_vtx.len() as u32;
            push_circle(pt.xz().length(), nsubdivs, ang_increment, pt.y, out_vtx);
            push_circle_outline_indices(out_idx, base..base + nsubdivs)
        }
    }
}
