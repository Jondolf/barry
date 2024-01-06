//! Axis Aligned Bounding Box.

use crate::bounding_volume::{BoundingSphere, BoundingVolume};
use crate::math::{Isometry, Real, UnitVector, Vector, DIM, TWO_DIM};
use crate::shape::{Cuboid, SupportMap};
use crate::utils::IsometryOps;
use arrayvec::ArrayVec;

#[cfg(feature = "rkyv")]
use rkyv::{bytecheck, CheckBytes};

/// An Axis Aligned Bounding Box.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(feature = "bytemuck", derive(bytemuck::Pod, bytemuck::Zeroable))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize, CheckBytes),
    archive(as = "Self")
)]
#[cfg_attr(feature = "cuda", derive(cust_core::DeviceCopy))]
#[derive(Debug, PartialEq, Copy, Clone)]
#[repr(C)]
pub struct Aabb {
    pub mins: Vector,
    pub maxs: Vector,
}

impl Aabb {
    /// The vertex indices of each edge of this `Aabb`.
    ///
    /// This gives, for each edge of this `Aabb`, the indices of its
    /// vertices when taken from the `self.vertices()` array.
    /// Here is how the faces are numbered, assuming
    /// a right-handed coordinate system:
    ///
    /// ```text
    ///    y             3 - 2
    ///    |           7 − 6 |
    ///    ___ x       |   | 1  (the zero is below 3 and on the left of 1,
    ///   /            4 - 5     hidden by the 4-5-6-7 face.)
    ///  z
    /// ```
    #[cfg(feature = "dim3")]
    pub const EDGES_VERTEX_IDS: [(usize, usize); 12] = [
        (0, 1),
        (1, 2),
        (3, 2),
        (0, 3),
        (4, 5),
        (5, 6),
        (7, 6),
        (4, 7),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
    ];

    /// The vertex indices of each face of this `Aabb`.
    ///
    /// This gives, for each face of this `Aabb`, the indices of its
    /// vertices when taken from the `self.vertices()` array.
    /// Here is how the faces are numbered, assuming
    /// a right-handed coordinate system:
    ///
    /// ```text
    ///    y             3 - 2
    ///    |           7 − 6 |
    ///    ___ x       |   | 1  (the zero is below 3 and on the left of 1,
    ///   /            4 - 5     hidden by the 4-5-6-7 face.)
    ///  z
    /// ```
    #[cfg(feature = "dim3")]
    pub const FACES_VERTEX_IDS: [(usize, usize, usize, usize); 6] = [
        (1, 2, 6, 5),
        (0, 3, 7, 4),
        (2, 3, 7, 6),
        (1, 0, 4, 5),
        (4, 5, 6, 7),
        (0, 1, 2, 3),
    ];

    /// Creates a new Aabb.
    ///
    /// # Arguments:
    ///   * `mins` - position of the point with the smallest coordinates.
    ///   * `maxs` - position of the point with the highest coordinates. Each component of `mins`
    ///   must be smaller than the related components of `maxs`.
    #[inline]
    pub fn new(mins: Vector, maxs: Vector) -> Aabb {
        Aabb { mins, maxs }
    }

    /// Creates an invalid `Aabb` with `mins` components set to `Real::max_values` and `maxs`components set to `-Real::max_values`.
    ///
    /// This is often used as the initial values of some `Aabb` merging algorithms.
    #[inline]
    pub fn new_invalid() -> Self {
        Self::new(Vector::MAX, -Vector::MAX)
    }

    /// Creates a new `Aabb` from its center and its half-extents.
    #[inline]
    pub fn from_half_extents(center: Vector, half_extents: Vector) -> Self {
        Self::new(center - half_extents, center + half_extents)
    }

    /// Creates a new `Aabb` from a set of points.
    pub fn from_points<'a, I>(pts: I) -> Self
    where
        I: IntoIterator<Item = &'a Vector>,
    {
        super::aabb_utils::local_point_cloud_aabb(pts)
    }

    /// The center of this `Aabb`.
    #[inline]
    pub fn center(&self) -> Vector {
        (self.mins + self.maxs) / 2.0
    }

    /// The half extents of this `Aabb`.
    #[inline]
    pub fn half_extents(&self) -> Vector {
        (self.maxs - self.mins) / 2.0
    }

    /// The volume of this `Aabb`.
    #[inline]
    pub fn volume(&self) -> Real {
        let extents = self.extents();
        #[cfg(feature = "dim2")]
        return extents.x * extents.y;
        #[cfg(feature = "dim3")]
        return extents.x * extents.y * extents.z;
    }

    /// The extents of this `Aabb`.
    #[inline]
    pub fn extents(&self) -> Vector {
        self.maxs - self.mins
    }

    /// Enlarges this `Aabb` so it also contains the point `pt`.
    pub fn take_point(&mut self, pt: Vector) {
        self.mins = self.mins.min(pt).into();
        self.maxs = self.maxs.max(pt).into();
    }

    /// Computes the `Aabb` bounding `self` transformed by `m`.
    #[inline]
    pub fn transform_by(&self, m: Isometry) -> Self {
        let ls_center = self.center();
        let center = m * ls_center;
        let ws_half_extents = m.absolute_transform_vector(self.half_extents());

        Aabb::new(center + (-ws_half_extents), center + ws_half_extents)
    }

    #[inline]
    pub fn scaled(self, scale: Vector) -> Self {
        let a = self.mins * scale;
        let b = self.maxs * scale;
        Self {
            mins: a.min(b).into(),
            maxs: a.max(b).into(),
        }
    }

    /// The smallest bounding sphere containing this `Aabb`.
    #[inline]
    pub fn bounding_sphere(&self) -> BoundingSphere {
        let center = self.center();
        let radius = self.mins.distance(self.maxs) / 2.0;
        BoundingSphere::new(center, radius)
    }

    #[inline]
    pub fn contains_local_point(&self, point: Vector) -> bool {
        for i in 0..DIM {
            if point[i] < self.mins[i] || point[i] > self.maxs[i] {
                return false;
            }
        }

        true
    }

    /// Computes the intersection of this `Aabb` and another one.
    pub fn intersection(&self, other: &Aabb) -> Option<Aabb> {
        let result = Aabb {
            mins: Vector::from(self.mins.max(other.mins)),
            maxs: Vector::from(self.maxs.min(other.maxs)),
        };

        for i in 0..DIM {
            if result.mins[i] > result.maxs[i] {
                return None;
            }
        }

        Some(result)
    }

    /// Returns the difference between this `Aabb` and `rhs`.
    ///
    /// Removing another `Aabb` from `self` will result in zero, one, or up to 4 (in 2D) or 8 (in 3D)
    /// new smaller Aabbs.
    pub fn difference(&self, rhs: &Aabb) -> ArrayVec<Self, TWO_DIM> {
        self.difference_with_cut_sequence(rhs).0
    }

    /// Returns the difference between this `Aabb` and `rhs`.
    ///
    /// Removing another `Aabb` from `self` will result in zero, one, or up to 4 (in 2D) or 8 (in 3D)
    /// new smaller Aabbs.
    ///
    /// # Return
    /// This returns a pair where the first item are the new Aabbs and the the second item is
    /// the sequance of cuts applied to `self` to obtain the new Aabbs. Each cut is performed
    /// along one axis identified by `-1, -2, -3` for `-X, -Y, -Z` and `1, 2, 3` for `+X, +Y, +Z`, and
    /// the plane’s bias.
    /// The cuts are applied sequancially. For example, if `result.1[0]` contains `1`, then it means
    /// that `result.0[0]` is equal to the piece of `self` lying in the negative half-space delimited
    /// by the plane with outward normal `+X`. Then, the other piece of `self` generated by this cut
    /// (i.e. the piece of `self` lying in the positive half-space delimited by the plane with outward
    /// normal `+X`) is the one that will be affected by the next cut.
    ///
    /// The returned cut sequence will be empty if the aabbs are disjoint.
    pub fn difference_with_cut_sequence(
        &self,
        rhs: &Aabb,
    ) -> (ArrayVec<Self, TWO_DIM>, ArrayVec<(i8, Real), TWO_DIM>) {
        let mut result = ArrayVec::new();
        let mut cut_sequence = ArrayVec::new();

        // NOTE: special case when the boxes are disjoint.
        //       This isn’t exactly the same as `!self.intersects(rhs)`
        //       because of the equality.
        for i in 0..DIM {
            if self.mins[i] >= rhs.maxs[i] || self.maxs[i] <= rhs.mins[i] {
                result.push(*self);
                return (result, cut_sequence);
            }
        }

        let mut rest = *self;

        for i in 0..DIM {
            if rhs.mins[i] > rest.mins[i] {
                let mut fragment = rest;
                fragment.maxs[i] = rhs.mins[i];
                rest.mins[i] = rhs.mins[i];
                result.push(fragment);
                cut_sequence.push((i as i8 + 1, rhs.mins[i]));
            }

            if rhs.maxs[i] < rest.maxs[i] {
                let mut fragment = rest;
                fragment.mins[i] = rhs.maxs[i];
                rest.maxs[i] = rhs.maxs[i];
                result.push(fragment);
                cut_sequence.push((-(i as i8 + 1), -rhs.maxs[i]));
            }
        }

        (result, cut_sequence)
    }

    /// Computes the vertices of this `Aabb`.
    #[inline]
    #[cfg(feature = "dim2")]
    pub fn vertices(&self) -> [Vector; 4] {
        [
            Vector::new(self.mins.x, self.mins.y),
            Vector::new(self.mins.x, self.maxs.y),
            Vector::new(self.maxs.x, self.mins.y),
            Vector::new(self.maxs.x, self.maxs.y),
        ]
    }

    /// Computes the vertices of this `Aabb`.
    #[inline]
    #[cfg(feature = "dim3")]
    pub fn vertices(&self) -> [Vector; 8] {
        [
            Vector::new(self.mins.x, self.mins.y, self.mins.z),
            Vector::new(self.maxs.x, self.mins.y, self.mins.z),
            Vector::new(self.maxs.x, self.maxs.y, self.mins.z),
            Vector::new(self.mins.x, self.maxs.y, self.mins.z),
            Vector::new(self.mins.x, self.mins.y, self.maxs.z),
            Vector::new(self.maxs.x, self.mins.y, self.maxs.z),
            Vector::new(self.maxs.x, self.maxs.y, self.maxs.z),
            Vector::new(self.mins.x, self.maxs.y, self.maxs.z),
        ]
    }

    /// Splits this `Aabb` at its center, into four parts (as in a quad-tree).
    #[inline]
    #[cfg(feature = "dim2")]
    pub fn split_at_center(&self) -> [Aabb; 4] {
        let center = self.center();

        [
            Aabb::new(self.mins, center),
            Aabb::new(
                Vector::new(center.x, self.mins.y),
                Vector::new(self.maxs.x, center.y),
            ),
            Aabb::new(center, self.maxs),
            Aabb::new(
                Vector::new(self.mins.x, center.y),
                Vector::new(center.x, self.maxs.y),
            ),
        ]
    }

    /// Splits this `Aabb` at its center, into eight parts (as in an octree).
    #[inline]
    #[cfg(feature = "dim3")]
    pub fn split_at_center(&self) -> [Aabb; 8] {
        let center = self.center();

        [
            Aabb::new(
                Vector::new(self.mins.x, self.mins.y, self.mins.z),
                Vector::new(center.x, center.y, center.z),
            ),
            Aabb::new(
                Vector::new(center.x, self.mins.y, self.mins.z),
                Vector::new(self.maxs.x, center.y, center.z),
            ),
            Aabb::new(
                Vector::new(center.x, center.y, self.mins.z),
                Vector::new(self.maxs.x, self.maxs.y, center.z),
            ),
            Aabb::new(
                Vector::new(self.mins.x, center.y, self.mins.z),
                Vector::new(center.x, self.maxs.y, center.z),
            ),
            Aabb::new(
                Vector::new(self.mins.x, self.mins.y, center.z),
                Vector::new(center.x, center.y, self.maxs.z),
            ),
            Aabb::new(
                Vector::new(center.x, self.mins.y, center.z),
                Vector::new(self.maxs.x, center.y, self.maxs.z),
            ),
            Aabb::new(
                Vector::new(center.x, center.y, center.z),
                Vector::new(self.maxs.x, self.maxs.y, self.maxs.z),
            ),
            Aabb::new(
                Vector::new(self.mins.x, center.y, center.z),
                Vector::new(center.x, self.maxs.y, self.maxs.z),
            ),
        ]
    }

    /// Projects every point of `Aabb` on an arbitrary axis.
    pub fn project_on_axis(&self, axis: UnitVector) -> (Real, Real) {
        let cuboid = Cuboid::new(self.half_extents());
        let shift = cuboid.local_support_point_toward(axis).dot(*axis).abs();
        let center = self.center().dot(*axis);
        (center - shift, center + shift)
    }
}

impl BoundingVolume for Aabb {
    #[inline]
    fn center(&self) -> Vector {
        self.center()
    }

    #[inline]
    fn intersects(&self, other: &Aabb) -> bool {
        self.mins.cmplt(other.maxs).all() && self.maxs.cmpgt(other.mins).all()
    }

    #[inline]
    fn contains(&self, other: &Aabb) -> bool {
        self.mins.cmplt(other.mins).all() && self.maxs.cmpgt(other.maxs).all()
    }

    #[inline]
    fn merge(&mut self, other: &Aabb) {
        self.mins = self.mins.min(other.mins);
        self.maxs = self.maxs.max(other.maxs);
    }

    #[inline]
    fn merged(&self, other: &Aabb) -> Aabb {
        Aabb {
            mins: self.mins.min(other.mins),
            maxs: self.maxs.max(other.maxs),
        }
    }

    #[inline]
    fn loosen(&mut self, amount: Real) {
        assert!(amount >= 0.0, "The loosening margin must be positive.");
        self.mins = self.mins + Vector::splat(-amount);
        self.maxs = self.maxs + Vector::splat(amount);
    }

    #[inline]
    fn loosened(&self, amount: Real) -> Aabb {
        assert!(amount >= 0.0, "The loosening margin must be positive.");
        Aabb {
            mins: self.mins + Vector::splat(-amount),
            maxs: self.maxs + Vector::splat(amount),
        }
    }

    #[inline]
    fn tighten(&mut self, amount: Real) {
        assert!(amount >= 0.0, "The tightening margin must be positive.");
        self.mins = self.mins + Vector::splat(amount);
        self.maxs = self.maxs + Vector::splat(-amount);
        assert!(
            self.mins.cmplt(self.maxs).all(),
            "The tightening margin is to large."
        );
    }

    #[inline]
    fn tightened(&self, amount: Real) -> Aabb {
        assert!(amount >= 0.0, "The tightening margin must be positive.");

        Aabb::new(
            self.mins + Vector::splat(amount),
            self.maxs + Vector::splat(-amount),
        )
    }
}
