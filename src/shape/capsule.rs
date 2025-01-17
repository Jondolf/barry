use crate::math::{Isometry, Real, Rotation, UnitVector, Vector};
use crate::shape::{Segment, SupportMap};

#[cfg(feature = "std")]
use either::Either;

#[cfg(feature = "rkyv")]
use rkyv::{bytecheck, CheckBytes};

#[derive(Copy, Clone, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(feature = "bytemuck", derive(bytemuck::Pod, bytemuck::Zeroable))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize, CheckBytes),
    archive(as = "Self")
)]
#[cfg_attr(feature = "cuda", derive(cust_core::DeviceCopy))]
#[repr(C)]
/// A capsule shape defined as a round segment.
pub struct Capsule {
    /// The endpoints of the capsule’s principal axis.
    pub segment: Segment,
    /// The radius of the capsule.
    pub radius: Real,
}

impl Capsule {
    /// Creates a new capsule aligned with the `x` axis and with the given half-height an radius.
    pub fn new_x(half_height: Real, radius: Real) -> Self {
        let b = Vector::X * half_height;
        Self::new(-b, b, radius)
    }

    /// Creates a new capsule aligned with the `y` axis and with the given half-height an radius.
    pub fn new_y(half_height: Real, radius: Real) -> Self {
        let b = Vector::Y * half_height;
        Self::new(-b, b, radius)
    }

    /// Creates a new capsule aligned with the `z` axis and with the given half-height an radius.
    #[cfg(feature = "dim3")]
    pub fn new_z(half_height: Real, radius: Real) -> Self {
        let b = Vector::Z * half_height;
        Self::new(-b, b, radius)
    }

    /// Creates a new capsule defined as the segment between `a` and `b` and with the given `radius`.
    pub fn new(a: Vector, b: Vector, radius: Real) -> Self {
        let segment = Segment::new(a, b);
        Self { segment, radius }
    }

    /// The height of this capsule.
    pub fn height(&self) -> Real {
        (self.segment.b - self.segment.a).length()
    }

    /// The half-height of this capsule.
    pub fn half_height(&self) -> Real {
        self.height() / 2.0
    }

    /// The center of this capsule.
    pub fn center(&self) -> Vector {
        (self.segment.a + self.segment.b) / 2.0
    }

    /// Creates a new capsule equal to `self` with all its endpoints transformed by `pos`.
    pub fn transform_by(&self, pos: Isometry) -> Self {
        Self::new(pos * self.segment.a, pos * self.segment.b, self.radius)
    }

    /// The transformation such that `t * Y` is collinear with `b - a` and `t * origin` equals
    /// the capsule's center.
    #[doc(alias = "transform_wrt_y")]
    pub fn canonical_transform(&self) -> Isometry {
        Isometry {
            translation: self.center(),
            rotation: self.rotation_wrt_y(),
        }
    }

    /// The rotation `r` such that `r * Y` is collinear with `b - a`.
    pub fn rotation_wrt_y(&self) -> Rotation {
        let mut dir = self.segment.b - self.segment.a;

        if dir == Vector::ZERO {
            return Rotation::IDENTITY;
        }

        if dir.y < 0.0 {
            dir = -dir;
        }

        #[cfg(feature = "dim2")]
        {
            Rotation::from_rotation_arc_colinear(Vector::Y, dir)
        }

        #[cfg(feature = "dim3")]
        {
            Rotation::from_rotation_arc_colinear(Vector::Y, dir).unwrap_or(Rotation::IDENTITY)
        }
    }

    /// Computes a scaled version of this capsule.
    ///
    /// If the scaling factor is non-uniform, then it can’t be represented as
    /// capsule. Instead, a convex polygon approximation (with `nsubdivs`
    /// subdivisions) is returned. Returns `None` if that approximation had degenerate
    /// normals (for example if the scaling factor along one axis is zero).
    #[cfg(all(feature = "dim2", feature = "std"))]
    pub fn scaled(
        self,
        scale: Vector,
        nsubdivs: u32,
    ) -> Option<Either<Self, super::ConvexPolygon>> {
        if scale.x != scale.y {
            // The scaled shape is not a capsule.
            let mut vtx = self.to_polyline(nsubdivs);
            vtx.iter_mut().for_each(|pt| *pt = *pt * scale);
            Some(Either::Right(super::ConvexPolygon::from_convex_polyline(
                vtx,
            )?))
        } else {
            let uniform_scale = scale.x;
            Some(Either::Left(Self::new(
                self.segment.a * uniform_scale,
                self.segment.b * uniform_scale,
                self.radius * uniform_scale.abs(),
            )))
        }
    }

    /// Computes a scaled version of this capsule.
    ///
    /// If the scaling factor is non-uniform, then it can’t be represented as
    /// capsule. Instead, a convex polygon approximation (with `nsubdivs`
    /// subdivisions) is returned. Returns `None` if that approximation had degenerate
    /// normals (for example if the scaling factor along one axis is zero).
    #[cfg(all(feature = "dim3", feature = "std"))]
    pub fn scaled(
        self,
        scale: Vector,
        nsubdivs: u32,
    ) -> Option<Either<Self, super::ConvexPolyhedron>> {
        if scale.x != scale.y || scale.x != scale.z || scale.y != scale.z {
            // The scaled shape is not a capsule.
            let (mut vtx, idx) = self.to_trimesh(nsubdivs, nsubdivs);
            vtx.iter_mut().for_each(|pt| *pt = *pt * scale);
            Some(Either::Right(super::ConvexPolyhedron::from_convex_mesh(
                vtx, &idx,
            )?))
        } else {
            let uniform_scale = scale.x;
            Some(Either::Left(Self::new(
                self.segment.a * uniform_scale,
                self.segment.b * uniform_scale,
                self.radius * uniform_scale.abs(),
            )))
        }
    }
}

impl SupportMap for Capsule {
    fn local_support_point(&self, dir: Vector) -> Vector {
        let dir = UnitVector::new(dir).unwrap_or(UnitVector::Y);
        self.local_support_point_toward(dir)
    }

    fn local_support_point_toward(&self, dir: UnitVector) -> Vector {
        if dir.dot(self.segment.a) > dir.dot(self.segment.b) {
            self.segment.a + *dir * self.radius
        } else {
            self.segment.b + *dir * self.radius
        }
    }
}
