/*!
barry
========

**barry** is a 2 and 3-dimensional geometric library written with
the rust programming language.

*/

#![deny(non_camel_case_types)]
#![deny(unused_parens)]
#![deny(non_upper_case_globals)]
#![deny(unused_results)]
#![warn(missing_docs)] // TODO: deny this
#![warn(unused_imports)]
#![allow(missing_copy_implementations)]
#![doc(html_root_url = "http://docs.rs/barry/0.1.1")]
#![cfg_attr(not(feature = "std"), no_std)]
#![cfg_attr(not(feature = "rkyv"), deny(unused_qualifications))] // TODO: deny that everytime

use core::ops::{Add, AddAssign, Div, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign};

use bevy_math::{DVec2, DVec3, EulerRot, Quat, Vec2, Vec3};
use math::{Matrix2, Matrix3, Real, UVector, UnitVector2, UnitVector3, Vector, Vector2, Vector3};

#[cfg(all(
    feature = "simd-is-enabled",
    not(feature = "simd-stable"),
    not(feature = "simd-nightly")
))]
std::compile_error!("The `simd-is-enabled` feature should not be enabled explicitly. Please enable the `simd-stable` or the `simd-nightly` feature instead.");
#[cfg(all(feature = "simd-is-enabled", feature = "enhanced-determinism"))]
std::compile_error!(
    "SIMD cannot be enabled when the `enhanced-determinism` feature is also enabled."
);

macro_rules! array(
    ($callback: expr; SIMD_WIDTH) => {
        {
            #[inline(always)]
            #[allow(dead_code)]
            fn create_arr<T>(mut callback: impl FnMut(usize) -> T) -> [T; SIMD_WIDTH] {
                [callback(0usize), callback(1usize), callback(2usize), callback(3usize)]
            }

            create_arr($callback)
        }
    }
);

#[cfg(all(feature = "alloc", not(feature = "std")))]
#[cfg_attr(test, macro_use)]
extern crate alloc;

#[cfg(not(feature = "std"))]
extern crate core as std;

#[cfg(feature = "serde")]
#[macro_use]
extern crate serde;
#[macro_use]
extern crate approx;
#[macro_use]
#[cfg(feature = "dim3")]
extern crate bitflags;
extern crate num_traits as num;

pub extern crate either;
pub extern crate simba;

pub mod bounding_volume;
pub mod eigen;
pub mod mass_properties;
pub mod partitioning;
pub mod query;
pub mod shape;
#[cfg(feature = "std")]
pub mod transformation;
pub mod utils;
pub trait UVectorExt {
    fn as_vector(self) -> Vector;
}

impl UVectorExt for UVector {
    #[cfg(all(feature = "dim2", feature = "f32"))]
    fn as_vector(self) -> Vector {
        self.as_vec2()
    }
    #[cfg(all(feature = "dim2", feature = "f64"))]
    fn as_vector(self) -> Vector {
        self.as_dvec2()
    }
    #[cfg(all(feature = "dim3", feature = "f32"))]
    fn as_vector(self) -> Vector {
        self.as_vec3()
    }
    #[cfg(all(feature = "dim3", feature = "f64"))]
    fn as_vector(self) -> Vector {
        self.as_dvec3()
    }
}

mod real {
    /// The scalar type used throughout this crate.
    #[cfg(feature = "f64")]
    pub type Real = f64;

    #[cfg(feature = "f64")]
    pub use std::f64::consts as real_consts;

    /// The scalar type used throughout this crate.
    #[cfg(feature = "f32")]
    pub type Real = f32;

    #[cfg(feature = "f32")]
    pub use std::f32::consts as real_consts;
}

pub trait Zero: PartialEq + Sized {
    const ZERO: Self;

    fn is_zero(&self) -> bool {
        *self == Self::ZERO
    }
}

pub trait AnyReal:
    Clone
    + Copy
    + PartialEq
    + PartialOrd
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
    + Zero
{
    const MAX: Self;
}

pub trait AnyVector:
    Clone
    + Copy
    + PartialEq
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
    + Index<usize, Output = Self::Real>
    + IndexMut<usize, Output = Self::Real>
    + Zero
{
    type Real: AnyReal;

    fn dot(self, other: Self) -> Self::Real;
    fn ith(i: usize, value: Self::Real) -> Self {
        let mut vector = Self::ZERO;
        vector[i] = value;
        vector
    }
}

impl AnyReal for f32 {
    const MAX: Self = f32::MAX;
}
impl Zero for f32 {
    const ZERO: Self = 0.0;
}

impl AnyReal for f64 {
    const MAX: Self = f64::MAX;
}
impl Zero for f64 {
    const ZERO: Self = 0.0;
}

impl AnyVector for Vec2 {
    type Real = f32;

    fn dot(self, other: Self) -> Self::Real {
        self.dot(other)
    }
}
impl Zero for Vec2 {
    const ZERO: Self = Self::ZERO;
}

impl AnyVector for Vec3 {
    type Real = f32;

    fn dot(self, other: Self) -> Self::Real {
        self.dot(other)
    }
}
impl Zero for Vec3 {
    const ZERO: Self = Self::ZERO;
}

impl AnyVector for DVec2 {
    type Real = f64;

    fn dot(self, other: Self) -> Self::Real {
        self.dot(other)
    }
}
impl Zero for DVec2 {
    const ZERO: Self = Self::ZERO;
}

impl AnyVector for DVec3 {
    type Real = f64;

    fn dot(self, other: Self) -> Self::Real {
        self.dot(other)
    }
}
impl Zero for DVec3 {
    const ZERO: Self = Self::ZERO;
}

pub trait MinMaxIndex {
    fn min_index(&self) -> usize;
    fn max_index(&self) -> usize;
}

impl MinMaxIndex for Vector2 {
    fn min_index(&self) -> usize {
        if self.x < self.y {
            0
        } else {
            1
        }
    }
    fn max_index(&self) -> usize {
        if self.x >= self.y {
            0
        } else {
            1
        }
    }
}

impl MinMaxIndex for Vector3 {
    fn min_index(&self) -> usize {
        if self.x < self.y && self.x < self.z {
            0
        } else if self.y < self.z {
            1
        } else {
            2
        }
    }
    fn max_index(&self) -> usize {
        if self.x >= self.y && self.x >= self.z {
            0
        } else if self.y >= self.z {
            1
        } else {
            2
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serialize", derive(serde::Serialize, serde::Deserialize))]
pub struct Rotation2 {
    /// The cosine of the rotation angle in radians.
    cos: Real,
    /// The sine of the rotation angle in radians.
    sin: Real,
}

impl Default for Rotation2 {
    fn default() -> Self {
        Self::IDENTITY
    }
}

impl Rotation2 {
    /// Zero rotation.
    pub const IDENTITY: Self = Self { cos: 1.0, sin: 0.0 };

    /// Returns the cosine of the rotation in radians.
    pub fn cos(&self) -> Real {
        self.cos
    }

    /// Returns the sine of the rotation in radians.
    pub fn sin(&self) -> Real {
        self.sin
    }

    /// Creates a [`Rotation2`] from radians.
    pub fn from_radians(radians: Real) -> Self {
        Self {
            cos: radians.cos(),
            sin: radians.sin(),
        }
    }

    /// Creates a [`Rotation2`] from the sine and cosine of an angle in radians.
    pub fn from_sin_cos(sin: Real, cos: Real) -> Self {
        Self { sin, cos }
    }

    /// Creates a [`Rotation2`] from degrees.
    pub fn from_degrees(degrees: Real) -> Self {
        Self::from_radians(degrees.to_radians())
    }

    pub fn from_rotation_arc_colinear(from: Vector2, to: Vector2) -> Self {
        let sang = from.perp_dot(to);
        let cang = from.dot(to);
        Self::from_radians(sang.atan2(cang))
    }

    /// The smallest rotation needed to make `from` and `to` collinear and point toward the same
    /// direction, raised to the power `scale`.
    pub fn from_scaled_rotation_arc_colinear(from: Vector2, to: Vector2, scale: Real) -> Self {
        let sang = from.perp_dot(to);
        let cang = from.dot(to);
        Self::from_radians(sang.atan2(cang) * scale)
    }

    /// Returns the rotation in radians.
    pub fn as_radians(&self) -> Real {
        Real::atan2(self.sin(), self.cos())
    }

    /// Returns the rotation in degrees.
    pub fn as_degrees(&self) -> Real {
        self.as_radians().to_degrees()
    }

    /// Inverts the rotation.
    pub fn inverse(&self) -> Self {
        Self {
            cos: self.cos,
            sin: -self.sin,
        }
    }
}

impl From<Rotation2> for Matrix2 {
    fn from(rot: Rotation2) -> Self {
        Matrix2::from_cols_array(&[rot.cos, -rot.sin, rot.sin, rot.cos])
    }
}

impl Add<Self> for Rotation2 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        self.mul(rhs)
    }
}

impl AddAssign<Self> for Rotation2 {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl Sub<Self> for Rotation2 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        self.mul(rhs.inverse())
    }
}

impl SubAssign<Self> for Rotation2 {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl Mul for Rotation2 {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        Self {
            cos: self.cos * rhs.cos() - self.sin * rhs.sin(),
            sin: self.sin * rhs.cos() + self.cos * rhs.sin(),
        }
    }
}

impl MulAssign for Rotation2 {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl Mul<Vector2> for Rotation2 {
    type Output = Vector2;
    fn mul(self, rhs: Vector2) -> Self::Output {
        Vector2::new(
            rhs.x * self.cos() - rhs.y * self.sin(),
            rhs.x * self.sin() + rhs.y * self.cos(),
        )
    }
}

impl Mul<UnitVector2> for Rotation2 {
    type Output = UnitVector2;
    fn mul(self, rhs: UnitVector2) -> Self::Output {
        UnitVector2::from_normalized(self * *rhs)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Rotation3(Quat);

impl Default for Rotation3 {
    fn default() -> Self {
        Rotation3(Quat::IDENTITY)
    }
}

impl Rotation3 {
    /// Zero rotation.
    pub const IDENTITY: Self = Self(Quat::IDENTITY);

    pub fn from_rotation_arc_colinear(from: Vector3, to: Vector3) -> Self {
        Self(Quat::from_rotation_arc_colinear(from, to))
    }

    fn from_scaled_rotation_arc_colinear(from: Vector3, to: Vector3, scale: Real) -> Option<Self> {
        let c = from.cross(to);
        let cos = from.dot(to);

        if let Ok(axis) = UnitVector3::new(c) {
            if cos <= -1.0 {
                None
            } else if cos >= 1.0 {
                // The cosinus may be out of [-1, 1] because of inaccuracies.
                Some(Self::IDENTITY)
            } else {
                Some(Self::from_axis_angle(*axis, cos.acos() * scale))
            }
        } else if cos < 0.0 {
            // PI. The rotation axis is undefined but the angle is not zero.
            None
        } else {
            // Zero
            Some(Self::IDENTITY)
        }
    }

    pub fn from_scaled_axis(axis: Vector3) -> Self {
        Self(Quat::from_scaled_axis(axis))
    }

    pub fn from_axis_angle(axis: Vector3, angle: Real) -> Self {
        Self(Quat::from_axis_angle(axis, angle))
    }

    pub fn from_euler(euler: EulerRot, a: Real, b: Real, c: Real) -> Self {
        Self(Quat::from_euler(euler, a, b, c))
    }

    pub fn normalize(self) -> Self {
        Self(self.0.normalize())
    }

    /// Inverts the rotation.
    pub fn inverse(self) -> Self {
        Self(self.0.inverse())
    }
}

impl From<Rotation3> for Matrix3 {
    fn from(rotation: Rotation3) -> Self {
        Matrix3::from_quat(rotation.0)
    }
}

impl Mul for Rotation3 {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        Self(self.0 * rhs.0)
    }
}

impl MulAssign for Rotation3 {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl Mul<Vector3> for Rotation3 {
    type Output = Vector3;
    fn mul(self, rhs: Vector3) -> Self::Output {
        self.0 * rhs
    }
}

impl Mul<UnitVector3> for Rotation3 {
    type Output = UnitVector3;
    fn mul(self, rhs: UnitVector3) -> Self::Output {
        UnitVector3::from_normalized(self.0 * *rhs)
    }
}

#[derive(Clone, Copy, Default, Debug, PartialEq)]
pub struct Iso2 {
    translation: Vector2,
    rotation: Rotation2,
}

impl Iso2 {
    /// An identity isometry.
    pub const IDENTITY: Self = Self {
        translation: Vector2::ZERO,
        rotation: Rotation2::IDENTITY,
    };

    pub fn new(translation: Vector2, angle: Real) -> Self {
        Self {
            translation,
            rotation: Rotation2::from_radians(angle),
        }
    }

    pub const fn from_translation(translation: Vector2) -> Self {
        Self {
            translation,
            rotation: Rotation2::IDENTITY,
        }
    }

    pub const fn from_rotation(rotation: Rotation2) -> Self {
        Self {
            translation: Vector2::ZERO,
            rotation,
        }
    }

    pub fn inverse(self) -> Self {
        let inv_rot = self.rotation.inverse();
        Self {
            translation: inv_rot * -self.translation,
            rotation: inv_rot,
        }
    }

    pub fn inv_mul(self, rhs: Iso2) -> Self {
        let inv_rot = self.rotation.inverse();
        let delta_translation = rhs.translation - self.translation;
        Self {
            translation: inv_rot * delta_translation,
            rotation: inv_rot * rhs.rotation,
        }
    }

    pub fn transform_point(self, point: Vector2) -> Vector2 {
        self.translation + self.rotation * point
    }

    pub fn inverse_transform_point(self, point: Vector2) -> Vector2 {
        self.rotation.inverse() * (point - self.translation)
    }
}

impl Add for Iso2 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        Self {
            translation: self.translation + rhs.translation,
            rotation: self.rotation * rhs.rotation,
        }
    }
}

impl Sub for Iso2 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            translation: self.translation - rhs.translation,
            rotation: self.rotation * rhs.rotation.inverse(),
        }
    }
}

impl Mul for Iso2 {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        Self {
            translation: self.translation + self.rotation * rhs.translation,
            rotation: self.rotation * rhs.rotation,
        }
    }
}

impl Mul<Vector2> for Iso2 {
    type Output = Vector2;
    fn mul(self, rhs: Vector2) -> Self::Output {
        self.rotation * rhs
    }
}

impl Mul<UnitVector2> for Iso2 {
    type Output = UnitVector2;
    fn mul(self, rhs: UnitVector2) -> Self::Output {
        UnitVector2::from_normalized(self.rotation * *rhs)
    }
}

#[derive(Clone, Copy, Default, Debug, PartialEq)]
pub struct Iso3 {
    translation: Vector3,
    rotation: Rotation3,
}

impl Iso3 {
    /// An identity isometry.
    pub const IDENTITY: Self = Self {
        translation: Vector3::ZERO,
        rotation: Rotation3::IDENTITY,
    };

    pub fn new(translation: Vector3, rotation: Vector3) -> Self {
        Self {
            translation,
            rotation: Rotation3(Quat::from_scaled_axis(rotation)),
        }
    }

    pub const fn from_translation(translation: Vector3) -> Self {
        Self {
            translation,
            rotation: Rotation3::IDENTITY,
        }
    }

    pub const fn from_rotation(rotation: Rotation3) -> Self {
        Self {
            translation: Vector3::ZERO,
            rotation,
        }
    }

    pub fn inverse(self) -> Self {
        let inv_rot = self.rotation.inverse();
        Self {
            translation: inv_rot * -self.translation,
            rotation: inv_rot,
        }
    }

    pub fn inv_mul(self, rhs: Iso3) -> Self {
        let inv_rot = self.rotation.inverse();
        let delta_translation = rhs.translation - self.translation;
        Self {
            translation: inv_rot * delta_translation,
            rotation: inv_rot * rhs.rotation,
        }
    }

    pub fn transform_point(self, point: Vector3) -> Vector3 {
        self.translation + self.rotation * point
    }

    pub fn inverse_transform_point(self, point: Vector3) -> Vector3 {
        self.rotation.inverse() * (point - self.translation)
    }
}

impl Add for Iso3 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        Self {
            translation: self.translation + rhs.translation,
            rotation: self.rotation * rhs.rotation,
        }
    }
}

impl Sub for Iso3 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            translation: self.translation - rhs.translation,
            rotation: self.rotation * rhs.rotation.inverse(),
        }
    }
}

impl Mul for Iso3 {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        Self {
            translation: self.translation + self.rotation * rhs.translation,
            rotation: self.rotation * rhs.rotation,
        }
    }
}

impl Mul<Vector3> for Iso3 {
    type Output = Vector3;
    fn mul(self, rhs: Vector3) -> Self::Output {
        self.rotation * rhs
    }
}

impl Mul<UnitVector3> for Iso3 {
    type Output = UnitVector3;
    fn mul(self, rhs: UnitVector3) -> Self::Output {
        UnitVector3::from_normalized(self.rotation * *rhs)
    }
}

/// Compilation flags dependent aliases for mathematical types.
#[cfg(feature = "dim3")]
pub mod math {
    use crate::eigen::SymmetricEigen3;
    use crate::Iso2;
    use crate::Iso3;
    use crate::Rotation3;

    pub use super::real::*;
    pub use super::simd::*;
    use bevy_math::prelude::*;
    use bevy_math::primitives::Direction2d;
    use bevy_math::primitives::Direction3d;

    /// The default tolerance used for geometric operations.
    pub const DEFAULT_EPSILON: Real = Real::EPSILON;

    /// The dimension of the space.
    pub const DIM: usize = 3;

    /// The dimension of the space multiplied by two.
    pub const TWO_DIM: usize = DIM * 2;

    /// The dimension of the ambient space.
    pub type Dim = u8;

    /// The dimension of a spatial vector.
    pub type SpatialDim = u8;

    /// The dimension of the rotations.
    pub type AngDim = u8;

    /// The angular vector type.
    pub type AngVector = Vec3;

    /// The vector type.
    pub type Vector = Vec3;

    /// The 2D vector type.
    pub type Vector2 = Vec2;

    /// The 3D vector type.
    pub type Vector3 = Vec3;

    /// The u32 vector type.
    pub type UVector = UVec3;

    /// The 2D u32 vector type.
    pub type UVector2 = UVec2;

    /// The 3D vector type.
    pub type UVector3 = UVec3;

    /// The unit vector type.
    pub type UnitVector = Direction3d;

    /// The 2D unit vector type.
    pub type UnitVector2 = Direction2d;

    /// The 3D unit vector type.
    pub type UnitVector3 = Direction3d;

    /// The matrix type.
    pub type Matrix = Mat3;

    /// The 2x2 matrix type.
    pub type Matrix2 = Mat2;

    /// The 3x3 matrix type.
    pub type Matrix3 = Mat3;

    /// The eigen decomposition of a symmetric 3x3 matrix.
    pub type SymmetricEigen = SymmetricEigen3;

    /// The vector type with dimension `SpatialDim × 1`.
    pub type SpatialVector = (Vec3, Vec3);

    /// The orientation type.
    pub type Orientation = Vec3;

    /// The transformation matrix type.
    pub type Isometry = Iso3;

    /// The 2D transformation matrix type.
    pub type Isometry2 = Iso2;

    /// The 3D transformation matrix type.
    pub type Isometry3 = Iso3;

    /// The rotation matrix type.
    pub type Rotation = Rotation3;

    /// The translation type.
    pub type Translation = Vec3;

    /// The angular inertia of a rigid body.
    pub type AngularInertia = crate::utils::SdpMatrix3;

    /// The principal angular inertia of a rigid body.
    pub type PrincipalAngularInertia = Vec3;

    /// A matrix that represent the cross product with a given vector.
    pub type CrossMatrix = Mat3;

    /// A vector with a dimension equal to the maximum number of degrees of freedom of a rigid body.
    pub type SpacialVector = (Vec3, Vec3);

    // A 3D symmetric-definite-positive matrix.
    pub type SdpMatrix = crate::utils::SdpMatrix3;
}

/// Compilation flags dependent aliases for mathematical types.
#[cfg(feature = "dim2")]
pub mod math {
    use crate::eigen::SymmetricEigen2;
    use crate::{Iso2, Iso3, Rotation2};

    pub use super::real::*;
    pub use super::simd::*;
    use bevy_math::prelude::*;
    use bevy_math::primitives::{Direction2d, Direction3d};

    /// The default tolerance used for geometric operations.
    pub const DEFAULT_EPSILON: Real = Real::EPSILON;

    /// The dimension of the space.
    pub const DIM: usize = 2;

    /// The dimension of the space multiplied by two.
    pub const TWO_DIM: usize = DIM * 2;

    /// The dimension of the ambient space.
    pub type Dim = u8;

    /// The dimension of the rotations.
    pub type AngDim = u8;

    /// The angular vector type.
    pub type AngVector = f32;

    /// The vector type.
    pub type Vector = Vec2;

    /// The 2D vector type.
    pub type Vector2 = Vec2;

    /// The 3D vector type.
    pub type Vector3 = Vec3;

    /// The u32 vector type.
    pub type UVector = UVec2;

    /// The 2D u32 vector type.
    pub type UVector2 = UVec2;

    /// The 3D vector type.
    pub type UVector3 = UVec3;

    /// The unit vector type.
    pub type UnitVector = Direction2d;

    /// The 2D unit vector type.
    pub type UnitVector2 = Direction2d;

    /// The 3D unit vector type.
    pub type UnitVector3 = Direction3d;

    /// The matrix type.
    pub type Matrix = Mat2;

    /// The 2x2 matrix type.
    pub type Matrix2 = Mat2;

    /// The 3x3 matrix type.
    pub type Matrix3 = Mat3;

    /// The eigen decomposition of a symmetric 2x2 matrix.
    pub type SymmetricEigen = SymmetricEigen2;

    /// The orientation type.
    pub type Orientation = f32;

    /// The transformation matrix type.
    pub type Isometry = Iso2;

    /// The 2D transformation matrix type.
    pub type Isometry2 = Iso2;

    /// The 3D transformation matrix type.
    pub type Isometry3 = Iso3;

    pub type Rotation = Rotation2;

    /// The translation type.
    pub type Translation = Vec2;

    /// The angular inertia of a rigid body.
    pub type AngularInertia = f32;

    /// The principal angular inertia of a rigid body.
    pub type PrincipalAngularInertia = f32;

    /// A matrix that represent the cross product with a given vector.
    pub type CrossMatrix = Vec2;

    /// A vector with a dimension equal to the maximum number of degrees of freedom of a rigid body.
    pub type SpacialVector = Vec3;

    // A 2D symmetric-definite-positive matrix.
    pub type SdpMatrix = crate::utils::SdpMatrix2;
}

#[cfg(not(feature = "simd-is-enabled"))]
mod simd {
    use simba::simd::AutoBoolx4;
    /// The number of lanes of a SIMD number.
    pub const SIMD_WIDTH: usize = 4;
    /// SIMD_WIDTH - 1
    pub const SIMD_LAST_INDEX: usize = 3;

    /// A SIMD float with SIMD_WIDTH lanes.
    #[cfg(feature = "f32")]
    pub type SimdReal = simba::simd::AutoF32x4;

    /// A SIMD float with SIMD_WIDTH lanes.
    #[cfg(feature = "f64")]
    pub type SimdReal = simba::simd::AutoF64x4;

    /// A SIMD bool with SIMD_WIDTH lanes.
    pub type SimdBool = AutoBoolx4;
}

#[cfg(feature = "simd-is-enabled")]
mod simd {
    #[cfg(all(feature = "simd-nightly", feature = "f32"))]
    pub use simba::simd::{f32x4 as SimdReal, m32x4 as SimdBool};
    #[cfg(all(feature = "simd-stable", feature = "f32"))]
    pub use simba::simd::{WideBoolF32x4 as SimdBool, WideF32x4 as SimdReal};

    #[cfg(all(feature = "simd-nightly", feature = "f64"))]
    pub use simba::simd::{f64x4 as SimdReal, m64x4 as SimdBool};
    #[cfg(all(feature = "simd-stable", feature = "f64"))]
    pub use simba::simd::{WideBoolF64x4 as SimdBool, WideF64x4 as SimdReal};

    /// The number of lanes of a SIMD number.
    pub const SIMD_WIDTH: usize = 4;
    /// SIMD_WIDTH - 1
    pub const SIMD_LAST_INDEX: usize = 3;
}
