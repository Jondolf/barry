use core::ops::{Add, Div, Index, IndexMut, Mul, Neg, Sub};

use bevy_math::{DVec2, DVec3, Vec2, Vec3};

mod eigen;
mod isometry;
mod rotation;
mod simd;

pub use eigen::*;
pub use isometry::*;
pub use real::*;
pub use rotation::*;
pub use simd::*;

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
{
    const ZERO: Self;
    const ONE: Self;
    const MAX: Self;

    fn sqrt(self) -> Self;
    fn from_real(number: Real) -> Self;
    fn to_real(self) -> Real;
    fn clamp(self, min: Self, max: Self) -> Self;
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
{
    type Real: AnyReal;

    const ZERO: Self;
    const ONE: Self;

    fn dot(self, other: Self) -> Self::Real;
    fn length(self) -> Self::Real {
        self.length_squared().sqrt()
    }
    fn length_squared(self) -> Self::Real;
    fn ith(i: usize, value: Self::Real) -> Self {
        let mut vector = Self::ZERO;
        vector[i] = value;
        vector
    }
}

impl AnyReal for f32 {
    const ZERO: Self = 0.0;
    const ONE: Self = 1.0;
    const MAX: Self = Self::MAX;

    fn sqrt(self) -> Self {
        self.sqrt()
    }
    fn from_real(number: Real) -> Self {
        number as Self
    }
    fn to_real(self) -> Real {
        self as Real
    }
    fn clamp(self, min: Self, max: Self) -> Self {
        self.clamp(min, max)
    }
}

impl AnyReal for f64 {
    const ZERO: Self = 0.0;
    const ONE: Self = 1.0;
    const MAX: Self = Self::MAX;

    fn sqrt(self) -> Self {
        self.sqrt()
    }
    fn from_real(number: Real) -> Self {
        number as Self
    }
    fn to_real(self) -> Real {
        self as Real
    }
    fn clamp(self, min: Self, max: Self) -> Self {
        self.clamp(min, max)
    }
}

impl AnyVector for Vec2 {
    const ZERO: Self = Self::ZERO;
    const ONE: Self = Self::ONE;

    type Real = f32;

    fn dot(self, other: Self) -> Self::Real {
        self.dot(other)
    }
    fn length_squared(self) -> Self::Real {
        self.length_squared()
    }
}

impl AnyVector for Vec3 {
    const ZERO: Self = Self::ZERO;
    const ONE: Self = Self::ONE;

    type Real = f32;

    fn dot(self, other: Self) -> Self::Real {
        self.dot(other)
    }
    fn length_squared(self) -> Self::Real {
        self.length_squared()
    }
}
impl Zero for Vec3 {
    const ZERO: Self = Self::ZERO;
}

impl AnyVector for DVec2 {
    const ZERO: Self = Self::ZERO;
    const ONE: Self = Self::ONE;

    type Real = f64;

    fn dot(self, other: Self) -> Self::Real {
        self.dot(other)
    }
    fn length_squared(self) -> Self::Real {
        self.length_squared()
    }
}

impl AnyVector for DVec3 {
    const ZERO: Self = Self::ZERO;
    const ONE: Self = Self::ONE;

    type Real = f64;

    fn dot(self, other: Self) -> Self::Real {
        self.dot(other)
    }
    fn length_squared(self) -> Self::Real {
        self.length_squared()
    }
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

#[cfg(feature = "dim2")]
pub use math2d::*;
#[cfg(feature = "dim3")]
pub use math3d::*;

use crate::MinMaxIndex;

/// Compilation flags dependent aliases for mathematical types.
#[cfg(feature = "dim2")]
mod math2d {
    use super::{eigen::SymmetricEigen2, Iso2, Iso3, Rotation2};

    pub use super::real::*;
    pub use crate::simd::*;
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

/// Compilation flags dependent aliases for mathematical types.
#[cfg(feature = "dim3")]
mod math3d {
    use super::{Iso2, Iso3, Real, Rotation3, SymmetricEigen3};

    pub use crate::simd::*;
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
