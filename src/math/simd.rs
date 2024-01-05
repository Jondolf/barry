use core::ops::{
    Add, AddAssign, Deref, DerefMut, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Rem,
    RemAssign, Sub, SubAssign,
};

use bevy_math::{EulerRot, Quat};
use simba::simd::{AutoSimd, SimdComplexField, SimdPartialOrd, SimdValue};

use crate::{
    math::{Isometry2, Isometry3, Real, Rotation2, Rotation3, Vector2, Vector3},
    simd::{SimdBool, SimdReal},
};

#[cfg(feature = "dim2")]
pub type SimdVector = SimdVec2;
#[cfg(feature = "dim3")]
pub type SimdVector = SimdVec3;

#[cfg(feature = "dim2")]
pub type SimdMatrix = SimdMat2;
#[cfg(feature = "dim3")]
pub type SimdMatrix = SimdMat3;

#[cfg(feature = "dim2")]
pub type SimdRotation = SimdRotation2;
#[cfg(feature = "dim3")]
pub type SimdRotation = SimdRotation3;

#[cfg(feature = "dim2")]
pub type SimdIsometry = SimdIso2;
#[cfg(feature = "dim3")]
pub type SimdIsometry = SimdIso3;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SimdVec2 {
    pub x: SimdReal,
    pub y: SimdReal,
}

impl Default for SimdVec2 {
    fn default() -> Self {
        Self::ZERO
    }
}

impl Index<usize> for SimdVec2 {
    type Output = SimdReal;
    fn index(&self, index: usize) -> &Self::Output {
        match index {
            0 => &self.x,
            1 => &self.y,
            _ => panic!("index {index} exceeds vector length 2"),
        }
    }
}

impl IndexMut<usize> for SimdVec2 {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            _ => panic!("index {index} exceeds vector length 2"),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SimdVec3 {
    pub x: SimdReal,
    pub y: SimdReal,
    pub z: SimdReal,
}

impl Default for SimdVec3 {
    fn default() -> Self {
        Self::ZERO
    }
}

impl Index<usize> for SimdVec3 {
    type Output = SimdReal;
    fn index(&self, index: usize) -> &Self::Output {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("index {index} exceeds vector length 3"),
        }
    }
}

impl IndexMut<usize> for SimdVec3 {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            _ => panic!("index {index} exceeds vector length 3"),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SimdQuat {
    pub x: SimdReal,
    pub y: SimdReal,
    pub z: SimdReal,
    pub w: SimdReal,
}

impl Default for SimdQuat {
    fn default() -> Self {
        Self::IDENTITY
    }
}

pub struct SimdMat2 {
    pub x_axis: SimdVec2,
    pub y_axis: SimdVec2,
}

pub struct SimdMat3 {
    pub x_axis: SimdVec3,
    pub y_axis: SimdVec3,
    pub z_axis: SimdVec3,
}

impl From<SimdRotation2> for SimdMat2 {
    fn from(rot: SimdRotation2) -> Self {
        Self {
            x_axis: SimdVec2 {
                x: rot.cos,
                y: -rot.sin,
            },
            y_axis: SimdVec2 {
                x: rot.sin,
                y: rot.cos,
            },
        }
    }
}

impl From<SimdRotation3> for SimdMat3 {
    #[inline]
    #[must_use]
    fn from(rotation: SimdRotation3) -> Self {
        let x2 = rotation.x + rotation.x;
        let y2 = rotation.y + rotation.y;
        let z2 = rotation.z + rotation.z;
        let xx = rotation.x * x2;
        let xy = rotation.x * y2;
        let xz = rotation.x * z2;
        let yy = rotation.y * y2;
        let yz = rotation.y * z2;
        let zz = rotation.z * z2;
        let wx = rotation.w * x2;
        let wy = rotation.w * y2;
        let wz = rotation.w * z2;

        let one = SimdReal::splat(1.0);

        Self {
            x_axis: SimdVec3 {
                x: one - (yy + zz),
                y: xy + wz,
                z: xz - wy,
            },
            y_axis: SimdVec3 {
                x: xy - wz,
                y: one - (xx + zz),
                z: yz + wx,
            },
            z_axis: SimdVec3 {
                x: xz + wy,
                y: yz - wx,
                z: one - (xx + yy),
            },
        }
    }
}

impl Mul<SimdVec2> for SimdMat2 {
    type Output = SimdVec2;
    fn mul(self, rhs: SimdVec2) -> Self::Output {
        self.x_axis * rhs.x + self.y_axis * rhs.y
    }
}

impl Mul<SimdVec3> for SimdMat3 {
    type Output = SimdVec3;
    fn mul(self, rhs: SimdVec3) -> Self::Output {
        self.x_axis * rhs.x + self.y_axis * rhs.y + self.z_axis * rhs.z
    }
}

impl SimdMat2 {
    pub fn abs(self) -> Self {
        Self {
            x_axis: self.x_axis.abs(),
            y_axis: self.y_axis.abs(),
        }
    }
}

impl SimdMat3 {
    pub fn abs(self) -> Self {
        Self {
            x_axis: self.x_axis.abs(),
            y_axis: self.y_axis.abs(),
            z_axis: self.z_axis.abs(),
        }
    }
}

impl SimdVec2 {
    pub const ZERO: Self = Self::splat_simd_real(AutoSimd([0.0; 4]));

    pub fn new(x: Real, y: Real) -> Self {
        Self {
            x: SimdReal::splat(x),
            y: SimdReal::splat(y),
        }
    }

    pub fn splat(vec: Vector2) -> Self {
        Self {
            x: SimdReal::splat(vec.x),
            y: SimdReal::splat(vec.y),
        }
    }

    pub const fn splat_simd_real(real: SimdReal) -> Self {
        Self { x: real, y: real }
    }

    pub fn from_vecs(vecs: [Vector2; 4]) -> Self {
        Self {
            x: SimdReal::new(vecs[0].x, vecs[1].x, vecs[2].x, vecs[3].x),
            y: SimdReal::new(vecs[0].y, vecs[1].y, vecs[2].y, vecs[3].y),
        }
    }

    pub fn extract(&self, lane: usize) -> Vector2 {
        Vector2 {
            x: self.x.extract(lane),
            y: self.y.extract(lane),
        }
    }

    pub fn map<F: Fn(SimdReal) -> f32>(self, f: F) -> Vector2 {
        Vector2 {
            x: f(self.x),
            y: f(self.y),
        }
    }

    pub fn replace(&mut self, lane: usize, new: Vector2) {
        self.x.replace(lane, new.x);
        self.y.replace(lane, new.y);
    }

    pub fn select(self, cond: SimdBool, other: Self) -> Self {
        Self {
            x: self.x.select(cond, other.x),
            y: self.y.select(cond, other.y),
        }
    }

    pub fn length(self) -> SimdReal {
        self.length_squared().simd_sqrt()
    }

    pub fn length_squared(self) -> SimdReal {
        self.dot(self)
    }

    pub fn dot(self, other: Self) -> SimdReal {
        self.x * other.x + self.y * other.y
    }

    pub fn min(self, other: Self) -> Self {
        Self {
            x: self.x.simd_min(other.x),
            y: self.y.simd_min(other.y),
        }
    }

    pub fn max(self, other: Self) -> Self {
        Self {
            x: self.x.simd_max(other.x),
            y: self.y.simd_max(other.y),
        }
    }

    pub fn abs(self) -> Self {
        Self {
            x: self.x.simd_abs(),
            y: self.y.simd_abs(),
        }
    }
}

impl Add<SimdVec2> for SimdVec2 {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self {
            x: self.x.add(rhs.x),
            y: self.y.add(rhs.y),
        }
    }
}

impl AddAssign<SimdVec2> for SimdVec2 {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.x.add_assign(rhs.x);
        self.y.add_assign(rhs.y);
    }
}

impl Add<SimdReal> for SimdVec2 {
    type Output = Self;
    #[inline]
    fn add(self, rhs: SimdReal) -> Self {
        Self {
            x: self.x.add(rhs),
            y: self.y.add(rhs),
        }
    }
}

impl AddAssign<SimdReal> for SimdVec2 {
    #[inline]
    fn add_assign(&mut self, rhs: SimdReal) {
        self.x.add_assign(rhs);
        self.y.add_assign(rhs);
    }
}

impl Add<SimdVec2> for SimdReal {
    type Output = SimdVec2;
    #[inline]
    fn add(self, rhs: SimdVec2) -> SimdVec2 {
        SimdVec2 {
            x: self.add(rhs.x),
            y: self.add(rhs.y),
        }
    }
}

impl Sub<SimdVec2> for SimdVec2 {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self {
            x: self.x.sub(rhs.x),
            y: self.y.sub(rhs.y),
        }
    }
}

impl SubAssign<SimdVec2> for SimdVec2 {
    #[inline]
    fn sub_assign(&mut self, rhs: SimdVec2) {
        self.x.sub_assign(rhs.x);
        self.y.sub_assign(rhs.y);
    }
}

impl Sub<SimdReal> for SimdVec2 {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: SimdReal) -> Self {
        Self {
            x: self.x.sub(rhs),
            y: self.y.sub(rhs),
        }
    }
}

impl SubAssign<SimdReal> for SimdVec2 {
    #[inline]
    fn sub_assign(&mut self, rhs: SimdReal) {
        self.x.sub_assign(rhs);
        self.y.sub_assign(rhs);
    }
}

impl Sub<SimdVec2> for SimdReal {
    type Output = SimdVec2;
    #[inline]
    fn sub(self, rhs: SimdVec2) -> SimdVec2 {
        SimdVec2 {
            x: self.sub(rhs.x),
            y: self.sub(rhs.y),
        }
    }
}

impl Mul<SimdVec2> for SimdVec2 {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        Self {
            x: self.x.mul(rhs.x),
            y: self.y.mul(rhs.y),
        }
    }
}

impl MulAssign<SimdVec2> for SimdVec2 {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        self.x.mul_assign(rhs.x);
        self.y.mul_assign(rhs.y);
    }
}

impl Mul<SimdReal> for SimdVec2 {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: SimdReal) -> Self {
        Self {
            x: self.x.mul(rhs),
            y: self.y.mul(rhs),
        }
    }
}

impl MulAssign<SimdReal> for SimdVec2 {
    #[inline]
    fn mul_assign(&mut self, rhs: SimdReal) {
        self.x.mul_assign(rhs);
        self.y.mul_assign(rhs);
    }
}

impl Mul<SimdVec2> for SimdReal {
    type Output = SimdVec2;
    #[inline]
    fn mul(self, rhs: SimdVec2) -> SimdVec2 {
        SimdVec2 {
            x: self.mul(rhs.x),
            y: self.mul(rhs.y),
        }
    }
}

impl Div<SimdVec2> for SimdVec2 {
    type Output = Self;
    #[inline]
    fn div(self, rhs: Self) -> Self {
        Self {
            x: self.x.div(rhs.x),
            y: self.y.div(rhs.y),
        }
    }
}

impl DivAssign<SimdVec2> for SimdVec2 {
    #[inline]
    fn div_assign(&mut self, rhs: Self) {
        self.x.div_assign(rhs.x);
        self.y.div_assign(rhs.y);
    }
}

impl Div<SimdReal> for SimdVec2 {
    type Output = Self;
    #[inline]
    fn div(self, rhs: SimdReal) -> Self {
        Self {
            x: self.x.div(rhs),
            y: self.y.div(rhs),
        }
    }
}

impl DivAssign<SimdReal> for SimdVec2 {
    #[inline]
    fn div_assign(&mut self, rhs: SimdReal) {
        self.x.div_assign(rhs);
        self.y.div_assign(rhs);
    }
}

impl Div<SimdVec2> for SimdReal {
    type Output = SimdVec2;
    #[inline]
    fn div(self, rhs: SimdVec2) -> SimdVec2 {
        SimdVec2 {
            x: self.div(rhs.x),
            y: self.div(rhs.y),
        }
    }
}

impl Rem<SimdVec2> for SimdVec2 {
    type Output = Self;
    #[inline]
    fn rem(self, rhs: Self) -> Self {
        Self {
            x: self.x.rem(rhs.x),
            y: self.y.rem(rhs.y),
        }
    }
}

impl RemAssign<SimdVec2> for SimdVec2 {
    #[inline]
    fn rem_assign(&mut self, rhs: Self) {
        self.x.rem_assign(rhs.x);
        self.y.rem_assign(rhs.y);
    }
}

impl Rem<SimdReal> for SimdVec2 {
    type Output = Self;
    #[inline]
    fn rem(self, rhs: SimdReal) -> Self {
        Self {
            x: self.x.rem(rhs),
            y: self.y.rem(rhs),
        }
    }
}

impl RemAssign<SimdReal> for SimdVec2 {
    #[inline]
    fn rem_assign(&mut self, rhs: SimdReal) {
        self.x.rem_assign(rhs);
        self.y.rem_assign(rhs);
    }
}

impl Rem<SimdVec2> for SimdReal {
    type Output = SimdVec2;
    #[inline]
    fn rem(self, rhs: SimdVec2) -> SimdVec2 {
        SimdVec2 {
            x: self.rem(rhs.x),
            y: self.rem(rhs.y),
        }
    }
}

impl Neg for SimdVec2 {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self {
            x: -self.x,
            y: -self.y,
        }
    }
}

impl SimdVec3 {
    pub const ZERO: Self = Self::splat_simd_real(AutoSimd([0.0; 4]));

    pub fn new(x: Real, y: Real, z: Real) -> Self {
        Self {
            x: SimdReal::splat(x),
            y: SimdReal::splat(y),
            z: SimdReal::splat(z),
        }
    }

    pub fn splat(vec: Vector3) -> Self {
        Self {
            x: SimdReal::splat(vec.x),
            y: SimdReal::splat(vec.y),
            z: SimdReal::splat(vec.z),
        }
    }

    pub const fn splat_simd_real(real: SimdReal) -> Self {
        Self {
            x: real,
            y: real,
            z: real,
        }
    }

    pub fn from_vecs(vecs: [Vector3; 4]) -> Self {
        Self {
            x: SimdReal::new(vecs[0].x, vecs[1].x, vecs[2].x, vecs[3].x),
            y: SimdReal::new(vecs[0].y, vecs[1].y, vecs[2].y, vecs[3].y),
            z: SimdReal::new(vecs[0].z, vecs[1].z, vecs[2].z, vecs[3].z),
        }
    }

    pub fn extract(&self, lane: usize) -> Vector3 {
        Vector3 {
            x: self.x.extract(lane),
            y: self.y.extract(lane),
            z: self.z.extract(lane),
        }
    }

    pub fn map<F: Fn(SimdReal) -> f32>(self, f: F) -> Vector3 {
        Vector3 {
            x: f(self.x),
            y: f(self.y),
            z: f(self.z),
        }
    }

    pub fn replace(&mut self, lane: usize, new: Vector3) {
        self.x.replace(lane, new.x);
        self.y.replace(lane, new.y);
        self.z.replace(lane, new.z);
    }

    pub fn select(self, cond: SimdBool, other: Self) -> Self {
        Self {
            x: self.x.select(cond, other.x),
            y: self.y.select(cond, other.y),
            z: self.z.select(cond, other.z),
        }
    }

    pub fn length(self) -> SimdReal {
        self.length_squared().simd_sqrt()
    }

    pub fn length_squared(self) -> SimdReal {
        self.dot(self)
    }

    pub fn dot(self, other: Self) -> SimdReal {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn cross(self, other: Self) -> SimdVec3 {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    pub fn min(self, other: Self) -> Self {
        Self {
            x: self.x.simd_min(other.x),
            y: self.y.simd_min(other.y),
            z: self.y.simd_min(other.z),
        }
    }

    pub fn max(self, other: Self) -> Self {
        Self {
            x: self.x.simd_max(other.x),
            y: self.y.simd_max(other.y),
            z: self.y.simd_max(other.z),
        }
    }

    pub fn abs(self) -> Self {
        Self {
            x: self.x.simd_abs(),
            y: self.y.simd_abs(),
            z: self.z.simd_abs(),
        }
    }
}

impl Add<SimdVec3> for SimdVec3 {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self {
            x: self.x.add(rhs.x),
            y: self.y.add(rhs.y),
            z: self.z.add(rhs.z),
        }
    }
}

impl AddAssign<SimdVec3> for SimdVec3 {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.x.add_assign(rhs.x);
        self.y.add_assign(rhs.y);
        self.z.add_assign(rhs.z);
    }
}

impl Add<SimdReal> for SimdVec3 {
    type Output = Self;
    #[inline]
    fn add(self, rhs: SimdReal) -> Self {
        Self {
            x: self.x.add(rhs),
            y: self.y.add(rhs),
            z: self.z.add(rhs),
        }
    }
}

impl AddAssign<SimdReal> for SimdVec3 {
    #[inline]
    fn add_assign(&mut self, rhs: SimdReal) {
        self.x.add_assign(rhs);
        self.y.add_assign(rhs);
        self.z.add_assign(rhs);
    }
}

impl Add<SimdVec3> for SimdReal {
    type Output = SimdVec3;
    #[inline]
    fn add(self, rhs: SimdVec3) -> SimdVec3 {
        SimdVec3 {
            x: self.add(rhs.x),
            y: self.add(rhs.y),
            z: self.add(rhs.z),
        }
    }
}

impl Sub<SimdVec3> for SimdVec3 {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self {
            x: self.x.sub(rhs.x),
            y: self.y.sub(rhs.y),
            z: self.z.sub(rhs.z),
        }
    }
}

impl SubAssign<SimdVec3> for SimdVec3 {
    #[inline]
    fn sub_assign(&mut self, rhs: SimdVec3) {
        self.x.sub_assign(rhs.x);
        self.y.sub_assign(rhs.y);
        self.z.sub_assign(rhs.z);
    }
}

impl Sub<SimdReal> for SimdVec3 {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: SimdReal) -> Self {
        Self {
            x: self.x.sub(rhs),
            y: self.y.sub(rhs),
            z: self.z.sub(rhs),
        }
    }
}

impl SubAssign<SimdReal> for SimdVec3 {
    #[inline]
    fn sub_assign(&mut self, rhs: SimdReal) {
        self.x.sub_assign(rhs);
        self.y.sub_assign(rhs);
        self.z.sub_assign(rhs);
    }
}

impl Sub<SimdVec3> for SimdReal {
    type Output = SimdVec3;
    #[inline]
    fn sub(self, rhs: SimdVec3) -> SimdVec3 {
        SimdVec3 {
            x: self.sub(rhs.x),
            y: self.sub(rhs.y),
            z: self.sub(rhs.z),
        }
    }
}

impl Mul<SimdVec3> for SimdVec3 {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        Self {
            x: self.x.mul(rhs.x),
            y: self.y.mul(rhs.y),
            z: self.z.mul(rhs.z),
        }
    }
}

impl MulAssign<SimdVec3> for SimdVec3 {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        self.x.mul_assign(rhs.x);
        self.y.mul_assign(rhs.y);
        self.z.mul_assign(rhs.z);
    }
}

impl Mul<SimdReal> for SimdVec3 {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: SimdReal) -> Self {
        Self {
            x: self.x.mul(rhs),
            y: self.y.mul(rhs),
            z: self.z.mul(rhs),
        }
    }
}

impl MulAssign<SimdReal> for SimdVec3 {
    #[inline]
    fn mul_assign(&mut self, rhs: SimdReal) {
        self.x.mul_assign(rhs);
        self.y.mul_assign(rhs);
        self.z.mul_assign(rhs);
    }
}

impl Mul<SimdVec3> for SimdReal {
    type Output = SimdVec3;
    #[inline]
    fn mul(self, rhs: SimdVec3) -> SimdVec3 {
        SimdVec3 {
            x: self.mul(rhs.x),
            y: self.mul(rhs.y),
            z: self.mul(rhs.z),
        }
    }
}

impl Div<SimdVec3> for SimdVec3 {
    type Output = Self;
    #[inline]
    fn div(self, rhs: Self) -> Self {
        Self {
            x: self.x.div(rhs.x),
            y: self.y.div(rhs.y),
            z: self.z.div(rhs.z),
        }
    }
}

impl DivAssign<SimdVec3> for SimdVec3 {
    #[inline]
    fn div_assign(&mut self, rhs: Self) {
        self.x.div_assign(rhs.x);
        self.y.div_assign(rhs.y);
        self.z.div_assign(rhs.z);
    }
}

impl Div<SimdReal> for SimdVec3 {
    type Output = Self;
    #[inline]
    fn div(self, rhs: SimdReal) -> Self {
        Self {
            x: self.x.div(rhs),
            y: self.y.div(rhs),
            z: self.z.div(rhs),
        }
    }
}

impl DivAssign<SimdReal> for SimdVec3 {
    #[inline]
    fn div_assign(&mut self, rhs: SimdReal) {
        self.x.div_assign(rhs);
        self.y.div_assign(rhs);
        self.z.div_assign(rhs);
    }
}

impl Div<SimdVec3> for SimdReal {
    type Output = SimdVec3;
    #[inline]
    fn div(self, rhs: SimdVec3) -> SimdVec3 {
        SimdVec3 {
            x: self.div(rhs.x),
            y: self.div(rhs.y),
            z: self.div(rhs.z),
        }
    }
}

impl Rem<SimdVec3> for SimdVec3 {
    type Output = Self;
    #[inline]
    fn rem(self, rhs: Self) -> Self {
        Self {
            x: self.x.rem(rhs.x),
            y: self.y.rem(rhs.y),
            z: self.z.rem(rhs.z),
        }
    }
}

impl RemAssign<SimdVec3> for SimdVec3 {
    #[inline]
    fn rem_assign(&mut self, rhs: Self) {
        self.x.rem_assign(rhs.x);
        self.y.rem_assign(rhs.y);
        self.z.rem_assign(rhs.z);
    }
}

impl Rem<SimdReal> for SimdVec3 {
    type Output = Self;
    #[inline]
    fn rem(self, rhs: SimdReal) -> Self {
        Self {
            x: self.x.rem(rhs),
            y: self.y.rem(rhs),
            z: self.z.rem(rhs),
        }
    }
}

impl RemAssign<SimdReal> for SimdVec3 {
    #[inline]
    fn rem_assign(&mut self, rhs: SimdReal) {
        self.x.rem_assign(rhs);
        self.y.rem_assign(rhs);
        self.z.rem_assign(rhs);
    }
}

impl Rem<SimdVec3> for SimdReal {
    type Output = SimdVec3;
    #[inline]
    fn rem(self, rhs: SimdVec3) -> SimdVec3 {
        SimdVec3 {
            x: self.rem(rhs.x),
            y: self.rem(rhs.y),
            z: self.rem(rhs.z),
        }
    }
}

impl Neg for SimdVec3 {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl SimdQuat {
    pub const IDENTITY: Self = Self {
        x: AutoSimd([0.0; 4]),
        y: AutoSimd([0.0; 4]),
        z: AutoSimd([0.0; 4]),
        w: AutoSimd([1.0; 4]),
    };

    pub fn splat(quat: Quat) -> Self {
        Self {
            x: SimdReal::splat(quat.x),
            y: SimdReal::splat(quat.y),
            z: SimdReal::splat(quat.z),
            w: SimdReal::splat(quat.w),
        }
    }

    pub fn extract(&self, lane: usize) -> Quat {
        Quat::from_xyzw(
            self.x.extract(lane),
            self.y.extract(lane),
            self.z.extract(lane),
            self.w.extract(lane),
        )
    }

    pub fn from_scaled_axis(axis: Vector3) -> Self {
        Self::splat(Quat::from_scaled_axis(axis))
    }

    pub fn from_axis_angle(axis: Vector3, angle: Real) -> Self {
        Self::splat(Quat::from_axis_angle(axis, angle))
    }

    pub fn from_euler(euler: EulerRot, a: Real, b: Real, c: Real) -> Self {
        Self::splat(Quat::from_euler(euler, a, b, c))
    }

    pub fn normalize(self) -> Self {
        let length =
            (self.x * self.x + self.y * self.y + self.y * self.y + self.z * self.z).simd_sqrt();
        Self {
            x: self.x / length,
            y: self.y / length,
            z: self.z / length,
            w: self.w / length,
        }
    }

    pub fn inverse(self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
            w: self.w,
        }
    }
}

impl Add<SimdQuat> for SimdQuat {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
            w: self.w + rhs.w,
        }
    }
}

impl Sub<SimdQuat> for SimdQuat {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
            w: self.w - rhs.w,
        }
    }
}

impl Mul<SimdReal> for SimdQuat {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: SimdReal) -> Self {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
            w: self.w * rhs,
        }
    }
}

impl Div<SimdReal> for SimdQuat {
    type Output = Self;
    #[inline]
    fn div(self, rhs: SimdReal) -> Self {
        Self {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
            w: self.w / rhs,
        }
    }
}

// TODO: Mul<Self>
impl Mul<SimdVec3> for SimdQuat {
    type Output = SimdVec3;
    #[inline]
    fn mul(self, v: SimdVec3) -> Self::Output {
        let u = SimdVec3 {
            x: self.x,
            y: self.y,
            z: self.z,
        };
        let two = SimdReal::splat(2.0);
        two * u.dot(v) * u + (self.w * self.w - u.dot(u)) * v + two * self.w * u.cross(v)
    }
}

impl Neg for SimdQuat {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        self * SimdReal::splat(-1.0)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serialize", derive(serde::Serialize, serde::Deserialize))]
pub struct SimdRotation2 {
    /// The cosine of the rotation angle in radians.
    pub cos: SimdReal,
    /// The sine of the rotation angle in radians.
    pub sin: SimdReal,
}

impl Default for SimdRotation2 {
    fn default() -> Self {
        Self::IDENTITY
    }
}

impl SimdRotation2 {
    /// Zero rotation.
    pub const IDENTITY: Self = Self {
        cos: AutoSimd([1.0; 4]),
        sin: AutoSimd([0.0; 4]),
    };

    /// Returns the cosine of the rotation in radians.
    pub fn cos(&self) -> SimdReal {
        self.cos
    }

    /// Returns the sine of the rotation in radians.
    pub fn sin(&self) -> SimdReal {
        self.sin
    }

    pub fn splat(rotation: Rotation2) -> Self {
        Self {
            cos: SimdReal::splat(rotation.cos),
            sin: SimdReal::splat(rotation.sin),
        }
    }

    pub fn extract(&self, lane: usize) -> Rotation2 {
        Rotation2 {
            cos: self.cos.extract(lane),
            sin: self.sin.extract(lane),
        }
    }

    /// Creates a [`SRotation2`] from radians.
    pub fn from_radians(radians: Real) -> Self {
        Self {
            cos: SimdReal::splat(radians.cos()),
            sin: SimdReal::splat(radians.sin()),
        }
    }

    /// Creates a [`SRotation2`] from the sine and cosine of an angle in radians.
    pub fn from_sin_cos(sin: Real, cos: Real) -> Self {
        Self {
            sin: SimdReal::splat(sin),
            cos: SimdReal::splat(cos),
        }
    }

    /// Creates a [`SRotation2`] from degrees.
    pub fn from_degrees(degrees: Real) -> Self {
        Self::from_radians(degrees.to_radians())
    }

    /// Inverts the rotation.
    pub fn inverse(&self) -> Self {
        Self {
            cos: self.cos,
            sin: -self.sin,
        }
    }
}

impl Add<Self> for SimdRotation2 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        self.mul(rhs)
    }
}

impl AddAssign<Self> for SimdRotation2 {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl Sub<Self> for SimdRotation2 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        self.mul(rhs.inverse())
    }
}

impl SubAssign<Self> for SimdRotation2 {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl Mul for SimdRotation2 {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        Self {
            cos: self.cos * rhs.cos() - self.sin * rhs.sin(),
            sin: self.sin * rhs.cos() + self.cos * rhs.sin(),
        }
    }
}

impl MulAssign for SimdRotation2 {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl Mul<SimdVec2> for SimdRotation2 {
    type Output = SimdVec2;
    fn mul(self, rhs: SimdVec2) -> Self::Output {
        SimdVec2 {
            x: rhs.x * self.cos() - rhs.y * self.sin(),
            y: rhs.x * self.sin() + rhs.y * self.cos(),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SimdRotation3(pub SimdQuat);

impl Default for SimdRotation3 {
    fn default() -> Self {
        SimdRotation3(SimdQuat::IDENTITY)
    }
}

impl Deref for SimdRotation3 {
    type Target = SimdQuat;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for SimdRotation3 {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl SimdRotation3 {
    /// Zero rotation.
    pub const IDENTITY: Self = Self(SimdQuat::IDENTITY);

    pub fn splat(rotation: Rotation3) -> Self {
        Self(SimdQuat::splat(rotation.0))
    }

    pub fn extract(&self, lane: usize) -> Rotation3 {
        Rotation3(self.0.extract(lane))
    }

    pub fn from_scaled_axis(axis: Vector3) -> Self {
        Self(SimdQuat::from_scaled_axis(axis))
    }

    pub fn from_axis_angle(axis: Vector3, angle: Real) -> Self {
        Self(SimdQuat::from_axis_angle(axis, angle))
    }

    pub fn from_euler(euler: EulerRot, a: Real, b: Real, c: Real) -> Self {
        Self(SimdQuat::from_euler(euler, a, b, c))
    }

    pub fn normalize(self) -> Self {
        Self(self.0.normalize())
    }

    /// Inverts the rotation.
    pub fn inverse(self) -> Self {
        Self(self.0.inverse())
    }
}

//TODO: Mul<Self>
impl Mul<SimdVec3> for SimdRotation3 {
    type Output = SimdVec3;
    fn mul(self, rhs: SimdVec3) -> Self::Output {
        self.0 * rhs
    }
}

#[derive(Clone, Copy, Default, Debug, PartialEq)]
pub struct SimdIso2 {
    pub translation: SimdVec2,
    pub rotation: SimdRotation2,
}

impl SimdIso2 {
    /// An identity isometry.
    pub const IDENTITY: Self = Self {
        translation: SimdVec2::ZERO,
        rotation: SimdRotation2::IDENTITY,
    };

    pub fn new(translation: Vector2, angle: Real) -> Self {
        Self {
            translation: SimdVec2::splat(translation),
            rotation: SimdRotation2::from_radians(angle),
        }
    }

    pub fn splat(iso: Isometry2) -> Self {
        Self {
            translation: SimdVec2::splat(iso.translation),
            rotation: SimdRotation2::splat(iso.rotation),
        }
    }

    pub fn extract(&self, lane: usize) -> Isometry2 {
        Isometry2 {
            translation: self.translation.extract(lane),
            rotation: self.rotation.extract(lane),
        }
    }

    pub fn from_translation(translation: Vector2) -> Self {
        Self {
            translation: SimdVec2::splat(translation),
            rotation: SimdRotation2::IDENTITY,
        }
    }

    pub fn from_rotation(rotation: Rotation2) -> Self {
        Self {
            translation: SimdVec2::ZERO,
            rotation: SimdRotation2::splat(rotation),
        }
    }

    pub fn inverse(self) -> Self {
        let inv_rot = self.rotation.inverse();
        Self {
            translation: inv_rot * -self.translation,
            rotation: inv_rot,
        }
    }

    pub fn inv_mul(self, rhs: SimdIso2) -> Self {
        let inv_rot = self.rotation.inverse();
        let delta_translation = rhs.translation - self.translation;
        Self {
            translation: inv_rot * delta_translation,
            rotation: inv_rot * rhs.rotation,
        }
    }

    pub fn transform_point(self, point: SimdVec2) -> SimdVec2 {
        self.translation + self.rotation * point
    }

    pub fn inverse_transform_point(self, point: SimdVec2) -> SimdVec2 {
        self.rotation.inverse() * (point - self.translation)
    }
}

impl Add for SimdIso2 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        Self {
            translation: self.translation + rhs.translation,
            rotation: self.rotation * rhs.rotation,
        }
    }
}

impl Sub for SimdIso2 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            translation: self.translation - rhs.translation,
            rotation: self.rotation * rhs.rotation.inverse(),
        }
    }
}

impl Mul for SimdIso2 {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        Self {
            translation: self.translation + self.rotation * rhs.translation,
            rotation: self.rotation * rhs.rotation,
        }
    }
}

impl Mul<SimdVec2> for SimdIso2 {
    type Output = SimdVec2;
    fn mul(self, rhs: SimdVec2) -> Self::Output {
        self.rotation * rhs
    }
}

#[derive(Clone, Copy, Default, Debug, PartialEq)]
pub struct SimdIso3 {
    pub translation: SimdVec3,
    pub rotation: SimdRotation3,
}

impl SimdIso3 {
    /// An identity isometry.
    pub const IDENTITY: Self = Self {
        translation: SimdVec3::ZERO,
        rotation: SimdRotation3::IDENTITY,
    };

    pub fn new(translation: Vector3, rotation: Vector3) -> Self {
        Self {
            translation: SimdVec3::splat(translation),
            rotation: SimdRotation3(SimdQuat::from_scaled_axis(rotation)),
        }
    }

    pub fn splat(iso: Isometry3) -> Self {
        Self {
            translation: SimdVec3::splat(iso.translation),
            rotation: SimdRotation3::splat(iso.rotation),
        }
    }

    pub fn extract(&self, lane: usize) -> Isometry3 {
        Isometry3 {
            translation: self.translation.extract(lane),
            rotation: self.rotation.extract(lane),
        }
    }

    pub fn from_translation(translation: Vector3) -> Self {
        Self {
            translation: SimdVec3::splat(translation),
            rotation: SimdRotation3::IDENTITY,
        }
    }

    pub fn from_rotation(rotation: Rotation3) -> Self {
        Self {
            translation: SimdVec3::ZERO,
            rotation: SimdRotation3::splat(rotation),
        }
    }

    pub fn inverse(self) -> Self {
        let inv_rot = self.rotation.inverse();
        Self {
            translation: inv_rot * -self.translation,
            rotation: inv_rot,
        }
    }

    pub fn transform_point(self, point: SimdVec3) -> SimdVec3 {
        self.translation + self.rotation * point
    }

    pub fn inverse_transform_point(self, point: SimdVec3) -> SimdVec3 {
        self.rotation.inverse() * (point - self.translation)
    }
}

impl Mul<SimdVec3> for SimdIso3 {
    type Output = SimdVec3;
    fn mul(self, rhs: SimdVec3) -> Self::Output {
        self.rotation * rhs
    }
}
