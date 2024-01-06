use crate::math::*;
use std::ops::{Add, Mul};

use num::Zero;
#[cfg(feature = "rkyv")]
use rkyv::{bytecheck, CheckBytes};

/// A 2x2 symmetric-definite-positive matrix.
#[derive(Copy, Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize, CheckBytes),
    archive(as = "Self", bound(archive = "N: rkyv::Archive<Archived = N>"))
)]
pub struct SdpMatrix2 {
    /// The component at the first row and first column of this matrix.
    pub m11: f32,
    /// The component at the first row and second column of this matrix.
    pub m12: f32,
    /// The component at the second row and second column of this matrix.
    pub m22: f32,
}

impl SdpMatrix2 {
    /// A new SDP 2x2 matrix with the given components.
    ///
    /// Because the matrix is symmetric, only the lower off-diagonal component is required.
    pub fn new(m11: f32, m12: f32, m22: f32) -> Self {
        Self { m11, m12, m22 }
    }

    /// Build an `SdpMatrix2` structure from a plain matrix, assuming it is SDP.
    ///
    /// No check is performed to ensure `mat` is actually SDP.
    pub fn from_sdp_matrix(mat: Matrix2) -> Self {
        Self {
            m11: mat.x_axis.x,
            m12: mat.x_axis.y,
            m22: mat.y_axis.y,
        }
    }

    /// Create a new SDP matrix filled with zeros.
    pub fn zero() -> Self {
        Self {
            m11: 0.0,
            m12: 0.0,
            m22: 0.0,
        }
    }

    /// Create a new SDP matrix with its diagonal filled with `val`, and its off-diagonal elements set to zero.
    pub fn diagonal(val: f32) -> Self {
        Self {
            m11: val,
            m12: 0.0,
            m22: val,
        }
    }

    /// Adds `val` to the diagonal components of `self`.
    pub fn add_diagonal(&mut self, elt: f32) -> Self {
        Self {
            m11: self.m11 + elt,
            m12: self.m12,
            m22: self.m22 + elt,
        }
    }

    /// Compute the inverse of this SDP matrix without performing any inversibility check.
    pub fn inverse_unchecked(&self) -> Self {
        let determinant = self.m11 * self.m22 - self.m12 * self.m12;
        let m11 = self.m22 / determinant;
        let m12 = -self.m12 / determinant;
        let m22 = self.m11 / determinant;

        Self { m11, m12, m22 }
    }

    /// Convert this SDP matrix to a regular matrix representation.
    pub fn into_matrix(self) -> Matrix2 {
        Matrix2::from_cols_array(&[self.m11, self.m12, self.m12, self.m22])
    }
}

impl Add<SdpMatrix2> for SdpMatrix2 {
    type Output = Self;

    fn add(self, rhs: SdpMatrix2) -> Self {
        Self::new(self.m11 + rhs.m11, self.m12 + rhs.m12, self.m22 + rhs.m22)
    }
}

impl Mul<Vector2> for SdpMatrix2 {
    type Output = Vector2;

    fn mul(self, rhs: Vector2) -> Self::Output {
        Vector2::new(
            self.m11 * rhs.x + self.m12 * rhs.y,
            self.m12 * rhs.x + self.m22 * rhs.y,
        )
    }
}

impl Mul<Real> for SdpMatrix2 {
    type Output = SdpMatrix2;

    fn mul(self, rhs: Real) -> Self::Output {
        SdpMatrix2::new(self.m11 * rhs, self.m12 * rhs, self.m22 * rhs)
    }
}

/// A 3x3 symmetric-definite-positive matrix.
#[derive(Copy, Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize, CheckBytes),
    archive(as = "Self", bound(archive = "N: rkyv::Archive<Archived = N>"))
)]
pub struct SdpMatrix3 {
    /// The component at the first row and first column of this matrix.
    pub m11: f32,
    /// The component at the first row and second column of this matrix.
    pub m12: f32,
    /// The component at the first row and third column of this matrix.
    pub m13: f32,
    /// The component at the second row and second column of this matrix.
    pub m22: f32,
    /// The component at the second row and third column of this matrix.
    pub m23: f32,
    /// The component at the third row and third column of this matrix.
    pub m33: f32,
}

impl SdpMatrix3 {
    /// A new SDP 3x3 matrix with the given components.
    ///
    /// Because the matrix is symmetric, only the lower off-diagonal components is required.
    pub fn new(m11: f32, m12: f32, m13: f32, m22: f32, m23: f32, m33: f32) -> Self {
        Self {
            m11,
            m12,
            m13,
            m22,
            m23,
            m33,
        }
    }

    /// Build an `SdpMatrix3` structure from a plain matrix, assuming it is SDP.
    ///
    /// No check is performed to ensure `mat` is actually SDP.
    pub fn from_sdp_matrix(mat: Matrix3) -> Self {
        Self {
            m11: mat.x_axis.x,
            m12: mat.x_axis.y,
            m13: mat.x_axis.z,
            m22: mat.y_axis.y,
            m23: mat.y_axis.z,
            m33: mat.z_axis.z,
        }
    }

    /// Create a new SDP matrix filled with zeros.
    pub fn zero() -> Self {
        Self {
            m11: 0.0,
            m12: 0.0,
            m13: 0.0,
            m22: 0.0,
            m23: 0.0,
            m33: 0.0,
        }
    }

    /// Create a new SDP matrix with its diagonal filled with `val`, and its off-diagonal elements set to zero.
    pub fn diagonal(val: f32) -> Self {
        Self {
            m11: val,
            m12: 0.0,
            m13: 0.0,
            m22: val,
            m23: 0.0,
            m33: val,
        }
    }

    /// Are all components of this matrix equal to zero?
    pub fn is_zero(&self) -> bool {
        self.m11.is_zero()
            && self.m12.is_zero()
            && self.m13.is_zero()
            && self.m22.is_zero()
            && self.m23.is_zero()
            && self.m33.is_zero()
    }

    /// Compute the inverse of this SDP matrix without performing any inversibility check.
    pub fn inverse_unchecked(&self) -> Self {
        let minor_m12_m23 = self.m22 * self.m33 - self.m23 * self.m23;
        let minor_m11_m23 = self.m12 * self.m33 - self.m13 * self.m23;
        let minor_m11_m22 = self.m12 * self.m23 - self.m13 * self.m22;

        let determinant =
            self.m11 * minor_m12_m23 - self.m12 * minor_m11_m23 + self.m13 * minor_m11_m22;
        let inv_det = 1.0 / determinant;

        SdpMatrix3 {
            m11: minor_m12_m23 * inv_det,
            m12: -minor_m11_m23 * inv_det,
            m13: minor_m11_m22 * inv_det,
            m22: (self.m11 * self.m33 - self.m13 * self.m13) * inv_det,
            m23: (self.m13 * self.m12 - self.m23 * self.m11) * inv_det,
            m33: (self.m11 * self.m22 - self.m12 * self.m12) * inv_det,
        }
    }

    /// Adds `elt` to the diagonal components of `self`.
    pub fn add_diagonal(&self, elt: f32) -> Self {
        Self {
            m11: self.m11 + elt,
            m12: self.m12,
            m13: self.m13,
            m22: self.m22 + elt,
            m23: self.m23,
            m33: self.m33 + elt,
        }
    }
}

impl Add<SdpMatrix3> for SdpMatrix3 {
    type Output = SdpMatrix3;

    fn add(self, rhs: SdpMatrix3) -> Self::Output {
        SdpMatrix3 {
            m11: self.m11 + rhs.m11,
            m12: self.m12 + rhs.m12,
            m13: self.m13 + rhs.m13,
            m22: self.m22 + rhs.m22,
            m23: self.m23 + rhs.m23,
            m33: self.m33 + rhs.m33,
        }
    }
}

impl Mul<Real> for SdpMatrix3 {
    type Output = SdpMatrix3;

    fn mul(self, rhs: Real) -> Self::Output {
        SdpMatrix3 {
            m11: self.m11 * rhs,
            m12: self.m12 * rhs,
            m13: self.m13 * rhs,
            m22: self.m22 * rhs,
            m23: self.m23 * rhs,
            m33: self.m33 * rhs,
        }
    }
}

impl Mul<Vector3> for SdpMatrix3 {
    type Output = Vector3;

    fn mul(self, rhs: Vector3) -> Self::Output {
        let x = self.m11 * rhs.x + self.m12 * rhs.y + self.m13 * rhs.z;
        let y = self.m12 * rhs.x + self.m22 * rhs.y + self.m23 * rhs.z;
        let z = self.m13 * rhs.x + self.m23 * rhs.y + self.m33 * rhs.z;
        Vector3::new(x, y, z)
    }
}

impl Mul<Matrix3> for SdpMatrix3 {
    type Output = Matrix3;

    fn mul(self, rhs: Matrix3) -> Self::Output {
        let x0 = self.m11 * rhs.x_axis.x + self.m12 * rhs.y_axis.x + self.m13 * rhs.z_axis.x;
        let y0 = self.m12 * rhs.x_axis.x + self.m22 * rhs.y_axis.x + self.m23 * rhs.z_axis.x;
        let z0 = self.m13 * rhs.x_axis.x + self.m23 * rhs.y_axis.x + self.m33 * rhs.z_axis.x;

        let x1 = self.m11 * rhs.x_axis.y + self.m12 * rhs.y_axis.y + self.m13 * rhs.z_axis.y;
        let y1 = self.m12 * rhs.x_axis.y + self.m22 * rhs.y_axis.y + self.m23 * rhs.z_axis.y;
        let z1 = self.m13 * rhs.x_axis.y + self.m23 * rhs.y_axis.y + self.m33 * rhs.z_axis.y;

        let x2 = self.m11 * rhs.x_axis.z + self.m12 * rhs.y_axis.z + self.m13 * rhs.z_axis.z;
        let y2 = self.m12 * rhs.x_axis.z + self.m22 * rhs.y_axis.z + self.m23 * rhs.z_axis.z;
        let z2 = self.m13 * rhs.x_axis.z + self.m23 * rhs.y_axis.z + self.m33 * rhs.z_axis.z;

        Matrix3::from_cols_array(&[x0, x1, x2, y0, y1, y2, z0, z1, z2])
    }
}

#[cfg(feature = "simd-nightly")]
impl From<[SdpMatrix3; 8]> for SdpMatrix3<simba::simd::f32x8> {
    fn from(data: [SdpMatrix3; 8]) -> Self {
        SdpMatrix3 {
            m11: simba::simd::f32x8::from([
                data[0].m11,
                data[1].m11,
                data[2].m11,
                data[3].m11,
                data[4].m11,
                data[5].m11,
                data[6].m11,
                data[7].m11,
            ]),
            m12: simba::simd::f32x8::from([
                data[0].m12,
                data[1].m12,
                data[2].m12,
                data[3].m12,
                data[4].m12,
                data[5].m12,
                data[6].m12,
                data[7].m12,
            ]),
            m13: simba::simd::f32x8::from([
                data[0].m13,
                data[1].m13,
                data[2].m13,
                data[3].m13,
                data[4].m13,
                data[5].m13,
                data[6].m13,
                data[7].m13,
            ]),
            m22: simba::simd::f32x8::from([
                data[0].m22,
                data[1].m22,
                data[2].m22,
                data[3].m22,
                data[4].m22,
                data[5].m22,
                data[6].m22,
                data[7].m22,
            ]),
            m23: simba::simd::f32x8::from([
                data[0].m23,
                data[1].m23,
                data[2].m23,
                data[3].m23,
                data[4].m23,
                data[5].m23,
                data[6].m23,
                data[7].m23,
            ]),
            m33: simba::simd::f32x8::from([
                data[0].m33,
                data[1].m33,
                data[2].m33,
                data[3].m33,
                data[4].m33,
                data[5].m33,
                data[6].m33,
                data[7].m33,
            ]),
        }
    }
}

#[cfg(feature = "simd-nightly")]
impl From<[SdpMatrix3; 16]> for SdpMatrix3<simba::simd::f32x16> {
    fn from(data: [SdpMatrix3; 16]) -> Self {
        SdpMatrix3 {
            m11: simba::simd::f32x16::from([
                data[0].m11,
                data[1].m11,
                data[2].m11,
                data[3].m11,
                data[4].m11,
                data[5].m11,
                data[6].m11,
                data[7].m11,
                data[8].m11,
                data[9].m11,
                data[10].m11,
                data[11].m11,
                data[12].m11,
                data[13].m11,
                data[14].m11,
                data[15].m11,
            ]),
            m12: simba::simd::f32x16::from([
                data[0].m12,
                data[1].m12,
                data[2].m12,
                data[3].m12,
                data[4].m12,
                data[5].m12,
                data[6].m12,
                data[7].m12,
                data[8].m12,
                data[9].m12,
                data[10].m12,
                data[11].m12,
                data[12].m12,
                data[13].m12,
                data[14].m12,
                data[15].m12,
            ]),
            m13: simba::simd::f32x16::from([
                data[0].m13,
                data[1].m13,
                data[2].m13,
                data[3].m13,
                data[4].m13,
                data[5].m13,
                data[6].m13,
                data[7].m13,
                data[8].m13,
                data[9].m13,
                data[10].m13,
                data[11].m13,
                data[12].m13,
                data[13].m13,
                data[14].m13,
                data[15].m13,
            ]),
            m22: simba::simd::f32x16::from([
                data[0].m22,
                data[1].m22,
                data[2].m22,
                data[3].m22,
                data[4].m22,
                data[5].m22,
                data[6].m22,
                data[7].m22,
                data[8].m22,
                data[9].m22,
                data[10].m22,
                data[11].m22,
                data[12].m22,
                data[13].m22,
                data[14].m22,
                data[15].m22,
            ]),
            m23: simba::simd::f32x16::from([
                data[0].m23,
                data[1].m23,
                data[2].m23,
                data[3].m23,
                data[4].m23,
                data[5].m23,
                data[6].m23,
                data[7].m23,
                data[8].m23,
                data[9].m23,
                data[10].m23,
                data[11].m23,
                data[12].m23,
                data[13].m23,
                data[14].m23,
                data[15].m23,
            ]),
            m33: simba::simd::f32x16::from([
                data[0].m33,
                data[1].m33,
                data[2].m33,
                data[3].m33,
                data[4].m33,
                data[5].m33,
                data[6].m33,
                data[7].m33,
                data[8].m33,
                data[9].m33,
                data[10].m33,
                data[11].m33,
                data[12].m33,
                data[13].m33,
                data[14].m33,
                data[15].m33,
            ]),
        }
    }
}
