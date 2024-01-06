use crate::math::{Isometry, Matrix, SimdIsometry, SimdMatrix, SimdVector, UnitVector, Vector};

/// Extra operations with isometries.
pub trait IsometryOps {
    /// Transform a vector by the absolute value of the homogeneous matrix
    /// equivalent to `self`.
    fn absolute_transform_vector(&self, v: Vector) -> Vector;
}

/// Extra operations with SIMD isometries.
pub trait SimdIsometryOps {
    /// Transform a vector by the absolute value of the homogeneous matrix
    /// equivalent to `self`.
    fn absolute_transform_vector(&self, v: SimdVector) -> SimdVector;
}

impl IsometryOps for Isometry {
    #[inline]
    fn absolute_transform_vector(&self, v: Vector) -> Vector {
        Matrix::from_cols_array(
            &Matrix::from(self.rotation)
                .to_cols_array()
                .map(|col| col.abs()),
        ) * v
    }
}

impl SimdIsometryOps for SimdIsometry {
    #[inline]
    fn absolute_transform_vector(&self, v: SimdVector) -> SimdVector {
        SimdMatrix::from(self.rotation).abs() * v
    }
}

/// Various operations usable with `Option<Isometry>` and `Option<Isometry>`
/// where `None` is assumed to be equivalent to the identity.
pub trait IsometryOpt {
    /// Computes `self.inverse() * rhs`.
    fn inv_mul(self, rhs: Isometry) -> Isometry;
    /// Computes `rhs * self`.
    fn prepend_to(self, rhs: Isometry) -> Isometry;
    /// Computes `self * p`.
    fn transform_point(self, p: Vector) -> Vector;
    /// Computes `self * v`.
    fn transform_vector(self, v: Vector) -> Vector;
    /// Computes `self * v`.
    fn transform_unit_vector(self, v: UnitVector) -> UnitVector;
    /// Computes `self.inverse() * p`.
    fn inverse_transform_point(self, p: Vector) -> Vector;
    /// Computes `self.inverse() * v`.
    fn inverse_transform_vector(self, v: Vector) -> Vector;
    /// Computes `self.inverse() * v`.
    fn inverse_transform_unit_vector(self, v: UnitVector) -> UnitVector;
}

impl IsometryOpt for Option<Isometry> {
    #[inline]
    fn inv_mul(self, rhs: Isometry) -> Isometry {
        if let Some(iso) = self {
            iso.inv_mul(rhs)
        } else {
            rhs
        }
    }

    #[inline]
    fn prepend_to(self, rhs: Isometry) -> Isometry {
        if let Some(iso) = self {
            rhs * iso
        } else {
            rhs
        }
    }

    #[inline]
    fn transform_point(self, p: Vector) -> Vector {
        if let Some(iso) = self {
            iso * p
        } else {
            p
        }
    }

    #[inline]
    fn transform_vector(self, v: Vector) -> Vector {
        if let Some(iso) = self {
            iso * v
        } else {
            v
        }
    }

    #[inline]
    fn transform_unit_vector(self, v: UnitVector) -> UnitVector {
        if let Some(iso) = self {
            iso * v
        } else {
            v
        }
    }

    #[inline]
    fn inverse_transform_point(self, p: Vector) -> Vector {
        if let Some(iso) = self {
            iso.inverse_transform_point(p)
        } else {
            p
        }
    }

    #[inline]
    fn inverse_transform_vector(self, v: Vector) -> Vector {
        if let Some(iso) = self {
            iso.rotation.inverse() * v
        } else {
            v
        }
    }

    #[inline]
    fn inverse_transform_unit_vector(self, v: UnitVector) -> UnitVector {
        if let Some(iso) = self {
            UnitVector::new(iso.rotation.inverse() * *v).unwrap()
        } else {
            v
        }
    }
}
