use crate::math::{Isometry, Real, UnitVector, Vector};
use crate::shape::SupportMap;

/// A support mapping that is a single point.
pub struct ConstantPoint(pub Vector);

impl SupportMap for ConstantPoint {
    #[inline]
    fn support_point(&self, m: Isometry, _: Vector) -> Vector {
        m * self.0
    }

    #[inline]
    fn support_point_toward(&self, m: Isometry, _: UnitVector) -> Vector {
        m * self.0
    }

    #[inline]
    fn local_support_point(&self, _: Vector) -> Vector {
        self.0
    }

    #[inline]
    fn local_support_point_toward(&self, _: UnitVector) -> Vector {
        self.0
    }
}

/// A support mapping that is the point at (0.0, 0.0, 0.0).
pub struct ConstantOrigin;

impl SupportMap for ConstantOrigin {
    #[inline]
    fn support_point(&self, m: Isometry, _: Vector) -> Vector {
        m.translation.into()
    }

    #[inline]
    fn support_point_toward(&self, m: Isometry, _: UnitVector) -> Vector {
        m.translation.into()
    }

    #[inline]
    fn local_support_point(&self, _: Vector) -> Vector {
        Vector::ZERO
    }

    #[inline]
    fn local_support_point_toward(&self, _: UnitVector) -> Vector {
        Vector::ZERO
    }
}

/// The Minkowski sum of a shape and a ball.
pub struct DilatedShape<'a, S: ?Sized + SupportMap> {
    /// The shape involved in the Minkowski sum.
    pub shape: &'a S,
    /// The radius of the ball involved in the Minkoski sum.
    pub radius: Real,
}

impl<'a, S: ?Sized + SupportMap> SupportMap for DilatedShape<'a, S> {
    #[inline]
    fn support_point(&self, m: Isometry, dir: Vector) -> Vector {
        self.support_point_toward(m, UnitVector::new(dir).unwrap())
    }

    #[inline]
    fn support_point_toward(&self, m: Isometry, dir: UnitVector) -> Vector {
        self.shape.support_point_toward(m, dir) + *dir * self.radius
    }

    #[inline]
    fn local_support_point(&self, dir: Vector) -> Vector {
        self.local_support_point_toward(UnitVector::new(dir).unwrap())
    }

    #[inline]
    fn local_support_point_toward(&self, dir: UnitVector) -> Vector {
        self.shape.local_support_point_toward(dir) + *dir * self.radius
    }
}
