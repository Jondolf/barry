use core::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

use bevy_math::{EulerRot, Quat};

use super::{Matrix2, Matrix3, Real, UnitVector2, UnitVector3, Vector2, Vector3};

#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serialize", derive(serde::Serialize, serde::Deserialize))]
pub struct Rotation2 {
    /// The cosine of the rotation angle in radians.
    pub cos: Real,
    /// The sine of the rotation angle in radians.
    pub sin: Real,
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
        Self::from_scaled_rotation_arc_colinear(from, to, 1.0)
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
pub struct Rotation3(pub Quat);

impl Default for Rotation3 {
    fn default() -> Self {
        Rotation3(Quat::IDENTITY)
    }
}

impl Rotation3 {
    /// Zero rotation.
    pub const IDENTITY: Self = Self(Quat::IDENTITY);

    pub fn from_rotation_arc_colinear(from: Vector3, to: Vector3) -> Option<Self> {
        Self::from_scaled_rotation_arc_colinear(from, to, 1.0)
    }

    pub fn from_scaled_rotation_arc_colinear(
        from: Vector3,
        to: Vector3,
        scale: Real,
    ) -> Option<Self> {
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

    pub fn is_finite(self) -> bool {
        self.0.is_finite()
    }
}

impl From<Rotation3> for Matrix3 {
    fn from(rotation: Rotation3) -> Self {
        Matrix3::from_quat(rotation.0)
    }
}

impl From<Matrix3> for Rotation3 {
    fn from(rotation: Matrix3) -> Self {
        Rotation3(Quat::from_mat3(&rotation))
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
