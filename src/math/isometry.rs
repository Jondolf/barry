use core::ops::{Add, Mul, Sub};

use bevy_math::Quat;

use super::{Real, Rotation2, Rotation3, UnitVector2, UnitVector3, Vector2, Vector3};

#[derive(Clone, Copy, Default, Debug, PartialEq)]
pub struct Iso2 {
    pub translation: Vector2,
    pub rotation: Rotation2,
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
    pub translation: Vector3,
    pub rotation: Rotation3,
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
