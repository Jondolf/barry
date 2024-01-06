use crate::mass_properties::MassProperties;
use crate::math::{self, PrincipalAngularInertia, Real, Vector};

impl MassProperties {
    pub(crate) fn ball_volume_unit_angular_inertia(
        radius: Real,
    ) -> (Real, PrincipalAngularInertia) {
        #[cfg(feature = "dim2")]
        {
            let volume = math::real_consts::PI * radius * radius;
            let i = radius * radius / 2.0;
            (volume, i)
        }
        #[cfg(feature = "dim3")]
        {
            let volume = math::real_consts::PI * radius * radius * radius * 4.0 / 3.0;
            let i = radius * radius * 2.0 / 5.0;

            (volume, Vector::splat(i))
        }
    }

    /// Computes the mass properties of a ball.
    pub fn from_ball(density: Real, radius: Real) -> Self {
        let (vol, unit_i) = Self::ball_volume_unit_angular_inertia(radius);
        let mass = vol * density;
        Self::new(Vector::ZERO, mass, unit_i * mass)
    }
}
