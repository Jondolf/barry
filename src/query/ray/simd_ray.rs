use crate::math::SimdVector;
use crate::query::Ray;

/// A structure representing 4 rays in an SIMD SoA fashion.
#[derive(Debug, Copy, Clone)]
pub struct SimdRay {
    /// The origin of the rays represented as a single SIMD point.
    pub origin: SimdVector,
    /// The direction of the rays represented as a single SIMD vector.
    pub dir: SimdVector,
}

impl SimdRay {
    /// Creates a new SIMD ray with all its lanes filled with the same ray.
    pub fn splat(ray: Ray) -> Self {
        Self {
            origin: SimdVector::splat(ray.origin),
            dir: SimdVector::splat(ray.dir),
        }
    }
}
