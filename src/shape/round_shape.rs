use crate::math::{Real, UnitVector, Vector};
use crate::shape::SupportMap;

#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize),
    archive(check_bytes)
)]
#[cfg_attr(feature = "cuda", derive(cust_core::DeviceCopy))]
#[derive(Copy, Clone, Debug)]
#[repr(C)]
/// A shape with rounded borders.
pub struct RoundShape<S> {
    /// The shape being rounded.
    pub inner_shape: S,
    /// The radius of the rounded border.
    pub border_radius: Real,
}

impl<S: SupportMap> SupportMap for RoundShape<S> {
    fn local_support_point(&self, dir: Vector) -> Vector {
        self.local_support_point_toward(UnitVector::new(dir).unwrap())
    }

    fn local_support_point_toward(&self, dir: UnitVector) -> Vector {
        self.inner_shape.local_support_point_toward(dir) + *dir * self.border_radius
    }
}
