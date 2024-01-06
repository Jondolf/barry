//! Support mapping based HalfSpace shape.
use crate::math::{UnitVector, Vector};

#[cfg(feature = "rkyv")]
use rkyv::{bytecheck, CheckBytes};

/// A half-space delimited by an infinite plane.
#[derive(PartialEq, Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize, CheckBytes),
    archive(as = "Self")
)]
#[cfg_attr(feature = "cuda", derive(cust_core::DeviceCopy))]
#[repr(C)]
pub struct HalfSpace {
    /// The halfspace planar boundary's outward normal.
    pub normal: UnitVector,
}

impl HalfSpace {
    /// Builds a new halfspace from its center and its normal.
    #[inline]
    pub fn new(normal: UnitVector) -> HalfSpace {
        HalfSpace { normal }
    }

    /// Computes a scaled version of this half-space.
    ///
    /// Returns `None` if `self.normal` scaled by `scale` is zero (the scaled half-space
    /// degenerates to a single point).
    pub fn scaled(self, scale: Vector) -> Option<Self> {
        UnitVector::new(*self.normal * scale)
            .map(|normal| Self { normal })
            .ok()
    }
}
