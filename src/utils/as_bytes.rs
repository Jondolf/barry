use std::mem;
use std::slice;

use crate::math::{Real, Vector2, Vector3};

/// Trait that transforms thing to a slice of u8.
pub trait AsBytes {
    /// Converts `self` to a slice of bytes.
    fn as_bytes(&self) -> &[u8];
}

macro_rules! generic_as_bytes_impl(
    ($t: ident, $dimension: expr) => (
        impl AsBytes for $t {
            #[inline(always)]
            fn as_bytes<'a>(&'a self) -> &'a [u8] {
                unsafe {
                    slice::from_raw_parts(mem::transmute(self), mem::size_of::<Real>() * $dimension)
                }
            }
        }
    )
);

generic_as_bytes_impl!(Vector2, 2);
generic_as_bytes_impl!(Vector3, 3);

// FIXME: implement for all `T: Copy` instead?
