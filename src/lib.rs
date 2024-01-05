/*!
barry
========

**barry** is a 2 and 3-dimensional geometric library written with
the rust programming language.

*/

#![deny(non_camel_case_types)]
#![deny(unused_parens)]
#![deny(non_upper_case_globals)]
#![deny(unused_results)]
#![warn(missing_docs)] // TODO: deny this
#![warn(unused_imports)]
#![allow(missing_copy_implementations)]
#![doc(html_root_url = "http://docs.rs/barry/0.1.1")]
#![cfg_attr(not(feature = "std"), no_std)]
#![cfg_attr(not(feature = "rkyv"), deny(unused_qualifications))] // TODO: deny that everytime

#[cfg(all(
    feature = "simd-is-enabled",
    not(feature = "simd-stable"),
    not(feature = "simd-nightly")
))]
std::compile_error!("The `simd-is-enabled` feature should not be enabled explicitly. Please enable the `simd-stable` or the `simd-nightly` feature instead.");
#[cfg(all(feature = "simd-is-enabled", feature = "enhanced-determinism"))]
std::compile_error!(
    "SIMD cannot be enabled when the `enhanced-determinism` feature is also enabled."
);

macro_rules! array(
    ($callback: expr; SIMD_WIDTH) => {
        {
            #[inline(always)]
            #[allow(dead_code)]
            fn create_arr<T>(mut callback: impl FnMut(usize) -> T) -> [T; SIMD_WIDTH] {
                [callback(0usize), callback(1usize), callback(2usize), callback(3usize)]
            }

            create_arr($callback)
        }
    }
);

#[cfg(all(feature = "alloc", not(feature = "std")))]
#[cfg_attr(test, macro_use)]
extern crate alloc;

#[cfg(not(feature = "std"))]
extern crate core as std;

#[cfg(feature = "serde")]
#[macro_use]
extern crate serde;
#[macro_use]
extern crate approx;
#[macro_use]
#[cfg(feature = "dim3")]
extern crate bitflags;
extern crate num_traits as num;

pub extern crate either;
pub extern crate simba;

pub mod bounding_volume;
pub mod mass_properties;
pub mod math;
pub mod partitioning;
pub mod query;
pub mod shape;
#[cfg(feature = "std")]
pub mod transformation;
pub mod utils;

/// A trait for getting the index of the smallest and largest value in a collection.
pub trait MinMaxIndex {
    fn min_index(&self) -> usize;
    fn max_index(&self) -> usize;
}

#[cfg(not(feature = "simd-is-enabled"))]
mod simd {
    use simba::simd::AutoBoolx4;
    /// The number of lanes of a SIMD number.
    pub const SIMD_WIDTH: usize = 4;
    /// SIMD_WIDTH - 1
    pub const SIMD_LAST_INDEX: usize = 3;

    /// A SIMD float with SIMD_WIDTH lanes.
    #[cfg(feature = "f32")]
    pub type SimdReal = simba::simd::AutoF32x4;

    /// A SIMD float with SIMD_WIDTH lanes.
    #[cfg(feature = "f64")]
    pub type SimdReal = simba::simd::AutoF64x4;

    /// A SIMD bool with SIMD_WIDTH lanes.
    pub type SimdBool = AutoBoolx4;
}

#[cfg(feature = "simd-is-enabled")]
mod simd {
    #[cfg(all(feature = "simd-nightly", feature = "f32"))]
    pub use simba::simd::{f32x4 as SimdReal, m32x4 as SimdBool};
    #[cfg(all(feature = "simd-stable", feature = "f32"))]
    pub use simba::simd::{WideBoolF32x4 as SimdBool, WideF32x4 as SimdReal};

    #[cfg(all(feature = "simd-nightly", feature = "f64"))]
    pub use simba::simd::{f64x4 as SimdReal, m64x4 as SimdBool};
    #[cfg(all(feature = "simd-stable", feature = "f64"))]
    pub use simba::simd::{WideBoolF64x4 as SimdBool, WideF64x4 as SimdReal};

    /// The number of lanes of a SIMD number.
    pub const SIMD_WIDTH: usize = 4;
    /// SIMD_WIDTH - 1
    pub const SIMD_LAST_INDEX: usize = 3;
}
