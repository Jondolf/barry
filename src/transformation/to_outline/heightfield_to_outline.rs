use crate::{
    math::Vector3,
    shape::{GenericHeightField, HeightFieldStorage},
};

impl<Storage: HeightFieldStorage> GenericHeightField<Storage> {
    /// Outlines this heightfield’s shape using polylines.
    pub fn to_outline(&self) -> (Vec<Vector3>, Vec<[u32; 2]>) {
        todo!()
    }
}
