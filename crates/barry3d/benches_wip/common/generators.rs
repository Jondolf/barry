use barry3d::math::{Vector2, Vector3};
use barry3d::shape::TriMesh;
use rand::Rng;

pub fn generate_trimesh_around_origin<R: Rng>(rng: &mut R) -> TriMesh {
    let pts = (0..3000).map(|_| rng.gen::<Vector3<f32>>() * 3.0).collect();
    let uvs = (0..3000).map(|_| rng.gen::<Vector2<f32>>() * 3.0).collect();
    let indices = (0..1000)
        .map(|i| Vector3::new(i * 3, i * 3 + 1, i * 3 + 2))
        .collect();

    TriMesh::new(pts, indices, Some(uvs))
}
