use barry3d::math::Vector3;
use barry3d::shape::TriMesh;

fn main() {
    let points = vec![
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(-1.0, -0.5, 0.0),
        Vector3::new(0.0, -0.5, -1.0),
        Vector3::new(1.0, -0.5, 0.0),
    ];

    let indices = vec![[0u32, 1, 2], [0, 2, 3], [0, 3, 1]];

    // Build the mesh.
    let mesh = TriMesh::new(points, indices);

    assert!(mesh.vertices().len() == 4);
}
