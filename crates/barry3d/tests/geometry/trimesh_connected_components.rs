use barry3d::math::Vector;
use barry3d::shape::{TriMesh, TriMeshFlags};

#[test]
// From https://github.com/dimforge/barry/issues/115
fn mesh_connected_components_grouped_faces() {
    let verts = vec![
        // Face 0
        Vector::new(15.82, 6.455, -0.15), // <- Vertex shared with face 1.
        Vector::new(9.915, 6.455, -0.15),
        Vector::new(9.915, 6.4, 0.0), // <- Vertex shared with face 1.
        // Face1
        Vector::new(15.82, 6.455, -0.15), // <- Vertex shared with face 0.
        Vector::new(9.915, 6.4, 0.0),     // <- Vertex shared with face 0.
        Vector::new(15.82, 6.4, 0.0),
    ];

    let mut roof = TriMesh::new(verts, vec![[0, 1, 2], [3, 4, 5]]);

    if let Err(e) =
        roof.set_flags(TriMeshFlags::MERGE_DUPLICATE_VERTICES | TriMeshFlags::CONNECTED_COMPONENTS)
    {
        dbg!(e);
        assert!(false);
    }

    let components = roof.connected_components().unwrap();
    println!("components: {:?}", components);
    assert_eq!(components.ranges.len(), 2); // Only one connected-component (two range values).
    assert_eq!(components.grouped_faces.len(), 2); // Only two faces in the connected-component.
}
