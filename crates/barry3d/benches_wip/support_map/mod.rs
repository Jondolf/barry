use crate::common::{generate, unref};
use barry3d::math::{Isometry3, Vector3};
use barry3d::shape::SupportMap;
use barry3d::shape::{Ball, Capsule, Cone, ConvexHull, Cuboid, Cylinder, Segment, Triangle};
use rand::SeedableRng;
use rand_isaac::IsaacRng;
use test::Bencher;

#[path = "../common/macros.rs"]
#[macro_use]
mod macros;

bench_method!(
    bench_ball_support_map,
    support_point,
    c: Ball,
    m: Isometry3<f32>,
    dir: Vector3<f32>
);
bench_method!(
    bench_cuboid_support_map,
    support_point,
    c: Cuboid,
    m: Isometry3<f32>,
    dir: Vector3<f32>
);
bench_method!(
    bench_capsule_support_map,
    support_point,
    c: Capsule,
    m: Isometry3<f32>,
    dir: Vector3<f32>
);
bench_method!(
    bench_cone_support_map,
    support_point,
    c: Cone,
    m: Isometry3<f32>,
    dir: Vector3<f32>
);
bench_method!(
    bench_cylinder_support_map,
    support_point,
    c: Cylinder,
    m: Isometry3<f32>,
    dir: Vector3<f32>
);
bench_method!(
    bench_segment_support_map,
    support_point,
    c: Segment,
    m: Isometry3<f32>,
    dir: Vector3<f32>
);
bench_method!(
    bench_triangle_support_map,
    support_point,
    c: Triangle,
    m: Isometry3<f32>,
    dir: Vector3<f32>
);
bench_method!(
    bench_convex_support_map,
    support_point,
    c: ConvexHull,
    m: Isometry3<f32>,
    dir: Vector3<f32>
);
