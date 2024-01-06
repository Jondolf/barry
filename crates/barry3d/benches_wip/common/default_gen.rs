use barry3d::bounding_volume::{Aabb, BoundingSphere};
use barry3d::math::{Isometry2, Isometry3, Matrix2, Matrix3, Vector2, Vector3};
use barry3d::math::{Real, Vector};
use barry3d::query::Ray;
use barry3d::shape::{Ball, Capsule, Cone, Cuboid, Cylinder, Segment, Triangle};
use rand::distributions::{Distribution, Standard};
use rand::Rng;

#[cfg(feature = "dim2")]
type ConvexHull = Vec<Vector>;
#[cfg(feature = "dim3")]
type ConvexHull = (Vec<Vector>, Vec<[u32; 3]>);

pub trait DefaultGen {
    fn generate<R: Rng>(rng: &mut R) -> Self;
}

pub fn generate<T: DefaultGen, R: Rng>(rng: &mut R) -> T {
    DefaultGen::generate(rng)
}

macro_rules! impl_rand_default_gen (
    ($t: ty) => {
        impl DefaultGen for $t {
            fn generate<R: Rng>(rng: &mut R) -> $t {
                rng.gen::<$t>()
            }
        }
    }
);

impl_rand_default_gen!(Vector2);
impl_rand_default_gen!(Vector3);
impl_rand_default_gen!(Matrix2);
impl_rand_default_gen!(Matrix3);
impl_rand_default_gen!(Isometry2);
impl_rand_default_gen!(Isometry3);
impl_rand_default_gen!(f32);
impl_rand_default_gen!(f64);
impl_rand_default_gen!(bool);

impl DefaultGen for Ball
where
    Standard: Distribution<Real>,
{
    fn generate<R: Rng>(rng: &mut R) -> Ball {
        Ball::new(rng.gen::<f32>().abs())
    }
}

impl DefaultGen for Cuboid
where
    Standard: Distribution<Real>,
{
    fn generate<R: Rng>(rng: &mut R) -> Cuboid {
        Cuboid::new(rng.gen::<Vector>().abs())
    }
}

impl DefaultGen for Capsule
where
    Standard: Distribution<Real>,
{
    fn generate<R: Rng>(rng: &mut R) -> Capsule {
        Capsule::new(
            rng.gen::<Vector>().abs(),
            rng.gen::<Vector>().abs(),
            rng.gen::<Real>().abs(),
        )
    }
}

impl DefaultGen for Cone
where
    Standard: Distribution<Real>,
{
    fn generate<R: Rng>(rng: &mut R) -> Cone {
        Cone::new(rng.gen::<Real>().abs(), rng.gen::<Real>().abs())
    }
}

impl DefaultGen for Cylinder
where
    Standard: Distribution<Real>,
{
    fn generate<R: Rng>(rng: &mut R) -> Cylinder {
        Cylinder::new(rng.gen::<Real>().abs(), rng.gen::<Real>().abs())
    }
}

impl DefaultGen for Segment
where
    Standard: Distribution<Real>,
{
    fn generate<R: Rng>(rng: &mut R) -> Segment {
        Segment::new(rng.gen(), rng.gen())
    }
}

impl DefaultGen for Triangle
where
    Standard: Distribution<Real>,
{
    fn generate<R: Rng>(rng: &mut R) -> Triangle {
        Triangle::new(rng.gen(), rng.gen(), rng.gen())
    }
}

impl DefaultGen for ConvexHull
where
    Standard: Distribution<Real>,
{
    fn generate<R: Rng>(rng: &mut R) -> ConvexHull {
        // It is recommended to have at most 100 points.
        // Otherwise, a smarter structure like the DK hierarchy would be needed.
        // let pts: Vec<_> = (0..100).map(|_| rng.gen()).collect();
        // ConvexHull::try_from_points(&pts).unwrap()
        unimplemented!()
    }
}

impl DefaultGen for Ray
where
    Standard: Distribution<Real>,
{
    fn generate<R: Rng>(rng: &mut R) -> Ray {
        // The generate ray will always point to the origin.
        let shift = rng.gen::<Vector>() * 10.0;
        Ray::new(Vector::ZERO + shift, -shift)
    }
}

impl DefaultGen for Aabb
where
    Standard: Distribution<Real>,
{
    fn generate<R: Rng>(rng: &mut R) -> Aabb {
        // an Aabb centered at the origin.
        let half_extents = rng.gen::<Vector>().abs();
        Aabb::new(Vector::ZERO + (-half_extents), Vector::ZERO + half_extents)
    }
}

impl DefaultGen for BoundingSphere
where
    Standard: Distribution<Real>,
{
    fn generate<R: Rng>(rng: &mut R) -> BoundingSphere {
        // a bounding sphere centered at the origin.
        BoundingSphere::new(Vector::ZERO, rng.gen::<Real>().abs())
    }
}
