use crate::math::{self, Real, Vector, Vector3, DIM};
use crate::shape::Ball;
use crate::transformation::utils;

impl Ball {
    /// Discretize the boundary of this ball as a triangle-mesh.
    pub fn to_trimesh(
        &self,
        ntheta_subdiv: u32,
        nphi_subdiv: u32,
    ) -> (Vec<Vector3>, Vec<[u32; 3]>) {
        let diameter = self.radius * 2.0;
        let (vtx, idx) = unit_sphere(ntheta_subdiv, nphi_subdiv);
        (utils::scaled(vtx, Vector::splat(diameter)), idx)
    }
}

fn unit_sphere(ntheta_subdiv: u32, nphi_subdiv: u32) -> (Vec<Vector3>, Vec<[u32; 3]>) {
    let dtheta = math::real_consts::TAU / (ntheta_subdiv as Real);
    let dphi = math::real_consts::PI / (nphi_subdiv as Real);

    let mut coords = Vec::new();
    let mut curr_phi: Real = -math::real_consts::FRAC_PI_2 + dphi;

    coords.push(Vector::new(0.0, -1.0, 0.0));

    for _ in 1..nphi_subdiv {
        utils::push_circle(
            curr_phi.cos(),
            ntheta_subdiv,
            dtheta,
            curr_phi.sin(),
            &mut coords,
        );
        curr_phi = curr_phi + dphi;
    }

    coords.push(Vector::new(0.0, 1.0, 0.0));

    let mut idx = Vec::new();

    utils::push_degenerate_top_ring_indices(1, 0, ntheta_subdiv, &mut idx);
    utils::reverse_clockwising(&mut idx);

    for i in 0..nphi_subdiv - 2 {
        utils::push_ring_indices(
            1 + i * ntheta_subdiv,
            1 + (i + 1) * ntheta_subdiv,
            ntheta_subdiv,
            &mut idx,
        );
    }

    utils::push_degenerate_top_ring_indices(
        coords.len() as u32 - 1 - ntheta_subdiv,
        coords.len() as u32 - 1,
        ntheta_subdiv,
        &mut idx,
    );

    (utils::scaled(coords, Vector::splat(0.5)), idx)
}

/// Creates an hemisphere with a diameter of 1.
pub(crate) fn unit_hemisphere(
    ntheta_subdiv: u32,
    nphi_subdiv: u32,
) -> (Vec<Vector>, Vec<[u32; DIM]>) {
    let dtheta = math::real_consts::TAU / (ntheta_subdiv as Real);
    let dphi = math::real_consts::FRAC_PI_2 / (nphi_subdiv as Real);

    let mut coords = Vec::new();
    let mut curr_phi: Real = 0.0;

    for _ in 0..nphi_subdiv {
        utils::push_circle(
            curr_phi.cos(),
            ntheta_subdiv,
            dtheta,
            curr_phi.sin(),
            &mut coords,
        );
        curr_phi = curr_phi + dphi;
    }

    coords.push(Vector::new(0.0, 1.0, 0.0));

    let mut idx = Vec::new();

    for i in 0..nphi_subdiv - 1 {
        utils::push_ring_indices(
            i * ntheta_subdiv,
            (i + 1) * ntheta_subdiv,
            ntheta_subdiv,
            &mut idx,
        );
    }

    utils::push_degenerate_top_ring_indices(
        (nphi_subdiv - 1) * ntheta_subdiv,
        coords.len() as u32 - 1,
        ntheta_subdiv,
        &mut idx,
    );

    (utils::scaled(coords, Vector::splat(0.5)), idx)
}
