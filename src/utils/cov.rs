use crate::math::{Matrix, Real, Vector};

/// Computes the covariance matrix of a set of points.
pub fn cov(pts: &[Vector]) -> Matrix {
    center_cov(pts).1
}

/// Computes the center and the covariance matrix of a set of points.
pub fn center_cov(pts: &[Vector]) -> (Vector, Matrix) {
    let center = crate::utils::center(pts);
    let mut cov: Matrix = Matrix::ZERO;
    let normalizer: Real = 1.0 / (pts.len() as Real);

    for p in pts.iter() {
        let cp = *p - center;
        let ncp = cp * normalizer;
        #[cfg(feature = "dim2")]
        {
            cov += Matrix::from_cols(cp * ncp.x, cp * ncp.y);
        }
        #[cfg(feature = "dim3")]
        {
            cov += Matrix::from_cols(cp * ncp.x, cp * ncp.y, cp * ncp.z);
        }
    }

    (center, cov)
}
