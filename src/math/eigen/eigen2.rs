use bevy_math::Vec2Swizzles;

use crate::math::{Matrix2, Real, Vector2};

/// The eigen decomposition of a symmetric 2x2 matrix.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serialize", derive(serde::Serialize, serde::Deserialize))]
pub struct SymmetricEigen2 {
    /// The eigenvalues of the symmetric 2x2 matrix.
    pub eigenvalues: Vector2,
    /// The two eigenvectors of the symmetric 2x2 matrix.
    pub eigenvectors: Matrix2,
}

impl SymmetricEigen2 {
    /// Computes the eigen decomposition of the given symmetric 2x2 matrix.
    ///
    /// The eigenvalues are returned in ascending order `eigen1 < eigen2 < eigen3`.
    /// This can be reversed with the [`reverse`](Self::reverse) method.
    pub fn new(mat: Matrix2) -> Self {
        let eigenvalues = Self::eigenvalues(mat);
        let eigenvector1 = Self::eigenvector(mat, eigenvalues.x);
        let eigenvector2 = Self::eigenvector(mat, eigenvalues.y);

        Self {
            eigenvalues,
            eigenvectors: Matrix2::from_cols(eigenvector1, eigenvector2),
        }
    }

    /// Reverses the order of the eigenvalues and their corresponding eigenvectors.
    pub fn reverse(&self) -> Self {
        Self {
            eigenvalues: self.eigenvalues.yx(),
            eigenvectors: Matrix2::from_cols(self.eigenvectors.y_axis, self.eigenvectors.x_axis),
        }
    }

    /// Computes the eigenvalues of a symmetric 2x2 matrix.
    ///
    /// Reference: <https://croninprojects.org/Vince/Geodesy/FindingEigenvectors.pdf>
    pub fn eigenvalues(mat: Matrix2) -> Vector2 {
        let [a, b, c] = [
            1.0,
            -(mat.x_axis.x + mat.y_axis.y),
            mat.x_axis.x * mat.y_axis.y - mat.x_axis.y * mat.y_axis.x,
        ];
        // The eigenvalues are the roots of the quadratic equation:
        // ax^2 + bx + c = 0
        // x = (-b ± sqrt(b^2 - 4ac)) / 2a
        let sqrt_part = (b.powi(2) - 4.0 * a * c).sqrt();
        let eigen1 = (-b + sqrt_part) / (2.0 * a);
        let eigen2 = (-b - sqrt_part) / (2.0 * a);
        Vector2::new(eigen1, eigen2)
    }

    /// Computes the unit-length eigenvector corresponding to the given `eigenvalue`
    /// of the symmetric 2x2 `mat`.
    ///
    /// Reference: <https://croninprojects.org/Vince/Geodesy/FindingEigenvectors.pdf>
    pub fn eigenvector(mat: Matrix2, eigenvalue: Real) -> Vector2 {
        Vector2::new(1.0, (eigenvalue - mat.x_axis.x) / mat.y_axis.x).normalize()
    }
}

#[cfg(test)]
mod test {
    use crate::math::{Matrix2, SymmetricEigen2, Vector2};

    #[test]
    fn eigen_2x2() {
        let mat = Matrix2::from_cols_array_2d(&[[6.0, 3.0], [3.0, 4.0]]);
        let eigen = SymmetricEigen2::new(mat);

        assert_relative_eq!(
            eigen.eigenvalues,
            Vector2::new(8.16228, 1.83772),
            epsilon = 0.001
        );
        assert_relative_eq!(
            Matrix2::from_cols(eigen.eigenvectors.x_axis, eigen.eigenvectors.y_axis,),
            Matrix2::from_cols(
                Vector2::new(0.811242, 0.58471),
                Vector2::new(0.58471, -0.811242),
            ),
            epsilon = 0.001
        );
    }
}
