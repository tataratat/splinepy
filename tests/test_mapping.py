try:
    from . import common as c
except BaseException:
    import common as c


class TestGeometryMapping(c.unittest.TestCase):
    """Test mapping of basis function derivatives"""

    def setUp(self) -> None:
        """Define geometries the function will be mapped to

        Create three simple geometries, two where results are known, and one
        with no zero entries in the jacobian.

        Functions will be 2D and 3D only
        """
        self.scaling_factors = 1 / c.np.array([2.0, 0.3, 1.5])
        self.scaling3D = c.splinepy.Bezier(
            degrees=[1, 1, 1],
            control_points=[
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [0.0, 0.3, 0.0],
                [2.0, 0.3, 0.0],
                [0.0, 0.0, 1.5],
                [2.0, 0.0, 1.5],
                [0.0, 0.3, 1.5],
                [2.0, 0.3, 1.5],
            ],
        )
        # Spline that rotates field by 0.17 radiants
        cc, ss = c.np.cos(0.17), c.np.sin(0.17)
        self.rotation_matrix = c.np.array(((cc, -ss), (ss, cc)))
        self.rotating2D = c.splinepy.Bezier(
            degrees=[1, 1],
            control_points=[
                [0, 0],
                [1, 0],
                [0, 1],
                [1, 1],
            ],
        )
        self.rotating2D.control_points = c.np.einsum(
            "ij,qj->qi", self.rotation_matrix, self.rotating2D.control_points
        )

        # Complicated spline
        self.askew_spline2D = c.splinepy.BSpline(
            degrees=[2, 2],
            control_points=[
                [0.0, 0.0],
                [1.0, 0.5],
                [2.0, 0.2],
                [0.5, 1.5],
                [1.0, 1.5],
                [1.5, 1.5],
                [0.0, 3.0],
                [1.0, 2.5],
                [2.0, 3.0],
            ],
            knot_vectors=[[0, 0, 0, 1, 1, 1], [0, 0, 0, 1, 1, 1]],
        )

        # Two solution fields
        self.solution_field_mono3D = c.splinepy.Bezier(
            degrees=[2, 1, 1], control_points=c.np.ones((12, 1)) * 0.2
        )
        self.solution_field_rando = c.splinepy.Bezier(
            degrees=[2, 1], control_points=c.np.random.rand(6, 1)
        )
        self.solution_field_rando2D = c.splinepy.Bezier(
            degrees=[2, 1], control_points=c.np.random.rand(6, 2)
        )

        self.query_points2D = c.np.random.rand(13, 2)
        self.query_points3D = c.np.random.rand(17, 3)

    def test_cross_evaluation_of_different_implementations(self):
        """Divergence and Laplacian are implemented differently when gradients
        are called at the same time
        """
        # Test Basis functions first
        mapper = self.solution_field_rando2D.geometry_mapper(
            self.askew_spline2D
        )
        bf_results = mapper.basis_function_derivatives(
            self.query_points2D, gradient=True, hessian=True, laplacian=True
        )
        bf_gradient, support_gradient = mapper.basis_gradient_and_support(
            self.query_points2D
        )
        bf_hessian, support_hessian = mapper.basis_hessian_and_support(
            self.query_points2D
        )
        bf_laplace, support_laplacian = mapper.basis_laplacian_and_support(
            self.query_points2D
        )

        # Laplacian is computed differently depending on function call
        self.assertTrue(c.np.allclose(bf_results["laplacian"], bf_laplace))
        self.assertTrue(
            c.np.allclose(bf_results["support"], support_laplacian)
        )

        # Check function calls
        self.assertTrue(c.np.allclose(bf_results["gradient"], bf_gradient))
        self.assertTrue(c.np.allclose(bf_results["support"], support_gradient))
        self.assertTrue(c.np.allclose(bf_results["hessian"], bf_hessian))
        self.assertTrue(c.np.allclose(bf_results["support"], support_hessian))

        # Try field values
        results = mapper.field_derivatives(
            self.query_points2D,
            gradient=True,
            divergence=True,
            hessian=True,
            laplacian=True,
            basis_function_values=True,
        )
        laplacian = mapper.laplacian(self.query_points2D)
        divergence = mapper.divergence(self.query_points2D)
        gradient = mapper.gradient(self.query_points2D)
        hessian = mapper.hessian(self.query_points2D)

        # Divergence and laplacian have different implementations
        self.assertTrue(c.np.allclose(results["laplacian"], laplacian))
        self.assertTrue(c.np.allclose(results["divergence"], divergence))

        # Check function calls
        self.assertTrue(c.np.allclose(results["hessian"], hessian))
        self.assertTrue(c.np.allclose(results["gradient"], gradient))

        # Check with previous test
        self.assertTrue(
            c.np.allclose(
                results["basis_function_values"]["gradient"],
                bf_results["gradient"],
            )
        )

    def check_assertions(self):
        mapper = self.solution_field_rando.geometry_mapper(self.askew_spline2D)
        self.assertRaises(mapper.divergence(self.query_points2D))

    def test_first_order_derivatives_analytical(self):
        mapper2D = self.solution_field_rando.geometry_mapper(self.rotating2D)
        mapper3D = self.solution_field_mono3D.geometry_mapper(self.scaling3D)
        bf_gradient, support = mapper2D.basis_gradient_and_support(
            self.query_points2D
        )
        (
            bf_reference0,
            supportb,
        ) = self.solution_field_rando.basis_derivative_and_support(
            self.query_points2D, [1, 0]
        )
        (
            bf_reference1,
            supportb,
        ) = self.solution_field_rando.basis_derivative_and_support(
            self.query_points2D, [0, 1]
        )
        bf_reference = c.np.dstack((bf_reference0, bf_reference1))

        self.assertTrue(c.np.allclose(support, supportb))
        # Rotate bf_reference
        bf_reference = c.np.einsum(
            "ij,qsj->qsi", self.rotation_matrix, bf_reference
        )
        self.assertTrue(c.np.allclose(bf_gradient, bf_reference))

        bf_gradient, support = mapper3D.basis_gradient_and_support(
            self.query_points3D
        )
        (
            bf_reference0,
            supportb,
        ) = self.solution_field_mono3D.basis_derivative_and_support(
            self.query_points3D, [1, 0, 0]
        )
        (
            bf_reference1,
            supportb,
        ) = self.solution_field_mono3D.basis_derivative_and_support(
            self.query_points3D, [0, 1, 0]
        )
        (
            bf_reference2,
            supportb,
        ) = self.solution_field_mono3D.basis_derivative_and_support(
            self.query_points3D, [0, 0, 1]
        )
        bf_reference = c.np.dstack(
            (bf_reference0, bf_reference1, bf_reference2)
        )
        bf_reference = c.np.einsum(
            "qsi,i->qsi", bf_reference, self.scaling_factors
        )
        self.assertTrue(c.np.allclose(bf_gradient, bf_reference))

    def test_second_order_analytical(self):
        pass

    def test_second_order_fd(self):
        pass


if __name__ == "__main__":
    c.unittest.main()
