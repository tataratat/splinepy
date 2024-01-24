try:
    from . import common as c
except BaseException:
    import common as c


class TestSplinepyEvaluation(c.SplineBasedTestCase):
    def setUp(self):
        self.bspline = self.bspline_2p2d()
        self.nurbs = self.nurbs_2p2d()
        self.rational = self.rational_bezier_2p2d()
        self.bezier = self.bezier_2p2d()

    def test_basis_and_support(self):
        """Test the correct calculation of the basis functions.
        (.basis_and_support())"""

        # reference solutions
        bspline_ref_basis_functions = c.np.array(
            [
                [
                    9.4128804e-01,
                    3.8615940e-02,
                    1.9602000e-04,
                    1.9015920e-02,
                    7.8012000e-04,
                    3.9600000e-06,
                    9.6040000e-05,
                    3.9400000e-06,
                    2.0000000e-08,
                ],
                [
                    2.4010000e-01,
                    9.8500000e-03,
                    5.0000000e-05,
                    4.8020000e-01,
                    1.9700000e-02,
                    1.0000000e-04,
                    2.4010000e-01,
                    9.8500000e-03,
                    5.0000000e-05,
                ],
                [
                    1.6200000e-02,
                    2.7540000e-01,
                    5.1840000e-01,
                    3.6000000e-03,
                    6.1200000e-02,
                    1.1520000e-01,
                    2.0000000e-04,
                    3.4000000e-03,
                    6.4000000e-03,
                ],
                [
                    7.2000000e-03,
                    5.0400000e-02,
                    3.2400000e-02,
                    3.3600000e-02,
                    2.3520000e-01,
                    1.5120000e-01,
                    3.9200000e-02,
                    2.7440000e-01,
                    1.7640000e-01,
                ],
                [
                    4.0000000e-06,
                    6.4000000e-05,
                    3.2000000e-05,
                    7.9200000e-04,
                    1.2672000e-02,
                    6.3360000e-03,
                    3.9204000e-02,
                    6.2726400e-01,
                    3.1363200e-01,
                ],
            ]
        )
        nurbs_ref_basis_functions = c.np.array(
            [
                [
                    9.75958864e-01,
                    1.39415582e-02,
                    9.95774782e-05,
                    9.85817035e-03,
                    1.40823820e-04,
                    1.00583311e-06,
                ],
                [
                    4.92908517e-01,
                    7.04119101e-03,
                    5.02916557e-05,
                    4.92908517e-01,
                    7.04119101e-03,
                    5.02916557e-05,
                ],
                [
                    9.50089457e-03,
                    1.20926646e-01,
                    7.69572460e-01,
                    1.05565495e-03,
                    1.34362940e-02,
                    8.55080511e-02,
                ],
                [
                    1.32410262e-02,
                    7.49025551e-02,
                    2.11856419e-01,
                    3.08957277e-02,
                    1.74772629e-01,
                    4.94331644e-01,
                ],
                [
                    4.18891419e-03,
                    3.94934617e-03,
                    1.86173964e-03,
                    4.14702505e-01,
                    3.90985271e-01,
                    1.84312224e-01,
                ],
            ]
        )

        # test basis functions
        self.assertTrue(
            c.np.allclose(
                self.bspline.basis_and_support(c.get_queries_2D())[0],
                bspline_ref_basis_functions,
            )
        )
        self.assertTrue(
            c.np.allclose(
                self.nurbs.basis_and_support(c.get_queries_2D())[0],
                nurbs_ref_basis_functions,
            )
        )

    def test_jacobian(self):
        """Test the correct evaluation of basis function derivatives"""
        # Test both for nurbs as for rational beziers for the integration of a
        # circle, just because we can
        # Get Legendre points
        interval = [0, 1]
        order = 7
        dim = 2
        positions, weights = c.np.polynomial.legendre.leggauss(order)
        positions = interval[0] + (interval[1] - interval[0]) / 2 * (
            positions + 1
        )
        weights = weights * (interval[1] - interval[0]) / 2
        positions = c.np.reshape(
            c.np.meshgrid(*[positions for _ in range(dim)]), (dim, -1)
        ).T
        weights = c.np.prod(
            c.np.reshape(
                c.np.meshgrid(*[weights for _ in range(dim)]), (dim, -1)
            ).T,
            axis=1,
        )
        inv_sqrt_2 = c.np.sqrt(0.5)
        circle_bezier = c.splinepy.RationalBezier(
            degrees=[2, 2],
            control_points=[
                [0, 0],
                [0.5, -0.5],
                [1, 0],
                [-0.5, 0.5],
                [0.5, 0.5],
                [1 + 0.5, 0.5],
                [0, 1],
                [0.5, 1 + 0.5],
                [1, 1],
            ],
            weights=[
                1,
                inv_sqrt_2,
                1,
                inv_sqrt_2,
                1,
                inv_sqrt_2,
                1,
                inv_sqrt_2,
                1,
            ],
        )
        circle_nurbs = c.splinepy.NURBS(
            degrees=circle_bezier.degrees,
            control_points=circle_bezier.control_points,
            knot_vectors=[[0, 0, 0, 1, 1, 1], [0, 0, 0, 1, 1, 1]],
            weights=circle_bezier.weights,
        )

        # Test the determinant
        area = c.np.sum(
            c.np.linalg.det(circle_bezier.jacobian(positions)) * weights
        )
        self.assertAlmostEqual(area, c.np.pi * 0.5)
        area = c.np.sum(
            c.np.linalg.det(circle_nurbs.jacobian(positions)) * weights
        )
        self.assertAlmostEqual(area, c.np.pi * 0.5)

        # Assuming that the derivative is correct
        test_spline = c.splinepy.Bezier(
            degrees=[3, 3], control_points=c.np.random.rand(16, 3)
        )
        query = c.np.random.rand(1, 2)
        jacs = test_spline.jacobian(query)[0]
        expected_jacs = c.np.vstack(
            (
                test_spline.derivative(query, [1, 0]),
                test_spline.derivative(query, [0, 1]),
            )
        ).T
        self.assertTrue(c.np.allclose(jacs, expected_jacs))

        # Test for scalar valued splines
        test_spline_scalar = c.splinepy.Bezier(
            degrees=[3, 3], control_points=c.np.random.rand(16, 1)
        )
        query = c.np.random.rand(1, 2)
        jacs = test_spline_scalar.jacobian(query)[0]
        expected_jacs = c.np.vstack(
            (
                test_spline_scalar.derivative(query, [1, 0]),
                test_spline_scalar.derivative(query, [0, 1]),
            )
        ).T
        self.assertTrue(c.np.allclose(jacs, expected_jacs))

    def test_basis_function_derivatives(self):
        """Test the correct evaluation of basis function derivatives"""
        # Cross-testing different libraries
        # use random query points
        q2D = c.np.random.rand(10, 2)

        # Rational Bezier and NURBS are equivalent but use different backends
        rational_c = self.rational_bezier_2p2d().copy()
        nurbs_c = self.nurbs.copy()
        # increase orders for derivatives
        for _ in range(2):
            rational_c.elevate_degrees([0, 1])
            nurbs_c.elevate_degrees([0, 1])

        # Test different derivatives - all have global supports
        self.assertTrue(
            c.np.allclose(
                rational_c.basis_derivative_and_support(q2D, [1, 0])[0],
                nurbs_c.basis_derivative_and_support(q2D, [1, 0])[0],
            )
        )
        self.assertTrue(
            c.np.allclose(
                rational_c.basis_derivative_and_support(q2D, [1, 0])[1],
                nurbs_c.basis_derivative_and_support(q2D, [1, 0])[1],
            )
        )
        self.assertTrue(
            c.np.allclose(
                rational_c.basis_derivative_and_support(q2D, [3, 2])[0],
                nurbs_c.basis_derivative_and_support(q2D, [3, 2])[0],
            )
        )
        self.assertTrue(
            c.np.allclose(
                rational_c.basis_derivative_and_support(q2D, [1, 3])[0],
                nurbs_c.basis_derivative_and_support(q2D, [1, 3])[0],
            )
        )

        # For polynomial splines
        bezier_c = self.bezier.copy()
        bspline_c = c.splinepy.BSpline(
            **self.bezier.todict(),
            knot_vectors=[
                [0] * (self.bezier.degrees[i] + 1)
                + [1] * (self.bezier.degrees[i] + 1)
                for i in range(self.bezier.para_dim)
            ],
        )

        # increase orders for derivatives
        for _i in range(2):
            bezier_c.elevate_degrees([0, 1])
            bspline_c.elevate_degrees([0, 1])

        # Test different derivatives - all have global supports
        self.assertTrue(
            c.np.allclose(
                bezier_c.basis_derivative_and_support(q2D, [0, 1])[0],
                bspline_c.basis_derivative_and_support(q2D, [0, 1])[0],
            )
        )
        self.assertTrue(
            c.np.allclose(
                bezier_c.basis_derivative_and_support(q2D, [0, 1])[1],
                bspline_c.basis_derivative_and_support(q2D, [0, 1])[1],
            )
        )
        self.assertTrue(
            c.np.allclose(
                bezier_c.basis_derivative_and_support(q2D, [2, 3])[0],
                bspline_c.basis_derivative_and_support(q2D, [2, 3])[0],
            )
        )
        self.assertTrue(
            c.np.allclose(
                bezier_c.basis_derivative_and_support(q2D, [3, 1])[0],
                bspline_c.basis_derivative_and_support(q2D, [3, 1])[0],
            )
        )

    def test_partition_of_unity(self):
        """Test the partition of unity of the calculated basis functions."""

        def basis_functions_sum(basis_functions):
            return basis_functions.sum(axis=1)

        # use random query points
        q2D = c.np.random.rand(10, 2)

        u_bspline = basis_functions_sum(self.bspline.basis_and_support(q2D)[0])
        self.assertTrue(
            c.np.allclose(u_bspline, c.np.ones(c.np.shape(u_bspline)))
        )
        u_nurbs = basis_functions_sum(self.nurbs.basis_and_support(q2D)[0])
        self.assertTrue(c.np.allclose(u_nurbs, c.np.ones(c.np.shape(u_nurbs))))

    def test_evaluate(self):
        """Test the correct spline evaluation in the physical space.
        (.evaluate())"""

        # reference solutions
        bspline_ref_evaluate = c.np.array(
            [
                [-0.019796, 0.04049403],
                [-0.9996, 0.069675],
                [2.256, 1.9691],
                [1.528, 4.0766],
                [-1.0264, 2.873488],
            ]
        )
        nurbs_ref_evaluate = c.np.array(
            [
                [-1.00989841, 0.01432479],
                [-1.49984913, 0.02127445],
                [-0.15941144, 1.0883878],
                [-0.49948029, 1.62496752],
                [-1.61951381, 1.15640608],
            ]
        )
        bezier_ref_evaluate = c.np.array(
            [
                [-1.009899, 0.020099],
                [-1.49985, 0.02985],
                [-0.209, 1.089],
                [-0.612, 1.632],
                [-1.6716, 1.2736],
            ]
        )
        rational_ref_evaluate = c.np.array(
            [
                [-1.00989841, 0.01432479],
                [-1.49984913, 0.02127445],
                [-0.15941144, 1.0883878],
                [-0.49948029, 1.62496752],
                [-1.61951381, 1.15640608],
            ]
        )

        # test evaluation
        self.assertTrue(
            c.np.allclose(
                self.bspline.evaluate(c.get_queries_2D()),
                bspline_ref_evaluate,
            )
        )
        self.assertTrue(
            c.np.allclose(
                self.nurbs.evaluate(c.get_queries_2D()),
                nurbs_ref_evaluate,
            )
        )
        self.assertTrue(
            c.np.allclose(
                self.bezier.evaluate(c.get_queries_2D()),
                bezier_ref_evaluate,
            )
        )
        self.assertTrue(
            c.np.allclose(
                self.rational.evaluate(c.get_queries_2D()),
                rational_ref_evaluate,
            )
        )

    def test_derivative(self):
        """Test the correct calculation of the first derivative.
        (.derivative())"""

        # reference solutions
        bspline_ref_derivative = c.np.array(
            [
                [0.08, 7.8812],
                [0.08, 4.02],
                [0.8, 1.16],
                [1.6, 1.84],
                [3.2, 3.232],
            ]
        )
        nurbs_ref_derivative = c.np.array(
            [
                [0.02017474, 1.42231975],
                [0.02017474, 1.42231975],
                [1.47716145, 0.21635343],
                [1.49159582, 0.45848468],
                [0.95624956, 1.33920031],
            ]
        )

        # order
        o1 = [1, 1]

        # test derivative evaluation
        self.assertTrue(
            c.np.allclose(
                self.bspline.derivative(c.get_queries_2D(), o1),
                bspline_ref_derivative,
            )
        )
        self.assertTrue(
            c.np.allclose(
                self.nurbs.derivative(c.get_queries_2D(), o1),
                nurbs_ref_derivative,
            )
        )

        # Test higher order derivatives against each other
        max_deg = 5
        dim = 3
        n_test_points = 10
        derivatives = c.np.random.randint(0, max_deg, size=(dim))
        queries = c.np.random.random((n_test_points, dim))
        randomized = {
            "degrees": [5] * dim,
            "control_points": c.np.random.random(((max_deg + 1) ** dim, dim)),
            "weights": c.np.random.random((max_deg + 1) ** dim),
            "knot_vectors": [[0] * (max_deg + 1) + [1] * (max_deg + 1)] * dim,
        }
        random_bezier = c.splinepy.Bezier(
            degrees=randomized["degrees"],
            control_points=randomized["control_points"],
        )
        random_rational = c.splinepy.RationalBezier(
            degrees=randomized["degrees"],
            control_points=randomized["control_points"],
            weights=randomized["weights"],
        )
        random_bspline = c.splinepy.BSpline(
            degrees=randomized["degrees"],
            control_points=randomized["control_points"],
            knot_vectors=randomized["knot_vectors"],
        )
        random_nurbs = c.splinepy.NURBS(
            degrees=randomized["degrees"],
            control_points=randomized["control_points"],
            knot_vectors=randomized["knot_vectors"],
            weights=randomized["weights"],
        )
        self.assertTrue(
            c.np.allclose(
                random_bezier.derivative(queries, derivatives),
                random_bspline.derivative(queries, derivatives),
            )
        )
        self.assertTrue(
            c.np.allclose(
                random_rational.derivative(queries, derivatives),
                random_nurbs.derivative(queries, derivatives),
            )
        )

    def test_assertions_evaluation(self):
        """Test if the right assertions are thrown"""
        # Test for evaluation
        for spline in [self.bezier, self.rational, self.nurbs, self.bspline]:
            # Check minimum with kwargs
            self.assertRaisesRegex(
                ValueError,
                r"Query request out of bounds in parametric dimension 0. "
                r"Detected query \[-0.1  0. \] at positions 0, which is out of"
                r" bounds with minimum values \[1. 1.\].",
                spline.evaluate,
                queries=[[-0.1, 0.0], [0.4, 0.6]],
            )
            # Check maximum with args
            self.assertRaisesRegex(
                ValueError,
                r"Query request out of bounds in parametric dimension 0. "
                r"Detected query \[1.1 0. \] at positions 0, which is out of "
                r"bounds with maximum values \[1. 1.\].",
                spline.evaluate,
                [[1.1, 0.0], [0.4, 0.6]],
            )
            # Check dimensions
            self.assertRaisesRegex(
                ValueError,
                r"Dimension mismatch between parametric dimension of spline "
                r"\(2\), and query-request \(3\)",
                spline.evaluate,
                [[1.0, 1.0, 1.0]],
            )

    def test_assertions_jacobian(self):
        """Test if the right assertions are thrown"""
        # Test for evaluation
        for spline in [self.bezier, self.rational, self.nurbs, self.bspline]:
            # Check minimum with kwargs
            self.assertRaisesRegex(
                ValueError,
                r"Query request out of bounds in parametric dimension 0. "
                r"Detected query \[-0.1  0. \] at positions 0, which is out of"
                r" bounds with minimum values \[1. 1.\].",
                spline.jacobian,
                queries=[[-0.1, 0.0], [0.4, 0.6]],
            )
            # Check maximum with args
            self.assertRaisesRegex(
                ValueError,
                r"Query request out of bounds in parametric dimension 0. "
                r"Detected query \[1.1 0. \] at positions 0, which is out of "
                r"bounds with maximum values \[1. 1.\].",
                spline.jacobian,
                [[1.1, 0.0], [0.4, 0.6]],
            )
            # Check dimensions
            self.assertRaisesRegex(
                ValueError,
                r"Dimension mismatch between parametric dimension of spline "
                r"\(2\), and query-request \(3\)",
                spline.jacobian,
                [[1.0, 1.0, 1.0]],
            )

    def test_basis_function_matrix(self):
        """Test the correct evaluation of basis function matrices"""
        make_matrix = c.splinepy.utils.data.make_matrix
        # Use
        q2D = c.np.random.rand(10, 2)

        # Test for both bezier and BSpline family  and rationals
        bezier = self.bezier_2p2d().copy()
        rational = self.rational_bezier_2p2d().copy()
        bspline = self.bspline_2p2d().copy()

        # Compute derivative matrix
        for spline in (bezier, rational, bspline):
            # Trivial evaluation
            matrix = make_matrix(
                *spline.basis_and_support(q2D), spline.cps.shape[0]
            )
            self.assertTrue(
                c.np.allclose(
                    matrix @ spline.cps,
                    spline.evaluate(q2D),
                )
            )

            # Test first order derivative
            matrix = make_matrix(
                *spline.basis_derivative_and_support(q2D, orders=[1, 0]),
                spline.cps.shape[0],
            )
            self.assertTrue(
                c.np.allclose(
                    matrix @ spline.cps,
                    spline.derivative(q2D, orders=[1, 0]),
                )
            )

            # Test seoncd order derivative as numpy (enforce)
            matrix = make_matrix(
                *spline.basis_derivative_and_support(q2D, orders=[1, 2]),
                spline.cps.shape[0],
                as_array=True,
            )
            self.assertTrue(
                c.np.allclose(
                    matrix @ spline.cps,
                    spline.derivative(q2D, orders=[1, 2]),
                )
            )

    def test_multiple_derivative_queries(self):
        """Test cartesian product queries of parametric coordinates and orders"""
        p_coord = {2: c.np.random.rand(10, 2), 3: c.np.random.rand(10, 3)}
        a = c.np.array  # shortcut
        orders = {
            2: c.splinepy.utils.data.cartesian_product(
                [a([0, 1, 2]), a([0, 1, 2])]
            ),
            3: c.splinepy.utils.data.cartesian_product(
                [a([0, 1, 2]), a([0, 1, 2]), a([0, 1, 2])]
            ),
        }
        jac_transpose_orders = {2: c.np.eye(2), 3: c.np.eye(3)}

        for spline in [*self.all_2p2d_splines(), *self.all_3p3d_splines()]:
            pd = spline.para_dim
            qs = p_coord[pd]
            jto = jac_transpose_orders[pd]
            os = orders[pd]

            # multiple queries
            multi0 = spline.derivative(qs, os)
            multi1 = spline.basis_derivative(qs, os)
            multi2, multi2_support = spline.basis_derivative_and_support(
                qs, os
            )

            # as well as jacobian transposed queries
            # however, we transpose results to enable direct comparison
            multi_jac = spline.derivative(qs, jto).transpose(0, 2, 1)
            single_jac = spline.jacobian(qs)

            # multiple single evaluations
            single0 = []
            single1 = []
            single2 = []
            single2_support = []
            # loop in order of the visit
            for q in qs:
                for o in os:
                    single0.append(spline.derivative([q], o))
                    single1.append(spline.basis_derivative([q], o))
                    s2, s2_s = spline.basis_derivative_and_support([q], o)
                    single2.append(s2)
                    single2_support.append(s2_s)

            # compare and see if they match
            # multi queries will return
            # (n_queries, n_orders, dim or n_support) shape
            # for comparison, reshape. we could also ravel both, alternatively.
            dim = multi0.shape[-1]
            assert c.np.allclose(multi0.reshape(-1, dim), c.np.vstack(single0))
            sup = multi1.shape[-1]
            assert c.np.allclose(multi1.reshape(-1, sup), c.np.vstack(single1))
            assert c.np.allclose(multi2.reshape(-1, sup), c.np.vstack(single2))
            assert (
                multi2_support
                == a(single2_support).reshape(-1, multi2_support.shape[-1])[
                    :: len(os)
                ]
            ).all()
            assert c.np.allclose(multi_jac, single_jac)


if __name__ == "__main__":
    c.unittest.main()
