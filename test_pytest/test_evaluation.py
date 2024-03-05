import numpy as np
import pytest

import splinepy

# frequently used fixtures
all_2p2d_splines = (
    "rational_bezier_2p2d",
    "bezier_2p2d",
    "bspline_2p2d",
    "nurbs_2p2d",
)


def setUp(request):
    bspline = request.getfixturevalue("bspline_2p2d")
    nurbs = request.getfixturevalue("nurbs_2p2d")
    rational = request.getfixturevalue("rational_bezier_2p2d")
    bezier = request.getfixturevalue("bezier_2p2d")
    return bspline, nurbs, rational, bezier


def test_basis_and_support(request):
    """Test the correct calculation of the basis functions.
    (.basis_and_support())"""

    bspline, nurbs, _, _ = setUp(request)

    # reference solutions
    bspline_ref_basis_functions = np.array(
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
    nurbs_ref_basis_functions = np.array(
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

    queries = request.getfixturevalue("get_queries_2D")

    # test basis functions
    assert np.allclose(
        bspline.basis_and_support(queries)[0],
        bspline_ref_basis_functions,
    )

    assert np.allclose(
        nurbs.basis_and_support(queries)[0],
        nurbs_ref_basis_functions,
    )


def test_jacobian(np_rng):
    """Test the correct evaluation of basis function derivatives"""
    # Test both for nurbs as for rational beziers for the integration of a
    # circle, just because we can
    # Get Legendre points
    interval = [0, 1]
    order = 7
    dim = 2
    positions, weights = np.polynomial.legendre.leggauss(order)
    positions = interval[0] + (interval[1] - interval[0]) / 2 * (positions + 1)
    weights = weights * (interval[1] - interval[0]) / 2
    positions = np.reshape(
        np.meshgrid(*[positions for _ in range(dim)]), (dim, -1)
    ).T
    weights = np.prod(
        np.reshape(np.meshgrid(*[weights for _ in range(dim)]), (dim, -1)).T,
        axis=1,
    )
    inv_sqrt_2 = np.sqrt(0.5)
    circle_bezier = splinepy.RationalBezier(
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
    circle_nurbs = splinepy.NURBS(
        degrees=circle_bezier.degrees,
        control_points=circle_bezier.control_points,
        knot_vectors=[[0, 0, 0, 1, 1, 1], [0, 0, 0, 1, 1, 1]],
        weights=circle_bezier.weights,
    )

    # Test the determinant
    area = np.sum(np.linalg.det(circle_bezier.jacobian(positions)) * weights)
    assert area == pytest.approx(np.pi * 0.5)

    area = np.sum(np.linalg.det(circle_nurbs.jacobian(positions)) * weights)
    assert area == pytest.approx(np.pi * 0.5)

    # Assuming that the derivative is correct
    test_spline = splinepy.Bezier(
        degrees=[3, 3], control_points=np_rng.random((16, 3))
    )
    query = np_rng.random((1, 2))
    jacs = test_spline.jacobian(query)[0]
    expected_jacs = np.vstack(
        (
            test_spline.derivative(query, [1, 0]),
            test_spline.derivative(query, [0, 1]),
        )
    ).T

    assert np.allclose(jacs, expected_jacs)

    # Test for scalar valued splines
    test_spline_scalar = splinepy.Bezier(
        degrees=[3, 3], control_points=np_rng.random((16, 1))
    )
    query = np_rng.random((1, 2))
    jacs = test_spline_scalar.jacobian(query)[0]
    expected_jacs = np.vstack(
        (
            test_spline_scalar.derivative(query, [1, 0]),
            test_spline_scalar.derivative(query, [0, 1]),
        )
    ).T

    assert np.allclose(jacs, expected_jacs)


def test_basis_function_derivatives(np_rng, request):
    """Test the correct evaluation of basis function derivatives"""
    # Cross-testing different libraries
    # use random query points
    q2D = np_rng.random((10, 2))

    # Rational Bezier and NURBS are equivalent but use different backends
    _, nurbs, rational, _ = setUp(request)

    # increase orders for derivatives
    for _ in range(2):
        rational.elevate_degrees([0, 1])
        nurbs.elevate_degrees([0, 1])

    # Test different derivatives - all have global supports
    assert np.allclose(
        rational.basis_derivative_and_support(q2D, [1, 0])[0],
        nurbs.basis_derivative_and_support(q2D, [1, 0])[0],
    )

    assert np.allclose(
        rational.basis_derivative_and_support(q2D, [1, 0])[1],
        nurbs.basis_derivative_and_support(q2D, [1, 0])[1],
    )

    assert np.allclose(
        rational.basis_derivative_and_support(q2D, [3, 2])[0],
        nurbs.basis_derivative_and_support(q2D, [3, 2])[0],
    )

    assert np.allclose(
        rational.basis_derivative_and_support(q2D, [1, 3])[0],
        nurbs.basis_derivative_and_support(q2D, [1, 3])[0],
    )

    # For polynomial splines
    bezier = request.getfixturevalue("bezier_2p2d")
    bspline_c = splinepy.BSpline(
        **bezier.todict(),
        knot_vectors=[
            [0] * (bezier.degrees[i] + 1) + [1] * (bezier.degrees[i] + 1)
            for i in range(bezier.para_dim)
        ],
    )

    # increase orders for derivatives
    for _i in range(2):
        bezier.elevate_degrees([0, 1])
        bspline_c.elevate_degrees([0, 1])

    # Test different derivatives - all have global supports
    assert np.allclose(
        bezier.basis_derivative_and_support(q2D, [0, 1])[0],
        bspline_c.basis_derivative_and_support(q2D, [0, 1])[0],
    )

    assert np.allclose(
        bezier.basis_derivative_and_support(q2D, [0, 1])[1],
        bspline_c.basis_derivative_and_support(q2D, [0, 1])[1],
    )

    assert np.allclose(
        bezier.basis_derivative_and_support(q2D, [2, 3])[0],
        bspline_c.basis_derivative_and_support(q2D, [2, 3])[0],
    )

    assert np.allclose(
        bezier.basis_derivative_and_support(q2D, [3, 1])[0],
        bspline_c.basis_derivative_and_support(q2D, [3, 1])[0],
    )


def test_partition_of_unity(request, np_rng):
    """Test the partition of unity of the calculated basis functions."""

    def basis_functions_sum(basis_functions):
        return basis_functions.sum(axis=1)

    bspline, nurbs, _, _ = setUp(request)

    # use random query points
    q2D = np_rng.random((10, 2))

    u_bspline = basis_functions_sum(bspline.basis_and_support(q2D)[0])
    assert np.allclose(u_bspline, np.ones(np.shape(u_bspline)))

    u_nurbs = basis_functions_sum(nurbs.basis_and_support(q2D)[0])
    assert np.allclose(u_nurbs, np.ones(np.shape(u_nurbs)))


def test_evaluate(request):
    """Test the correct spline evaluation in the physical space.
    (.evaluate())"""

    bspline, nurbs, rational, bezier = setUp(request)

    # reference solutions
    bspline_ref_evaluate = np.array(
        [
            [-0.019796, 0.04049403],
            [-0.9996, 0.069675],
            [2.256, 1.9691],
            [1.528, 4.0766],
            [-1.0264, 2.873488],
        ]
    )
    nurbs_ref_evaluate = np.array(
        [
            [-1.00989841, 0.01432479],
            [-1.49984913, 0.02127445],
            [-0.15941144, 1.0883878],
            [-0.49948029, 1.62496752],
            [-1.61951381, 1.15640608],
        ]
    )
    bezier_ref_evaluate = np.array(
        [
            [-1.009899, 0.020099],
            [-1.49985, 0.02985],
            [-0.209, 1.089],
            [-0.612, 1.632],
            [-1.6716, 1.2736],
        ]
    )
    rational_ref_evaluate = np.array(
        [
            [-1.00989841, 0.01432479],
            [-1.49984913, 0.02127445],
            [-0.15941144, 1.0883878],
            [-0.49948029, 1.62496752],
            [-1.61951381, 1.15640608],
        ]
    )

    queries = request.getfixturevalue("get_queries_2D")
    # test evaluation
    assert np.allclose(
        bspline.evaluate(queries),
        bspline_ref_evaluate,
    )

    assert np.allclose(
        nurbs.evaluate(queries),
        nurbs_ref_evaluate,
    )

    assert np.allclose(
        bezier.evaluate(queries),
        bezier_ref_evaluate,
    )

    assert np.allclose(
        rational.evaluate(queries),
        rational_ref_evaluate,
    )
