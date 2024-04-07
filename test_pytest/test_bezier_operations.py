"""
Test composition sensitivity, i.e., the derivative concerning the
deformation function's control points.
"""

import numpy as np

import splinepy


def test_composition():
    # Here we only try 2D compositions (surface-line is tested in bezman)
    # Init Splines to be tested
    surface_polynomial = splinepy.Bezier(
        degrees=[1, 2],
        control_points=[[1, 1], [2, 1], [1, 2], [2, 2], [1, 3], [2, 3]],
    )
    surface_rational = splinepy.RationalBezier(
        degrees=[2, 1],
        control_points=[[1, 1], [2, 1], [3, 2], [1, 2], [2, 2], [3, 3]],
        weights=[1.0, 0.8, 1.0, 1.0, 0.8, 1.0],
    )
    inner_polynomial = splinepy.Bezier(
        degrees=[1, 1],
        control_points=[[0.2, 0.0], [1.0, 0.2], [0.0, 0.8], [0.8, 1.0]],
    )
    inner_rational = splinepy.RationalBezier(
        degrees=[2, 1],
        control_points=[
            [0.2, 0.0],
            [0.5, 0.5],
            [1.0, 0.2],
            [0.0, 0.8],
            [0.5, 0.5],
            [0.8, 1.0],
        ],
        weights=[1, 0.5, 1, 1, 0.5, 1],
    )

    dx = 1e-4
    tolerance = 1e-4

    def test_splines(outer, inner):
        # Determine base line
        composed, composed_sensitivities = outer.compose(
            inner, compute_sensitivities=True
        )

        # Check against FD approximation
        for i in range(outer.cps.shape[0]):
            for j in range(outer.cps.shape[1]):
                outer.cps[i, j] += dx

                # FD values
                composed_dx = outer.compose(inner)

                # FD
                derivative_fd = (composed_dx.cps - composed.cps) / dx
                expected = np.zeros(composed.cps.shape)
                expected[:, j] = composed_sensitivities[i].cps.flatten()

                assert np.allclose(derivative_fd, expected, atol=tolerance)

                outer.cps[i, j] -= dx

    test_splines(surface_polynomial, inner_polynomial)
    test_splines(surface_polynomial, inner_rational)
    test_splines(surface_rational, inner_polynomial)
    test_splines(surface_rational, inner_rational)


def test_composition_sensitivities_on_bsplines(bspline_2p2d):
    """Combine Composition sensitivities with BSpline extraction"""

    # Initialize outer function
    bspline = bspline_2p2d

    # Initialize inner function
    inner_function = splinepy.Bezier(
        degrees=[1, 1],
        control_points=[
            [0.5, 0.0],
            [1.0, 0.5],
            [0.0, 0.5],
            [0.5, 1.0],
        ],
    )

    # Initialize values for finite differences
    bspline_dx = bspline.copy()
    bspline.insert_knots(1, 0.5)
    bspline_dx.insert_knots(1, 0.5)
    dx = 1e-5

    # Perform computations
    extract_beziers = bspline.extract_bezier_patches()
    extraction_matrices = bspline.knot_insertion_matrix(beziers=True)

    composed_der_ctps = []
    beziers = []
    for bez, mat in zip(extract_beziers, extraction_matrices):
        # Composition
        composed, derivatives = bez.compose(
            inner_function, compute_sensitivities=True
        )
        # Form matrix
        deriv_points = np.hstack([der.cps for der in derivatives])

        # Compute derivatives
        composed_der_ctps.append(deriv_points @ mat)
        beziers.append(composed)

    for cps in range(bspline.cps.shape[0]):
        for dim in range(bspline.cps.shape[1]):
            # Modify dx spline
            bspline_dx.cps[cps, dim] += dx

            # Extract Beziers
            extract_beziers_dx = bspline_dx.extract_bezier_patches()
            for bez, bez_dx, comps in zip(
                beziers, extract_beziers_dx, composed_der_ctps
            ):
                # Compose finite differences spline
                composed_dx = bez_dx.compose(inner_function)
                fd_ctps = (composed_dx.cps - bez.cps) / dx

                assert np.allclose(fd_ctps[:, dim], comps[:, cps])

            # Reset dx spline
            bspline_dx.cps[cps, dim] -= dx


def test_sum(np_rng, bezier_2p2d):
    # Create two splines
    bezier1 = bezier_2p2d
    bezier2 = splinepy.Bezier(
        degrees=[1, 1],
        control_points=[
            [0.5, 0.0],
            [1.0, 0.5],
            [0.0, 0.5],
            [0.5, 1.0],
        ],
    )
    bezier_sum = bezier1 + bezier2

    # Create queries
    n_test_points = 10
    para_dim = 2
    queries = np_rng.random((n_test_points, para_dim))

    assert np.allclose(
        bezier2.evaluate(queries) + bezier1.evaluate(queries),
        bezier_sum.evaluate(queries),
    )


def test_multiply(np_rng, bezier_2p2d):
    # Create two splines
    bezier1 = bezier_2p2d
    bezier2 = splinepy.Bezier(
        degrees=[1, 1],
        control_points=[
            [0.0],
            [0.5],
            [0.5],
            [1.0],
        ],
    )
    bezier_prod = bezier1 * bezier2

    # Create queries
    n_test_points = 10
    para_dim = 2
    queries = np_rng.random((n_test_points, para_dim))

    assert np.allclose(
        bezier2.evaluate(queries) * bezier1.evaluate(queries),
        bezier_prod.evaluate(queries),
    )


def test_close_form_derivative(np_rng):
    # Construct a random spline
    degrees = [5, 5, 5]
    orders = [3, 2, 4]
    bezier = splinepy.Bezier(
        degrees=degrees,
        control_points=np_rng.random((np.prod(np.array(degrees) + 1), 3)),
    )
    close_form = bezier.derivative_spline(orders)

    # Create queries
    n_test_points = 100
    para_dim = bezier.para_dim
    queries = np_rng.random((n_test_points, para_dim))

    assert np.allclose(
        bezier.derivative(queries, orders),
        close_form.evaluate(queries),
    )

    # Repeat for rationals (low orders required, because orders square with
    # every derivative)
    degrees = [2, 1, 2]
    orders = [1, 0, 1]
    n_cps = np.prod(np.array(degrees) + 1)
    rationals = splinepy.RationalBezier(
        degrees=degrees,
        control_points=np_rng.random((n_cps, 3)),
        weights=np.abs(np_rng.random(n_cps)),
    )
    close_form = rationals.derivative_spline(orders)

    # Create queries
    n_test_points = 100
    para_dim = rationals.para_dim
    queries = np_rng.random((n_test_points, para_dim))

    assert np.allclose(
        rationals.derivative(queries, orders),
        close_form.evaluate(queries),
    )
