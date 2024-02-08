import numpy as np
import pytest

import splinepy

# frequently used fixtures
all_splinetypes = ("bspline_3p3d", "nurbs_3p3d", "bspline_2p2d", "nurbs_2p2d")


def test_extraction():
    """
    test the extraction of beziers
    """
    # Define some splines
    b = splinepy.BSpline(
        knot_vectors=[[0, 0, 0, 1, 2, 2, 3, 4, 4, 4]],
        degrees=[2],
        control_points=np.random.rand(7, 2),
    )
    n = splinepy.NURBS(
        knot_vectors=[[0, 0, 0, 0, 1, 1, 2, 3, 4, 4, 4, 4]],
        degrees=[3],
        control_points=np.random.rand(8, 2),
        weights=np.random.rand(8, 1),
    )

    # Extract Beziers
    b_beziers = b.extract_bezier_patches()
    n_beziers = n.extract_bezier_patches()

    # Loop over knot_spans and test at random points
    for offset in range(4):
        queries = np.random.rand(20, 1)
        assert np.allclose(
            b.evaluate(queries + offset),
            b_beziers[offset].evaluate(queries),
        )

        assert np.allclose(
            n.evaluate(queries + offset),
            n_beziers[offset].evaluate(queries),
        )


@pytest.mark.parametrize("splinetype", all_splinetypes)
def test_extraction_matrices(splinetype, request):
    spline = request.getfixturevalue(splinetype)

    if splinetype == "bspline_3p3d":
        """Create matrices to extract splines"""
        # Init b-splines
        bspline = spline
        bspline.elevate_degrees([0, 1, 2])
        bspline.insert_knots(0, np.random.rand(3))
        bspline.insert_knots(1, np.random.rand(3))
        bspline.insert_knots(2, np.random.rand(3))

        # Extract splines
        # BSpline
        beziers_b = bspline.extract_bezier_patches()
        b_matrices = bspline.knot_insertion_matrix(beziers=True)
        for m, b in zip(b_matrices, beziers_b):
            # Test matrices m against spline ctps
            assert np.allclose(b.control_points, m @ bspline.control_points)

    elif splinetype == "nurbs_3p3d":
        """Create matrices to extract splines"""
        # Init nurbs
        nurbs = spline
        nurbs.elevate_degrees([0, 1, 2])
        nurbs.insert_knots(0, np.random.rand(3))
        nurbs.insert_knots(1, np.random.rand(3))
        nurbs.insert_knots(2, np.random.rand(3))
        # NURBS
        n_matrices = nurbs.knot_insertion_matrix(beziers=True)
        beziers_n = nurbs.extract_bezier_patches()
        for m, b in zip(n_matrices, beziers_n):
            # Test matrices m against spline ctps
            assert np.allclose(b.weights, m @ nurbs.weights)
            assert np.allclose(
                b.control_points,
                (m @ (nurbs.control_points * nurbs.weights)) / b.weights,
            )

    elif splinetype == "bspline_2p2d":
        """Create matrices to extract splines"""
        # Init nurbs
        bspline = spline
        bspline.elevate_degrees([0, 1])
        bspline.insert_knots(0, np.random.rand(3))
        bspline.insert_knots(1, np.random.rand(3))
        # NURBS
        n_matrices = bspline.knot_insertion_matrix(beziers=True)
        beziers_n = bspline.extract_bezier_patches()
        for m, b in zip(n_matrices, beziers_n):
            # Test matrices m against spline ctps
            assert np.allclose(b.control_points, m @ bspline.control_points)

    elif splinetype == "nurbs_2p2d":
        """Create matrices to extract splines"""
        # Init nurbs
        nurbs = spline
        nurbs.elevate_degrees([0, 1])
        nurbs.insert_knots(0, np.random.rand(3))
        nurbs.insert_knots(1, np.random.rand(3))
        # NURBS
        n_matrices = nurbs.knot_insertion_matrix(beziers=True)
        beziers_n = nurbs.extract_bezier_patches()
        for m, b in zip(n_matrices, beziers_n):
            # Test matrices m against spline ctps
            assert np.allclose(b.weights, m @ nurbs.weights)
            assert np.allclose(
                b.control_points,
                (m @ (nurbs.control_points * nurbs.weights)) / b.weights,
            )

    else:
        pytest.fail()
