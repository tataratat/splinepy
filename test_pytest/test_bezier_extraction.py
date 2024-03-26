import numpy as np
import pytest

import splinepy

# frequently used fixtures
all_splinetypes = ("bspline_3p3d", "nurbs_3p3d", "bspline_2p2d", "nurbs_2p2d")


def test_extraction(np_rng):
    """
    test the extraction of beziers
    """
    # Define some splines
    b = splinepy.BSpline(
        knot_vectors=[[0, 0, 0, 1, 2, 2, 3, 4, 4, 4]],
        degrees=[2],
        control_points=np_rng.random((7, 2)),
    )
    n = splinepy.NURBS(
        knot_vectors=[[0, 0, 0, 0, 1, 1, 2, 3, 4, 4, 4, 4]],
        degrees=[3],
        control_points=np_rng.random((8, 2)),
        weights=np_rng.random((8, 1)),
    )

    # Extract Beziers
    b_beziers = b.extract_bezier_patches()
    n_beziers = n.extract_bezier_patches()

    # Loop over knot_spans and test at random points
    for offset in range(4):
        queries = np_rng.random((20, 1))
        assert np.allclose(
            b.evaluate(queries + offset),
            b_beziers[offset].evaluate(queries),
        )

        assert np.allclose(
            n.evaluate(queries + offset),
            n_beziers[offset].evaluate(queries),
        )


@pytest.mark.parametrize("splinetype", all_splinetypes)
def test_extraction_matrices(splinetype, np_rng, request):
    spline = request.getfixturevalue(splinetype)

    if "3d" in splinetype:
        spline.elevate_degrees([0, 1, 2])
        spline.insert_knots(0, np_rng.random(3))
        spline.insert_knots(1, np_rng.random(3))
        spline.insert_knots(2, np_rng.random(3))
    elif "2d" in splinetype:
        spline.elevate_degrees([0, 1])
        spline.insert_knots(0, np_rng.random(3))
        spline.insert_knots(1, np_rng.random(3))
    else:
        pytest.fail()

    n_matrices = spline.knot_insertion_matrix(beziers=True)
    beziers_n = spline.extract_bezier_patches()
    for m, b in zip(n_matrices, beziers_n):
        # Test matrices m against spline ctps
        if "nurbs" in splinetype:
            assert np.allclose(b.weights, m @ spline.weights)
            assert np.allclose(
                b.control_points,
                (m @ (spline.control_points * spline.weights)) / b.weights,
            )
        elif "bspline" in splinetype:
            assert np.allclose(b.control_points, m @ spline.control_points)
        else:
            pytest.fail()
