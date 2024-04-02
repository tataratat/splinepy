import numpy as np
import pytest

# fixtures for test_elevate_degree
all_2p2d_splines = (
    "rational_bezier_2p2d",
    "bezier_2p2d",
    "bspline_2p2d",
    "nurbs_2p2d",
)

# fixtures for test_reduce_degree
nurbs_and_bspline = (
    "nurbs_2p2d",
    "bspline_2p2d",
)


@pytest.mark.parametrize("splinetype", all_2p2d_splines)
def test_elevate_degree(splinetype, request, np_rng, are_items_close):
    """Test the order elevation function (.elevate_degrees())."""

    # reference solution
    bspline_ref_kv = [
        [0, 0, 0, 0, 0.5, 0.5, 1, 1, 1, 1],
        [0, 0, 0, 1, 1, 1],
    ]
    nurbs_ref_kv = [
        [0, 0, 0, 0, 1, 1, 1, 1],
        [0, 0, 1, 1],
    ]

    # define the spline and a copy as reference
    spline = request.getfixturevalue(splinetype)
    ref_spline = spline.copy()

    # elevate order
    spline.elevate_degrees(0)

    # test degrees and knot-vectors
    if splinetype == "bspline_2p2d":
        assert np.allclose(spline.degrees, [3, 2])
        assert are_items_close(spline.knot_vectors, bspline_ref_kv)
    elif splinetype == "nurbs_2p2d":
        assert np.allclose(spline.degrees, [3, 1])
        assert are_items_close(spline.knot_vectors, nurbs_ref_kv)
    else:
        assert np.allclose(spline.degrees, [3, 1])

    # use random query points
    q2D = np_rng.random((10, 2))

    # test evaluation
    assert np.allclose(spline.evaluate(q2D), ref_spline.evaluate(q2D))


@pytest.mark.parametrize("splinetype", nurbs_and_bspline)
def test_reduce_degree(splinetype, request, are_items_close, np_rng):
    """Test the function .reduce_degrees.
    This test also depends on the function .elevate_degrees!"""

    # define the spline and a copy as reference
    spline = request.getfixturevalue(splinetype)
    ref_spline = spline.copy()

    # elevate and reduce order
    spline.elevate_degrees(0)
    spline.reduce_degrees(0)

    # test degrees
    assert np.allclose(spline.degrees, ref_spline.degrees)

    # test knot_vectors
    assert are_items_close(spline.knot_vectors, ref_spline.knot_vectors)

    # use random query points
    q2D = np_rng.random((10, 2))

    # test evaluation
    assert np.allclose(spline.evaluate(q2D), ref_spline.evaluate(q2D))
