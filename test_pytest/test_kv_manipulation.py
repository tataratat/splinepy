import numpy as np

# def setUp():
#   bspline = bspline_2p2d()
#   nurbs = nurbs_2p2d()


def test_insert_knot(np_rng, are_items_close, bspline_2p2d, nurbs_2p2d):
    """Test the knot insertion function (.insert_knot())."""

    # creating a reference copy of the splines
    bspline_copy = bspline_2p2d
    nurbs_copy = nurbs_2p2d

    # reference solutions
    bspline_ref_kv = [
        [0, 0, 0, 0.2, 0.5, 0.7, 1, 1, 1],
        [0, 0, 0, 0.2, 0.5, 0.7, 1, 1, 1],
    ]
    nurbs_ref_kv = [
        [0, 0, 0, 0.2, 0.5, 0.7, 1, 1, 1],
        [0, 0, 0.2, 0.5, 0.7, 1, 1],
    ]

    # insert knots
    bspline_2p2d.insert_knots(0, [0.2, 0.7])
    bspline_2p2d.insert_knots(
        1,
        [
            0.2,
            0.5,
            0.7,
        ],
    )

    nurbs_2p2d.insert_knots(
        0,
        [
            0.2,
            0.5,
            0.7,
        ],
    )
    nurbs_2p2d.insert_knots(
        1,
        [
            0.2,
            0.5,
            0.7,
        ],
    )

    # test knot_vectors
    assert are_items_close(bspline_2p2d.knot_vectors, bspline_ref_kv)

    assert are_items_close(nurbs_2p2d.knot_vectors, nurbs_ref_kv)

    # use random query points
    q2D = np_rng.random((10, 2))

    # test evaluation
    assert np.allclose(bspline_2p2d.evaluate(q2D), bspline_copy.evaluate(q2D))

    assert np.allclose(nurbs_2p2d.evaluate(q2D), nurbs_copy.evaluate(q2D))


def test_insert_knot_with_matrix(np_rng, bspline_2p2d, nurbs_2p2d):
    """Test the knot insertion function (.insert_knot())."""

    # creating a reference copy of the splines
    bspline_copy = bspline_2p2d
    nurbs_copy = nurbs_2p2d

    # BSpline Data
    b_spline_knots_0 = np_rng.random(8)
    b_spline_knots_1 = np_rng.random(9)
    matrix_bspline = bspline_2p2d.knot_insertion_matrix(0, b_spline_knots_0)
    bspline_2p2d.insert_knots(
        0,
        b_spline_knots_0,
    )
    matrix_bspline = (
        bspline_2p2d.knot_insertion_matrix(1, b_spline_knots_1)
        @ matrix_bspline
    )
    bspline_2p2d.insert_knots(1, b_spline_knots_1)

    # NURBS Data
    nurbs_knots_0 = np_rng.random(10)
    nurbs_knots_1 = np_rng.random(11)
    matrix_nurbs = nurbs_2p2d.knot_insertion_matrix(0, nurbs_knots_0)
    nurbs_2p2d.insert_knots(0, nurbs_knots_0)
    matrix_nurbs = (
        nurbs_2p2d.knot_insertion_matrix(1, nurbs_knots_1) @ matrix_nurbs
    )
    nurbs_2p2d.insert_knots(1, nurbs_knots_1)

    # use random query points
    q2D = np_rng.random((50, 2))

    # test evaluation
    assert np.allclose(bspline_2p2d.evaluate(q2D), bspline_copy.evaluate(q2D))

    assert np.allclose(nurbs_2p2d.evaluate(q2D), nurbs_copy.evaluate(q2D))

    # Test control points and weights
    assert np.allclose(
        bspline_2p2d.control_points,
        matrix_bspline @ bspline_copy.control_points,
    )

    assert np.allclose(nurbs_2p2d.weights, matrix_nurbs @ nurbs_copy.weights)

    assert np.allclose(
        nurbs_2p2d.control_points,
        matrix_nurbs
        @ (nurbs_copy.weights * nurbs_copy.control_points)
        / nurbs_2p2d.weights,
    )


def test_remove_knot(np_rng, are_items_close, bspline_2p2d, nurbs_2p2d):
    """Test the function .remove_knots.
    This test also depends on the function .insert_knots!"""

    # creating a reference copy of the splines
    bspline_copy = bspline_2p2d
    nurbs_copy = nurbs_2p2d

    # insert and remove knots
    bspline_2p2d.insert_knots(
        0,
        [
            0.2,
            0.7,
        ],
    )
    bspline_2p2d.insert_knots(
        1,
        [
            0.2,
            0.5,
            0.7,
        ],
    )
    bspline_2p2d.remove_knots(0, [0.2, 0.7])
    bspline_2p2d.remove_knots(
        1,
        [
            0.2,
            0.5,
            0.7,
        ],
    )

    nurbs_2p2d.insert_knots(
        0,
        [
            0.2,
            0.5,
            0.7,
        ],
    )
    nurbs_2p2d.insert_knots(
        1,
        [
            0.2,
            0.5,
            0.7,
        ],
    )
    nurbs_2p2d.remove_knots(
        0,
        [
            0.2,
            0.5,
            0.7,
        ],
    )
    nurbs_2p2d.remove_knots(
        1,
        [
            0.2,
            0.5,
            0.7,
        ],
    )

    # test knot_vectors
    assert are_items_close(
        bspline_2p2d.knot_vectors, bspline_copy.knot_vectors
    )

    assert are_items_close(nurbs_2p2d.knot_vectors, nurbs_copy.knot_vectors)

    # use random query points
    q2D = np_rng.random((10, 2))

    # test evaluation
    assert np.allclose(bspline_2p2d.evaluate(q2D), bspline_copy.evaluate(q2D))

    assert np.allclose(nurbs_2p2d.evaluate(q2D), nurbs_copy.evaluate(q2D))
