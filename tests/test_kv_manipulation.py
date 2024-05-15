import numpy as np

import splinepy


def test_insert_knot(np_rng, are_items_close, bspline_2p2d, nurbs_2p2d):
    """Test the knot insertion function (.insert_knot())."""

    # creating a reference copy of the splines
    bspline_original = bspline_2p2d.copy()
    nurbs_original = nurbs_2p2d.copy()

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
    bspline_2p2d.insert_knots(1, [0.2, 0.5, 0.7])

    nurbs_2p2d.insert_knots(0, [0.2, 0.5, 0.7])
    nurbs_2p2d.insert_knots(1, [0.2, 0.5, 0.7])

    # test knot_vectors
    assert are_items_close(bspline_2p2d.knot_vectors, bspline_ref_kv)

    assert are_items_close(nurbs_2p2d.knot_vectors, nurbs_ref_kv)

    # use random query points
    q2D = np_rng.random((10, 2))

    # test evaluation
    assert np.allclose(
        bspline_2p2d.evaluate(q2D), bspline_original.evaluate(q2D)
    )

    assert np.allclose(nurbs_2p2d.evaluate(q2D), nurbs_original.evaluate(q2D))


def test_insert_knot_status_return(np_rng):
    # create 3d box and loop for bspline families
    box = splinepy.helpme.create.box(1, 1, 1)
    for s in [box.bspline, box.nurbs]:
        s.elevate_degrees([0, 1, 2])
        degrees = s.degrees

        # create random knots to insert
        knots = np_rng.random(5)
        for i in range(s.para_dim):
            # repeat knots for c^-1 repetition - they should all be possible
            repeat = degrees[i] + 1
            repeated = np.repeat(knots, repeat)

            # we create reference
            # last n_knots should fail
            ref = np.ones(repeated.size + knots.size, dtype=bool)
            ref[-knots.size :] = False

            # schuffle
            np_rng.shuffle(repeated)
            repeated = np.append(repeated, knots)

            successful = s.insert_knots(i, repeated)
            assert all(successful == ref)


def test_insert_knot_with_matrix(np_rng, bspline_2p2d, nurbs_2p2d):
    """Test the knot insertion function (.insert_knot())."""

    # creating a reference copy of the splines
    bspline_original = bspline_2p2d.copy()
    nurbs_original = nurbs_2p2d.copy()

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
    assert np.allclose(
        bspline_2p2d.evaluate(q2D), bspline_original.evaluate(q2D)
    )

    assert np.allclose(nurbs_2p2d.evaluate(q2D), nurbs_original.evaluate(q2D))

    # Test control points and weights
    assert np.allclose(
        bspline_2p2d.control_points,
        matrix_bspline @ bspline_original.control_points,
    )

    assert np.allclose(
        nurbs_2p2d.weights, matrix_nurbs @ nurbs_original.weights
    )

    assert np.allclose(
        nurbs_2p2d.control_points,
        matrix_nurbs
        @ (nurbs_original.weights * nurbs_original.control_points)
        / nurbs_2p2d.weights,
    )


def test_remove_knot(np_rng, are_items_close, bspline_2p2d, nurbs_2p2d):
    """Test the function .remove_knots.
    This test also depends on the function .insert_knots!"""

    # creating a reference copy of the splines
    bspline_original = bspline_2p2d.copy()
    nurbs_original = nurbs_2p2d.copy()

    # insert and remove knots
    bspline_2p2d.insert_knots(0, [0.2, 0.7])
    bspline_2p2d.insert_knots(1, [0.2, 0.5, 0.7])
    bspline_2p2d.remove_knots(0, [0.2, 0.7])
    bspline_2p2d.remove_knots(1, [0.2, 0.5, 0.7])

    nurbs_2p2d.insert_knots(0, [0.2, 0.5, 0.7])
    nurbs_2p2d.insert_knots(1, [0.2, 0.5, 0.7])
    nurbs_2p2d.remove_knots(0, [0.2, 0.5, 0.7])
    nurbs_2p2d.remove_knots(1, [0.2, 0.5, 0.7])

    # test knot_vectors
    assert are_items_close(
        bspline_2p2d.knot_vectors, bspline_original.knot_vectors
    )

    assert are_items_close(
        nurbs_2p2d.knot_vectors, nurbs_original.knot_vectors
    )

    # use random query points
    q2D = np_rng.random((10, 2))

    # test evaluation
    assert np.allclose(
        bspline_2p2d.evaluate(q2D), bspline_original.evaluate(q2D)
    )

    assert np.allclose(nurbs_2p2d.evaluate(q2D), nurbs_original.evaluate(q2D))


def test_uniform_refine():
    """test uniform refine from bspline families by
    checking resulting number of elements"""

    # create single element bspline
    b = splinepy.helpme.create.box(1, 2, 3).bspline
    n = b.nurbs

    for spline in [b, n]:
        # check default - one subdivision
        spline.uniform_refine()
        unique_knots1 = spline.unique_knots
        n_elem1 = 1
        for uk in unique_knots1:
            n_elem1 *= len(uk) - 1

        assert n_elem1 == 2**3

        # add different numbers
        n_knots = [1, 2, 3]
        spline.uniform_refine(n_knots=n_knots)
        unique_knots2 = spline.unique_knots
        n_elem2 = 1
        for uk in unique_knots2:
            n_elem2 *= len(uk) - 1

        # compute what it should be
        n_elem2_ref = 1
        for uk, nk in zip(unique_knots1, n_knots):
            n_elem2_ref *= (len(uk) - 1) * (nk + 1)

        assert n_elem2_ref == n_elem2

        # try only one axis
        n_knots = [2]
        p_dim = [0]
        spline.uniform_refine(p_dim, n_knots)
        unique_knots3 = spline.unique_knots
        n_elem3 = 1
        for uk in unique_knots3:
            n_elem3 *= len(uk) - 1
        assert n_elem2 * (n_knots[0] + 1) == n_elem3
