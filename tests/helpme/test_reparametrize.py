import numpy as np
import pytest

import splinepy

# fixtures used
all_3p3d_splines = (
    "rational_bezier_3p3d",
    "bezier_3p3d",
    "bspline_3p3d",
    "nurbs_3p3d",
)


@pytest.mark.parametrize("splinetype", all_3p3d_splines)
def test_permute_parametric_axes(request, splinetype, queries_3D):
    """
    test permute
    """
    # Define spline
    spline = request.getfixturevalue(splinetype)

    # define permutation
    permutation = [2, 0, 1]

    # return permuted
    # make more work
    spline.elevate_degrees(1)
    spline.elevate_degrees(2)
    spline.elevate_degrees(2)
    if "knot_vectors" in spline.required_properties:
        spline.insert_knots(0, [0.4, 0.7, 0.8])
        spline.insert_knots(1, [0.1, 0.2])
        spline.insert_knots(2, [0.3, 0.5, 0.6, 0.9])

    perm = splinepy.helpme.reparametrize.permute_parametric_axes(
        spline, permutation, inplace=False
    )
    queries = np.asarray(queries_3D)

    assert np.allclose(
        spline.evaluate(queries),
        perm.evaluate(queries[:, permutation]),
    ), f"{perm.whatami} failed to permute."

    # inplace
    perm = spline.copy()
    splinepy.helpme.reparametrize.permute_parametric_axes(
        perm, permutation, inplace=True
    )
    queries = np.asarray(queries_3D)

    assert np.allclose(
        spline.evaluate(queries),
        perm.evaluate(queries[:, permutation]),
    ), f"{perm.whatami} failed to permute inplace."


@pytest.mark.parametrize("splinetype", all_3p3d_splines)
def test_flip_axes(splinetype, request, queries_3D):
    """
    test invert axes
    """
    # Define spline
    spline = request.getfixturevalue(splinetype)

    flipped_axes = [1, 2]

    queries = np.asarray(queries_3D)
    flipped_queries = queries.copy()
    for axis in flipped_axes:
        # Queries are all in unit cube
        flipped_queries[:, axis] = 1 - flipped_queries[:, axis]

    # With copy
    new_spline = splinepy.helpme.reparametrize.flip_axes(
        spline, flipped_axes, inplace=False
    )

    assert np.allclose(
        spline.evaluate(queries),
        new_spline.evaluate(flipped_queries),
    ), f"{spline.whatami} failed to flip axis."

    # in place
    spline_copy = spline.copy()
    new_spline2 = splinepy.helpme.reparametrize.flip_axes(
        spline_copy, flipped_axes, inplace=True
    )

    # Check if truly in place
    assert (
        new_spline2 is spline_copy
    ), "Operation has not been performed in place"

    assert np.allclose(
        new_spline2.evaluate(queries),
        new_spline.evaluate(queries),
    ), f"{spline.whatami} failed to flip axis."
