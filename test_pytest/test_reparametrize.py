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
def test_permute_parametric_axes(request, splinetype, get_queries_3D):
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
    queries = np.asarray(get_queries_3D)

    assert np.allclose(
        spline.evaluate(queries),
        perm.evaluate(queries[:, permutation]),
    ), f"{perm.whatami} failed to permute."

    # inplace
    perm = spline.copy()
    splinepy.helpme.reparametrize.permute_parametric_axes(
        perm, permutation, inplace=True
    )
    queries = np.asarray(get_queries_3D)

    assert np.allclose(
        spline.evaluate(queries),
        perm.evaluate(queries[:, permutation]),
    ), f"{perm.whatami} failed to permute inplace."
