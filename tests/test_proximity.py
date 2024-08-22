import numpy as np
import pytest

# fixtures used
all_splines = (
    "rational_bezier_2p2d",
    "bezier_2p2d",
    "bspline_2p2d",
    "nurbs_2p2d",
    "rational_bezier_3p3d",
    "bezier_3p3d",
    "bspline_3p3d",
    "nurbs_3p3d",
)


@pytest.mark.parametrize("splinetype", all_splines)
def test_queries_inside_spline_initial_guess_with_kdt(
    request, splinetype, np_rng
):
    """
    Initial guess made with kdt.
    """
    spline = request.getfixturevalue(splinetype)
    # form parametric queries
    # -> parametric space is [0, 1], so no further manipulation
    para_q = np_rng.random((10, spline.para_dim))

    # form physical queries
    phys_q = spline.evaluate(para_q)

    # proximity - single thread exe
    prox_r = spline.proximities(
        queries=phys_q,
        initial_guess_sample_resolutions=[10] * spline.para_dim,
        nthreads=1,
        return_verbose=True,
    )

    assert np.allclose(
        para_q, prox_r[0]
    ), f"WRONG proximity query for {spline.whatami}"
