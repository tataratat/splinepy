import numpy as np
import pytest

# frequently used fixtures
all_splinetypes = ("bspline_3p3d", "nurbs_3p3d", "bspline_2p2d", "nurbs_2p2d")


@pytest.mark.parametrize("splinetype", all_splinetypes)
def test_knot_vectors(splinetype, request):
    spline = request.getfixturevalue(splinetype)
    # test unique_knots
    copy_knot_vectors = spline.knot_vectors[:]
    unique_knots = [np.unique(ckvs) for ckvs in copy_knot_vectors]

    for uk, uk_fct in zip(unique_knots, spline.unique_knots):
        assert np.allclose(uk, uk_fct)

    # test knot_multiplicities
    multiplicity = [
        np.unique(ckvs, return_counts=True) for ckvs in copy_knot_vectors
    ]

    for m_u, m_fct in zip(multiplicity, spline.knot_multiplicities):
        assert np.allclose(m_u[1], m_fct)

    # test knot_vector creation
    for u_kv, kn_m, spl_kv in zip(
        spline.unique_knots,
        spline.knot_multiplicities,
        spline.knot_vectors,
    ):
        assert np.allclose(np.array(spl_kv), np.repeat(u_kv, kn_m))
