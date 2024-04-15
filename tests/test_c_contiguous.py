import numpy as np
import pytest

# used fixtures for parametrization
all_2D_splinetypes = (
    "bspline_2p2d",
    "nurbs_2p2d",
    "bezier_2p2d",
    "rational_bezier_2p2d",
)


@pytest.mark.parametrize("splinetype", all_2D_splinetypes)
def test_c_contiguous_array_input(splinetype, request):
    """Test whether cpp sides receives contiguous array"""

    # define spline and its properties
    spl = request.getfixturevalue(splinetype)
    properties = spl.todict()

    # make f contiguous
    f_contig_orig = np.asarray(properties["control_points"], order="F")
    assert not f_contig_orig.flags["C_CONTIGUOUS"]
    assert f_contig_orig.flags["F_CONTIGUOUS"]

    # set non contiguous array as prop
    properties["control_points"] = f_contig_orig

    # check round_trip cps
    round_trip_cps = spl.cps

    # they should be the same
    assert np.allclose(round_trip_cps, properties["control_points"])

    # test setter
    # make sure this array is still not c contiguous
    assert not f_contig_orig.flags["C_CONTIGUOUS"]

    spl.cps = f_contig_orig  # setter

    assert spl.cps.flags["C_CONTIGUOUS"]
    assert np.allclose(spl.cps, properties["control_points"])
