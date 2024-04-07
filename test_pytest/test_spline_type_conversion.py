import numpy as np
import pytest

# for parametrization used fixtures
all_2p2d_splines = (
    "bezier_2p2d",
    "rational_bezier_2p2d",
    "bspline_2p2d",
    "nurbs_2p2d",
)


def convert_and_compare_samples(spline, conversion_list):
    for convert_to in conversion_list:
        converted = getattr(spline, convert_to)

        sample_res = [3] * spline.para_dim
        ref_sample = spline.sample(sample_res)

        # check same sample result
        assert np.allclose(ref_sample, converted.sample(sample_res))

        # check same type
        assert type(converted).__qualname__.lower().startswith(convert_to)


@pytest.mark.parametrize("splinetype", all_2p2d_splines)
def test_conversion(splinetype, request):
    spline = request.getfixturevalue(splinetype)

    if splinetype == "bezier_2p2d":
        test_conversions = ["bezier", "rationalbezier", "bspline", "nurbs"]
    elif splinetype == "rational_bezier_2p2d":
        test_conversions = ["rationalbezier", "nurbs"]
    elif splinetype == "bspline_2p2d":
        test_conversions = ["bspline", "nurbs"]
    elif splinetype == "nurbs_2p2d":
        test_conversions = ["nurbs"]
    else:
        pytest.fail()

    convert_and_compare_samples(spline, test_conversions)
