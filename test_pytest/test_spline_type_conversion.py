import numpy as np

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


def test_bezier_conversion(bezier_2p2d):
    test_conversions = ["bezier", "rationalbezier", "bspline", "nurbs"]
    convert_and_compare_samples(bezier_2p2d, test_conversions)


def test_rational_bezier_conversion(rational_bezier_2p2d):
    test_conversions = ["rationalbezier", "nurbs"]
    convert_and_compare_samples(rational_bezier_2p2d, test_conversions)


def test_bspline_conversion(bspline_2p2d):
    test_conversions = ["bspline", "nurbs"]
    convert_and_compare_samples(bspline_2p2d, test_conversions)


def test_nurbs_conversion(nurbs_2p2d):
    test_conversions = ["nurbs"]
    convert_and_compare_samples(nurbs_2p2d, test_conversions)
