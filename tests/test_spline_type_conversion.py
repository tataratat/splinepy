try:
    from . import common as c
except BaseException:
    import common as c


def convert_and_compare_samples(spline, conversion_list):
    for convert_to in conversion_list:
        converted = getattr(spline, convert_to)

        sample_res = [3] * spline.para_dim
        ref_sample = spline.sample(sample_res)

        # check same sample result
        assert c.np.allclose(ref_sample, converted.sample(sample_res))

        # check same type
        assert type(converted).__qualname__.lower().startswith(convert_to)


class SplineTypeConversionTest(c.SplineBasedTestCase):
    def setUp(self):
        splines = c.spline_types_as_list()
        self.b = splines[0]
        self.n = splines[1]
        self.z = splines[2]
        self.r = splines[3]

    def test_bezier_conversions(self):
        test_conversions = ["bezier", "rationalbezier", "bspline", "nurbs"]

        convert_and_compare_samples(self.z, test_conversions)

    def test_rational_bezier_conversions(self):
        test_conversions = ["rationalbezier", "nurbs"]

        convert_and_compare_samples(self.r, test_conversions)

    def test_bspline_conversions(self):
        test_conversions = ["bspline", "nurbs"]

        convert_and_compare_samples(self.b, test_conversions)

    def test_nurbs_conversions(self):
        test_conversions = ["nurbs"]

        convert_and_compare_samples(self.n, test_conversions)


if __name__ == "__main__":
    c.unittest.main()
