import tempfile

try:
    from . import common as c
except BaseException:
    import common as c


class IOTestIRIT(c.SplineBasedTestCase):
    def test_export_import(self):
        """Test import export routine for irit files"""
        # Create some splines
        bspline = c.bspline_2p2d()
        bezier = c.bezier_3p3d()
        nurbs = c.nurbs_2p2d()
        # Create a line
        rational = c.nurbs_2p2d().extract.spline(
            splitting_plane=0, interval=0.3
        )

        # Make it more tricky
        bspline.elevate_degrees([0, 1, 0])
        nurbs.elevate_degrees([1])
        bezier.elevate_degrees([0, 1, 2])
        nurbs.insert_knots(1, [0.5])
        bspline.insert_knots(1, [0.5])

        list_of_splines = [bezier, nurbs, bspline, rational]

        with tempfile.TemporaryDirectory() as tmpd:
            tmpf = c.to_tmpf(tmpd)
            c.splinepy.io.irit.export(tmpf, list_of_splines)
            list_of_splines_loaded = c.splinepy.io.irit.load(tmpf)
            self.assertTrue(
                all(
                    c.are_splines_equal(a, b)
                    for a, b in zip(list_of_splines, list_of_splines_loaded)
                )
            )


if __name__ == "__main__":
    c.unittest.main()
