import tempfile

try:
    from . import common as c
except BaseException:
    import common as c


class npzExportTest(c.unittest.TestCase):
    def test_gismo_import(self):
        """
        Test npz export-import routine
        """
        bez_el0 = c.splinepy.Bezier(
            degrees=[1, 1], control_points=[[0, 0], [1, 0], [0, 1], [1, 1]]
        )
        rbz_el1 = c.splinepy.RationalBezier(
            degrees=[1, 1],
            control_points=[[1, 0], [2, 0], [1, 1], [2, 1]],
            weights=[1, 1, 1, 1],
        )
        bsp_el2 = c.splinepy.BSpline(
            degrees=[1, 1],
            control_points=[[0, 1], [1, 1], [0, 2], [1, 2]],
            knot_vectors=[[0, 0, 1, 1], [0, 0, 1, 1]],
        )
        nur_el3 = c.splinepy.NURBS(
            degrees=[1, 1],
            control_points=[[1, 1], [2, 1], [1, 2], [2, 2]],
            weights=[1, 1, 1, 1],
            knot_vectors=[[0, 0, 1, 1], [0, 0, 1, 1]],
        )

        # Make it more tricky
        bsp_el2.elevate_degrees(0)
        bsp_el2.elevate_degrees(1)
        nur_el3.elevate_degrees(0)
        nur_el3.elevate_degrees(1)
        bsp_el2.insert_knots(1, [0.5])
        nur_el3.insert_knots(1, [0.5])

        # Test output against input
        list_of_splines = [bez_el0, rbz_el1, bsp_el2, nur_el3]
        for spline_in_list in list_of_splines:
            with tempfile.TemporaryDirectory() as tmpd:
                tmpf = c.to_tmpf(tmpd)
                c.splinepy.io.npz.export(
                    tmpf, spline_in_list
                )
                spline_in_list_loaded = c.splinepy.io.npz.load(tmpf + ".npz")[0]
                self.assertTrue(
                        c.are_splines_equal(spline_in_list, spline_in_list_loaded)
                )

if __name__ == "__main__":
    c.unittest.main()
