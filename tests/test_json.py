import tempfile
from sys import version as python_version

try:
    from . import common as c
except BaseException:
    import common as c

class jsonExportTest(c.unittest.TestCase):
    def test_gismo_import(self):
        """
        Test json export-import routine
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
        bsp_el2.elevate_degree(0)
        bsp_el2.elevate_degree(1)
        nur_el3.elevate_degree(0)
        nur_el3.elevate_degree(1)
        bsp_el2.insert_knots(1, [0.5])
        nur_el3.insert_knots(1, [0.5])

        # Test Output against input
        list_of_splines = [bez_el0,rbz_el1,bsp_el2, nur_el3]
        with tempfile.NamedTemporaryFile() as tmpf:
            c.splinepy.io.json.export(
    tmpf.name, list_of_splines, base64encoding=True)
            list_of_splines_loaded = c.splinepy.io.json.load(tmpf.name)
            self.assertTrue(
                all(
                    [
                        c.are_splines_equal(a, b)
                        for a, b in zip(
                            list_of_splines, list_of_splines_loaded
                        )
                    ]
                )
            )
        with tempfile.NamedTemporaryFile() as tmpf:
            c.splinepy.io.json.export(
    tmpf.name, list_of_splines, base64encoding=False)
            list_of_splines_loaded = c.splinepy.io.json.load(tmpf.name)
            self.assertTrue(
                all(
                    [
                        c.are_splines_equal(a, b)
                        for a, b in zip(
                            list_of_splines, list_of_splines_loaded
                        )
                    ]
                )
            )


if __name__ == "__main__":
    c.unittest.main()
