import tempfile

import splinepy
import numpy as np
try:
    from . import common as c
except BaseException:
    import common as c

_gismo_export_ref_2d = [
]

_gismo_export_ref_3d = [
]


class gismoExportTest(c.unittest.TestCase):

    def test_gismo_export(self):
        """
        Test gismo export routine
        """
        # Define some splines
        bez_el0 = splinepy.Bezier(
            degrees=[1, 1], control_points=[[0, 0], [1, 0], [0, 1], [1, 1]]
        )
        rbz_el1 = splinepy.RationalBezier(
            degrees=[1, 1],
            control_points=[[1, 0], [2, 0], [1, 1], [2, 1]],
            weights=[1, 1, 1, 1]
        )
        bsp_el2 = splinepy.BSpline(
            degrees=[1, 1],
            control_points=[[0, 1], [1, 1], [0, 2], [1, 2]],
            knot_vectors=[[0, 0, 1, 1], [0, 0, 1, 1]]
        )
        nur_el3 = splinepy.NURBS(
            degrees=[1, 1],
            control_points=[[1, 1], [2, 1], [1, 2], [2, 2]],
            weights=[1, 1, 1, 1],
            knot_vectors=[[0, 0, 1, 1], [0, 0, 1, 1]]
        )

        # Make it more tricky
        bez_el0.elevate_degree(0)
        bsp_el2.elevate_degree(0)

        bsp_el2.insert_knots(1, [0.5])
        nur_el3.insert_knots(1, [0.5])

        # Test Output
        # @todo
        with tempfile.NamedTemporaryFile() as tmpf:
            splinepy.io.gismo.export(
                "test2D.xml", [bez_el0, bsp_el2, nur_el3, rbz_el1]
            )

            with open(tmpf.name, "r") as tmp_read:
                self.assertTrue(_gismo_export_ref_2d == tmp_read.readlines())

        # Test Also 3D Meshes
        bez_el0 = splinepy.Bezier(
            degrees=[1, 1, 1],
            control_points=[
                [0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [0, 0, 1],
                [1, 0, 1], [0, 1, 1], [1, 1, 1]
            ]
        )
        rbz_el1 = splinepy.RationalBezier(
            degrees=[1, 1, 1],
            control_points=[
                [1, 0, 0], [2, 0, 0], [1, 1, 0], [2, 1, 0], [1, 0, 1],
                [2, 0, 1], [1, 1, 1], [2, 1, 1]
            ],
            weights=[1] * 8
        )
        bsp_el2 = splinepy.BSpline(
            degrees=[1, 1, 1],
            control_points=[
                [0, 1, 0], [1, 1, 0], [0, 2, 0], [1, 2, 0], [0, 1, 1],
                [1, 1, 1], [0, 2, 1], [1, 2, 1]
            ],
            knot_vectors=[[0, 0, 1, 1]] * 3
        )
        nur_el3 = splinepy.NURBS(
            degrees=[1, 1, 1],
            control_points=[
                [1, 1, 0], [2, 1, 0], [1, 2, 0], [2, 2, 0], [1, 1, 1],
                [2, 1, 1], [1, 2, 1], [2, 2, 1]
            ],
            weights=[1] * 8,
            knot_vectors=[[0, 0, 1, 1]] * 3
        )

        # Make it more tricky
        bez_el0.elevate_degree(0)
        bsp_el2.elevate_degree(0)
        bsp_el2.insert_knots(1, [0.5])
        nur_el3.insert_knots(1, [0.5])
        for s in [bez_el0, bsp_el2, nur_el3, rbz_el1]:
            s.elevate_degree(2)

        # Test output
        with tempfile.NamedTemporaryFile() as tmpf:
            splinepy.io.gismo.export(
                "test.xml", [bez_el0, bsp_el2, nur_el3, rbz_el1]
            )

            with open(tmpf.name, "r") as tmp_read:
                self.assertTrue(_gismo_export_ref_3d == tmp_read.readlines())


if __name__ == "__main__":
    c.unittest.main()
