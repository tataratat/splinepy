import splinepy
import numpy as np
try:
    from . import common as c
except:
    import common as c


class ProximityTest(c.unittest.TestCase):

    def test_mfem_export(self):
        """
        Test MFEM export routine
        """
        # Define some splines
        bez_el0 = splinepy.Bezier(
            degrees=[1, 1],
            control_points=[
                [0, 0],
                [1, 0],
                [0, 1],
                [1, 1]
            ]
        )
        rbz_el1 = splinepy.RationalBezier(
            degrees=[1, 1],
            control_points=[
                [1, 0],
                [2, 0],
                [1, 1],
                [2, 1]
            ],
            weights=[1, 1, 1, 1]
        )
        bsp_el2 = splinepy.BSpline(
            degrees=[1, 1],
            control_points=[
                [0, 1],
                [1, 1],
                [0, 2],
                [1, 2]
            ],
            knot_vectors=[[0, 0, 1, 1], [0, 0, 1, 1]]
        )
        nur_el3 = splinepy.NURBS(
            degrees=[1, 1],
            control_points=[
                [1, 1],
                [2, 1],
                [1, 2],
                [2, 2]
            ],
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
        splinepy.io.mfem.export_cartesian(
            "test_mfem_r0.mesh", [bez_el0, bsp_el2, nur_el3, rbz_el1])


if __name__ == "__main__":
    c.unittest.main()
