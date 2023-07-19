import tempfile

try:
    from . import common as c
except BaseException:
    import common as c

_mfem_export_ref_2d = [
    "MFEM NURBS mesh v1.0\n",
    "\n",
    "dimension\n",
    "2\n",
    "\n",
    "elements\n",
    "4\n",
    "1 3 0 1 2 3\n",
    "1 3 3 2 4 5\n",
    "1 3 2 6 7 4\n",
    "1 3 1 8 6 2\n",
    "\n",
    "boundary\n",
    "8\n",
    "1 1 0 3\n",
    "1 1 0 1\n",
    "1 1 3 5\n",
    "1 1 5 4\n",
    "1 1 6 7\n",
    "1 1 4 7\n",
    "1 1 8 6\n",
    "1 1 1 8\n",
    "\n",
    "vertices\n",
    "9\n",
    "\n",
    "patches\n",
    "\n",
    "knotvectors\n",
    "2\n",
    "2 3 0.0 0.0 0.0 1.0 1.0 1.0 \n",
    "1 2 0.0 0.0 1.0 1.0 \n",
    "dimension\n",
    "2\n",
    "controlpoints_cartesian\n",
    "0.0 0.0 1.0\n",
    "0.5 0.0 1.0\n",
    "1.0 0.0 1.0\n",
    "0.0 1.0 1.0\n",
    "0.5 1.0 1.0\n",
    "1.0 1.0 1.0\n",
    "\n",
    "knotvectors\n",
    "2\n",
    "2 3 0.0 0.0 0.0 1.0 1.0 1.0\n",
    "1 3 0.0 0.0 0.5 1.0 1.0\n",
    "dimension\n",
    "2\n",
    "controlpoints_cartesian\n",
    "0.0 1.0 1.0\n",
    "0.5 1.0 1.0\n",
    "1.0 1.0 1.0\n",
    "0.0 1.5 1.0\n",
    "0.5 1.5 1.0\n",
    "1.0 1.5 1.0\n",
    "0.0 2.0 1.0\n",
    "0.5 2.0 1.0\n",
    "1.0 2.0 1.0\n",
    "\n",
    "knotvectors\n",
    "2\n",
    "1 2 0.0 0.0 1.0 1.0\n",
    "1 3 0.0 0.0 0.5 1.0 1.0\n",
    "dimension\n",
    "2\n",
    "controlpoints_cartesian\n",
    "1.0 1.0 1.0\n",
    "2.0 1.0 1.0\n",
    "1.0 1.5 1.0\n",
    "2.0 1.5 1.0\n",
    "1.0 2.0 1.0\n",
    "2.0 2.0 1.0\n",
    "\n",
    "knotvectors\n",
    "2\n",
    "1 2 0.0 0.0 1.0 1.0 \n",
    "1 2 0.0 0.0 1.0 1.0 \n",
    "dimension\n",
    "2\n",
    "controlpoints_cartesian\n",
    "1.0 0.0 1.0\n",
    "2.0 0.0 1.0\n",
    "1.0 1.0 1.0\n",
    "2.0 1.0 1.0\n",
    "\n",
]

_mfem_export_ref_3d = [
    "MFEM NURBS mesh v1.0\n",
    "\n",
    "dimension\n",
    "3\n",
    "\n",
    "elements\n",
    "4\n",
    "1 5 0 1 2 3 4 5 6 7\n",
    "1 5 3 2 8 9 7 6 10 11\n",
    "1 5 2 12 13 8 6 14 15 10\n",
    "1 5 1 16 12 2 5 17 14 6\n",
    "\n",
    "boundary\n",
    "16\n",
    "1 3 0 3 7 4\n",
    "1 3 0 1 5 4\n",
    "1 3 0 1 2 3\n",
    "1 3 4 5 6 7\n",
    "1 3 3 9 11 7\n",
    "1 3 9 8 10 11\n",
    "1 3 3 2 8 9\n",
    "1 3 7 6 10 11\n",
    "1 3 12 13 15 14\n",
    "1 3 8 13 15 10\n",
    "1 3 2 12 13 8\n",
    "1 3 6 14 15 10\n",
    "1 3 16 12 14 17\n",
    "1 3 1 16 17 5\n",
    "1 3 1 16 12 2\n",
    "1 3 5 17 14 6\n",
    "\n",
    "vertices\n",
    "18\n",
    "\n",
    "patches\n",
    "\n",
    "knotvectors\n",
    "3\n",
    "2 3 0.0 0.0 0.0 1.0 1.0 1.0 \n",
    "1 2 0.0 0.0 1.0 1.0 \n",
    "2 3 0.0 0.0 0.0 1.0 1.0 1.0 \n",
    "dimension\n",
    "3\n",
    "controlpoints_cartesian\n",
    "0.0 0.0 0.0 1.0\n",
    "0.5 0.0 0.0 1.0\n",
    "1.0 0.0 0.0 1.0\n",
    "0.0 1.0 0.0 1.0\n",
    "0.5 1.0 0.0 1.0\n",
    "1.0 1.0 0.0 1.0\n",
    "0.0 0.0 0.5 1.0\n",
    "0.5 0.0 0.5 1.0\n",
    "1.0 0.0 0.5 1.0\n",
    "0.0 1.0 0.5 1.0\n",
    "0.5 1.0 0.5 1.0\n",
    "1.0 1.0 0.5 1.0\n",
    "0.0 0.0 1.0 1.0\n",
    "0.5 0.0 1.0 1.0\n",
    "1.0 0.0 1.0 1.0\n",
    "0.0 1.0 1.0 1.0\n",
    "0.5 1.0 1.0 1.0\n",
    "1.0 1.0 1.0 1.0\n",
    "\n",
    "knotvectors\n",
    "3\n",
    "2 3 0.0 0.0 0.0 1.0 1.0 1.0\n",
    "1 3 0.0 0.0 0.5 1.0 1.0\n",
    "2 3 0.0 0.0 0.0 1.0 1.0 1.0\n",
    "dimension\n",
    "3\n",
    "controlpoints_cartesian\n",
    "0.0 1.0 0.0 1.0\n",
    "0.5 1.0 0.0 1.0\n",
    "1.0 1.0 0.0 1.0\n",
    "0.0 1.5 0.0 1.0\n",
    "0.5 1.5 0.0 1.0\n",
    "1.0 1.5 0.0 1.0\n",
    "0.0 2.0 0.0 1.0\n",
    "0.5 2.0 0.0 1.0\n",
    "1.0 2.0 0.0 1.0\n",
    "0.0 1.0 0.5 1.0\n",
    "0.5 1.0 0.5 1.0\n",
    "1.0 1.0 0.5 1.0\n",
    "0.0 1.5 0.5 1.0\n",
    "0.5 1.5 0.5 1.0\n",
    "1.0 1.5 0.5 1.0\n",
    "0.0 2.0 0.5 1.0\n",
    "0.5 2.0 0.5 1.0\n",
    "1.0 2.0 0.5 1.0\n",
    "0.0 1.0 1.0 1.0\n",
    "0.5 1.0 1.0 1.0\n",
    "1.0 1.0 1.0 1.0\n",
    "0.0 1.5 1.0 1.0\n",
    "0.5 1.5 1.0 1.0\n",
    "1.0 1.5 1.0 1.0\n",
    "0.0 2.0 1.0 1.0\n",
    "0.5 2.0 1.0 1.0\n",
    "1.0 2.0 1.0 1.0\n",
    "\n",
    "knotvectors\n",
    "3\n",
    "1 2 0.0 0.0 1.0 1.0\n",
    "1 3 0.0 0.0 0.5 1.0 1.0\n",
    "2 3 0.0 0.0 0.0 1.0 1.0 1.0\n",
    "dimension\n",
    "3\n",
    "controlpoints_cartesian\n",
    "1.0 1.0 0.0 1.0\n",
    "2.0 1.0 0.0 1.0\n",
    "1.0 1.5 0.0 1.0\n",
    "2.0 1.5 0.0 1.0\n",
    "1.0 2.0 0.0 1.0\n",
    "2.0 2.0 0.0 1.0\n",
    "1.0 1.0 0.5 1.0\n",
    "2.0 1.0 0.5 1.0\n",
    "1.0 1.5 0.5 1.0\n",
    "2.0 1.5 0.5 1.0\n",
    "1.0 2.0 0.5 1.0\n",
    "2.0 2.0 0.5 1.0\n",
    "1.0 1.0 1.0 1.0\n",
    "2.0 1.0 1.0 1.0\n",
    "1.0 1.5 1.0 1.0\n",
    "2.0 1.5 1.0 1.0\n",
    "1.0 2.0 1.0 1.0\n",
    "2.0 2.0 1.0 1.0\n",
    "\n",
    "knotvectors\n",
    "3\n",
    "1 2 0.0 0.0 1.0 1.0 \n",
    "1 2 0.0 0.0 1.0 1.0 \n",
    "2 3 0.0 0.0 0.0 1.0 1.0 1.0 \n",
    "dimension\n",
    "3\n",
    "controlpoints_cartesian\n",
    "1.0 0.0 0.0 1.0\n",
    "2.0 0.0 0.0 1.0\n",
    "1.0 1.0 0.0 1.0\n",
    "2.0 1.0 0.0 1.0\n",
    "1.0 0.0 0.5 1.0\n",
    "2.0 0.0 0.5 1.0\n",
    "1.0 1.0 0.5 1.0\n",
    "2.0 1.0 0.5 1.0\n",
    "1.0 0.0 1.0 1.0\n",
    "2.0 0.0 1.0 1.0\n",
    "1.0 1.0 1.0 1.0\n",
    "2.0 1.0 1.0 1.0\n",
    "\n",
]


class MFEMExportTest(c.unittest.TestCase):
    def test_mfem_export(self):
        """
        Test MFEM export routine
        """
        # Define some splines
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
        bez_el0.elevate_degrees(0)
        bsp_el2.elevate_degrees(0)

        bsp_el2.insert_knots(1, [0.5])
        nur_el3.insert_knots(1, [0.5])

        # Test Output
        with tempfile.TemporaryDirectory() as tmpd:
            tmpf = c.to_tmpf(tmpd)
            c.splinepy.io.mfem.export_cartesian(
                tmpf, [bez_el0, bsp_el2, nur_el3, rbz_el1]
            )

            with open(tmpf) as tmp_read:
                self.assertTrue(
                    c.are_items_same(_mfem_export_ref_2d, tmp_read.readlines())
                )

        # Test Also 3D Meshes
        bez_el0 = c.splinepy.Bezier(
            degrees=[1, 1, 1],
            control_points=[
                [0, 0, 0],
                [1, 0, 0],
                [0, 1, 0],
                [1, 1, 0],
                [0, 0, 1],
                [1, 0, 1],
                [0, 1, 1],
                [1, 1, 1],
            ],
        )
        rbz_el1 = c.splinepy.RationalBezier(
            degrees=[1, 1, 1],
            control_points=[
                [1, 0, 0],
                [2, 0, 0],
                [1, 1, 0],
                [2, 1, 0],
                [1, 0, 1],
                [2, 0, 1],
                [1, 1, 1],
                [2, 1, 1],
            ],
            weights=[1] * 8,
        )
        bsp_el2 = c.splinepy.BSpline(
            degrees=[1, 1, 1],
            control_points=[
                [0, 1, 0],
                [1, 1, 0],
                [0, 2, 0],
                [1, 2, 0],
                [0, 1, 1],
                [1, 1, 1],
                [0, 2, 1],
                [1, 2, 1],
            ],
            knot_vectors=[[0, 0, 1, 1]] * 3,
        )
        nur_el3 = c.splinepy.NURBS(
            degrees=[1, 1, 1],
            control_points=[
                [1, 1, 0],
                [2, 1, 0],
                [1, 2, 0],
                [2, 2, 0],
                [1, 1, 1],
                [2, 1, 1],
                [1, 2, 1],
                [2, 2, 1],
            ],
            weights=[1] * 8,
            knot_vectors=[[0, 0, 1, 1]] * 3,
        )

        # Make it more tricky
        bez_el0.elevate_degrees(0)
        bsp_el2.elevate_degrees(0)
        bsp_el2.insert_knots(1, [0.5])
        nur_el3.insert_knots(1, [0.5])
        for s in [bez_el0, bsp_el2, nur_el3, rbz_el1]:
            s.elevate_degrees(2)

        # Test output
        with tempfile.TemporaryDirectory() as tmpd:
            tmpf = c.to_tmpf(tmpd)
            c.splinepy.io.mfem.export_cartesian(
                tmpf, [bez_el0, bsp_el2, nur_el3, rbz_el1]
            )

            with open(tmpf) as tmp_read:
                self.assertTrue(
                    c.are_items_same(_mfem_export_ref_3d, tmp_read.readlines())
                )


if __name__ == "__main__":
    c.unittest.main()
