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
    "1 3 0 1 3 2\n",
    "1 3 2 3 6 4\n",
    "1 3 3 7 8 6\n",
    "1 3 1 5 7 3\n",
    "\n",
    "boundary\n",
    "8\n",
    "1 1 0 1\n",
    "1 1 0 2\n",
    "1 1 4 6\n",
    "1 1 2 4\n",
    "1 1 7 8\n",
    "1 1 6 8\n",
    "1 1 1 5\n",
    "1 1 5 7\n",
    "\n",
    "edges\n",
    "12\n",
    "0 0 1\n",
    "1 0 2\n",
    "1 1 3\n",
    "3 1 5\n",
    "0 2 3\n",
    "2 2 4\n",
    "2 3 6\n",
    "3 3 7\n",
    "0 4 6\n",
    "1 5 7\n",
    "3 6 8\n",
    "2 7 8\n",
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
    "1 5 0 2 6 3 1 4 9 5\n",
    "1 5 3 6 12 7 5 9 14 10\n",
    "1 5 6 13 16 12 9 15 17 14\n",
    "1 5 2 8 13 6 4 11 15 9\n",
    "\n",
    "boundary\n",
    "16\n",
    "1 3 0 2 6 3\n",
    "1 3 0 3 5 1\n",
    "1 3 0 2 4 1\n",
    "1 3 1 4 9 5\n",
    "1 3 3 6 12 7\n",
    "1 3 7 12 14 10\n",
    "1 3 3 7 10 5\n",
    "1 3 5 9 14 10\n",
    "1 3 6 13 16 12\n",
    "1 3 13 16 17 15\n",
    "1 3 12 16 17 14\n",
    "1 3 9 15 17 14\n",
    "1 3 2 8 13 6\n",
    "1 3 8 13 15 11\n",
    "1 3 2 8 11 4\n",
    "1 3 4 11 15 9\n",
    "\n",
    "edges\n",
    "33\n",
    "2 0 1\n",
    "0 0 2\n",
    "1 0 3\n",
    "0 1 4\n",
    "1 1 5\n",
    "2 2 4\n",
    "1 2 6\n",
    "4 2 8\n",
    "2 3 5\n",
    "0 3 6\n",
    "3 3 7\n",
    "1 4 9\n",
    "4 4 11\n",
    "0 5 9\n",
    "3 5 10\n",
    "2 6 9\n",
    "3 6 12\n",
    "4 6 13\n",
    "2 7 10\n",
    "0 7 12\n",
    "2 8 11\n",
    "1 8 13\n",
    "3 9 14\n",
    "4 9 15\n",
    "0 10 14\n",
    "1 11 15\n",
    "2 12 14\n",
    "4 12 16\n",
    "2 13 15\n",
    "3 13 16\n",
    "4 14 17\n",
    "3 15 17\n",
    "2 16 17\n",
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
        bez_el0.elevate_degree(0)
        bsp_el2.elevate_degree(0)

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
        bez_el0.elevate_degree(0)
        bsp_el2.elevate_degree(0)
        bsp_el2.insert_knots(1, [0.5])
        nur_el3.insert_knots(1, [0.5])
        for s in [bez_el0, bsp_el2, nur_el3, rbz_el1]:
            s.elevate_degree(2)

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
