import tempfile

try:
    from . import common as c
except BaseException:
    import common as c


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

            with open(tmpf) as tmp_read, open(
                c.os.path.dirname(__file__) + "/data/mfem_cartesian_2d.mesh"
            ) as base_file:
                self.assertTrue(
                    c.are_stripped_lines_same(
                        base_file.readlines(), tmp_read.readlines(), True
                    )
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

            with open(tmpf) as tmp_read, open(
                c.os.path.dirname(__file__) + "/data/mfem_cartesian_3d.mesh"
            ) as base_file:
                self.assertTrue(
                    c.are_stripped_lines_same(
                        base_file.readlines(), tmp_read.readlines(), True
                    )
                )


if __name__ == "__main__":
    c.unittest.main()
