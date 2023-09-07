import tempfile
from sys import version as python_version

try:
    from . import common as c
except BaseException:
    import common as c


class gismoExportTest(c.unittest.TestCase):
    def test_gismo_export(self):
        """
        Test gismo export routine
        """
        if int(python_version.split(".")[1]) < 8:
            c.splinepy.utils.log.info(
                "gismo export is only tested here from python3.8+. "
                "Skipping test, because current version is: "
                f"{python_version}"
            )
            return True

        ###########
        # 2D Mesh #
        ###########
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

        # Init multipatch
        multipatch = c.splinepy.Multipatch(
            splines=[bez_el0, rbz_el1, bsp_el2, nur_el3]
        )

        # Define some functions for boundary identification
        def is_bottom(x):
            return x[:, 0] < 0.01

        def is_top(x):
            return x[:, 0] > 1.99

        # Add boundary
        multipatch.boundary_from_function(is_bottom)
        multipatch.boundary_from_function(is_top)

        # Test Output
        with tempfile.TemporaryDirectory() as tmpd:
            tmpf = c.to_tmpf(tmpd)
            c.splinepy.io.gismo.export(
                tmpf,
                multipatch=multipatch,
                indent=False,
                labeled_boundaries=False,
            )

            with open(tmpf) as tmp_read, open(
                c.os.path.dirname(__file__)
                + "/data/gismo_noindent_nolabels_ascii_2d.xml"
            ) as base_file:
                self.assertTrue(
                    c.are_stripped_lines_same(
                        base_file.readlines(), tmp_read.readlines(), True
                    )
                )

        # for python version > 3.9, test indented version
        if int(python_version.split(".")[1]) >= 9:
            with tempfile.TemporaryDirectory() as tmpd:
                tmpf = c.to_tmpf(tmpd)
                c.splinepy.io.gismo.export(
                    tmpf,
                    multipatch=multipatch,
                    indent=True,
                    labeled_boundaries=False,
                )

                with open(tmpf) as tmp_read, open(
                    c.os.path.dirname(__file__)
                    + "/data/gismo_indent_nolabels_ascii_2d.xml"
                ) as base_file:
                    self.assertTrue(
                        c.are_stripped_lines_same(
                            base_file.readlines(), tmp_read.readlines(), True
                        )
                    )

        ########################
        # 2D Mesh - new format #
        ########################
        with tempfile.TemporaryDirectory() as tmpd:
            tmpf = c.to_tmpf(tmpd)
            c.splinepy.io.gismo.export(
                tmpf,
                multipatch=multipatch,
                indent=False,
                labeled_boundaries=True,
            )

            with open(tmpf) as tmp_read, open(
                c.os.path.dirname(__file__)
                + "/data/gismo_noindent_labels_ascii_2d.xml"
            ) as base_file:
                self.assertTrue(
                    c.are_stripped_lines_same(
                        base_file.readlines(), tmp_read.readlines(), True
                    )
                )

        # for python version > 3.9, test indented version
        if int(python_version.split(".")[1]) >= 9:
            with tempfile.TemporaryDirectory() as tmpd:
                tmpf = c.to_tmpf(tmpd)
                c.splinepy.io.gismo.export(
                    tmpf,
                    multipatch=multipatch,
                    indent=True,
                    labeled_boundaries=True,
                )

                with open(tmpf) as tmp_read, open(
                    c.os.path.dirname(__file__)
                    + "/data/gismo_indent_labels_ascii_2d.xml"
                ) as base_file:
                    self.assertTrue(
                        c.are_stripped_lines_same(
                            base_file.readlines(), tmp_read.readlines(), True
                        )
                    )

        ###########
        # 3D Mesh #
        ###########

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

        # Define some functions for boundary identification
        def is_bottom(x):
            return x[:, 0] < 0.01

        def is_top(x):
            return x[:, 0] > 1.99

        # Add boundary
        multipatch = c.splinepy.Multipatch(
            splines=[bez_el0, bsp_el2, nur_el3, rbz_el1]
        )
        multipatch.boundary_from_function(is_bottom)
        multipatch.boundary_from_function(is_top)

        # Test output
        with tempfile.TemporaryDirectory() as tmpd:
            tmpf = c.to_tmpf(tmpd)
            c.splinepy.io.gismo.export(
                tmpf,
                multipatch=multipatch,
                indent=False,
                labeled_boundaries=False,
            )

            with open(tmpf) as tmp_read, open(
                c.os.path.dirname(__file__)
                + "/data/gismo_noindent_nolabels_ascii_3d.xml"
            ) as base_file:
                self.assertTrue(
                    c.are_stripped_lines_same(
                        base_file.readlines(), tmp_read.readlines(), True
                    )
                )

        # for python version > 3.9, test indented version
        if int(python_version.split(".")[1]) >= 9:
            with tempfile.TemporaryDirectory() as tmpd:
                tmpf = c.to_tmpf(tmpd)
                c.splinepy.io.gismo.export(
                    tmpf,
                    multipatch=multipatch,
                    indent=True,
                    labeled_boundaries=False,
                )

                with open(tmpf) as tmp_read, open(
                    c.os.path.dirname(__file__)
                    + "/data/gismo_indent_nolabels_ascii_3d.xml"
                ) as base_file:
                    self.assertTrue(
                        c.are_stripped_lines_same(
                            base_file.readlines(), tmp_read.readlines(), True
                        )
                    )

    def test_gismo_import(self):
        """
        Test gismo export routine
        """
        if int(python_version.split(".")[1]) < 8:
            c.splinepy.utils.log.info(
                "gismo export is only tested here from python3.8+. "
                "Skipping test, because current version is: "
                f"{python_version}"
            )
            return True

        # Define some splines
        bsp_el2 = c.splinepy.BSpline(
            degrees=[1, 1],
            control_points=[[0, 0], [1, 0], [0, 1], [1, 1]],
            knot_vectors=[[0, 0, 1, 1], [0, 0, 1, 1]],
        )
        nur_el3 = c.splinepy.NURBS(
            degrees=[1, 1],
            control_points=[[1, 0], [2, 0], [1, 1], [2, 1]],
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

        # Test Output against input
        multipatch_geometry = c.splinepy.Multipatch([bsp_el2, nur_el3])
        with tempfile.TemporaryDirectory() as tmpd:
            tmpf = c.to_tmpf(tmpd)
            c.splinepy.io.gismo.export(
                tmpf,
                multipatch=multipatch_geometry,
                indent=False,
                labeled_boundaries=False,
            )
            multipatch_geometry_loaded = c.splinepy.io.gismo.load(
                tmpf, load_options=False
            )
            self.assertTrue(
                all(
                    c.are_splines_equal(a, b)
                    for a, b in zip(
                        multipatch_geometry.patches,
                        multipatch_geometry_loaded.patches,
                    )
                )
            )
            self.assertTrue(
                c.np.allclose(
                    multipatch_geometry.interfaces,
                    multipatch_geometry_loaded.interfaces,
                )
            )

            # Now with modified boundaries
            multipatch_geometry.boundaries_from_continuity()
        with tempfile.TemporaryDirectory() as tmpd:
            tmpf = c.to_tmpf(tmpd)
            c.splinepy.io.gismo.export(
                tmpf,
                multipatch=multipatch_geometry,
                indent=False,
                labeled_boundaries=True,
            )
            multipatch_geometry_loaded = c.splinepy.io.gismo.load(
                tmpf, load_options=False
            )
            self.assertTrue(
                all(
                    c.are_splines_equal(a, b)
                    for a, b in zip(
                        multipatch_geometry.patches,
                        multipatch_geometry_loaded.patches,
                    )
                )
            )
            self.assertTrue(
                c.np.allclose(
                    multipatch_geometry.interfaces,
                    multipatch_geometry_loaded.interfaces,
                )
            )

    def test_gismo_import_with_options(self):
        """
        Test gismo export routine
        """
        if int(python_version.split(".")[1]) < 8:
            c.splinepy.utils.log.info(
                "gismo export is only tested here from python3.8+. "
                "Skipping test, because current version is: "
                f"{python_version}"
            )
            return True

        # Define some splines
        bsp_el2 = c.splinepy.BSpline(
            degrees=[1, 1],
            control_points=[[0, 0], [1, 0], [0, 1], [1, 1]],
            knot_vectors=[[0, 0, 1, 1], [0, 0, 1, 1]],
        )
        nur_el3 = c.splinepy.NURBS(
            degrees=[1, 1],
            control_points=[[1, 0], [2, 0], [1, 1], [2, 1]],
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

        # Test Output against input
        multipatch_geometry = c.splinepy.Multipatch([bsp_el2, nur_el3])

        # Set some options
        gismo_options = [
            {
                "tag": "OptionNumber1",
                "attributes": {"Mambo": "No. 5", "id": "5"},
                "text": "One, two, three, four, five\nEverybody in the car, "
                "so come on, let's ride",
                "children": [
                    {
                        "tag": "AnotherOne",
                        "text": "0 ,0",
                        "attributes": {},
                        "children": [],
                    }
                ],
            }
        ]
        with tempfile.TemporaryDirectory() as tmpd:
            tmpf = c.to_tmpf(tmpd)
            c.splinepy.io.gismo.export(
                tmpf,
                multipatch=multipatch_geometry,
                indent=False,
                labeled_boundaries=False,
                options=gismo_options,
            )
            (
                multipatch_geometry_loaded,
                gismo_options_loaded,
            ) = c.splinepy.io.gismo.load(tmpf, load_options=True)
            self.assertTrue(
                all(
                    c.are_splines_equal(a, b)
                    for a, b in zip(
                        multipatch_geometry.patches,
                        multipatch_geometry_loaded.patches,
                    )
                )
            )
            self.assertTrue(
                c.np.allclose(
                    multipatch_geometry.interfaces,
                    multipatch_geometry_loaded.interfaces,
                )
            )
            self.assertEqual(gismo_options_loaded, gismo_options)

    def test_gismo_io_binary(self):
        """Test the base64 io-routines"""
        # We test this with just one (big, 3D) spline
        nurbs = c.splinepy.NURBS(
            degrees=[1, 1, 1],
            knot_vectors=[[0, 0, 1, 1], [0, 0, 1, 1], [0, 0, 1, 1]],
            control_points=c.np.ones((8, 3)),
            weights=c.np.ones((8, 1)),
        )
        nurbs.elevate_degrees([0, 1, 2, 2])
        c.np.random.seed(19284918)
        for i in range(3):
            nurbs.insert_knots(i, c.np.random.random(4))

        # Randomize points
        nurbs.cps = c.np.random.random(nurbs.cps.shape)
        nurbs.weights = c.np.random.random(nurbs.weights.shape)

        # Create a multipatch geometry
        multipatch_geometry = c.splinepy.Multipatch([nurbs])

        # Export
        with tempfile.TemporaryDirectory() as tmpd:
            tmpf = c.to_tmpf(tmpd)
            c.splinepy.io.gismo.export(
                tmpf,
                multipatch=multipatch_geometry,
                indent=False,
                labeled_boundaries=False,
                as_base64=True,
            )
            (
                multipatch_geometry_loaded,
                gismo_options_loaded,
            ) = c.splinepy.io.gismo.load(tmpf, load_options=True)

            with open(tmpf) as tmp_read, open(
                c.os.path.dirname(__file__)
                + "/data/gismo_noindent_nolabels_b64_3d.xml"
            ) as base_file:
                self.assertTrue(
                    c.are_stripped_lines_same(
                        base_file.readlines(), tmp_read.readlines(), True
                    )
                )
            self.assertTrue(
                all(
                    c.are_splines_equal(a, b)
                    for a, b in zip(
                        multipatch_geometry.patches,
                        multipatch_geometry_loaded.patches,
                    )
                )
            )
            self.assertTrue(
                c.np.allclose(
                    multipatch_geometry.interfaces,
                    multipatch_geometry_loaded.interfaces,
                )
            )


if __name__ == "__main__":
    c.unittest.main()
