import os
import tempfile
from sys import version as python_version

import gustaf as gus
import numpy as np

import splinepy

try:
    import vedo  # noqa F401

    has_vedo = True
except ImportError:
    has_vedo = False


def test_svg_export(to_tmpf, are_stripped_lines_same):
    """
    Test export routine of svg module
    """
    if not has_vedo:
        return

    bez_el0 = splinepy.Bezier(
        degrees=[1, 1], control_points=[[0, 0], [1, 0], [0, 1], [1, 1]]
    )
    rbz_el1 = splinepy.RationalBezier(
        degrees=[1, 1],
        control_points=[[1, 0], [2, 0], [1, 1], [2, 1]],
        weights=[1, 1, 1, 1],
    )
    bsp_el2 = splinepy.BSpline(
        degrees=[1, 1],
        control_points=[[0, 1], [1, 1], [0, 2], [1, 2]],
        knot_vectors=[[0, 0, 1, 1], [0, 0, 1, 1]],
    )
    nur_el3 = splinepy.NURBS(
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

    # Set some show_options
    nur_el3.show_options["control_points"] = False
    bez_el0.show_options["c"] = "grey"
    bez_el0.show_options["alpha"] = 0.5
    bez_el0.show_options["control_mesh_lw"] = 0.2
    bez_el0.show_options["control_mesh_c"] = "orange"
    bez_el0.show_options["knot_c"] = "blue"
    bsp_el2.show_options["control_point_ids"] = False
    nur_el3.show_options["knots"] = False

    # Test Output
    with tempfile.TemporaryDirectory() as tmpd:
        tmpf = to_tmpf(tmpd)

        # With a list
        splinepy.io.svg.export(
            tmpf,
            [bez_el0, rbz_el1, bsp_el2, nur_el3],
            indent=False,
            box_margins=0.7,
            font_family="sans-serif",
        )

        with open(tmpf) as tmp_read, open(
            os.path.dirname(__file__) + "/../data/svg_noindent.svg"
        ) as base_file:
            assert are_stripped_lines_same(
                base_file.readlines(), tmp_read.readlines(), True
            )

        # Without a list
        splinepy.io.svg.export(
            tmpf,
            bez_el0,
            rbz_el1,
            bsp_el2,
            nur_el3,
            indent=False,
            box_margins=0.7,
            font_family="sans-serif",
        )

        with open(tmpf) as tmp_read, open(
            os.path.dirname(__file__) + "/../data/svg_noindent.svg"
        ) as base_file:
            assert are_stripped_lines_same(
                base_file.readlines(), tmp_read.readlines(), True
            )

    # for python version > 3.9, test indented version
    if int(python_version.split(".")[1]) >= 9:
        with tempfile.TemporaryDirectory() as tmpd:
            tmpf = to_tmpf(tmpd)
            # With a list and also more options
            splinepy.io.svg.export(
                tmpf,
                [bez_el0, rbz_el1, bsp_el2, nur_el3],
                indent=True,
                box_margins=0.5,
                # Set some show_options
                c="blue",
                alpha=0.2,
                control_point_r=0.1,
                knots=False,
                control_points=True,
                control_mesh_c="purple",
                control_point_c="black",
                control_mesh_lw=0.05,
                font_family="sans-serif",
            )

            with open(tmpf) as tmp_read, open(
                os.path.dirname(__file__) + "/../data/svg_indent.svg"
            ) as base_file:
                assert are_stripped_lines_same(
                    base_file.readlines(), tmp_read.readlines(), True
                )

    # Check some spline_color and arrows
    askew_spline2D = splinepy.Bezier(
        degrees=[2, 2],
        control_points=[
            [0.0, 0.0],
            [1.0, 0.5],
            [2.0, 0.2],
            [0.5, 1.5],
            [1.0, 1.5],
            [1.5, 1.5],
            [0.0, 3.0],
            [1.0, 2.5],
            [2.0, 3.0],
        ],
    )

    def color_function(data, on):  # noqa ARG001
        """Some random function for colors"""
        return np.hstack(
            (
                np.cos(2 * np.pi * on[:, 0]).reshape(-1, 1),
                np.sin(2 * np.pi * on[:, 1]).reshape(-1, 1),
            )
        )

    plot_func_data = splinepy.SplineDataAdaptor(
        askew_spline2D, function=color_function
    )
    askew_spline2D.spline_data["field_function"] = plot_func_data
    askew_spline2D.show_options["data"] = "field_function"
    askew_spline2D.show_options["knot_c"] = "black"
    askew_spline2D.show_options["control_points"] = False
    askew_spline2D.show_options["resolutions"] = 100
    with tempfile.TemporaryDirectory() as tmpd:
        tmpf = to_tmpf(tmpd)

        # Without a list
        splinepy.io.svg.export(
            tmpf,
            askew_spline2D,
            indent=False,
            box_margins=0.7,
        )

        with open(tmpf) as tmp_read, open(
            os.path.dirname(__file__) + "/../data/svg_data_field.svg"
        ) as base_file:
            assert are_stripped_lines_same(
                base_file.readlines(), tmp_read.readlines(), True
            )

    # Plot arrows
    askew_spline2D.show_options["arrow_data_scale"] = 0.15
    askew_spline2D.show_options.pop("data")
    askew_spline2D.show_options["arrow_data"] = "field_function"
    askew_spline2D.show_options["knot_lw"] = 0.4
    askew_spline2D.show_options["c"] = (100, 100, 100)
    askew_spline2D.show_options["arrow_data_on"] = (
        splinepy.utils.data.cartesian_product(
            [np.linspace(0, 1, 20) for _ in range(2)]
        )
    )
    with tempfile.TemporaryDirectory() as tmpd:
        tmpf = to_tmpf(tmpd)

        # Without a list
        splinepy.io.svg.export(
            tmpf,
            askew_spline2D,
            indent=False,
            box_margins=0.7,
        )

        with open(tmpf) as tmp_read, open(
            os.path.dirname(__file__) + "/../data/svg_arrow_data.svg"
        ) as base_file:
            assert are_stripped_lines_same(
                base_file.readlines(), tmp_read.readlines(), True
            )

    # One last example minimal to reuse file for documentation

    spline = splinepy.Bezier(
        degrees=[2, 2],
        control_points=[
            [1, 0],
            [3, 1],
            [5, 0],  # First row
            [0, 2],
            [2, 3],
            [4, 2],  # Second row
            [1, 4],
            [3, 5],
            [5, 4],  # Third row
        ],
    )
    spline.show_options["c"] = (157, 157, 156)
    spline.show_options["knot_c"] = "red"
    spline.show_options["control_mesh"] = True
    spline.show_options["control_mesh_lw"] = 0.05
    spline.show_options["control_point_c"] = (0, 102, 157)
    spline.show_options["control_mesh_c"] = (0, 102, 157)
    spline.show_options["control_point_ids"] = False
    spline.show_options["control_point_r"] = 0.2

    with tempfile.TemporaryDirectory() as tmpd:
        tmpf = to_tmpf(tmpd)
        splinepy.io.svg.export(
            tmpf, spline, indent=False, box_margins=0.2, background=None
        )

        with open(tmpf) as tmp_read, open(
            os.path.dirname(__file__) + "/../data/svg_mini_example.svg"
        ) as base_file:
            assert are_stripped_lines_same(
                base_file.readlines(), tmp_read.readlines(), True
            )

    # Test Export of gus vertices
    vertices = gus.Vertices(spline.cps)
    vertices.vertex_data["field"] = np.linspace(0, 1, spline.cps.shape[0])

    with tempfile.TemporaryDirectory() as tmpd:
        tmpf = to_tmpf(tmpd)
        splinepy.io.svg.export(
            tmpf,
            vertices,
            indent=False,
            box_margins=0.4,
            background="black",
        )

        with open(tmpf) as tmp_read, open(
            os.path.dirname(__file__) + "/../data/svg_vertices.svg"
        ) as base_file:
            assert are_stripped_lines_same(
                base_file.readlines(), tmp_read.readlines(), True
            )
