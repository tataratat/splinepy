import os
import tempfile
from sys import version as python_version

import numpy as np
import pytest

import splinepy


@pytest.fixture
def gismo_2D_multipatch():
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
    bez_el0.elevate_degrees(0)
    bsp_el2.elevate_degrees(0)

    bsp_el2.insert_knots(1, [0.5])
    nur_el3.insert_knots(1, [0.5])

    # Init multipatch
    multipatch = splinepy.Multipatch(
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

    return multipatch


@pytest.fixture
def gismo_multipatch_3D():
    bez_el0 = splinepy.Bezier(
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
    rbz_el1 = splinepy.RationalBezier(
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
    bsp_el2 = splinepy.BSpline(
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
    nur_el3 = splinepy.NURBS(
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
    multipatch = splinepy.Multipatch(
        splines=[bez_el0, bsp_el2, nur_el3, rbz_el1]
    )
    multipatch.boundary_from_function(is_bottom)
    multipatch.boundary_from_function(is_top)
    return multipatch


def test_gismo_export_2D(
    to_tmpf, gismo_2D_multipatch, are_stripped_lines_same
):
    """
    Test gismo export routine for 2D
    """
    ###########
    # 2D Mesh #
    ###########
    # Define some splines

    # Test Output
    with tempfile.TemporaryDirectory() as tmpd:
        tmpf = to_tmpf(tmpd)
        splinepy.io.gismo.export(
            tmpf,
            multipatch=gismo_2D_multipatch,
            indent=False,
            labeled_boundaries=False,
        )

        with open(tmpf) as tmp_read, open(
            os.path.dirname(__file__)
            + "/../data/gismo_noindent_nolabels_ascii_2d.xml"
        ) as base_file:
            assert are_stripped_lines_same(
                base_file.readlines(), tmp_read.readlines(), True
            )


def test_gismo_export_2D_indented(
    gismo_2D_multipatch, to_tmpf, are_stripped_lines_same
):
    # for python version > 3.9, test indented version
    if int(python_version.split(".")[1]) >= 9:
        with tempfile.TemporaryDirectory() as tmpd:
            tmpf = to_tmpf(tmpd)
            splinepy.io.gismo.export(
                tmpf,
                multipatch=gismo_2D_multipatch,
                indent=True,
                labeled_boundaries=False,
            )

            with open(tmpf) as tmp_read, open(
                os.path.dirname(__file__)
                + "/../data/gismo_indent_nolabels_ascii_2d.xml"
            ) as base_file:
                assert are_stripped_lines_same(
                    base_file.readlines(), tmp_read.readlines(), True
                )


def test_gismo_export_2D_labels(
    gismo_2D_multipatch, to_tmpf, are_stripped_lines_same
):
    ########################
    # 2D Mesh - new format #
    ########################
    with tempfile.TemporaryDirectory() as tmpd:
        tmpf = to_tmpf(tmpd)
        splinepy.io.gismo.export(
            tmpf,
            multipatch=gismo_2D_multipatch,
            indent=False,
            labeled_boundaries=True,
        )

        with open(tmpf) as tmp_read, open(
            os.path.dirname(__file__)
            + "/../data/gismo_noindent_labels_ascii_2d.xml"
        ) as base_file:
            assert are_stripped_lines_same(
                base_file.readlines(), tmp_read.readlines(), True
            )


def test_gismo_export_2D_labels_indented(
    gismo_2D_multipatch, to_tmpf, are_stripped_lines_same
):
    # for python version > 3.9, test indented version
    if int(python_version.split(".")[1]) >= 9:
        with tempfile.TemporaryDirectory() as tmpd:
            tmpf = to_tmpf(tmpd)
            splinepy.io.gismo.export(
                tmpf,
                multipatch=gismo_2D_multipatch,
                indent=True,
                labeled_boundaries=True,
            )

            with open(tmpf) as tmp_read, open(
                os.path.dirname(__file__)
                + "/../data/gismo_indent_labels_ascii_2d.xml"
            ) as base_file:
                assert are_stripped_lines_same(
                    base_file.readlines(), tmp_read.readlines(), True
                )


def test_gismo_export_3D(
    gismo_multipatch_3D, to_tmpf, are_stripped_lines_same
):
    # Test Also 3D Meshes

    # Test output
    with tempfile.TemporaryDirectory() as tmpd:
        tmpf = to_tmpf(tmpd)
        splinepy.io.gismo.export(
            tmpf,
            multipatch=gismo_multipatch_3D,
            indent=False,
            labeled_boundaries=False,
        )

        with open(tmpf) as tmp_read, open(
            os.path.dirname(__file__)
            + "/../data/gismo_noindent_nolabels_ascii_3d.xml"
        ) as base_file:
            assert are_stripped_lines_same(
                base_file.readlines(), tmp_read.readlines(), True
            )


def test_gismo_export_additional_blocks(
    gismo_2D_multipatch, to_tmpf, are_stripped_lines_same
):
    """Test if the export of additional xml-blocks are correctly handled by splinepy's
    g+smo export function. This inlucdes scalar and vector-valued function and
    boundary condition blocks and also the default assembly options.
    """
    # Test only runs for Python version 3.9 and higher
    if int(python_version.split(".")[1]) >= 9:
        # Initialize AdditionalBlocks to manage various blocks for export
        additional_blocks = splinepy.io.gismo.AdditionalBlocks()

        # Function block for scalar variable
        additional_blocks.add_function(
            dim=2,
            block_id=1,
            function_string="2*pi^2*sin(pi*x)*sin(pi*y)",
            comment="Scalar-valued function",
        )
        # Function block for vector-valued variable
        additional_blocks.add_function(
            dim=2,
            block_id=2,
            function_string=("x", "y"),
            comment="Vector-valued function",
        )

        # Create block for boundary conditions
        additional_blocks.add_boundary_conditions(
            block_id=3,
            dim=2,
            function_list=[
                "sin(pi*x) * sin(pi*y)",
                ("pi*cos(pi*x) * sin(pi*y)", "pi*sin(pi*x) * cos(pi*y)"),
                "0",
            ],
            bc_list=[("BID2", "Dirichlet", 0), ("BID1", "Neumann", 1)],
            cv_list=[
                ("0", "0", "1", "0"),
                ("0", "0", "2", "sin(x)"),
            ],  # unknown, patch, corner, function
            unknown_id=0,
            multipatch_id=0,
            comment="Boundary conditions",
        )

        # Create dictionary for assembly options
        additional_blocks.add_assembly_options(
            block_id=4, comment="Assembler options"
        )

        with tempfile.TemporaryDirectory() as tmpd:
            tmpf = to_tmpf(tmpd)
            splinepy.io.gismo.export(
                tmpf,
                multipatch=gismo_2D_multipatch,
                indent=True,
                labeled_boundaries=True,
                additional_blocks=additional_blocks.to_list(),
            )

            with open(tmpf) as tmp_read, open(
                os.path.dirname(__file__)
                + "/../data/gismo_additional_blocks.xml"
            ) as base_file:
                assert are_stripped_lines_same(
                    base_file.readlines(), tmp_read.readlines(), True
                )


def test_gismo_export_3D_indented(
    gismo_multipatch_3D, to_tmpf, are_stripped_lines_same
):
    # for python version > 3.9, test indented version
    if int(python_version.split(".")[1]) >= 9:
        with tempfile.TemporaryDirectory() as tmpd:
            tmpf = to_tmpf(tmpd)
            splinepy.io.gismo.export(
                tmpf,
                multipatch=gismo_multipatch_3D,
                indent=True,
                labeled_boundaries=False,
            )

            with open(tmpf) as tmp_read, open(
                os.path.dirname(__file__)
                + "/../data/gismo_indent_nolabels_ascii_3d.xml"
            ) as base_file:
                assert are_stripped_lines_same(
                    base_file.readlines(), tmp_read.readlines(), True
                )


def test_gismo_import(to_tmpf, are_splines_equal):
    """
    Test gismo export routine
    """
    # Define some splines
    bsp_el2 = splinepy.BSpline(
        degrees=[1, 1],
        control_points=[[0, 0], [1, 0], [0, 1], [1, 1]],
        knot_vectors=[[0, 0, 1, 1], [0, 0, 1, 1]],
    )
    nur_el3 = splinepy.NURBS(
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
    multipatch_geometry = splinepy.Multipatch([bsp_el2, nur_el3])
    with tempfile.TemporaryDirectory() as tmpd:
        tmpf = to_tmpf(tmpd)
        splinepy.io.gismo.export(
            tmpf,
            multipatch=multipatch_geometry,
            indent=False,
            labeled_boundaries=False,
        )
        multipatch_geometry_loaded = splinepy.io.gismo.load(
            tmpf, load_options=False
        )
        assert all(
            are_splines_equal(a, b)
            for a, b in zip(
                multipatch_geometry.patches,
                multipatch_geometry_loaded.patches,
            )
        )

        assert np.allclose(
            multipatch_geometry.interfaces,
            multipatch_geometry_loaded.interfaces,
        )

        # Now with modified boundaries
        multipatch_geometry.boundaries_from_continuity()
    with tempfile.TemporaryDirectory() as tmpd:
        tmpf = to_tmpf(tmpd)
        splinepy.io.gismo.export(
            tmpf,
            multipatch=multipatch_geometry,
            indent=False,
            labeled_boundaries=True,
        )
        multipatch_geometry_loaded = splinepy.io.gismo.load(
            tmpf, load_options=False
        )
        assert all(
            are_splines_equal(a, b)
            for a, b in zip(
                multipatch_geometry.patches,
                multipatch_geometry_loaded.patches,
            )
        )

        assert np.allclose(
            multipatch_geometry.interfaces,
            multipatch_geometry_loaded.interfaces,
        )


def test_gismo_import_with_options(to_tmpf, are_splines_equal):
    """
    Test gismo export routine
    """
    if int(python_version.split(".")[1]) < 8:
        splinepy.utils.log.info(
            "gismo export is only tested here from python3.8+. "
            "Skipping test, because current version is: "
            f"{python_version}"
        )
        return True

    # Define some splines
    bsp_el2 = splinepy.BSpline(
        degrees=[1, 1],
        control_points=[[0, 0], [1, 0], [0, 1], [1, 1]],
        knot_vectors=[[0, 0, 1, 1], [0, 0, 1, 1]],
    )
    nur_el3 = splinepy.NURBS(
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
    multipatch_geometry = splinepy.Multipatch([bsp_el2, nur_el3])

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
        tmpf = to_tmpf(tmpd)
        splinepy.io.gismo.export(
            tmpf,
            multipatch=multipatch_geometry,
            indent=False,
            labeled_boundaries=False,
            additional_blocks=gismo_options,
        )
        (
            multipatch_geometry_loaded,
            gismo_options_loaded,
        ) = splinepy.io.gismo.load(tmpf, load_options=True)
        assert all(
            are_splines_equal(a, b)
            for a, b in zip(
                multipatch_geometry.patches,
                multipatch_geometry_loaded.patches,
            )
        )

        assert np.allclose(
            multipatch_geometry.interfaces,
            multipatch_geometry_loaded.interfaces,
        )

        assert gismo_options_loaded == gismo_options


def test_gismo_io_binary(to_tmpf, are_stripped_lines_same, are_splines_equal):
    """Test the base64 io-routines"""
    # We test this with just one (big, 3D) spline
    nurbs = splinepy.NURBS(
        degrees=[1, 1, 1],
        knot_vectors=[[0, 0, 1, 1], [0, 0, 1, 1], [0, 0, 1, 1]],
        control_points=np.ones((8, 3)),
        weights=np.ones((8, 1)),
    )
    nurbs.elevate_degrees([0, 1, 2, 2])
    np_rng = np.random.RandomState(19284918)
    for i in range(3):
        nurbs.insert_knots(i, np_rng.random(4))

    # Randomize points
    nurbs.cps = np_rng.random(nurbs.cps.shape)
    nurbs.weights = np_rng.random(nurbs.weights.shape)

    # Create a multipatch geometry
    multipatch_geometry = splinepy.Multipatch([nurbs])

    # Export
    with tempfile.TemporaryDirectory() as tmpd:
        tmpf = to_tmpf(tmpd)
        splinepy.io.gismo.export(
            tmpf,
            multipatch=multipatch_geometry,
            indent=False,
            labeled_boundaries=False,
            as_base64=True,
        )
        (
            multipatch_geometry_loaded,
            gismo_options_loaded,
        ) = splinepy.io.gismo.load(tmpf, load_options=True)

        with open(tmpf) as tmp_read, open(
            os.path.dirname(os.path.dirname(__file__))
            + "/data/gismo_noindent_nolabels_b64_3d.xml"
        ) as base_file:
            assert are_stripped_lines_same(
                base_file.readlines(), tmp_read.readlines(), True
            )

        assert all(
            are_splines_equal(a, b)
            for a, b in zip(
                multipatch_geometry.patches,
                multipatch_geometry_loaded.patches,
            )
        )

        assert np.allclose(
            multipatch_geometry.interfaces,
            multipatch_geometry_loaded.interfaces,
        )
