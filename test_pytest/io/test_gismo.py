import os
import tempfile
from sys import version as python_version

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
