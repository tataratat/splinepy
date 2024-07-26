import os
import tempfile

import splinepy


def test_mfem_export(to_tmpf, are_stripped_lines_same):
    """
    Test MFEM export routine
    """
    # Define some splines
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

    # Test Output
    with tempfile.TemporaryDirectory() as tmpd:
        tmpf = to_tmpf(tmpd)
        splinepy.io.mfem.export_cartesian(
            tmpf, [bez_el0, bsp_el2, nur_el3, rbz_el1]
        )

        with open(tmpf) as tmp_read, open(
            os.path.dirname(__file__) + "/../data/mfem_cartesian_2d.mesh"
        ) as base_file:
            assert are_stripped_lines_same(
                base_file.readlines(), tmp_read.readlines(), True
            )

    # Test Also 3D Meshes
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

    # Test output
    with tempfile.TemporaryDirectory() as tmpd:
        tmpf = to_tmpf(tmpd)
        splinepy.io.mfem.export_cartesian(
            tmpf, [bez_el0, bsp_el2, nur_el3, rbz_el1]
        )

        with open(tmpf) as tmp_read, open(
            os.path.dirname(__file__) + "/../data/mfem_cartesian_3d.mesh"
        ) as base_file:
            assert are_stripped_lines_same(
                base_file.readlines(), tmp_read.readlines(), True
            )


def test_mfem_single_patch_export(to_tmpf, are_splines_equal, np_rng):
    """
    Test MFEM single patch export routine
    """
    # create 2d/3d splines of all supporting types
    box = splinepy.helpme.create.box(1, 2)
    boxes = []
    boxes.extend([box, box.rationalbezier, box.bspline, box.nurbs])
    boxes.extend([b.create.extruded([0, 0, 3]) for b in boxes[:4]])
    for spline in boxes:
        dim = spline.dim

        # make it a bit complex
        for i in range(dim):
            if not spline.has_knot_vectors:
                break
            spline.insert_knots(i, np_rng.random([5]))

        spline.elevate_degrees(list(range(dim)) * 2)

        for i in range(dim):
            if not spline.has_knot_vectors:
                break
            spline.insert_knots(i, np_rng.random([5]))

        # Test Output
        with tempfile.TemporaryDirectory() as tmpd:
            tmpf = to_tmpf(tmpd)
            splinepy.io.mfem.export(tmpf, spline)

            loaded_spline = splinepy.io.mfem.load(tmpf)

            # all mfem export is nurbs
            assert are_splines_equal(spline.nurbs, loaded_spline)
