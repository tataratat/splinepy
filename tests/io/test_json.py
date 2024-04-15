import tempfile

import splinepy


def test_json_import(to_tmpf, are_splines_equal):
    """
    Test json export-import routine
    """
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

    # Test Output against input using base64 encoding
    list_of_splines = [bez_el0, rbz_el1, bsp_el2, nur_el3]
    with tempfile.TemporaryDirectory() as tmpd:
        tmpf = to_tmpf(tmpd)
        splinepy.io.json.export(tmpf, list_of_splines, base64encoding=True)
        list_of_splines_loaded = splinepy.io.json.load(tmpf)
        assert all(
            are_splines_equal(a, b)
            for a, b in zip(list_of_splines, list_of_splines_loaded)
        )

    # Test Import export with non base64 encoding
    with tempfile.TemporaryDirectory() as tmpd:
        tmpf = to_tmpf(tmpd)
        splinepy.io.json.export(tmpf, list_of_splines, base64encoding=False)
        list_of_splines_loaded = splinepy.io.json.load(tmpf)
        assert all(
            are_splines_equal(a, b)
            for a, b in zip(list_of_splines, list_of_splines_loaded)
        )
