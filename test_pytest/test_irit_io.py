import tempfile

import splinepy


def test_export_import(
    bspline_2p2d, nurbs_2p2d, bezier_3p3d, to_tmpf, are_splines_equal
):
    """Test import export routine for irit files"""
    # Create some splines
    bspline = bspline_2p2d
    bezier = bezier_3p3d
    nurbs = nurbs_2p2d
    # Create a line
    rational = nurbs_2p2d.extract.spline(splitting_plane=0, interval=0.3)

    # Make it more tricky
    bspline.elevate_degrees([0, 1, 0])
    nurbs.elevate_degrees([1])
    bezier.elevate_degrees([0, 1, 2])
    nurbs.insert_knots(1, [0.5])
    bspline.insert_knots(1, [0.5])

    list_of_splines = [bezier, nurbs, bspline, rational]

    with tempfile.TemporaryDirectory() as tmpd:
        tmpf = to_tmpf(tmpd)
        splinepy.io.irit.export(tmpf, list_of_splines)
        list_of_splines_loaded = splinepy.io.irit.load(tmpf)
        assert all(
            are_splines_equal(a, b)
            for a, b in zip(list_of_splines, list_of_splines_loaded)
        )
