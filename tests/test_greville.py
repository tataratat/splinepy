import numpy as np

import splinepy


def test_greville_points_bezier(bezier_3p3d):
    """
    test greville point calculation for bezier
    """
    greville_points = splinepy.utils.data.cartesian_product(
        [
            np.linspace(0, 1, bezier_3p3d.degrees[i] + 1)
            for i in range(bezier_3p3d.para_dim)
        ]
    )

    assert np.allclose(greville_points, bezier_3p3d.greville_abscissae())


def test_greville_points_rational_bezier(rational_bezier_3p3d):
    """
    test greville point calculation for rational bezier
    """
    greville_points = splinepy.utils.data.cartesian_product(
        [
            np.linspace(0, 1, rational_bezier_3p3d.degrees[i] + 1)
            for i in range(rational_bezier_3p3d.para_dim)
        ]
    )
    assert np.allclose(
        greville_points, rational_bezier_3p3d.greville_abscissae()
    )


def test_greville_points_bspline(bspline_3p3d):
    """
    test greville point calculation for bspline
    """
    # Modify some splines
    bspline_3p3d.insert_knots(0, [0.2, 0.7])
    bspline_3p3d.insert_knots(1, [0.2, 0.5, 0.7])
    bspline_3p3d.insert_knots(2, [0.5])

    # Test non-unifrom Splines
    greville_points = splinepy.utils.data.cartesian_product(
        [
            np.convolve(
                bspline_3p3d.knot_vectors[i][1:-1],
                np.ones(bspline_3p3d.degrees[i]),
            )
            for i in range(bspline_3p3d.para_dim)
        ]
    )
    assert np.allclose(greville_points, bspline_3p3d.greville_abscissae())


def test_greville_points_nurbs(nurbs_3p3d):
    """
    test greville point calculation for nurbs
    """
    nurbs_3p3d.elevate_degrees([0, 0, 1, 2])
    nurbs_3p3d.insert_knots(0, [0.25])
    # Test Uniform Types
    greville_points = splinepy.utils.data.cartesian_product(
        [
            np.convolve(
                nurbs_3p3d.knot_vectors[i][1:-1],
                np.ones(nurbs_3p3d.degrees[i]) / nurbs_3p3d.degrees[i],
                mode="valid",
            )
            for i in range(nurbs_3p3d.para_dim)
        ]
    )
    assert np.allclose(greville_points, nurbs_3p3d.greville_abscissae())


def test_greville_with_duplicate_points():
    """Tests if duplicate points are filtered out, as is required for the
    construction of c^(-1) spline"""
    # Construct c^(-1) spline
    a = splinepy.BSpline(
        degrees=[2],
        knot_vectors=[[0, 0, 0, 0.25, 0.5, 0.5, 0.5, 0.75, 1, 1, 1]],
        control_points=np.ones((8, 2)),
    )

    assert np.allclose(
        a.greville_abscissae(duplicate_tolerance=-1.0).ravel(),
        [0, 0.125, 0.375, 0.5, 0.5, 0.625, 0.875, 1.0],
    )

    assert np.allclose(
        a.greville_abscissae(
            duplicate_tolerance=splinepy.settings.TOLERANCE
        ).ravel(),
        [0, 0.125, 0.375, 0.4375, 0.5625, 0.625, 0.875, 1.0],
    )
