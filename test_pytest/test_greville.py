import numpy as np

import splinepy


def test_greville_points(
    bezier_3p3d, rational_bezier_3p3d, bspline_3p3d, nurbs_3p3d
):
    """
    test permute
    """
    # Define some splines
    z = bezier_3p3d
    r = rational_bezier_3p3d
    b = bspline_3p3d
    n = nurbs_3p3d

    # Modify some splines
    b.insert_knots(0, [0.2, 0.7])
    b.insert_knots(1, [0.2, 0.5, 0.7])
    b.insert_knots(2, [0.5])
    n.elevate_degrees([0, 0, 1, 2])
    n.insert_knots(0, [0.25])
    # Test Uniform Types
    greville_points = splinepy.utils.data.cartesian_product(
        [np.linspace(0, 1, z.degrees[i] + 1) for i in range(z.para_dim)]
    )

    assert np.allclose(greville_points, z.greville_abscissae())

    greville_points = splinepy.utils.data.cartesian_product(
        [np.linspace(0, 1, r.degrees[i] + 1) for i in range(r.para_dim)]
    )
    assert np.allclose(greville_points, r.greville_abscissae())

    # Test non-unifrom Splines
    greville_points = splinepy.utils.data.cartesian_product(
        [
            np.convolve(b.knot_vectors[i][1:-1], np.ones(b.degrees[i]))
            for i in range(b.para_dim)
        ]
    )
    assert np.allclose(greville_points, b.greville_abscissae())

    greville_points = splinepy.utils.data.cartesian_product(
        [
            np.convolve(
                n.knot_vectors[i][1:-1],
                np.ones(n.degrees[i]) / n.degrees[i],
                mode="valid",
            )
            for i in range(n.para_dim)
        ]
    )
    assert np.allclose(greville_points, n.greville_abscissae())


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
