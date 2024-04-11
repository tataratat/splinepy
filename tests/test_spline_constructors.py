import numpy as np
import pytest

import splinepy


def test_BSplines(np_rng):
    """
    Perform checks on BSpline class
    """
    # No assertion query
    splinepy.BSpline(
        degrees=[1, 1],
        knot_vectors=[[0, 0, 1, 1], [0, 0, 0.4, 0.5, 0.7, 1, 1]],
        control_points=np_rng.random((10, 2)),
    )
    # Check control points
    with pytest.raises(
        RuntimeError,
        match=r"SPLINEPY ERROR - Invalid number of control points. "
        r"2 expected, but 9 were given.",
    ):
        splinepy.BSpline(
            degrees=[1],
            knot_vectors=[[0, 0, 1, 1]],
            control_points=np.ones((9, 3)),
        )

    with pytest.raises(
        RuntimeError,
        match=r"SPLINEPY ERROR - Zero dimensional data points "
        r"can not constitute a spline.",
    ):
        splinepy.BSpline(
            degrees=[1],
            knot_vectors=[[0, 0, 1, 1]],
            control_points=np.ones((2, 0)),
        )

    with pytest.raises(
        ValueError,
        match=r"len\(knot_vectors\) \(2\) should match len\(self.degrees\) \(1\).",
    ):
        splinepy.BSpline(
            degrees=[1],
            knot_vectors=[[0, 0, 0.6, 0.3, 1, 1], [0, 0, 1, 1]],
            control_points=np.ones((8, 2)),
        )

    # Check degrees
    with pytest.raises(
        RuntimeError,
        match=r"SPLINEPY ERROR - Invalid degree, degrees need to be positive. "
        r"Detected degree -4 along parametric dimension: 0",
    ):
        splinepy.BSpline(
            degrees=[-4],
            knot_vectors=[[0, 0, 1, 1]],
            control_points=np.ones((1, 1)),
        )


def test_NURBS(np_rng):
    """
    Perform checks on Nurbs class
    """
    # No assertion querz
    splinepy.NURBS(
        degrees=[1, 1],
        knot_vectors=[[0, 0, 1, 1], [0, 0, 0.4, 0.5, 0.7, 1, 1]],
        control_points=np_rng.random((10, 2)),
        weights=np.ones(10),
    )
    # Check control points
    with pytest.raises(
        ValueError,
        match=r"len\(weights\) \(2\) should match len\(control_points\) \(9\).",
    ):
        splinepy.NURBS(
            degrees=[1],
            knot_vectors=[[0, 0, 1, 1]],
            control_points=np.ones((9, 3)),
            weights=np.ones(2),
        )

    with pytest.raises(
        RuntimeError,
        match=r"SPLINEPY ERROR - Zero dimensional data points "
        r"can not constitute a spline.",
    ):
        splinepy.NURBS(
            degrees=[1],
            knot_vectors=[[0, 0, 1, 1]],
            control_points=np.ones((2, 0)),
            weights=np.ones(2),
        )

    with pytest.raises(
        RuntimeError,
        match=r"SPLINEPY ERROR - 0.6 0.3 Knots of parametric dimension \( 0 \) "
        r"are not in increasing order.",
    ):
        splinepy.NURBS(
            degrees=[1, 1],
            knot_vectors=[[0, 0, 0.6, 0.3, 1, 1], [0, 0, 1, 1]],
            control_points=np.ones((8, 2)),
            weights=np.ones(8),
        )

    # Check degrees
    with pytest.raises(
        RuntimeError,
        match=r"SPLINEPY ERROR - Invalid degree, degrees need to be positive. "
        r"Detected degree -4 along parametric dimension: 0",
    ):
        splinepy.NURBS(
            degrees=[-4],
            knot_vectors=[[0, 0, 1, 1]],
            control_points=np.ones((1, 1)),
            weights=np.ones(1),
        )

    # Check weights
    with pytest.raises(
        ValueError,
        match=r"len\(weights\) \(19\) should match "
        r"len\(control_points\) \(2\).",
    ):
        splinepy.NURBS(
            degrees=[1],
            knot_vectors=[[0, 0, 1, 1]],
            control_points=np.ones((2, 1)),
            weights=np.ones(19),
        )


def test_Bezier():
    """Feed Bezier constructor with random data and see if the errors are
    meaningful"""
    # First create a reasonable spline
    splinepy.Bezier(degrees=[1, 1], control_points=np.ones((4, 1)))

    # Check assertions
    with pytest.raises(
        RuntimeError,
        match=r"SPLINEPY ERROR - Invalid degree, degrees need to be positive. "
        r"Detected degree -4 along parametric dimension: 0",
    ):
        splinepy.Bezier(
            degrees=[-4],
            control_points=np.ones((1, 1)),
        )

    with pytest.raises(
        RuntimeError,
        match=r"SPLINEPY ERROR - Invalid number of control points. "
        r"9 expected, but 4 were given.",
    ):
        splinepy.Bezier(
            degrees=[2, 2],
            control_points=np.ones((4, 2)),
        )

    with pytest.raises(
        RuntimeError,
        match=r"SPLINEPY ERROR - Invalid number of control points. "
        r"24 expected, but 9 were given.",
    ):
        splinepy.Bezier(
            degrees=[2, 3, 1],
            control_points=np.ones((9, 3)),
        )


def test_RationalBezier(np_rng):
    """Feed Bezier constructor with random data and see if the errors are
    meaningful"""
    # First create a reasonable spline
    splinepy.RationalBezier(
        degrees=[1, 1],
        control_points=np.ones((4, 1)),
        weights=np_rng.random(4),
    )

    # Check assertions
    with pytest.raises(
        RuntimeError,
        match=r"SPLINEPY ERROR - Invalid degree, degrees need to be positive. "
        r"Detected degree -4 along parametric dimension: 0",
    ):
        splinepy.RationalBezier(
            degrees=[-4],
            control_points=np.ones((2, 4)),
            weights=np_rng.random(2),
        )

    with pytest.raises(
        RuntimeError,
        match=r"SPLINEPY ERROR - Invalid number of control points. "
        r"9 expected, but 4 were given.",
    ):
        splinepy.RationalBezier(
            degrees=[2, 2],
            control_points=np.ones((4, 1)),
            weights=np_rng.random(4),
        )

    with pytest.raises(
        RuntimeError,
        match=r"SPLINEPY ERROR - Invalid number of control points. "
        r"24 expected, but 9 were given.",
    ):
        splinepy.RationalBezier(
            degrees=[2, 3, 1],
            control_points=np.ones((9, 2)),
            weights=np_rng.random(9),
        )

    with pytest.raises(
        RuntimeError,
        match=r"SPLINEPY ERROR - Zero dimensional data "
        r"points can not constitute a spline.",
    ):
        splinepy.RationalBezier(
            degrees=[2, 3, 1],
            control_points=np.ones((24, 0)),
            weights=np_rng.random(24),
        )

    with pytest.raises(
        ValueError,
        match=r"len\(weights\) \(9\) should match "
        r"len\(control_points\) \(24\).",
    ):
        splinepy.RationalBezier(
            degrees=[2, 3, 1],
            control_points=np.ones((24, 2)),
            weights=np_rng.random(9),
        )
