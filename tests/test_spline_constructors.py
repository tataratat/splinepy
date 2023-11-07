try:
    from . import common as c
except BaseException:
    import common as c


class TestSplineConstructors(c.unittest.TestCase):
    """Perform checks and see if spline data is coherent"""

    def testBSplines(self):
        """
        Perform checks on BSpline class
        """
        # No assertion querz
        c.splinepy.BSpline(
            degrees=[1, 1],
            knot_vectors=[[0, 0, 1, 1], [0, 0, 0.4, 0.5, 0.7, 1, 1]],
            control_points=c.np.random.rand(10, 2),
        )
        # Check control points
        self.assertRaisesRegex(
            RuntimeError,
            "SPLINEPY ERROR - Invalid number of control points. 2 expected, "
            "but 9 were given.",
            c.splinepy.BSpline,
            degrees=[1],
            knot_vectors=[[0, 0, 1, 1]],
            control_points=c.np.ones((9, 3)),
        )
        self.assertRaisesRegex(
            RuntimeError,
            "SPLINEPY ERROR - Zero dimensional data points can not constitute "
            "a spline.",
            c.splinepy.BSpline,
            degrees=[1],
            knot_vectors=[[0, 0, 1, 1]],
            control_points=c.np.ones((2, 0)),
        )
        self.assertRaisesRegex(
            ValueError,
            r"len\(knot_vectors\) \(2\) should match len\(self.degrees\) "
            r"\(1\).",
            c.splinepy.BSpline,
            degrees=[1],
            knot_vectors=[[0, 0, 0.6, 0.3, 1, 1], [0, 0, 1, 1]],
            control_points=c.np.ones((8, 2)),
        )

        # Check degrees
        self.assertRaisesRegex(
            RuntimeError,
            "SPLINEPY ERROR - Invalid degree, degrees need to be positive. "
            "Detected degree -4 along parametric dimension: 0",
            c.splinepy.BSpline,
            degrees=[-4],
            knot_vectors=[[0, 0, 1, 1]],
            control_points=c.np.ones((1, 1)),
        )

    def testNURBS(self):
        """
        Perform checks on BSpline class
        """
        # No assertion querz
        c.splinepy.NURBS(
            degrees=[1, 1],
            knot_vectors=[[0, 0, 1, 1], [0, 0, 0.4, 0.5, 0.7, 1, 1]],
            control_points=c.np.random.rand(10, 2),
            weights=c.np.ones(10),
        )
        # Check control points
        self.assertRaisesRegex(
            ValueError,
            r"len\(weights\) \(2\) should match len\(control_points\) \(9\).",
            c.splinepy.NURBS,
            degrees=[1],
            knot_vectors=[[0, 0, 1, 1]],
            control_points=c.np.ones((9, 3)),
            weights=c.np.ones(2),
        )
        self.assertRaisesRegex(
            RuntimeError,
            "SPLINEPY ERROR - Zero dimensional data points can not constitute "
            "a spline.",
            c.splinepy.NURBS,
            degrees=[1],
            knot_vectors=[[0, 0, 1, 1]],
            control_points=c.np.ones((2, 0)),
            weights=c.np.ones(2),
        )
        self.assertRaisesRegex(
            RuntimeError,
            r"SPLINEPY ERROR - 0.6 0.3 Knots of parametric dimension \( 0 \) "
            "are not in increasing order.",
            c.splinepy.NURBS,
            degrees=[1, 1],
            knot_vectors=[[0, 0, 0.6, 0.3, 1, 1], [0, 0, 1, 1]],
            control_points=c.np.ones((8, 2)),
            weights=c.np.ones(8),
        )

        # Check degrees
        self.assertRaisesRegex(
            RuntimeError,
            "SPLINEPY ERROR - Invalid degree, degrees need to be positive. "
            "Detected degree -4 along parametric dimension: 0",
            c.splinepy.NURBS,
            degrees=[-4],
            knot_vectors=[[0, 0, 1, 1]],
            control_points=c.np.ones((1, 1)),
            weights=c.np.ones(1),
        )

        # Check weights
        self.assertRaisesRegex(
            ValueError,
            r"len\(weights\) \(19\) should match len\(control_points\) \(2\).",
            c.splinepy.NURBS,
            degrees=[1],
            knot_vectors=[[0, 0, 1, 1]],
            control_points=c.np.ones((2, 1)),
            weights=c.np.ones(19),
        )

    def testBezier(self):
        """Feed Bezier constructor with random data and see if the errors are
        meaningful"""
        # First create a reasonable spline
        c.splinepy.Bezier(degrees=[1, 1], control_points=c.np.ones((4, 1)))

        # Check assertions
        self.assertRaisesRegex(
            RuntimeError,
            "SPLINEPY ERROR - Invalid degree, degrees need to be positive. "
            "Detected degree -4 along parametric dimension: 0",
            c.splinepy.Bezier,
            degrees=[-4],
            control_points=c.np.ones((1, 1)),
        )
        self.assertRaisesRegex(
            RuntimeError,
            "SPLINEPY ERROR - Invalid number of control points. 9 expected, "
            "but 4 were given.",
            c.splinepy.Bezier,
            degrees=[2, 2],
            control_points=c.np.ones((4, 2)),
        )
        self.assertRaisesRegex(
            RuntimeError,
            "SPLINEPY ERROR - Invalid number of control points. 24 expected, "
            "but 9 were given.",
            c.splinepy.Bezier,
            degrees=[2, 3, 1],
            control_points=c.np.ones((9, 3)),
        )

    def testRationalBezier(self):
        """Feed Bezier constructor with random data and see if the errors are
        meaningful"""
        # First create a reasonable spline
        c.splinepy.RationalBezier(
            degrees=[1, 1],
            control_points=c.np.ones((4, 1)),
            weights=c.np.random.rand(4),
        )

        # Check assertions
        self.assertRaisesRegex(
            RuntimeError,
            "SPLINEPY ERROR - Invalid degree, degrees need to be positive. "
            "Detected degree -4 along parametric dimension: 0",
            c.splinepy.RationalBezier,
            degrees=[-4],
            control_points=c.np.ones((2, 4)),
            weights=c.np.random.rand(2),
        )
        self.assertRaisesRegex(
            RuntimeError,
            "SPLINEPY ERROR - Invalid number of control points. 9 expected, "
            "but 4 were given.",
            c.splinepy.RationalBezier,
            degrees=[2, 2],
            control_points=c.np.ones((4, 1)),
            weights=c.np.random.rand(4),
        )
        self.assertRaisesRegex(
            RuntimeError,
            "SPLINEPY ERROR - Invalid number of control points. 24 expected, "
            "but 9 were given.",
            c.splinepy.RationalBezier,
            degrees=[2, 3, 1],
            control_points=c.np.ones((9, 2)),
            weights=c.np.random.rand(9),
        )
        self.assertRaisesRegex(
            RuntimeError,
            "SPLINEPY ERROR - Zero dimensional data points can not constitute "
            "a spline.",
            c.splinepy.RationalBezier,
            degrees=[2, 3, 1],
            control_points=c.np.ones((24, 0)),
            weights=c.np.random.rand(24),
        )
        self.assertRaisesRegex(
            ValueError,
            r"len\(weights\) \(9\) should match len\(control_points\) \(24\).",
            c.splinepy.RationalBezier,
            degrees=[2, 3, 1],
            control_points=c.np.ones((24, 2)),
            weights=c.np.random.rand(9),
        )


if __name__ == "__main__":
    c.unittest.main()
