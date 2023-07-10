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
        self.assertRaisesRegex(
            RuntimeError,
            "SPLINEPY ERROR - Invalid number of control points. 2 exepcted, "
            "but 9 were given.",
            c.splinepy.BSpline,
            degrees=[1],
            knot_vectors=[[0, 0, 1, 1]],
            control_points=c.np.ones((9, 3)),
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
            "Detected degree -4 along parametric dimension: 1",
            c.splinepy.Bezier,
            degrees=[-4],
            control_points=c.np.ones((1, 1)),
        )
        self.assertRaisesRegex(
            RuntimeError,
            "SPLINEPY ERROR - Invalid number of control points. 9 exepcted, "
            "but 4 were given.",
            c.splinepy.Bezier,
            degrees=[2, 2],
            control_points=c.np.ones((4, 2)),
        )
        self.assertRaisesRegex(
            RuntimeError,
            "SPLINEPY ERROR - Invalid number of control points. 24 exepcted, "
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
            weights=c.np.random.rand(4)
            )

        # Check assertions
        self.assertRaisesRegex(
            RuntimeError,
            "SPLINEPY ERROR - Invalid degree, degrees need to be positive. "
            "Detected degree -4 along parametric dimension: 1",
            c.splinepy.RationalBezier,
            degrees=[-4],
            control_points=c.np.ones((2, 4)),
            weights=c.np.random.rand(0)
        )
        self.assertRaisesRegex(
            RuntimeError,
            "SPLINEPY ERROR - Invalid number of control points. 9 exepcted, "
            "but 4 were given.",
            c.splinepy.RationalBezier,
            degrees=[2, 2],
            control_points=c.np.ones((4, 1)),
            weights=c.np.random.rand(0)
        )
        self.assertRaisesRegex(
            RuntimeError,
            "SPLINEPY ERROR - Invalid number of control points. 24 exepcted, "
            "but 9 were given.",
            c.splinepy.RationalBezier,
            degrees=[2, 3, 1],
            control_points=c.np.ones((9, 2)),
            weights=c.np.random.rand(0)
        )
        self.assertRaisesRegex(
            RuntimeError,
            "SPLINEPY ERROR - Zero dimensional data points can not constitute "
            "a spline.",
            c.splinepy.RationalBezier,
            degrees=[2, 3, 1],
            control_points=c.np.ones((24, 0)),
            weights=c.np.random.rand(9)
        )
        self.assertRaisesRegex(
            RuntimeError,
            "SPLINEPY ERROR - Number of weights \( 9 \) does not match number "
            "of control points \( 24 \).",
            c.splinepy.RationalBezier,
            degrees=[2, 3, 1],
            control_points=c.np.ones((24, 2)),
            weights=c.np.random.rand(9)
        )


if __name__ == "__main__":
    c.unittest.main()
