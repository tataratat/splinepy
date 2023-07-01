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
            "but 0 were given.",
            c.splinepy.BSpline,
            degrees=[1],
            knot_vectors=[[0, 0, 1, 1]],
            control_points=c.np.ones((0, 0)),
        )


if __name__ == "__main__":
    c.unittest.main()
