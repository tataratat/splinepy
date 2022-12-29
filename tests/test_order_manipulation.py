try:
    from . import common as c
except BaseException:
    import common as c


class TestSplinepyOrderManipulation(c.unittest.TestCase):
    def setUp(self):
        self.b2P2D = c.b2P2D.copy()
        self.n2P2D = c.n2P2D.copy()
        self.z2P2D = c.z2P2D.copy()
        self.r2P2D = c.r2P2D.copy()
        self.bspline = c.splinepy.BSpline(**self.b2P2D)
        self.nurbs = c.splinepy.NURBS(**self.n2P2D)
        self.bezier = c.splinepy.Bezier(**self.z2P2D)
        self.rational = c.splinepy.RationalBezier(**self.r2P2D)
        self.ref_bspline = c.splinepy.BSpline(**c.b2P2D)
        self.ref_nurbs = c.splinepy.NURBS(**c.n2P2D)
        self.ref_bezier = c.splinepy.Bezier(**c.z2P2D)
        self.ref_rational = c.splinepy.RationalBezier(**c.r2P2D)

    def test_elevate_degree(self):
        """Test the order elevation function (.elevate_degree())."""

        # reference solution
        bspline_ref_kv = [
            [0, 0, 0, 0, 0.5, 0.5, 1, 1, 1, 1],
            [0, 0, 0, 1, 1, 1],
        ]
        nurbs_ref_kv = [
            [0, 0, 0, 0, 1, 1, 1, 1],
            [0, 0, 1, 1],
        ]

        # elevate order
        self.bspline.elevate_degrees(0)

        self.nurbs.elevate_degrees(0)

        self.bezier.elevate_degrees(0)

        self.rational.elevate_degrees(0)

        # test degrees
        self.assertTrue(c.np.allclose(self.bspline.degrees, [3, 2]))
        self.assertTrue(c.np.allclose(self.nurbs.degrees, [3, 1]))
        self.assertTrue(c.np.allclose(self.bezier.degrees, [3, 1]))
        self.assertTrue(c.np.allclose(self.rational.degrees, [3, 1]))

        # test knot_vectors
        self.assertTrue(
            c.are_items_close(self.bspline.knot_vectors, bspline_ref_kv)
        )
        self.assertTrue(
            c.are_items_close(self.nurbs.knot_vectors, nurbs_ref_kv)
        )

        # use random query points
        q2D = c.np.random.rand(10, 2)

        # test evaluation
        self.assertTrue(
            c.np.allclose(
                self.bspline.evaluate(q2D), self.ref_bspline.evaluate(q2D)
            )
        )
        self.assertTrue(
            c.np.allclose(
                self.nurbs.evaluate(q2D), self.ref_nurbs.evaluate(q2D)
            )
        )
        self.assertTrue(
            c.np.allclose(
                self.bezier.evaluate(q2D), self.ref_bezier.evaluate(q2D)
            )
        )
        self.assertTrue(
            c.np.allclose(
                self.rational.evaluate(q2D), self.ref_rational.evaluate(q2D)
            )
        )

    def test_reduce_degree(self):
        """Test the function .reduce_degree.
        This test also depends on the function .elevate_degree!"""

        # elevate and reduce order
        self.bspline.elevate_degrees(0)
        self.bspline.reduce_degrees(0)

        self.nurbs.elevate_degrees(0)
        self.nurbs.reduce_degrees(0)

        # test degrees
        self.assertTrue(
            c.np.allclose(self.bspline.degrees, self.ref_bspline.degrees)
        )
        self.assertTrue(
            c.np.allclose(self.nurbs.degrees, self.ref_nurbs.degrees)
        )

        # test knot_vectors
        self.assertTrue(
            c.are_items_close(
                self.bspline.knot_vectors, self.ref_bspline.knot_vectors
            )
        )
        self.assertTrue(
            c.are_items_close(
                self.nurbs.knot_vectors, self.ref_nurbs.knot_vectors
            )
        )

        # use random query points
        q2D = c.np.random.rand(10, 2)

        # test evaluation
        self.assertTrue(
            c.np.allclose(
                self.bspline.evaluate(q2D), self.ref_bspline.evaluate(q2D)
            )
        )
        self.assertTrue(
            c.np.allclose(
                self.nurbs.evaluate(q2D), self.ref_nurbs.evaluate(q2D)
            )
        )


if __name__ == "__main__":
    c.unittest.main()
