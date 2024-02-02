try:
    from . import common as c
except BaseException:
    import common as c


class TestSplinepyOrderManipulation(c.SplineBasedTestCase):
    def test_elevate_degree(self):
        """Test the order elevation function (.elevate_degrees())."""

        # reference solution
        bspline_ref_kv = [
            [0, 0, 0, 0, 0.5, 0.5, 1, 1, 1, 1],
            [0, 0, 0, 1, 1, 1],
        ]
        nurbs_ref_kv = [
            [0, 0, 0, 0, 1, 1, 1, 1],
            [0, 0, 1, 1],
        ]

        bspline = self.bspline_2p2d()
        nurbs = self.nurbs_2p2d()
        bezier = self.bezier_2p2d()
        rational = self.rational_bezier_2p2d()

        # elevate order
        bspline.elevate_degrees(0)

        nurbs.elevate_degrees(0)

        bezier.elevate_degrees(0)

        rational.elevate_degrees(0)

        # test degrees
        self.assertTrue(c.np.allclose(bspline.degrees, [3, 2]))
        self.assertTrue(c.np.allclose(nurbs.degrees, [3, 1]))
        self.assertTrue(c.np.allclose(bezier.degrees, [3, 1]))
        self.assertTrue(c.np.allclose(rational.degrees, [3, 1]))

        # test knot_vectors
        self.assertTrue(
            c.are_items_close(bspline.knot_vectors, bspline_ref_kv)
        )
        self.assertTrue(c.are_items_close(nurbs.knot_vectors, nurbs_ref_kv))

        # use random query points
        q2D = c.np.random.rand(10, 2)

        # test evaluation
        self.assertTrue(
            c.np.allclose(
                bspline.evaluate(q2D), self.bspline_2p2d().evaluate(q2D)
            )
        )
        self.assertTrue(
            c.np.allclose(nurbs.evaluate(q2D), self.nurbs_2p2d().evaluate(q2D))
        )
        self.assertTrue(
            c.np.allclose(
                bezier.evaluate(q2D), self.bezier_2p2d().evaluate(q2D)
            )
        )
        self.assertTrue(
            c.np.allclose(
                rational.evaluate(q2D),
                self.rational_bezier_2p2d().evaluate(q2D),
            )
        )

    def test_reduce_degree(self):
        """Test the function .reduce_degrees.
        This test also depends on the function .elevate_degrees!"""

        # elevate and reduce order
        bspline = self.bspline_2p2d()
        nurbs = self.nurbs_2p2d()

        bspline.elevate_degrees(0)
        bspline.reduce_degrees(0)

        nurbs.elevate_degrees(0)
        nurbs.reduce_degrees(0)

        # test degrees
        self.assertTrue(
            c.np.allclose(bspline.degrees, self.bspline_2p2d().degrees)
        )
        self.assertTrue(
            c.np.allclose(nurbs.degrees, self.nurbs_2p2d().degrees)
        )

        # test knot_vectors
        self.assertTrue(
            c.are_items_close(
                bspline.knot_vectors, self.bspline_2p2d().knot_vectors
            )
        )
        self.assertTrue(
            c.are_items_close(
                nurbs.knot_vectors, self.nurbs_2p2d().knot_vectors
            )
        )

        # use random query points
        q2D = c.np.random.rand(10, 2)

        # test evaluation
        self.assertTrue(
            c.np.allclose(
                bspline.evaluate(q2D), self.bspline_2p2d().evaluate(q2D)
            )
        )
        self.assertTrue(
            c.np.allclose(nurbs.evaluate(q2D), self.nurbs_2p2d().evaluate(q2D))
        )


if __name__ == "__main__":
    c.unittest.main()
