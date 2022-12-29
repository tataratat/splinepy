try:
    from . import common as c
except BaseException:
    import common as c


class TestSplinepyKnotVectorManipulation(c.unittest.TestCase):
    def setUp(self):
        self.b2P2D = c.b2P2D.copy()
        self.n2P2D = c.n2P2D.copy()
        self.bspline = c.splinepy.BSpline(**self.b2P2D)
        self.nurbs = c.splinepy.NURBS(**self.n2P2D)
        self.ref_bspline = c.splinepy.BSpline(**c.b2P2D)
        self.ref_nurbs = c.splinepy.NURBS(**c.n2P2D)

    def test_insert_knot(self):
        """Test the knot insertion function (.insert_knot())."""

        # reference solutions
        bspline_ref_kv = [
            [0, 0, 0, 0.2, 0.5, 0.7, 1, 1, 1],
            [0, 0, 0, 0.2, 0.5, 0.7, 1, 1, 1],
        ]
        nurbs_ref_kv = [
            [0, 0, 0, 0.2, 0.5, 0.7, 1, 1, 1],
            [0, 0, 0.2, 0.5, 0.7, 1, 1],
        ]

        # insert knots
        self.bspline.insert_knots(0, [0.2, 0.7])
        self.bspline.insert_knots(
            1,
            [
                0.2,
                0.5,
                0.7,
            ],
        )

        self.nurbs.insert_knots(
            0,
            [
                0.2,
                0.5,
                0.7,
            ],
        )
        self.nurbs.insert_knots(
            1,
            [
                0.2,
                0.5,
                0.7,
            ],
        )

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

    def test_remove_knot(self):
        """Test the function .remove_knots.
        This test also depends on the function .insert_knots!"""

        # insert and remove knots
        self.bspline.insert_knots(
            0,
            [
                0.2,
                0.7,
            ],
        )
        self.bspline.insert_knots(
            1,
            [
                0.2,
                0.5,
                0.7,
            ],
        )
        self.bspline.remove_knots(0, [0.2, 0.7])
        self.bspline.remove_knots(
            1,
            [
                0.2,
                0.5,
                0.7,
            ],
        )

        self.nurbs.insert_knots(
            0,
            [
                0.2,
                0.5,
                0.7,
            ],
        )
        self.nurbs.insert_knots(
            1,
            [
                0.2,
                0.5,
                0.7,
            ],
        )
        self.nurbs.remove_knots(
            0,
            [
                0.2,
                0.5,
                0.7,
            ],
        )
        self.nurbs.remove_knots(
            1,
            [
                0.2,
                0.5,
                0.7,
            ],
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
