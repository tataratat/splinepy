import splinepy
import numpy as np
import unittest
try:
    from . import common as c
except:
    import common as c



class TestSplinepyKnotVectorManipulation(unittest.TestCase):

    def setUp(self):
        self.b2P2D = c.b2P2D.copy()
        self.n2P2D = c.n2P2D.copy()
        self.z2P2D = c.z2P2D.copy()
        self.r2P2D = c.r2P2D.copy()
        self.bspline = splinepy.BSpline(**self.b2P2D)
        self.nurbs = splinepy.NURBS(**self.n2P2D)
        self.bezier = splinepy.Bezier(**self.z2P2D)
        self.rational = splinepy.RationalBezier(**self.r2P2D)
        self.ref_bspline = splinepy.BSpline(**c.b2P2D)
        self.ref_nurbs = splinepy.NURBS(**c.n2P2D)

    def test_insert_knot(self):
        """ Test the knot insertion function (.insert_knot()). """ 

        # reference solutions
        bspline_ref_kv = [
            [0, 0, 0, .2, 0.5, .7, 1, 1, 1],
            [0, 0, 0, .2, 0.5, .7, 1, 1, 1],
        ]
        nurbs_ref_kv=[
            [0, 0, 0, .2, .5, .7, 1, 1, 1],
            [0, 0, .2, .5, .7, 1, 1],
        ]

        # insert knots
        self.bspline.insert_knots(0,[.2, .7])
        self.bspline.insert_knots(1,[.2, .5, .7,])

        self.nurbs.insert_knots(0,[.2, .5, .7,])
        self.nurbs.insert_knots(1,[.2, .5, .7,])

        # test knot_vectors
        self.assertEqual(self.bspline.knot_vectors, bspline_ref_kv)
        self.assertEqual(self.nurbs.knot_vectors, nurbs_ref_kv)

        # test evaluation
        self.assertTrue(np.allclose(
            self.bspline.evaluate(c.q2D), 
            self.ref_bspline.evaluate(c.q2D)
            )
        )
        self.assertTrue(np.allclose(
            self.nurbs.evaluate(c.q2D), 
            self.ref_nurbs.evaluate(c.q2D)
            )
        )

    def test_remove_knot(self):
        """ Test the function .remove_knots. 
        This test also depends on the function .insert_knots! """

        # reference solutions

        # insert and remove knots
        self.bspline.insert_knots(0,[.2, .7,])
        self.bspline.insert_knots(1,[.2, .5, .7,])
        self.bspline.remove_knots(0,[.2, .7])
        self.bspline.remove_knots(1,[.2, .5, .7,])

        self.nurbs.insert_knots(0,[.2, .5, .7,])
        self.nurbs.insert_knots(1,[.2, .5, .7,])
        self.nurbs.remove_knots(0,[.2, .5, .7,])
        self.nurbs.remove_knots(1,[.2, .5, .7,])

        # test knot_vectors
        self.assertEqual(
            self.bspline.knot_vectors, 
            self.ref_bspline.knot_vectors
        )
        self.assertEqual(
            self.nurbs.knot_vectors, 
            self.ref_nurbs.knot_vectors
            )

        # test evaluation
        self.assertTrue(np.allclose(
            self.bspline.evaluate(c.q2D), 
            self.ref_bspline.evaluate(c.q2D)
            )
        )
        self.assertTrue(np.allclose(
            self.nurbs.evaluate(c.q2D), 
            self.ref_nurbs.evaluate(c.q2D)
            )
        )

if __name__ == "__main__":
    unittest.main()