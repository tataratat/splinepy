import splinepy
import numpy as np
import unittest
try:
    from . import common as c
except:
    import common as c

class TestSplinepyInitialization(unittest.TestCase):

    def setUp(self):
        self.b2P2D = c.b2P2D.copy()
        self.n2P2D = c.n2P2D.copy()
        self.z2P2D = c.z2P2D.copy()
        self.r2P2D = c.r2P2D.copy()

    def test_constructor_raise_nonmatching_degree(self):
        """ Test the raise Error for mismatch of degree and par_dim. """

        self.b2P2D["degrees"] = [2,]
        self.n2P2D["degrees"] = [2,]
        self.z2P2D["degrees"] = [2,]
        self.r2P2D["degrees"] = [2,]

        with self.assertRaises(splinepy._spline.InputDimensionError) as context:
            splinepy.BSpline(**self.b2P2D)
        with self.assertRaises(splinepy._spline.InputDimensionError) as context:
            splinepy.NURBS(**self.n2P2D)
        with self.assertRaises(splinepy._spline.InputDimensionError) as context:
            splinepy.Bezier(**self.z2P2D)
        with self.assertRaises(splinepy._spline.InputDimensionError) as context:
            splinepy.RationalBezier(**self.r2P2D)

    def test_constructor_raise_nonmatching_kv_number(self):
        """ Test the raise Error for missing knot vectors. """

        self.b2P2D["knot_vectors"] = [[0, 0, 0, 0.5, 1, 1, 1]]
        self.n2P2D["knot_vectors"] = [[0, 0, 0, 1, 1, 1]]

        with self.assertRaises(splinepy._spline.InputDimensionError) as context:
            splinepy.BSpline(**self.b2P2D)
        with self.assertRaises(splinepy._spline.InputDimensionError) as context:
            splinepy.NURBS(**self.n2P2D)

    def test_constructor_raise_nonmatching_kv_length(self):
        """ Test the raise Error for length of knot vectors. """

        self.b2P2D["knot_vectors"] = [
            [0, 0, 0.5, 1, 1], 
            [0, 0, 1, 1]
        ]
        self.n2P2D["knot_vectors"] = [
            [0, 0, 1, 1], 
            [0, 0, 1, 1]
        ]

        with self.assertRaises(splinepy._spline.InputDimensionError) as context:
            splinepy.BSpline(**self.b2P2D)
        with self.assertRaises(splinepy._spline.InputDimensionError) as context:
            splinepy.NURBS(**self.n2P2D)

    def test_constructor_raise_nonmatching_cp_number(self):
        """ Test the raise Error for invalid number of control points. """

        self.b2P2D["control_points"] = [
            [0, 0],
            [0, 1],
            [1, 1.5],
            [3, 1.5],
            [-2, 0],
            [-2, 2],
            [1, 5],
            [3, 5],
        ]
        self.n2P2D["control_points"] = [
            [-1.,  0.],
            [ 0.,  1.],
            [-2.,  0.],
            [0.,  2.],
        ]
        self.n2P2D["weights"] = [
            [1.],
            [1.],
            [1.],
            [1.],
        ]
        self.z2P2D["control_points"] = [
            [-1.,  0.],
            [ 0.,  1.],
            [-2.,  0.],
            [-2.,  2.], 
            [0.,   2.]
        ]         
        self.r2P2D["control_points"] = [
            [-1.,  0.],
            [ 0.,  1.],
            [-2.,  0.],
            [-2.,  2.], 
            [0.,   2.]
        ]
        self.r2P2D["weights"] = [
            [1.],
            [1.],
            [1.],
            [(2** -.5)],
            [1.],
        ]
        
        with self.assertRaises(splinepy._spline.InputDimensionError) as context:
            splinepy.BSpline(**self.b2P2D)
        with self.assertRaises(splinepy._spline.InputDimensionError) as context:
            splinepy.NURBS(**self.n2P2D)
        with self.assertRaises(splinepy._spline.InputDimensionError) as context:
            splinepy.Bezier(**self.z2P2D)
        with self.assertRaises(splinepy._spline.InputDimensionError) as context:
            splinepy.RationalBezier(**self.r2P2D)

if __name__ == "__main__":
    unittest.main()