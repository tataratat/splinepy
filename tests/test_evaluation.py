import splinepy
import numpy as np
import unittest
try:
    from . import common as c
except:
    import common as c


class TestSplinepyEvaluation(unittest.TestCase):

    def setUp(self):
        self.b2P2D = c.b2P2D.copy()
        self.n2P2D = c.n2P2D.copy()
        self.z2P2D = c.z2P2D.copy()
        self.r2P2D = c.r2P2D.copy()
        self.bspline = splinepy.BSpline(**self.b2P2D)
        self.nurbs = splinepy.NURBS(**self.n2P2D)
        self.bezier = splinepy.Bezier(**self.z2P2D)
        self.rational = splinepy.RationalBezier(**self.r2P2D)

    def test_basis_functions(self):
        """ Test the correct calculation of the basis functions. 
        (.basis_functions()) """

        # reference solutions
        bspline_ref_basis_functions = np.array(
            [[9.4128804e-01, 3.8615940e-02, 1.9602000e-04, 1.9015920e-02, 
              7.8012000e-04, 3.9600000e-06, 9.6040000e-05, 3.9400000e-06, 
              2.0000000e-08],
             [2.4010000e-01, 9.8500000e-03, 5.0000000e-05, 4.8020000e-01, 
              1.9700000e-02, 1.0000000e-04, 2.4010000e-01, 9.8500000e-03, 
              5.0000000e-05],
             [1.6200000e-02, 2.7540000e-01, 5.1840000e-01, 3.6000000e-03, 
              6.1200000e-02, 1.1520000e-01, 2.0000000e-04, 3.4000000e-03, 
              6.4000000e-03],
             [7.2000000e-03, 5.0400000e-02, 3.2400000e-02, 3.3600000e-02, 
              2.3520000e-01, 1.5120000e-01, 3.9200000e-02, 2.7440000e-01, 
              1.7640000e-01],
             [4.0000000e-06, 6.4000000e-05, 3.2000000e-05, 7.9200000e-04, 
              1.2672000e-02, 6.3360000e-03, 3.9204000e-02, 6.2726400e-01, 
              3.1363200e-01]] 
        )
        nurbs_ref_basis_functions = np.array(
            [[9.75958864e-01, 1.39415582e-02, 9.95774782e-05, 9.85817035e-03,
                1.40823820e-04, 1.00583311e-06],
              [4.92908517e-01, 7.04119101e-03, 5.02916557e-05, 4.92908517e-01,
                7.04119101e-03, 5.02916557e-05],
              [9.50089457e-03, 1.20926646e-01, 7.69572460e-01, 1.05565495e-03,
                1.34362940e-02, 8.55080511e-02],
              [1.32410262e-02, 7.49025551e-02, 2.11856419e-01, 3.08957277e-02,
                1.74772629e-01, 4.94331644e-01],
              [4.18891419e-03, 3.94934617e-03, 1.86173964e-03, 4.14702505e-01,
                3.90985271e-01, 1.84312224e-01]]
        )

        # test basis functions
        self.assertTrue(np.allclose(
            self.bspline.basis_functions(c.q2D)[0], 
            bspline_ref_basis_functions
            )
        )
        self.assertTrue(np.allclose(
            self.nurbs.basis_functions(c.q2D)[0], 
            nurbs_ref_basis_functions
            )
        )
    
    def test_partition_of_unity(self):
        """ Test the partition of unity of the calculated basis functions. """

        def basis_functions_sum(basis_functions):
            bf = np.array(basis_functions[0])
            u = np.zeros(np.shape(bf)[0])
            for i in range(0,np.shape(bf)[0]):
                u[i] = np.sum(bf[i])
            return u

        # use random query points
        q2D = np.random.rand(10,2)

        u_bspline = basis_functions_sum(self.bspline.basis_functions(q2D))
        self.assertTrue(np.allclose(
            u_bspline, 
            np.ones(np.shape(u_bspline))
            )
        )
        u_nurbs = basis_functions_sum(self.nurbs.basis_functions(q2D))
        self.assertTrue(np.allclose(
            u_nurbs, 
            np.ones(np.shape(u_nurbs))
            )
        )
    
    def test_evaluate(self):
        """ Test the correct spline evaluation in the physical space.
        (.evaluate()) """

        # reference solutions
        bspline_ref_evaluate = np.array(
            [[-0.019796,    0.04049403],
             [-0.9996,      0.069675  ],
             [ 2.256,       1.9691    ],
             [ 1.528,       4.0766    ],
             [-1.0264,      2.873488  ]]
        )
        nurbs_ref_evaluate = np.array(
            [[-1.00989841,  0.01432479],
             [-1.49984913,  0.02127445],
             [-0.15941144,  1.0883878 ],
             [-0.49948029,  1.62496752],
             [-1.61951381,  1.15640608]] 
        )
        bezier_ref_evaluate = np.array(
            [[-1.009899,  0.020099],
             [-1.49985,   0.02985 ],
             [-0.209,     1.089   ],
             [-0.612,     1.632   ],
             [-1.6716,    1.2736  ]] 
        )
        rational_ref_evaluate = np.array(
            [[-1.00989841,  0.01432479],
             [-1.49984913,  0.02127445],
             [-0.15941144,  1.0883878 ],
             [-0.49948029,  1.62496752],
             [-1.61951381,  1.15640608]]
        )

        # test evaluation
        self.assertTrue(np.allclose(
            np.array(self.bspline.evaluate(c.q2D)), 
            bspline_ref_evaluate
            )
        )
        self.assertTrue(np.allclose(
            np.array(self.nurbs.evaluate(c.q2D)), 
            nurbs_ref_evaluate
            )
        )
        self.assertTrue(np.allclose(
            np.array(self.bezier.evaluate(c.q2D)), 
            bezier_ref_evaluate
            )
        )
        self.assertTrue(np.allclose(
            np.array(self.rational.evaluate(c.q2D)), 
            rational_ref_evaluate
            )
        )

    def test_derivative(self):
        """ Test the correct calculation of the first derivative.  
        (.derivative()) """

        # reference solutions
        bspline_ref_derivative = np.array(
            [[0.0408,   4.019206],
             [0.08,     6.935   ],
             [6.88,     0.318   ],
             [6.72,     1.884   ],
             [4.768,    6.36784 ]] 
        )
        nurbs_ref_derivative = np.array(
            [[0.02037649, 1.43654295],
             [0.03026211, 2.13347963],
             [1.62487759, 0.23798877],
             [2.5357129,  0.77942396],
             [1.90293663, 2.66500861]] 
        )

        # order
        o1 = [[1,1]]

        # test derivative evaluation
        self.assertTrue(np.allclose(
            np.array(self.bspline.derivative(c.q2D,o1)), 
            bspline_ref_derivative
            )
        )
        self.assertTrue(np.allclose(
            np.array(self.nurbs.derivative(c.q2D,o1)), 
            nurbs_ref_derivative
            )
        )

if __name__ == "__main__":
    unittest.main()