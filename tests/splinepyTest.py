import splinepy
import numpy as np
import unittest

## abbreviations:
# z: bezier
# b: bspline
# n: nurbs


class TestSplinepy(unittest.TestCase):

    def setUp(self): 
        """ sets up the spline definitions that should match the 
        reference solutions hardcoded in the testfunctions. """

        # set up BSpline
        self.b2P2D = dict(
            degrees=[2, 2],
            knot_vectors=[
                [0, 0, 0, 0.5, 1, 1, 1],
                [0, 0, 0, 1, 1, 1],
            ],
            control_points=[
                [0, 0],
                [0, 1],
                [1, 1.5],
                [3, 1.5],
                [-1, 0],
                [-1, 2],
                [1, 4],
                [3, 4],
                [-2, 0],
                [-2, 2],
                [1, 5],
                [3, 5],
            ],
        )
        self.b = splinepy.BSpline(**self.b2P2D)

        # set up NURBS
        self.n2P2D = dict(
            degrees=[2, 1],
            knot_vectors=[
                [0, 0, 0, 1, 1, 1],
                [0, 0, 1, 1],
            ],
            control_points=[
                [-1.,  0.],
                [-1.,  1.],
                [ 0.,  1.],
                [-2.,  0.],
                [-2.,  2.],
                [ 0.,  2.],
            ],
           weights=[
               [1.],
               [(2 ** .5) / 2],
               [1.],
               [1.],
               [(2 ** .5) / 2],
               [1.],
            ],
        )
        self.n = splinepy.NURBS(**self.n2P2D)

        # set up Bezier spline
        self.z2P2D = dict(
            degrees=[2, 1],
            control_points=[
                [-1.,  0.],
                [-1.,  1.],
                [ 0.,  1.],
                [-2.,  0.],
                [-2.,  2.],
                [ 0.,  2.],
            ],
        )
        self.z = splinepy.Bezier(**self.z2P2D)

        # set up queries
        self.q2D = [
            [.01, .01], 
            [.01, .5], 
            [.9, .1], 
            [.8, .7], 
            [.4, .99],
        ]

    def test_constructor_degree_check(self):
        """ Test the raise Error for mismatch of degree and par_dim. """
        # not for Bezier

        self.b2P2D.update(degrees = [2,])
        self.n2P2D.update(degrees = [2,])

        with self.assertRaises(splinepy._spline.InputDimensionError) as context:
            splinepy.BSpline(**self.b2P2D)
        with self.assertRaises(splinepy._spline.InputDimensionError) as context:
            splinepy.NURBS(**self.n2P2D)

    def test_constructor_kv_nr_check(self):
        """ Test the raise Error for missing knot vectors. """
        # not for Bezier

        self.b2P2D.update(knot_vectors = [[0, 0, 0, 0.5, 1, 1, 1]])
        self.n2P2D.update(knot_vectors = [[0, 0, 0, 1, 1, 1]])

        with self.assertRaises(splinepy._spline.InputDimensionError) as context:
            splinepy.BSpline(**self.b2P2D)
        with self.assertRaises(splinepy._spline.InputDimensionError) as context:
            splinepy.NURBS(**self.n2P2D)

    def test_constructor_kv_length_check(self):
        """ Test the raise Error for length of knot vectors. """
        # not for Bezier

        self.b2P2D.update(knot_vectors = [
            [0, 0, 0.5, 1, 1], 
            [0, 0, 1, 1]
            ]
        )
        self.n2P2D.update(knot_vectors = [
            [0, 0, 1, 1], 
            [0, 0, 1, 1]
            ]
        )

        with self.assertRaises(splinepy._spline.InputDimensionError) as context:
            splinepy.BSpline(**self.b2P2D)
        with self.assertRaises(splinepy._spline.InputDimensionError) as context:
            splinepy.NURBS(**self.n2P2D)

    def test_constructor_cp_nr_check(self):
        """ Test the raise Error for invalid number of control points. """
        # not for Bezier

        self.b2P2D.update(
            control_points = [
                [0, 0],
                [0, 1],
                [1, 1.5],
                [3, 1.5],
                [-2, 0],
                [-2, 2],
                [1, 5],
                [3, 5],
                ]
        )
        self.n2P2D.update(
            control_points = [
                [-1.,  0.],
                [ 0.,  1.],
                [-2.,  0.],
                [0.,  2.],
                ],
            weights = [
                [1.],
                [1.],
                [1.],
                [1.],
                ],         
        )
        
        with self.assertRaises(splinepy._spline.InputDimensionError) as context:
            splinepy.BSpline(**self.b2P2D)
        with self.assertRaises(splinepy._spline.InputDimensionError) as context:
            splinepy.NURBS(**self.n2P2D)

    def test_basis_functions(self):
        """ Test the correct calculation of the basis functions. 
        (.basis_functions()) """
        # not for Bezier

        # reference solutions
        b_ref_basis_functions = np.array(
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
        n_ref_basis_functions = np.array(
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

        self.assertTrue(np.allclose(
                            b_ref_basis_functions, 
                            self.b.basis_functions(self.q2D)[0], 
                            )
        )
        self.assertTrue(np.allclose(
                            n_ref_basis_functions, 
                            self.n.basis_functions(self.q2D)[0], 
                            )
        )

    def test_evaluate(self):
        """ Test the correct spline evaluation in the physical space.
        (.evaluate()) """

        # reference solutions
        b_ref_evaluate = np.array(
            [[-0.019796,    0.04049403],
             [-0.9996,      0.069675  ],
             [ 2.256,       1.9691    ],
             [ 1.528,       4.0766    ],
             [-1.0264,      2.873488  ]]
        )
        n_ref_evaluate = np.array(
            [[-1.00989841,  0.01432479],
             [-1.49984913,  0.02127445],
             [-0.15941144,  1.0883878 ],
             [-0.49948029,  1.62496752],
             [-1.61951381,  1.15640608]] 
        )
        z_ref_evaluate = np.array(
            [[-1.009899,  0.020099],
             [-1.49985,   0.02985 ],
             [-0.209,     1.089   ],
             [-0.612,     1.632   ],
             [-1.6716,    1.2736  ]] 
        )

        self.assertTrue(np.allclose(
                            b_ref_evaluate, 
                            np.array(self.b.evaluate(self.q2D)), 
                            rtol=1.e-5, 
                            atol=1.e-8
                            )
        )
        self.assertTrue(np.allclose(
                            n_ref_evaluate, 
                            np.array(self.n.evaluate(self.q2D)), 
                            rtol=1.e-5, 
                            atol=1.e-8
                            )
        )
        self.assertTrue(np.allclose(
                            z_ref_evaluate, 
                            np.array(self.z.evaluate(self.q2D)), 
                            rtol=1.e-5, 
                            atol=1.e-8
                            )
        )

    def test_derivative(self):
        """ Test the correct calculation of the first derivative.  
        (.derivative()) """
        # not for Bezier

        # reference solutions
        b_ref_derivative = np.array(
            [[0.0408,   4.019206],
             [0.08,     6.935   ],
             [6.88,     0.318   ],
             [6.72,     1.884   ],
             [4.768,    6.36784 ]] 
        )
        n_ref_derivative = np.array(
            [[0.02037649, 1.43654295],
             [0.03026211, 2.13347963],
             [1.62487759, 0.23798877],
             [2.5357129,  0.77942396],
             [1.90293663, 2.66500861]] 
        )

        # order
        o1 = [[1,1]]

        self.assertTrue(np.allclose(
                            b_ref_derivative,  
                            np.array(self.b.derivative(self.q2D,o1)), 
                            rtol=1.e-5, 
                            atol=1.e-8
                            )
        )
        self.assertTrue(np.allclose(
                            n_ref_derivative,  
                            np.array(self.n.derivative(self.q2D,o1)), 
                            rtol=1.e-5, 
                            atol=1.e-8
                            )
        )

    def test_insert_knot(self):
        """ Test the knot insertion function (.insert_knot()). """ 
        # not for Bezier

        self.b.insert_knots(0,[.2, .7])
        self.b.insert_knots(1,[.2, .5, .7,])

        self.n.insert_knots(0,[.2, .5, .7,])
        self.n.insert_knots(1,[.2, .5, .7,])

        # reference solutions
        b_ref_kv = [
            [0, 0, 0, .2, 0.5, .7, 1, 1, 1],
            [0, 0, 0, .2, 0.5, .7, 1, 1, 1],
        ]
        b_ref_basis_functions = np.array(
            [[8.14506250e-01, 8.70912500e-02, 9.02500000e-04, 8.70912500e-02,
              9.31225000e-03, 9.65000000e-05, 9.02500000e-04, 9.65000000e-05,
              1.00000000e-06],
             [3.61000000e-01, 3.86000000e-02, 4.00000000e-04, 5.41500000e-01,
              5.79000000e-02, 6.00000000e-04, 0.00000000e+00, 0.00000000e+00,
              0.00000000e+00],
             [1.66666667e-02, 1.22222222e-01, 1.11111111e-01, 4.33333333e-02,
              3.17777778e-01, 2.88888889e-01, 6.66666667e-03, 4.88888889e-02,
              4.44444444e-02],
             [1.60000000e-01, 3.73333333e-01, 6.66666667e-02, 1.06666667e-01,
              2.48888889e-01, 4.44444444e-02, 0.00000000e+00, 0.00000000e+00,
              0.00000000e+00],
             [4.44444444e-05, 4.44444444e-04, 1.77777778e-04, 4.32592593e-03,
              4.32592593e-02, 1.73037037e-02, 6.22962963e-02, 6.22962963e-01,
              2.49185185e-01]] 
        )
        n_ref_kv=[
            [0, 0, 0, .2, .5, .7, 1, 1, 1],
            [0, 0, .2, .5, .7, 1, 1],
        ]
        n_ref_basis_functions = np.array(
            [[8.62376166e-01, 8.68082286e-02, 8.15605651e-04, 4.53882193e-02,
               4.56885413e-03, 4.29266132e-05],
              [9.07764385e-01, 9.13770827e-02, 8.58532265e-04, 0.00000000e+00,
               0.00000000e+00, 0.00000000e+00],
              [3.00352621e-02, 2.35374749e-01, 2.34589989e-01, 3.00352621e-02,
               2.35374749e-01, 2.34589989e-01],
              [2.51153840e-01, 6.26244066e-01, 1.22602094e-01, 0.00000000e+00,
               0.00000000e+00, 0.00000000e+00],
              [2.43427980e-03, 2.20707525e-02, 8.82830101e-03, 7.05941142e-02,
               6.40051823e-01, 2.56020729e-01]]
        )

        # test knot_vectors
        self.assertEqual(self.b.knot_vectors, b_ref_kv)
        self.assertEqual(self.n.knot_vectors, n_ref_kv)

        # test basis_functions
        self.assertTrue(np.allclose(
            self.b.basis_functions(self.q2D)[0], 
            b_ref_basis_functions
            )
        )
        self.assertTrue(np.allclose(
            self.n.basis_functions(self.q2D)[0], 
            n_ref_basis_functions
            )
        )

    def test_remove_knot(self):
        """ Test the function .remove_knots. 
        This test also depends on the function .insert_knots! """

        # insert knots
        self.b.insert_knots(0,[.2, .7,])
        self.b.insert_knots(1,[.2, .5, .7,])

        self.n.insert_knots(0,[.2, .5, .7,])
        self.n.insert_knots(1,[.2, .5, .7,])

        # remove knots
        self.b.remove_knots(0,[.2, .7])
        self.b.remove_knots(1,[.2, .5, .7,])

        self.n.remove_knots(0,[.2, .5, .7,])
        self.n.remove_knots(1,[.2, .5, .7,])

        # reference solutions
        b_ref_kv = [
            [0, 0, 0, 0.5, 1, 1, 1], 
            [0, 0, 0, 1, 1, 1]
        ]
        n_ref_kv = [
            [0, 0, 0, 1, 1, 1], 
            [0, 0, 1, 1]
        ]
        b_ref_basis_functions = np.array(
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
        n_ref_basis_functions = np.array(
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

        # test knot_vectors
        self.assertEqual(self.b.knot_vectors, b_ref_kv)
        self.assertEqual(self.n.knot_vectors, n_ref_kv)
        
        # test basis functions
        self.assertTrue(np.allclose(
                            b_ref_basis_functions, 
                            self.b.basis_functions(self.q2D)[0], 
                            )
        )
        self.assertTrue(np.allclose(
                            n_ref_basis_functions, 
                            self.n.basis_functions(self.q2D)[0], 
                            )
        )

    def test_elevate_degree(self):
        """ Test the order elevation function (.elevate_degree()). """
        # not for Bezier

        self.b.elevate_degree(0)
        self.n.elevate_degree(0)

        # reference solutions
        b_ref_kv = [
            [0, 0, 0, 0, 0.5, 0.5, 1, 1, 1, 1], 
            [0, 0, 0, 1, 1, 1]
        ]
        b_ref_basis_functions = np.array(
            [[9.22462279e-01, 5.64772824e-02, 1.15651800e-03, 3.92040000e-06,
              1.86356016e-02, 1.14095520e-03, 2.33640000e-05, 7.92000000e-08,
              9.41192000e-05, 5.76240000e-06, 1.18000000e-07, 4.00000000e-10],
             [2.35298000e-01, 1.44060000e-02, 2.95000000e-04, 1.00000000e-06,
              4.70596000e-01, 2.88120000e-02, 5.90000000e-04, 2.00000000e-06,
              2.35298000e-01, 1.44060000e-02, 2.95000000e-04, 1.00000000e-06],
             [3.24000000e-03, 8.10000000e-02, 3.11040000e-01, 4.14720000e-01,
              7.20000000e-04, 1.80000000e-02, 6.91200000e-02, 9.21600000e-02,
              4.00000000e-05, 1.00000000e-03, 3.84000000e-03, 5.12000000e-03],
             [2.88000000e-03, 2.88000000e-02, 3.88800000e-02, 1.94400000e-02,
              1.34400000e-02, 1.34400000e-01, 1.81440000e-01, 9.07200000e-02,
              1.56800000e-02, 1.56800000e-01, 2.11680000e-01, 1.05840000e-01],
             [8.00000000e-07, 9.60000000e-06, 6.40000000e-05, 2.56000000e-05,
              1.58400000e-04, 1.90080000e-03, 1.26720000e-02, 5.06880000e-03,
              7.84080000e-03, 9.40896000e-02, 6.27264000e-01, 2.50905600e-01]] 
        )
        n_ref_kv =  [
            [0, 0, 0, 0, 1, 1, 1, 1],
            [0, 0, 1, 1],
        ]
        n_ref_basis_functions = np.array(
            [[9.66199276e-01, 2.35617313e-02, 2.37997286e-04, 9.95774782e-07,
               9.75958864e-03, 2.37997286e-04, 2.40401298e-06, 1.00583311e-08],
              [4.87979432e-01, 1.18998643e-02, 1.20200649e-04, 5.02916557e-07,
               4.87979432e-01, 1.18998643e-02, 1.20200649e-04, 5.02916557e-07],
              [9.50089457e-04, 2.06434697e-02, 1.85791227e-01, 6.92615214e-01,
               1.05565495e-04, 2.29371885e-03, 2.06434697e-02, 7.69572460e-02],
              [2.64820523e-03, 2.55733320e-02, 1.02293328e-01, 1.69485135e-01,
               6.17914555e-03, 5.96711079e-02, 2.38684432e-01, 3.95465315e-01],
              [2.51334851e-03, 4.04517338e-03, 2.69678225e-03, 7.44695856e-04,
               2.48821503e-01, 4.00472164e-01, 2.66981443e-01, 7.37248897e-02]]
        )

        # test degree
        self.assertTrue(np.allclose(self.b.degrees, [3,2]))
        self.assertTrue(np.allclose(self.n.degrees, [3,1]))

        # test knot_vectors
        self.assertEqual(self.b.knot_vectors, b_ref_kv)
        self.assertEqual(self.n.knot_vectors, n_ref_kv)

        # test basis functions
        self.assertTrue(np.allclose(
            self.b.basis_functions(self.q2D)[0], 
            b_ref_basis_functions
            )
        )
        self.assertTrue(np.allclose(
            self.n.basis_functions(self.q2D)[0], 
            n_ref_basis_functions
            )
        )

    def test_reduce_degree(self):
        """ Test the function .reduce_degree. 
        This test also depends on the function .elevate_degree! """

        # elevate order
        self.b.elevate_degree(0)
        self.n.elevate_degree(0)

        # reduce order
        self.b.reduce_degree(0)    
        self.n.reduce_degree(0)

        # reference solutions
        b_ref_kv = [
            [0, 0, 0, 0.5, 1, 1, 1], 
            [0, 0, 0, 1, 1, 1]
        ]
        n_ref_kv = [
            [0, 0, 0, 1, 1, 1], 
            [0, 0, 1, 1]
        ]
        b_ref_basis_functions = np.array(
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
        n_ref_basis_functions = np.array(
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

        # test degrees
        self.assertTrue(np.allclose(self.b.degrees, [2, 2]))
        self.assertTrue(np.allclose(self.n.degrees, [2, 1]))

        # test knot_vectors
        self.assertTrue(self.b.knot_vectors, b_ref_kv)
        self.assertTrue(self.n.knot_vectors, n_ref_kv)

        # test basis functions
        self.assertTrue(np.allclose(
                            b_ref_basis_functions, 
                            self.b.basis_functions(self.q2D)[0], 
                            )
        )
        self.assertTrue(np.allclose(
                            n_ref_basis_functions, 
                            self.n.basis_functions(self.q2D)[0], 
                            )
        )

if __name__ == "__main__":
    unittest.main()