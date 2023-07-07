try:
    from . import common as c
except BaseException:
    import common as c


class VolumeIntegrationTest(c.unittest.TestCase):
    def test_volume_integration_1D(self):
        """
        Test volume integration for splines using numerical integration of the
        Jacobi-Determinant
        """
        # Test 1D
        bezier = c.splinepy.Bezier(degrees=[1], control_points=[[0], [2]])

        self.assertTrue(c.np.allclose(bezier.integrate.volume(), 2.0))

        # test for other types same spline
        self.assertTrue(c.np.allclose(bezier.bspline.integrate.volume(), 2.0))
        self.assertTrue(
            c.np.allclose(bezier.rationalbezier.integrate.volume(), 2.0)
        )
        self.assertTrue(c.np.allclose(bezier.nurbs.integrate.volume(), 2.0))

        # Check if equal after refinement
        bezier.elevate_degrees([0, 0, 0])

        self.assertTrue(c.np.allclose(bezier.integrate.volume(), 2.0))

        # Test knot insertion
        bspline = bezier.bspline
        bspline.insert_knots(0, c.np.random.rand(10))
        self.assertTrue(c.np.allclose(bezier.integrate.volume(), 2.0))

    def test_volume_integration_2D(self):
        """
        Test volume integration for splines using numerical integration of the
        Jacobi-Determinant
        """
        # Test 1D
        bezier = c.splinepy.Bezier(
            degrees=[1, 1], control_points=[[0, 0], [2, 0], [0, 1], [2, 1]]
        )

        self.assertTrue(c.np.allclose(bezier.integrate.volume(), 2.0))

        # test for other types same spline
        self.assertTrue(c.np.allclose(bezier.bspline.integrate.volume(), 2.0))
        self.assertTrue(
            c.np.allclose(bezier.rationalbezier.integrate.volume(), 2.0)
        )
        self.assertTrue(c.np.allclose(bezier.nurbs.integrate.volume(), 2.0))

        # Check if equal after refinement
        bezier.elevate_degrees([0, 0, 1])

        self.assertTrue(c.np.allclose(bezier.integrate.volume(), 2.0))

        # Test knot insertion
        bspline = bezier.bspline
        bspline.insert_knots(0, c.np.random.rand(10))
        self.assertTrue(c.np.allclose(bezier.integrate.volume(), 2.0))

        # Move control points along axis does not change volume
        mi = bspline.multi_index
        bspline.cps[mi[3, :]][:, 1] += 10
        self.assertTrue(c.np.allclose(bezier.integrate.volume(), 2.0))

    def test_volume_integration_3D(self):
        """
        Test volume integration for splines using numerical integration of the
        Jacobi-Determinant
        """
        # Test 1D
        bezier = c.splinepy.Bezier(
            degrees=[1, 1, 1],
            control_points=[
                [0, 0, 0],
                [1, 0, 0],
                [0, 1, 0],
                [1, 1, 0],
                [0, 0, 3],
                [1, 0, 3],
                [0, 1, 3],
                [1, 1, 3],
            ],
        )

        self.assertTrue(c.np.allclose(bezier.integrate.volume(), 3.0))

        # Check if equal after refinement
        bezier.elevate_degrees([0, 0, 1, 1, 2, 2])

        self.assertTrue(c.np.allclose(bezier.integrate.volume(), 3.0))

        # Test knot insertion
        bspline = bezier.bspline
        bspline.insert_knots(0, c.np.random.rand(10))
        self.assertTrue(c.np.allclose(bezier.integrate.volume(), 3.0))

        # Move control points along axis does not change volume
        mi = bspline.multi_index
        bspline.cps[mi[3, :, :]][:, 1] += 10
        self.assertTrue(c.np.allclose(bezier.integrate.volume(), 3.0))

    def test_complex_geometry(self):
        """Test on a more complex geometry"""
        bezier = c.splinepy.Bezier(
            degrees=[1, 1, 1],
            control_points=[
                [0, 0, 0],
                [1, 0, 0],
                [0, 1, 0],
                [1, 1, 0],
                [0, 0, 3],
                [1, 0, 3],
                [0, 1, 3],
                [1, 1, 3],
            ],
        )
        bezier.elevate_degrees([0, 1, 2, 0, 1, 2])
        bezier.cps += c.np.random.rand(*bezier.cps.shape) * 0.1

        bezier_c = bezier.copy()

        # Move around internal control points

        mi = bezier.multi_index
        bezier.cps[mi[1:3, 1:3, 1:3]][:] += 0.05 * c.np.random.rand(
            *bezier.cps[mi[1:3, 1:3, 1:3]][:].shape
        )
        self.assertTrue(
            c.np.allclose(
                bezier.integrate.volume(), bezier_c.integrate.volume()
            )
        )

    def test_assertions(self):
        """Test the assertions in volume function"""
        bezier = c.splinepy.Bezier(
            degrees=[1, 2], control_points=c.np.random.rand(6, 3)
        )
        self.assertRaisesRegex(
            ValueError,
            "`Volume` of embedded spline depends on projection, "
            "integration is aborted",
            bezier.integrate.volume,
        )
        bezier = c.splinepy.Bezier(
            degrees=[2, 2], control_points=c.np.random.rand(9, 2)
        )
        self.assertRaisesRegex(
            ValueError,
            "Integration order must be array of size para_dim",
            bezier.integrate.volume,
            orders=[2, 4, 5],
        )


if __name__ == "__main__":
    c.unittest.main()
