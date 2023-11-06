try:
    from . import common as c
except BaseException:
    import common as c


class GrevilleAbscissaeTest(c.SplineBasedTestCase):
    def test_greville_points(self):
        """
        test permute
        """
        # Define some splines
        z = self.bezier_3p3d()
        r = self.rational_bezier_3p3d()
        b = self.bspline_3p3d()
        n = self.nurbs_3p3d()

        # Modify some splines
        b.insert_knots(0, [0.2, 0.7])
        b.insert_knots(1, [0.2, 0.5, 0.7])
        b.insert_knots(2, [0.5])
        n.elevate_degrees([0, 0, 1, 2])
        n.insert_knots(0, [0.25])
        # Test Uniform Types
        greville_points = c.splinepy.utils.data.cartesian_product(
            [c.np.linspace(0, 1, z.degrees[i] + 1) for i in range(z.para_dim)]
        )

        self.assertTrue(c.np.allclose(greville_points, z.greville_abscissae()))
        greville_points = c.splinepy.utils.data.cartesian_product(
            [c.np.linspace(0, 1, r.degrees[i] + 1) for i in range(r.para_dim)]
        )
        self.assertTrue(c.np.allclose(greville_points, r.greville_abscissae()))

        # Test non-unifrom Splines
        greville_points = c.splinepy.utils.data.cartesian_product(
            [
                c.np.convolve(b.knot_vectors[i][1:-1], c.np.ones(b.degrees[i]))
                for i in range(b.para_dim)
            ]
        )
        self.assertTrue(c.np.allclose(greville_points, b.greville_abscissae()))
        greville_points = c.splinepy.utils.data.cartesian_product(
            [
                c.np.convolve(
                    n.knot_vectors[i][1:-1],
                    c.np.ones(n.degrees[i]) / n.degrees[i],
                    mode="valid",
                )
                for i in range(n.para_dim)
            ]
        )
        self.assertTrue(c.np.allclose(greville_points, n.greville_abscissae()))

    def test_greville_with_duplicate_points(self):
        """Tests if duplicate points are filtered out, as is required for the
        construction of c^(-1) spline"""
        # Construct c^(-1) spline
        a = c.splinepy.BSpline(
            degrees=[2],
            knot_vectors=[[0, 0, 0, 0.25, 0.5, 0.5, 0.5, 0.75, 1, 1, 1]],
            control_points=c.np.ones((8, 2)),
        )

        self.assertTrue(
            c.np.allclose(
                a.greville_abscissae(duplicate_tolerance=-1.0).ravel(),
                [0, 0.125, 0.375, 0.5, 0.5, 0.625, 0.875, 1.0],
            )
        )
        self.assertTrue(
            c.np.allclose(
                a.greville_abscissae(
                    duplicate_tolerance=c.splinepy.settings.TOLERANCE
                ).ravel(),
                [0, 0.125, 0.375, 0.4375, 0.5625, 0.625, 0.875, 1.0],
            )
        )


if __name__ == "__main__":
    c.unittest.main()
