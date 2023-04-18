try:
    from . import common as c
except BaseException:
    import common as c


class GrevilleAbscissaeTest(c.unittest.TestCase):
    def test_greville_points(self):
        """
        test permute
        """
        # Define some splines
        z = c.splinepy.Bezier(**c.z3P3D)
        r = c.splinepy.RationalBezier(**c.r3P3D)
        b = c.splinepy.BSpline(**c.b3P3D)
        n = c.splinepy.NURBS(**c.n3P3D)

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

        self.assertTrue(c.np.allclose(greville_points, z.greville_abscissae))
        greville_points = c.splinepy.utils.data.cartesian_product(
            [c.np.linspace(0, 1, r.degrees[i] + 1) for i in range(r.para_dim)]
        )
        self.assertTrue(c.np.allclose(greville_points, r.greville_abscissae))

        # Test non-unifrom Splines
        greville_points = c.splinepy.utils.data.cartesian_product(
            [
                c.np.convolve(b.knot_vectors[i][1:-1], c.np.ones(b.degrees[i]))
                for i in range(b.para_dim)
            ]
        )
        self.assertTrue(c.np.allclose(greville_points, b.greville_abscissae))
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
        self.assertTrue(c.np.allclose(greville_points, n.greville_abscissae))


if __name__ == "__main__":
    c.unittest.main()
