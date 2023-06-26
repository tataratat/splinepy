try:
    from . import common as c
except BaseException:
    import common as c


class BezierExtractionTest(c.unittest.TestCase):
    def test_extraction(self):
        """
        test the extraction of beziers
        """
        # Define some splines
        b = c.splinepy.BSpline(
            knot_vectors=[[0, 0, 0, 1, 2, 2, 3, 4, 4, 4]],
            degrees=[2],
            control_points=c.np.random.rand(7, 2),
        )
        n = c.splinepy.BSpline(
            knot_vectors=[[0, 0, 0, 0, 1, 1, 2, 3, 4, 4, 4, 4]],
            degrees=[3],
            control_points=c.np.random.rand(8, 2),
        )

        # Extract Beziers
        b_beziers = b.extract_bezier_patches()
        n_beziers = n.extract_bezier_patches()

        # Loop over knot_spans and test at random points
        for offset in range(4):
            queries = c.np.random.rand(20, 1)
            self.assertTrue(
                c.np.allclose(
                    b.evaluate(queries + offset),
                    b_beziers[offset].evaluate(queries),
                )
            )
            self.assertTrue(
                c.np.allclose(
                    n.evaluate(queries + offset),
                    n_beziers[offset].evaluate(queries),
                )
            )

    def test_extraction_matrices(self):
        """Create matrices to extract splines"""
        # Init splines
        bspline = c.splinepy.BSpline(**c.b3P3D)
        bspline.elevate_degrees([0, 1, 2])
        bspline.insert_knots(0, c.np.random.rand(3))
        bspline.insert_knots(1, c.np.random.rand(3))
        bspline.insert_knots(2, c.np.random.rand(3))
        # nurbs = c.n2P2D.copy()

        # Extract splines
        beziers_b = bspline.extract_bezier_patches()
        b_matrices = bspline.knot_insertion_matrix(beziers=True)
        for m, b in zip(b_matrices, beziers_b):
            # Test matrices m against spline ctps
            self.assertTrue(
                c.np.allclose(b.control_points, m @ bspline.control_points)
            )


if __name__ == "__main__":
    c.unittest.main()
