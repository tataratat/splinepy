try:
    from . import common as c
except BaseException:
    import common as c


class BezierExtractionTest(c.unittest.TestCase):
    def test_greville_points(self):
        """
        test permute
        """
        # Define some splines
        b = c.splinepy.BSpline(
          knot_vectors=[[0,0,0,1,2,2,3,4,4,4]], 
          degrees=[2],
          control_points=c.np.random.rand(7,2))
        n = c.splinepy.BSpline(
          knot_vectors=[[0,0,0,0,1,1,2,3,4,4,4,4]], 
          degrees=[3],
          control_points=c.np.random.rand(8,2))

        # Extract Beziers
        b_beziers=b.extract_bezier_patches()
        n_beziers=n.extract_bezier_patches()

        # Loop over knot_spans and test at random points
        for offset in range(4):
            queries = c.np.random.rand(20,1)
            self.assertTrue(
                c.np.allclose(
                  b.evaluate(queries + offset),
                  b_beziers[offset].evaluate(queries)
                )
            )
            self.assertTrue(
                c.np.allclose(
                  n.evaluate(queries + offset),
                  n_beziers[offset].evaluate(queries)
                )
            )


if __name__ == "__main__":
    c.unittest.main()
