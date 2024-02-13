try:
    from . import common as c
except BaseException:
    import common as c


class reparametrizeTest(c.SplineBasedTestCase):
    def test_permute_parametric_axes(self):
        """
        test permute
        """
        # Define some splines
        z = self.bezier_3p3d()
        r = self.rational_bezier_3p3d()
        b = self.bspline_3p3d()
        n = self.nurbs_3p3d()
        originals = (z, r, b, n)

        # define permutation
        permutation = [2, 0, 1]

        # return permuted
        for orig in originals:
            # make more work
            orig.elevate_degrees(1)
            orig.elevate_degrees(2)
            orig.elevate_degrees(2)
            if "knot_vectors" in orig.required_properties:
                orig.insert_knots(0, [0.4, 0.7, 0.8])
                orig.insert_knots(1, [0.1, 0.2])
                orig.insert_knots(2, [0.3, 0.5, 0.6, 0.9])

            perm = c.splinepy.helpme.reparametrize.permute_parametric_axes(
                orig, permutation, inplace=False
            )
            queries = c.np.asarray(c.get_queries_3D())

            self.assertTrue(
                c.np.allclose(
                    orig.evaluate(queries),
                    perm.evaluate(queries[:, permutation]),
                ),
                f"{perm.whatami} failed to permute.",
            )

        # ic.nplace
        for orig in originals:
            perm = orig.copy()
            c.splinepy.helpme.reparametrize.permute_parametric_axes(
                perm, permutation, inplace=True
            )
            queries = c.np.asarray(c.get_queries_3D())

            self.assertTrue(
                c.np.allclose(
                    orig.evaluate(queries),
                    perm.evaluate(queries[:, permutation]),
                ),
                f"{perm.whatami} failed to permute inplace.",
            )

    def test_flip_axes(self):
        """
        test invert axes
        """
        # Define some splines
        z = self.bezier_3p3d()
        r = self.rational_bezier_3p3d()
        b = self.bspline_3p3d()
        n = self.nurbs_3p3d()

        flipped_axes = [1, 2]

        queries = c.np.asarray(c.get_queries_3D())
        flipped_queries = queries.copy()
        for axis in flipped_axes:
            # Queries are all in unit cube
            flipped_queries[:, axis] = 1 - flipped_queries[:, axis]

        for spline in [z, r, b, n]:
            # With copy
            new_spline = c.splinepy.helpme.reparametrize.flip_axes(
                spline, flipped_axes, inplace=False
            )

            self.assertTrue(
                c.np.allclose(
                    spline.evaluate(queries),
                    new_spline.evaluate(flipped_queries),
                ),
                f"{spline.whatami} failed to flip axis.",
            )

            # in place
            spline_copy = spline.copy()
            new_spline2 = c.splinepy.helpme.reparametrize.flip_axes(
                spline_copy, flipped_axes, inplace=True
            )

            # Check if truly in place
            self.assertTrue(
                new_spline2 is spline_copy,
                "Operation has not been performed in place",
            )

            self.assertTrue(
                c.np.allclose(
                    new_spline2.evaluate(queries),
                    new_spline.evaluate(queries),
                ),
                f"{spline.whatami} failed to flip axis.",
            )


if __name__ == "__main__":
    c.unittest.main()
