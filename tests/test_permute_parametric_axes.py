import tempfile

import splinepy
import numpy as np
try:
    from . import common as c
except:
    import common as c


class PermuteParametricAxesTest(c.unittest.TestCase):

    def test_permute_parametric_axes(self):
        """
        test permute
        """
        # Define some splines
        z = splinepy.Bezier(**c.z3P3D)
        r = splinepy.RationalBezier(**c.r3P3D)
        b = splinepy.BSpline(**c.b3P3D)
        n = splinepy.NURBS(**c.n3P3D)
        originals = (z, r, b, n)

        # define permutation
        permutation = [2, 0, 1]

        # return permuted
        for orig in originals:
            # make more work
            orig.elevate_degree(1)
            orig.elevate_degree(2)
            orig.elevate_degree(2)
            if "knot_vectors" in orig.required_properties:
                orig.insert_knots(0, [.4, .7, .8])
                orig.insert_knots(1, [.1, .2])
                orig.insert_knots(2, [.3, .5, .6, .9])

            perm = orig.permute_parametric_axes(permutation, inplace=False)
            queries = np.asarray(c.q3D)

            self.assertTrue(
                    np.allclose(
                        orig.evaluate(queries),
                        perm.evaluate(queries[:, permutation]),
                    ),
                    f"{perm.whatami} failed to permute.",
            )

        # inplace
        for orig in originals:
            perm = orig.copy()
            perm.permute_parametric_axes(permutation, inplace=True)
            queries = np.asarray(c.q3D)

            self.assertTrue(
                    np.allclose(
                        orig.evaluate(queries),
                        perm.evaluate(queries[:, permutation])
                    ),
                    f"{perm.whatami} failed to permute inplace.",
            )
        

if __name__ == "__main__":
    c.unittest.main()
