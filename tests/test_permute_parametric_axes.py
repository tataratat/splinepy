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
        z = splinepy.Bezier(**c.z2P2D)
        r = splinepy.RationalBezier(**c.r2P2D)
        b = splinepy.BSpline(**c.b2P2D)
        n = splinepy.NURBS(**c.n2P2D)
        tobepermuted = (z, r, b, n)
        permutation = [1, 0]

        # save original
        zo, ro, bo, no = z.copy(), r.copy(), b.copy(), n.copy()
        originals = (zo, ro, bo, no)

        for orig, tbp in zip(originals, tobepermuted):
            perm = tbp.permute_parametric_axes(permutation)
            queries = np.asarray(c.q2D)

            self.assertTrue(
                np.allclose(
                    orig.evaluate(queries),
                    perm.evaluate(queries[:, permutation]),
                ),
                f"{tbp.whatami} failed to permute.",
            )


if __name__ == "__main__":
    c.unittest.main()
