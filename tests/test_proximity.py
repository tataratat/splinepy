import splinepy
import numpy as np
try:
    from . import common as c
except:
    import common as c


class ProximityTest(c.unittest.TestCase):

    def test_queries_inside_bspline_initial_guess_with_kdt(self):
        """
        Initial guess made with kdt. Mid-point as initial guess tends to fail,
        so excluded from test.
        """
        # make "iga book" bspline
        bspline = splinepy.BSpline(**c.b2P2D)

        # form parametric queries
        # -> parametric space is [0, 1], so no further manipulation
        para_q = np.random.random((10, bspline.para_dim))

        # form physical queries
        phys_q = bspline.evaluate(para_q)

        # proximity - single thread exe
        prox_r = bspline.nearest_pcoord(
            queries=phys_q,
            kdt_resolutions=[10] * bspline.para_dim,
            n_threads=1,
        )

        assert np.allclose(para_q, prox_r), "WRONG `Spline.nearest_pcoord`"

    def test_queries_inside_nurbs_initial_guess_with_kdt(self):
        # make half-half circle. also known as quarter circle
        nurbs = splinepy.NURBS(**c.n2P2D)

        # form parametric queries
        # -> parametric space is [0, 1], so no further manipulation
        para_q = np.random.random((10, nurbs.para_dim))

        # form physical queries
        phys_q = nurbs.evaluate(para_q)

        # proximity - single thread exe
        prox_r = nurbs.nearest_pcoord(
            queries=phys_q,
            kdt_resolutions=[10] * nurbs.para_dim,
            n_threads=1,
        )

        assert np.allclose(para_q, prox_r), "WRONG `Spline.nearest_pcoord`"


if __name__ == "__main__":
    c.unittest.main()
