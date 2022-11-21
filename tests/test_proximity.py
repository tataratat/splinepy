try:
    from . import common as c
except BaseException:
    import common as c


class ProximityTest(c.unittest.TestCase):

    def test_queries_inside_bspline_initial_guess_with_kdt(self):
        """
        Initial guess made with kdt. Mid-point as initial guess tends to fail,
        so excluded from test.
        """
        # make "iga book" bspline
        bspline = c.splinepy.BSpline(**c.b2P2D)

        # form parametric queries
        # -> parametric space is [0, 1], so no further manipulation
        para_q = c.np.random.random((10, bspline.para_dim))

        # form physical queries
        phys_q = bspline.evaluate(para_q)

        # proximity - single thread exe
        prox_r = bspline.proximities(
                queries=phys_q,
                initial_guess_sample_resolutions=[10] * bspline.para_dim,
                nthreads=1,
        )

        assert c.np.allclose(para_q, prox_r[0]), "WRONG `Spline.nearest_pcoord`"

    def test_queries_inside_nurbs_initial_guess_with_kdt(self):
        # make half-half circle. also known as quarter circle
        nurbs = c.splinepy.NURBS(**c.n2P2D)

        # form parametric queries
        # -> parametric space is [0, 1], so no further manipulation
        para_q = c.np.random.random((10, nurbs.para_dim))

        # form physical queries
        phys_q = nurbs.evaluate(para_q)

        # proximity - single thread exe
        prox_r = nurbs.proximities(
                queries=phys_q,
                initial_guess_sample_resolutions=[10] * nurbs.para_dim,
                nthreads=1,
        )

        assert c.np.allclose(para_q, prox_r[0]), "WRONG `Spline.nearest_pcoord`"


if __name__ == "__main__":
    c.unittest.main()
