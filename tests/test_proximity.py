try:
    from . import common as c
except BaseException:
    import common as c


class ProximityTest(c.unittest.TestCase):
    def test_queries_inside_spline_initial_guess_with_kdt(self):
        """
        Initial guess made with kdt. Mid-point as initial guess tends to fail,
        so excluded from test.
        """
        for spline in c.spline_types_as_list():
            # form parametric queries
            # -> parametric space is [0, 1], so no further manipulation
            para_q = c.np.random.random((10, spline.para_dim))

            # form physical queries
            phys_q = spline.evaluate(para_q)

            # proximity - single thread exe
            prox_r = spline.proximities(
                queries=phys_q,
                initial_guess_sample_resolutions=[10] * spline.para_dim,
                nthreads=1,
                return_verbose=True,
            )

            assert c.np.allclose(
                para_q, prox_r[0]
            ), f"WRONG proximity query for {spline.whatami}"


if __name__ == "__main__":
    c.unittest.main()
