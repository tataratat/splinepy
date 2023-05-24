try:
    from . import common as c
except BaseException:
    import common as c


class compositionSensitivitiesTest(c.unittest.TestCase):
    """
    Test composition sensitivity, i.e., the derivative concerning the
    deformation function's control points.
    """

    def test_composition(self):
        # Here we only try 2D compositions (surface-line is tested in bezman)
        # Init Splines to be tested
        surface_polynomial = c.splinepy.Bezier(
            degrees=[1, 2],
            control_points=[[1, 1], [2, 1], [1, 2], [2, 2], [1, 3], [2, 3]],
        )
        surface_rational = c.splinepy.RationalBezier(
            degrees=[2, 1],
            control_points=[[1, 1], [2, 1], [3, 2], [1, 2], [2, 2], [3, 3]],
            weights=[1.0, 0.8, 1.0, 1.0, 0.8, 1.0],
        )
        inner_polynomial = c.splinepy.Bezier(
            degrees=[1, 1],
            control_points=[[0.2, 0.0], [0.0, 0.2], [0.0, 0.8], [0.8, 0]],
        )
        inner_rational = c.splinepy.RationalBezier(
            degrees=[2, 1],
            control_points=[
                [0.2, 0.0],
                [0.5, 0.5],
                [0.0, 0.2],
                [0.0, 0.8],
                [0.5, 0.5],
                [0.8, 0],
            ],
            weights=[1, 0.5, 1, 1, 0.5, 1],
        )

        dx = 1e-4
        tolerance = 1e-4

        def test_splines(outer, inner):
            # Determine base line
            composed = outer.compose(inner)
            composed_sensitivities = outer.compose_sensitivities(inner)

            #
            for i in range(outer.cps.shape[0]):
                for j in range(outer.cps.shape[1]):
                    outer.cps[i, j] += dx

                    # FD values
                    composed_dx = outer.compose(inner)

                    # FD
                    derivative_fd = (composed_dx.cps - composed.cps) / dx
                    expected = c.np.zeros(composed.cps.shape)
                    expected[:, j] = composed_sensitivities[i].cps.flatten()

                    self.assertTrue(
                        c.np.allclose(derivative_fd, expected, atol=tolerance)
                    )

                    outer.cps[i, j] -= dx

        test_splines(surface_polynomial, inner_polynomial)
        test_splines(surface_polynomial, inner_rational)
        test_splines(surface_rational, inner_polynomial)
        test_splines(surface_rational, inner_rational)


if __name__ == "__main__":
    c.unittest.main()
