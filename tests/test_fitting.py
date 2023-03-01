try:
    from . import common as c
except BaseException:
    import common as c


class TestSplinepyFitting(c.unittest.TestCase):
    """The available fitting functions are
    - interpolate_curve
    - approximate_curve
    - interpolate_surface
    - approximate_surface
    These functions determine the parametrization and control point
    positions. As the parametrization can be arbitrary, the test cases
    only consider if the spline is close in physical space to the dataset,
    which it shall approximate
    The test cases generate a dataset from a known function, which will
    be used for validation.
    """

    def test_interpolate_curve_2d(self):
        n_sample = 20
        resolution = [3 * n_sample]
        degree = 2

        x = c.np.linspace(-2, 2, n_sample)

        def f(x):
            return x**2

        target_points = c.np.vstack((x, f(x))).T.reshape(-1, 2)

        interpolation_spline = c.splinepy.BSpline.interpolate_curve(
            target_points, degree=degree
        )
        interpolated_points = interpolation_spline.sample(
            resolutions=resolution
        )

        self.assertTrue(
            c.np.allclose(
                f(interpolated_points[:, 0]),
                interpolated_points[:, 1],
                atol=1e-3,
            )
        )

    def test_approximate_curve_2d(self):
        n_sample = 20
        resolution = [3 * n_sample]
        degree = 2

        x = c.np.linspace(-2, 2, n_sample)

        def f(x):
            return x**2

        target_points = c.np.vstack((x, f(x))).T.reshape(-1, 2)

        approximated_spline = c.splinepy.BSpline.approximate_curve(
            target_points, degree=degree, num_control_points=n_sample - 5
        )
        approximated_points = approximated_spline.sample(resolution)

        self.assertTrue(
            c.np.allclose(
                approximated_points[:, 1],
                f(approximated_points[:, 0]),
                atol=1e-2,
            )
        )

    def test_interpolate_surface_3d(self):
        sample_size = [10, 10]
        x = [c.np.linspace(-2, 2, n) for n in sample_size]
        xx, yy = c.np.meshgrid(x[0], x[1])

        def f(x, y):
            return x**2 + y**2

        target_points = c.np.vstack(
            (xx.flatten(), yy.flatten(), f(xx, yy).flatten())
        ).T

        bspline = c.splinepy.BSpline.interpolate_surface(
            target_points,
            size_u=sample_size[0],
            size_v=sample_size[1],
            degree_u=2,
            degree_v=2,
            reorganize=False,
        )
        approximated_points = bspline.sample([sample_size])

        self.assertTrue(
            c.np.allclose(
                f(approximated_points[:, 0], approximated_points[:, 1]),
                approximated_points[:, -1],
                atol=1e-2,
            )
        )

    def test_approximate_surface_3d(self):
        sample_size = [10, 10]
        x = [c.np.linspace(-1, 1, n) for n in sample_size]
        xx, yy = c.np.meshgrid(x[0], x[1])

        def f(x, y):
            return 0.5 * (x**2 + y**2)

        target_points = c.np.vstack(
            (xx.flatten(), yy.flatten(), f(xx, yy).flatten())
        ).T

        bspline = c.splinepy.BSpline.approximate_surface(
            target_points,
            size_u=sample_size[0] - 2,
            size_v=sample_size[1] - 2,
            num_points_u=sample_size[0],
            num_points_v=sample_size[1],
            degree_u=3,
            degree_v=3,
            reorganize=False,
        )
        approximated_points = bspline.sample(sample_size)

        self.assertTrue(
            c.np.allclose(
                f(approximated_points[:, 0], approximated_points[:, 1]),
                approximated_points[:, -1],
                atol=1e-3,
            )
        )


if __name__ == "__main__":
    c.unittest.main()
