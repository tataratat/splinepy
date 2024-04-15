import numpy as np

import splinepy


def test_fit_curve_2d_interpolation():
    n_sample = 20
    resolution = [3 * n_sample]
    degree = 2

    x = np.linspace(-2, 2, n_sample)

    def f(x):
        return x**2

    fitting_points = np.vstack((x, f(x))).T.reshape(-1, 2)

    interpolation_spline, _ = splinepy.helpme.fit.curve(
        fitting_points=fitting_points, degree=degree
    )
    interpolated_points = interpolation_spline.sample(resolutions=resolution)

    assert np.allclose(
        f(interpolated_points[:, 0]),
        interpolated_points[:, 1],
        atol=1e-3,
    )


def test_fit_curve_2d_approximation():
    n_sample = 20
    resolution = [3 * n_sample]
    degree = 2

    x = np.linspace(-2, 2, n_sample)

    def f(x):
        return x**2

    fitting_points = np.vstack((x, f(x))).T.reshape(-1, 2)

    approximated_spline, _ = splinepy.helpme.fit.curve(
        fitting_points=fitting_points,
        degree=degree,
        n_control_points=n_sample - 5,
    )
    approximated_points = approximated_spline.sample(resolution)

    assert np.allclose(
        approximated_points[:, 1],
        f(approximated_points[:, 0]),
        atol=1e-2,
    )


def test_fit_surface_3d_interpolation():
    sample_size = [15, 10]
    x = [np.linspace(-2, 2, n) for n in sample_size]
    xx, yy = np.meshgrid(x[0], x[1])

    def f(x, y):
        return x**2 + y**2

    fitting_points = np.vstack(
        (xx.flatten(), yy.flatten(), f(xx, yy).flatten())
    ).T

    bspline, _ = splinepy.helpme.fit.surface(
        fitting_points=fitting_points,
        size=[sample_size[0], sample_size[1]],
        n_control_points=[sample_size[0], sample_size[1]],
        degrees=[2, 2],
    )
    approximated_points = bspline.sample(sample_size)

    assert np.allclose(
        f(approximated_points[:, 0], approximated_points[:, 1]),
        approximated_points[:, -1],
        atol=1e-2,
    )


def test_fit_surface_3d_approximation():
    sample_size = [10, 15]
    x = [np.linspace(-1, 1, n) for n in sample_size]
    xx, yy = np.meshgrid(x[0], x[1])

    def f(x, y):
        return 0.5 * (x**2 + y**2)

    fitting_points = np.vstack(
        (xx.flatten(), yy.flatten(), f(xx, yy).flatten())
    ).T

    bspline, _ = splinepy.helpme.fit.surface(
        fitting_points=fitting_points,
        size=[sample_size[0], sample_size[1]],
        n_control_points=[sample_size[0] - 2, sample_size[1] - 2],
        degrees=[3, 3],
    )
    approximated_points = bspline.sample(sample_size)

    assert np.allclose(
        f(approximated_points[:, 0], approximated_points[:, 1]),
        approximated_points[:, -1],
        atol=1e-3,
    )
