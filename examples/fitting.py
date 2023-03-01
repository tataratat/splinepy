import numpy as np

import splinepy

if __name__ == "__main__":
    try:
        import gustaf as gus

        gustaf_available = True
    except ImportError:
        gustaf_available = False

    # 1. Curve fitting

    # Create a function that creates the points for approximation
    def f(x):
        return x**4 - 2 * x**2 + x + 0.2 * np.sin(10 * x)

    def evaluate_approximation(spline, points, n_initial_guess=50):
        # Caluclate proximities between the spline and the given points
        prox_r = spline.proximities(
            queries=points,
            initial_guess_sample_resolutions=[n_initial_guess]
            * spline.para_dim,
            nthreads=1,
        )
        # Just evaluate the distance
        distances = prox_r[3]

        return np.max(distances), np.mean(distances)

    n_sample = 25
    n_resolution = 3 * n_sample
    x = np.linspace(-1, 1, n_sample)
    target_points = np.vstack([x, f(x)]).T.reshape(-1, 2)

    # Curve interpolation creates Bspline where the number of control points
    # equals the sampling points
    interpolated_curve = splinepy.BSpline.interpolate_curve(
        target_points, degree=2
    )

    # In curve approximation the number of control points must be choosen
    n_cps = 10
    approximated_curve = splinepy.BSpline.approximate_curve(
        target_points, degree=3, num_control_points=n_cps
    )
    interp_errs = evaluate_approximation(interpolated_curve, target_points)
    approx_errs = evaluate_approximation(approximated_curve, target_points)

    if gustaf_available:
        verts = gus.Vertices(target_points)
        gus.show.show_vedo(
            [
                f"Interpolated curve \nMax. distance {interp_errs[0]:.2f}, "
                f"mean distance {interp_errs[1]:.2f}",
                verts,
                gus.BSpline(**interpolated_curve.todict()),
            ],
            [
                f"Approximated curve with {n_cps} CPs \nMax. distance "
                f"{approx_errs[0]:.2f}, mean distance {approx_errs[1]:.2f}",
                gus.BSpline(**approximated_curve.todict()),
            ],
        )

    approximated_curves = []

    # The number of control points affects the approximation error
    for n_cps in range(n_sample, 5, -5):
        approximated_curve = splinepy.BSpline.approximate_curve(
            target_points, degree=3, num_control_points=n_cps
        )

        n_cps = approximated_curve.control_mesh_resolutions[0]
        max_error, mean_error = evaluate_approximation(
            approximated_curve, target_points
        )
        approximated_curves.append(
            [n_cps, approximated_curve, [max_error, mean_error]]
        )

    if gustaf_available:
        gus_list = [
            [
                f"Approximated curve with {cps} CPs \n Max. "
                f"distance {errs[0]:.2f}, mean distance {errs[1]:.2f}",
                gus.BSpline(**curve.todict()),
                verts,
            ]
            for (cps, curve, errs) in approximated_curves
        ]
        gus.show.show_vedo(*gus_list)

    # 2. Surface fitting
    n_sample = [20, 20]
    resolution = [3 * n for n in n_sample]

    # Create dataset
    x = [np.linspace(-1, 1, n) for n in n_sample]
    xx, yy = np.meshgrid(x[0], x[1])

    def h(x, y):
        return -0.9 * (x**3) + 0.7 * y**2 + x * y + 0.1 * np.sin(5 * x)

    target_points = np.vstack(
        (xx.flatten(), yy.flatten(), h(xx, yy).flatten())
    ).T

    interpolated_surface = splinepy.BSpline.interpolate_surface(
        target_points,
        size_u=n_sample[0],
        size_v=n_sample[1],
        degree_u=3,
        degree_v=2,
    )

    n_cps = [int(0.75 * n_sample[0]), int(0.6 * n_sample[1])]
    approximated_surface = splinepy.BSpline.approximate_surface(
        target_points,
        num_points_u=n_sample[0],
        num_points_v=n_sample[1],
        size_u=n_cps[0],
        size_v=n_cps[1],
        degree_u=3,
        degree_v=2,
    )

    if gustaf_available:
        verts = gus.Vertices(target_points)
        bspline = gus.BSpline(**interpolated_surface.todict())
        gus.show.show_vedo(
            ["Target points", verts],
            [
                "Interpolated surface",
                gus.BSpline(**interpolated_surface.todict()),
                verts,
            ],
            [
                f"Approximated surface with ({n_cps}) CPs",
                gus.BSpline(**approximated_surface.todict()),
                verts,
            ],
        )

    approximated_surfaces = []
    # The approximation error changes with the number of control points
    for r in range(0, min(n_sample) - 2, 5):
        approximated_surface = splinepy.BSpline.approximate_surface(
            target_points,
            num_points_u=n_sample[0],
            num_points_v=n_sample[1],
            size_u=n_sample[0] - int(r / 2),
            size_v=n_sample[1] - r,
            degree_u=3,
            degree_v=2,
        )

        n_cps = approximated_surface.control_mesh_resolutions
        max_error, mean_error = evaluate_approximation(
            approximated_surface, target_points
        )

        approximated_surfaces.append(
            [n_cps, approximated_surface, [max_error, mean_error]]
        )

    if gustaf_available:
        verts = gus.Vertices(target_points)
        gus_list = [
            [
                f"Approximated surface with {cps} CPs\n Max. distance "
                f"{errs[0]:.2f}, mean distance {errs[1]:.2f}",
                gus.BSpline(**surface.todict()),
                verts,
            ]
            for (cps, surface, errs) in approximated_surfaces
        ]
        gus.show.show_vedo(*gus_list)
