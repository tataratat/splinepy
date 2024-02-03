import gustaf as gus
import numpy as np

import splinepy

if __name__ == "__main__":
    # 1. Curve fitting
    # Create a function that creates the points for approximation
    def f(x):
        return x**4 - 2 * x**2 + x + 0.2 * np.sin(10 * x)

    def evaluate_approximation(spline, points, n_initial_guess=500):
        # Calculate proximities between the spline and the given points
        prox_r = spline.proximities(
            queries=points,
            initial_guess_sample_resolutions=[n_initial_guess]
            * spline.para_dim,
            return_verbose=True,
            nthreads=1,
        )
        # Just evaluate the distance
        distances = prox_r[3]
        return np.max(distances), np.mean(distances)

    n_sample = 25
    n_resolution = 3 * n_sample
    x = np.linspace(-1, 1, n_sample)
    fitting_points = np.vstack([x, f(x)]).T.reshape(-1, 2)
    pts = gus.create.vertices.Vertices(fitting_points)
    pts.show_options["c"] = "blue"
    pts.show_options["r"] = 10

    # fit_curve creates Bspline where the n_control_points = n_fitting_points
    # if no further information is given
    interpolated_curve, _ = splinepy.helpme.fit.curve(
        fitting_points=fitting_points, degree=2
    )

    # In curve approximation the number of control points must be chosen
    n_cps = 10
    approximated_curve, _ = splinepy.helpme.fit.curve(
        fitting_points=fitting_points, degree=3, n_control_points=n_cps
    )

    interp_errs = evaluate_approximation(interpolated_curve, fitting_points)
    approx_errs = evaluate_approximation(approximated_curve, fitting_points)

    gus.show(
        [
            f"Interpolated curve \nMax. distance {interp_errs[0]:.3f}, "
            f"mean distance {interp_errs[1]:.3f}",
            interpolated_curve,
            pts,
        ],
        [
            f"Approximated curve with {n_cps} CPs \nMax. distance "
            f"{approx_errs[0]:.3f}, mean distance {approx_errs[1]:.3f}",
            approximated_curve,
            pts,
        ],
    )

    approximated_curves = []

    # The number of control points affects the approximation error
    for n_cps in range(n_sample, 5, -5):
        approximated_curve, _ = splinepy.helpme.fit.curve(
            fitting_points, degree=3, n_control_points=n_cps
        )

        max_error, mean_error = evaluate_approximation(
            approximated_curve, fitting_points
        )
        approximated_curves.append(
            [n_cps, approximated_curve, [max_error, mean_error]]
        )

    gus_list = [
        [
            f"PY: Approximated curve with {cps} CPs \n Max. "
            f"distance {errs[0]:.2f}, mean distance {errs[1]:.2f}",
            curve,
            pts,
        ]
        for (cps, curve, errs) in approximated_curves
    ]
    gus.show(*gus_list)

    # 2. Surface fitting
    n_sample = [20, 20]
    resolution = [3 * n for n in n_sample]

    # Create dataset
    x = [np.linspace(-1, 1, n) for n in n_sample]
    xx, yy = np.meshgrid(x[0], x[1])

    def h(x, y):
        return -0.9 * (x**3) + 0.7 * y**2 + x * y + 0.1 * np.sin(5 * x)

    fitting_points = np.vstack(
        (xx.flatten(), yy.flatten(), h(xx, yy).flatten())
    ).T

    interpolated_surface, _ = splinepy.helpme.fit.surface(
        fitting_points=fitting_points,
        size=n_sample,
        n_control_points=n_sample,
        degrees=[3, 2],
    )

    n_cps = [int(0.75 * n_sample[0]), int(0.6 * n_sample[1])]

    approximated_surface, _ = splinepy.helpme.fit.surface(
        fitting_points=fitting_points,
        size=n_sample,
        n_control_points=n_cps,
        degrees=[3, 2],
        centripetal=True,
    )

    n_cps_int = interpolated_surface.control_mesh_resolutions
    max_error_int, mean_error_int = evaluate_approximation(
        interpolated_surface, fitting_points
    )

    n_cps = approximated_surface.control_mesh_resolutions
    max_error, mean_error = evaluate_approximation(
        approximated_surface, fitting_points
    )

    pts = gus.create.vertices.Vertices(fitting_points)
    pts.show_options["c"] = "blue"
    pts.show_options["r"] = 10

    gus.show(
        [
            f"Interpolated surface with ({n_cps_int}) CPs\n"
            f"Max. Distance {max_error_int:.3f},"
            f"mean distance: {mean_error_int:.3f}",
            interpolated_surface,
            pts,
        ],
        [
            f"Approximated surface with ({n_cps}) CPs\n"
            f"Max. Distance {max_error:.3f},"
            f"mean distance: {mean_error:.3f}",
            approximated_surface,
            pts,
        ],
    )

    approximated_surfaces = []
    # The approximation error changes with the number of control points
    for r in range(0, min(n_sample) - 2, 5):
        approximated_surface, _ = splinepy.helpme.fit.surface(
            fitting_points=fitting_points,
            size=n_sample,
            n_control_points=[n_sample[0] - int(r / 2), n_sample[1] - r],
            degrees=[3, 2],
            knot_vectors=None,
        )

        n_cps = approximated_surface.control_mesh_resolutions

        max_error, mean_error = evaluate_approximation(
            approximated_surface, fitting_points
        )

        approximated_surfaces.append(
            [n_cps, approximated_surface, [max_error, mean_error]]
        )

    gus_list = [
        [
            f"Approximated surface with {cps} CPs\n Max. distance "
            f"{errs[0]:.3f}, mean distance {errs[1]:.3f}",
            surface,
            pts,
        ]
        for (cps, surface, errs) in approximated_surfaces
    ]
    gus.show(*gus_list)
