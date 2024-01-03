import time

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
    target_points = np.vstack([x, f(x)]).T.reshape(-1, 2)

    # fit_curve creates Bspline where the n_control_points = n_target_points
    # if no further information is given
    interpolated_curve_py, _ = splinepy.helpme.fit.fit_curve(
        points=target_points, degree=2
    )
    interpolated_curve = splinepy.BSpline.interpolate_curve(target_points, 2)

    # In curve approximation the number of control points must be chosen
    n_cps = 10
    approximated_curve_py, _ = splinepy.helpme.fit.fit_curve(
        points=target_points, degree=3, n_control_points=n_cps
    )
    approximated_curve = splinepy.BSpline.approximate_curve(
        target_points, 3, n_cps
    )
    interp_errs = evaluate_approximation(interpolated_curve, target_points)
    approx_errs = evaluate_approximation(approximated_curve, target_points)
    interp_errs_py = evaluate_approximation(
        interpolated_curve_py, target_points
    )
    approx_errs_py = evaluate_approximation(
        approximated_curve_py, target_points
    )

    gus.show(
        [
            f"C++: Interpolated curve \nMax. distance {interp_errs[0]:.5f}, "
            f"mean distance {interp_errs[1]:.5f}",
            interpolated_curve,
        ],
        [
            f"C++: Approximated curve with {n_cps} CPs \nMax. distance "
            f"{approx_errs[0]:.5f}, mean distance {approx_errs[1]:.5f}",
            approximated_curve,
        ],
        [
            f"PY: Interpolated curve \nMax. distance {interp_errs_py[0]:.5f},"
            f" mean distance {interp_errs_py[1]:.5f}",
            interpolated_curve_py,
        ],
        [
            f"PY: Approximated curve with {n_cps} CPs \nMax. distance "
            f"{approx_errs_py[0]:.5f}, mean distance {approx_errs_py[1]:.5f}",
            approximated_curve_py,
        ],
    )

    approximated_curves = []
    approximated_curves_cpp = []

    # The number of control points affects the approximation error
    for n_cps in range(n_sample, 5, -5):
        approximated_curve, _ = splinepy.helpme.fit.fit_curve(
            target_points, degree=3, n_control_points=n_cps
        )

        max_error, mean_error = evaluate_approximation(
            approximated_curve, target_points
        )
        approximated_curves.append(
            [n_cps, approximated_curve, [max_error, mean_error]]
        )

        approximated_curve_cpp = splinepy.BSpline.approximate_curve(
            target_points, 3, n_cps
        )

        max_error_cpp, mean_error_cpp = evaluate_approximation(
            approximated_curve_cpp, target_points
        )
        approximated_curves_cpp.append(
            [n_cps, approximated_curve_cpp, [max_error_cpp, mean_error_cpp]]
        )
    gus_list_cpp = [
        [
            f"C++: Approximated curve with {cps} CPs \n Max. "
            f"distance {errs[0]:.5f}, mean distance {errs[1]:.5f}",
            curve,
        ]
        for (cps, curve, errs) in approximated_curves_cpp
    ]
    gus_list = [
        [
            f"PY: Approximated curve with {cps} CPs \n Max. "
            f"distance {errs[0]:.5f}, mean distance {errs[1]:.5f}",
            curve,
        ]
        for (cps, curve, errs) in approximated_curves
    ]
    gus.show(*gus_list_cpp, *gus_list)

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
    t = time.time()
    interpolated_surface = splinepy.BSpline.interpolate_surface(
        target_points,
        size_u=n_sample[0],
        size_v=n_sample[1],
        degree_u=3,
        degree_v=2,
    )
    elapsed = time.time() - t
    print("C++: ", elapsed)

    t = time.time()
    interpolated_surface_py, _ = splinepy.helpme.fit.fit_surface(
        points=target_points,
        size=[n_sample[0], n_sample[1]],
        knot_vectors=interpolated_surface.knot_vectors,
        n_control_points=[n_sample[0], n_sample[1]],
        degrees=[3, 2],
    )
    elapsed = time.time() - t
    print("Python: ", elapsed)
    print(np.allclose(interpolated_surface.cps, interpolated_surface_py.cps))

    n_cps = [int(0.75 * n_sample[0]), int(0.6 * n_sample[1])]

    t = time.time()
    approximated_surface = splinepy.BSpline.approximate_surface(
        target_points,
        num_points_u=n_sample[0],
        num_points_v=n_sample[1],
        size_u=n_cps[0],
        size_v=n_cps[1],
        degree_u=3,
        degree_v=2,
    )
    elapsed = time.time() - t
    print("C++: ", elapsed)

    t = time.time()
    approximated_surface_py, _ = splinepy.helpme.fit.fit_surface(
        points=target_points,
        size=[n_sample[0], n_sample[1]],
        n_control_points=[n_cps[0], n_cps[1]],
        degrees=[3, 2],
        centripetal=True,
    )
    elapsed = time.time() - t
    print("Python: ", elapsed)

    n_cps_int = interpolated_surface.control_mesh_resolutions
    max_error_int, mean_error_int = evaluate_approximation(
        interpolated_surface, target_points
    )
    n_cps_int_py = interpolated_surface_py.control_mesh_resolutions
    max_error_int_py, mean_error_int_py = evaluate_approximation(
        interpolated_surface_py, target_points
    )
    n_cps = approximated_surface.control_mesh_resolutions
    max_error, mean_error = evaluate_approximation(
        approximated_surface, target_points
    )
    n_cps_py = approximated_surface_py.control_mesh_resolutions
    max_error_py, mean_error_py = evaluate_approximation(
        approximated_surface_py, target_points
    )

    pts = gus.create.vertices.Vertices(target_points)
    pts.show_options["c"] = "blue"
    pts.show_options["r"] = 10
    interpolated_surface_py.show_options["c"] = "red"
    approximated_surface_py.show_options["c"] = "red"
    gus.show(
        [
            f"C++: Interpolated surface with ({n_cps_int}) CPs\n"
            f"Max. Distance {max_error_int:.5f},\n"
            f"Mean Distance: {mean_error_int:.5f}",
            interpolated_surface,
            pts,
        ],
        [
            f"C++: Approximated surface with ({n_cps}) CPs\n"
            f"Max. Distance {max_error:.5f},\n"
            f"Mean Distance: {mean_error:.5f}",
            approximated_surface,
            pts,
        ],
        [
            f"PY: Interpolated surface with ({n_cps_int_py}) CPs\n"
            f"Max. Distance {max_error_int_py:.5f},\n"
            f"Mean Distance: {mean_error_int_py:.5f}",
            interpolated_surface,
            interpolated_surface_py,
            pts,
        ],
        [
            f"PY: Approximated surface with ({n_cps_py}) CPs\n"
            f"Max. Distance {max_error_py:.5f},\n"
            f"Mean Distance: {mean_error_py:.5f}",
            approximated_surface,
            approximated_surface_py,
            pts,
        ],
    )
    print(
        "same cps?\n",
        np.allclose(approximated_surface.cps, approximated_surface_py.cps),
    )

    approximated_surfaces = []
    approximated_surfaces_py = []
    # The approximation error changes with the number of control points
    for r in range(0, min(n_sample) - 2, 5):
        t = time.time()
        approximated_surface = splinepy.BSpline.approximate_surface(
            target_points,
            num_points_u=n_sample[0],
            num_points_v=n_sample[1],
            size_u=n_sample[0] - int(r / 2),
            size_v=n_sample[1] - r,
            degree_u=3,
            degree_v=2,
        )
        elapsed = time.time() - t
        print("C++: ", elapsed)

        n_cps = approximated_surface.control_mesh_resolutions
        max_error, mean_error = evaluate_approximation(
            approximated_surface, target_points
        )
        t = time.time()
        approximated_surface_py, _ = splinepy.helpme.fit.fit_surface(
            points=target_points,
            size=[n_sample[0], n_sample[1]],
            n_control_points=[n_sample[0] - int(r / 2), n_sample[1] - r],
            degrees=[3, 2],
            knot_vectors=[None, None],
        )
        elapsed = time.time() - t
        print("Python: ", elapsed)
        n_cps_py = approximated_surface_py.control_mesh_resolutions
        max_error_py, mean_error_py = evaluate_approximation(
            approximated_surface_py, target_points
        )

        approximated_surfaces.append(
            [n_cps, approximated_surface, [max_error, mean_error]]
        )
        approximated_surfaces_py.append(
            [n_cps_py, approximated_surface_py, [max_error_py, mean_error_py]]
        )

        print(
            "same cps?\n",
            np.allclose(approximated_surface.cps, approximated_surface_py.cps),
        )

    gus_list = [
        [
            f"C++: Approximated surface with {cps} CPs\n Max. distance "
            f"{errs[0]:.5f}, mean distance {errs[1]:.5f}",
            surface,
        ]
        for (cps, surface, errs) in approximated_surfaces
    ]
    gus_list_py = [
        [
            f"PY: Approximated surface with {cps} CPs\n Max. distance "
            f"{errs[0]:.5f}, mean distance {errs[1]:.5f}",
            surface,
        ]
        for (cps, surface, errs) in approximated_surfaces_py
    ]
    gus.show(*gus_list, *gus_list_py)
