import gustaf as gus
import numpy as np

import splinepy

if __name__ == "__main__":
    # trajectory
    # define degrees
    ds1 = [2]
    # define knot vectors
    kvs1 = [[0.0, 0.0, 0.0, 1.0, 1.0, 1.0]]
    # define control points
    cps1 = np.array(
        [
            [0.0, 0.0, 0.0],
            [3.0, 1, 0.0],
            [4.0, 3, 0.0],
        ]
    )
    # init trajectory as bspline
    trajectory = splinepy.BSpline(
        degrees=ds1,
        knot_vectors=kvs1,
        control_points=cps1,
    )

    # define sections along trajectory
    nsect = len(kvs1[0]) - ds1[0] - 1

    # cross section
    # define degrees
    ds_cs = [3]

    # define knot vectors
    kvs_cs = [[0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]]

    # define control points
    cps_cs = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 1.0, 0.0],
        ]
    )

    # init cross section as bspline
    cross_section = splinepy.BSpline(
        degrees=ds_cs,
        knot_vectors=kvs_cs,
        control_points=cps_cs,
    )

    # user should define the normal vector of the cross section
    cs_nv = np.array([0, 0, 1])

    # create interpolated surface
    interpolated_surface = splinepy.helpme.create.swept(
        trajectory=trajectory,
        cross_section=cross_section,
        nsections=nsect,
        cross_section_normal=cs_nv,
    )

    # testing derivative and evaluation
    # derivative = trajectory.derivative([[0.5], [1]], [1])
    # evals = trajectory.evaluate([[0.5]])

    # show
    def der(data, on):
        return data.derivative(on, [1])

    trajectory.spline_data["der"] = splinepy.SplineDataAdaptor(
        trajectory, function=der
    )
    trajectory.show_options["arrow_data"] = "der"
    trajectory.show_options["arrow_data_on"] = np.linspace(0, 1, 20).reshape(
        -1, 1
    )
    trajectory.show_options["control_mesh"] = False
    cross_section.show_options["control_mesh"] = False
    interpolated_surface.show_options["control_mesh"] = False

    gus.show(
        ["Trajectory", trajectory],
        ["Cross section", cross_section],
        resolution=50,
    )

    gus.show(
        ["Surface", interpolated_surface],
        resolution=50,
    )
