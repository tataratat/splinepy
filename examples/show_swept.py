import gustaf as gus
import numpy as np

import splinepy

if __name__ == "__main__":

    ### TRAJECTORY ###

    # define a questionmark-trajectory
    dict_trajectory = {
        "degrees": [3],
        "knot_vectors": [
            [0.0, 0.0, 0.0, 0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0, 1.0, 1.0, 1.0]
        ],
        "control_points": np.array(
            [
                [0.5, 0], 
                [0.5, 2],
                [1.0, 3],
                [2.0, 4],
                [2.15, 5],
                [1.8, 5.9],
                [1.0, 6.2],
                [-0.25, 6],
                [-0.5, 5],
            ]
        ),
    }

    # create spline of trajectory dict
    trajectory = splinepy.BSpline(**dict_trajectory)

    # refine trajectory by inserting knots and control points
    trajectory.uniform_refine(0, 1)

    ### CROSS SECTION ###

    #  helpme to create a circular cross section
    cross_section = splinepy.helpme.create.surface_circle(0.5).nurbs

    # user can  define the normal vector of the cross section, in case
    # the cross section is not planar in the x-y plane (default)
    cs_nv = np.array([0, 0, 1])

    ### SWEEP ###

    # create swept surface
    swept_surface = splinepy.helpme.create.swept(
        trajectory=trajectory,
        cross_section=cross_section,
        cross_section_normal=cs_nv,
        set_on_trajectory=False,
    )

    ### VISUALIZATION ###

    trajectory.show_options["control_mesh"] = False
    cross_section.show_options["control_mesh"] = False
    swept_surface.show_options["control_mesh"] = False

    gus.show(
        ["Trajectory", trajectory],
        ["Cross Section", cross_section],
        ["Swept Surface", swept_surface],
        resolution=101,
    )
