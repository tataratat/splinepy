import sys

import gustaf as gus
import numpy as np

import splinepy

if __name__ == "__main__":

    ### TRAJECTORY ###

    # arbitrary trajectory
    # dict_trajectory = {
    #     "degrees": [2],
    #     "knot_vectors": [[0.0, 0.0, 0.0, 0.333, 0.666, 1.0, 1.0, 1.0]],
    #     "control_points": np.array(
    #         [
    #             [0.0, 0.0, 0.0],
    #             [0.0, 0.0, 5.0],
    #             [10.0, 5.0, 0.0],
    #             [15.0, 0.0, -5.0],
    #             [20.0, 0.0, 0.0],
    #         ]
    #     ),
    # }

    # 2D questionmark
    # dict_trajectory = {
    #     "degrees": [3],
    #     "knot_vectors": [[0.0, 0.0, 0.0, 0.0, 0.2, 0.4,
    #                   0.6, 0.8, 0.9, 1.0, 1.0, 1.0, 1.0]],
    #     "control_points": np.array([
    #             [0.5, 0],  # Startpunkt
    #             [0.5, 2],
    #             [1.0, 3],
    #             [2.0, 4],
    #             [2.15, 5],
    #             [1.8, 5.9],
    #             [1.0, 6.2],
    #             [-0.25, 6],
    #             [-0.5, 5],])
    # }
    # init trajectory as bspline

    # closed 3D questionmark
    dict_trajectory = {
        "degrees": [3],
        "knot_vectors": [
            [
                0.0,
                0.0,
                0.0,
                0.0,
                0.1,
                0.2,
                0.3,
                0.4,
                0.5,
                0.6,
                0.8,
                0.9,
                1.0,
                1.0,
                1.0,
                1.0,
            ]
        ],
        "control_points": np.array(
            [
                [0.5, 0, 0],
                [0.5, 2, 0.3],
                [1.0, 3, 0.1],
                [2.0, 4, -0.1],
                [2.15, 5, -0.2],
                [1.8, 5.9, -0.4],
                [1.0, 6.2, -0.3],
                [-0.25, 6, -0.1],
                [-0.5, 5.0, 0.1],
                [-2.0, 4.0, 0.2],
                [-1, 3.0, 0.1],
                [0.5, 0.0, 0.0],
            ]
        ),
    }
    trajectory = splinepy.BSpline(**dict_trajectory)

    # alternatively, use helpme to create a trajectory
    # trajectory = splinepy.helpme.create.circle(10)

    # insert knots and control points
    trajectory.uniform_refine(0, 1)

    ### CROSS SECTION ###
    dict_cross_section = {
        "degrees": [3],
        "knot_vectors": [[0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0]],
        "control_points": np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 2.0, 0.0],
                [2.0, 0.0, 0.0],
                [3.0, -2.0, 0.0],
                [4.0, 0.0, 0.0],
            ]
        ),
    }

    # init cross section as bspline
    cross_section = splinepy.BSpline(**dict_cross_section)

    # alternatively, use helpme to create a cross section
    # cross_section = splinepy.helpme.create.surface_circle(0.5).nurbs

    # user can  define the normal vector of the cross section, in case
    # the cross section is not planar in the x-y plane (default)
    cs_nv = np.array([0, 0, 1])

    ### SWEEP ###
    swept_surface = splinepy.helpme.create.swept(
        trajectory=trajectory,
        cross_section=cross_section,
        cross_section_normal=cs_nv,
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

    sys.exit()

    ### EXPORT A SWEPT SPLINE ###
    dict_export_cs = {
        "degrees": [1],
        "knot_vectors": [[0.0, 0.0, 1.0, 1.0]],
        "control_points": np.array(
            [
                [0.0, 0.0],
                [1.0, 0.0],
            ]
        ),
    }
    export_cs = splinepy.BSpline(**dict_export_cs)

    dict_export_traj = {
        "degrees": [3],
        "knot_vectors": [[0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0]],
        "control_points": np.array(
            [
                [0.0, 0.0],
                [1.0, 2.0],
                [2.0, 0.0],
                [3.0, -2.0],
                [4.0, 0.0],
            ]
        ),
    }
    export_traj = splinepy.BSpline(**dict_export_traj)

    swept_surface = splinepy.helpme.create.swept(
        trajectory=export_traj,
        cross_section=export_cs,
        cross_section_normal=[-1, 0, 0],
    )

    gus.show(
        ["Trajectory", export_traj],
        ["Cross Section", export_cs],
        ["Swept Surface", swept_surface],
        resolution=101,
    )

    projection = swept_surface.create.embedded(2)
    gus.show(["Projection", projection], resolution=101)
    splinepy.io.mfem.export("testmeshmesh.mesh", projection)
