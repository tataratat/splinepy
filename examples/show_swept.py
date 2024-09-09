import gustaf as gus
import numpy as np

import splinepy

if __name__ == "__main__":

    ### TRAJECTORY ###

    # define a hook-trajectory
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

    ### CROSS SECTIONS ###

    # 1. create a circular line-cross-section
    cross_section_circle1D = splinepy.helpme.create.circle(0.5).nurbs

    # 2. create a circular surface-cross-section
    cross_section_circle2D = splinepy.helpme.create.surface_circle(0.5).nurbs

    # 3. create a rectangular surface-cross-section
    cross_section_rect2D = splinepy.helpme.create.box(1, 1).nurbs

    # user can  define the normal vector of the cross section, in case
    # the cross section is not planar in the x-y plane (or the user
    # wants crooked sweeping)
    cs_nv = np.array([1, 0, 1])

    ### SWEEP ###

    # create swept surface
    swept_surface_1D_circle = splinepy.helpme.create.swept(
        trajectory=trajectory,
        cross_section=cross_section_circle1D,
        cross_section_normal=None,
        set_on_trajectory=False,
        rotation_adaption=None,
    )

    # create swept solid (circular nr. 1)
    # the cross-sections are set on the trajectory's control points (default)
    swept_surface_2D_circle_1 = splinepy.helpme.create.swept(
        trajectory=trajectory,
        cross_section=cross_section_circle2D,
        cross_section_normal=None,
        set_on_trajectory=False,
        rotation_adaption=None,
    )

    # create crooked swept solid (circular nr. 2)
    # the cross-sections are set on the trajectory's evaluation points
    swept_surface_2D_circle_2 = splinepy.helpme.create.swept(
        trajectory=trajectory,
        cross_section=cross_section_circle2D,
        cross_section_normal=None,
        set_on_trajectory=True,
        rotation_adaption=None,
    )

    # create swept solid set on trajectory (circular nr. 3)
    # the cross-section's normal vector is not default; crooked sweeping
    swept_surface_2D_circle_3 = splinepy.helpme.create.swept(
        trajectory=trajectory,
        cross_section=cross_section_circle2D,
        cross_section_normal=cs_nv,
        set_on_trajectory=False,
        rotation_adaption=None,
    )

    # create swept solid (rectangular nr. 1)
    swept_surface_2D_rect_1 = splinepy.helpme.create.swept(
        trajectory=trajectory,
        cross_section=cross_section_rect2D,
        cross_section_normal=None,
        set_on_trajectory=False,
        rotation_adaption=None,
    )

    # create swept solid (rectangular nr. 2)
    # rotation adaption with 45 degrees
    swept_surface_2D_rect_2 = splinepy.helpme.create.swept(
        trajectory=trajectory,
        cross_section=cross_section_rect2D,
        cross_section_normal=None,
        set_on_trajectory=False,
        rotation_adaption=45 * np.pi / 180,
    )

    ### VISUALIZATION ###

    gus.show(
        ["Trajectory", trajectory],
        ["1D Cross Section", cross_section_circle1D],
        ["Swept Surface", swept_surface_1D_circle],
        resolution=101,
        control_mesh=False,
    )

    gus.show(
        ["2D Cross Section", cross_section_circle2D],
        ["Swept Solid - Set on Control Points", swept_surface_2D_circle_1],
        ["Swept Solid - Set on Evaluation Points", swept_surface_2D_circle_2],
        ["Swept Solid - Crooked Sweeping", swept_surface_2D_circle_3],
        resolution=101,
        control_mesh=False,
    )

    gus.show(
        ["New Cross Section", cross_section_rect2D],
        ["Swept Solid without Rotation", swept_surface_2D_rect_1, trajectory],
        ["Swept Solid with 45° Rotation", swept_surface_2D_rect_2, trajectory],
        resolution=101,
        control_mesh=False,
    )
