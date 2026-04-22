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

    # 1. create a circular 1D-line-cross-section
    cross_section_circle = splinepy.helpme.create.circle(0.5).nurbs

    # 2. create a circular 2D-surface-cross-section
    cross_section_disk = splinepy.helpme.create.surface_circle(0.5).nurbs

    # 3. create a rectangular 2D-surface-cross-section
    cross_section_plate = splinepy.helpme.create.box(1, 1).nurbs

    # 4. create a more complex 1D-line-cross-section
    # use an asymmetric rational curve so control-box and geometry-box anchors
    # visibly differ
    cross_section_rectangle_line = splinepy.helpme.create.circle(0.55).nurbs
    cross_section_rectangle_line.control_points[:, 0] *= 1.2
    cross_section_rectangle_line.control_points[:, 1] *= 0.8
    cross_section_rectangle_line.control_points[1] += np.array([0.9, -0.15])
    cross_section_rectangle_line.control_points[2] += np.array([0.45, 0.55])
    cross_section_rectangle_line.control_points[5] += np.array([-0.35, -0.45])
    cross_section_rectangle_line.weights[1] *= 0.25
    cross_section_rectangle_line.weights[2] *= 0.45
    cross_section_rectangle_line.weights[5] *= 0.5

    # Define a custom normal vector for the cross-section when:
    # a) The cross-section is not planar in the x-y plane, or
    # b) Crooked sweeping is desired
    cs_nv = np.array([1, 0, 1])

    ### SWEEP ###

    # create swept surface
    swept_surface_circle = splinepy.helpme.create.swept(
        trajectory=trajectory,
        cross_section=cross_section_circle,
        cross_section_normal=None,
        anchor="auto",
        set_on_trajectory=False,
        rotation_adaption=None,
    )

    # create crooked swept solid (circular nr. 1)
    # the cross-section's normal vector is not default; crooked sweeping
    swept_solid_disk_1 = splinepy.helpme.create.swept(
        trajectory=trajectory,
        cross_section=cross_section_disk,
        cross_section_normal=cs_nv,
        anchor="auto",
        set_on_trajectory=False,
        rotation_adaption=None,
    )

    # create swept solid (circular nr. 2)
    # the cross-sections are set on the trajectory's control points (default)
    swept_solid_disk_2 = splinepy.helpme.create.swept(
        trajectory=trajectory,
        cross_section=cross_section_disk,
        cross_section_normal=None,
        anchor="auto",
        set_on_trajectory=False,
        rotation_adaption=None,
    )

    # create swept solid (circular nr. 3)
    # the cross-sections are set on the trajectory's evaluation points
    swept_solid_disk_3 = splinepy.helpme.create.swept(
        trajectory=trajectory,
        cross_section=cross_section_disk,
        cross_section_normal=None,
        anchor="auto",
        set_on_trajectory=True,
        rotation_adaption=None,
    )

    # create swept solid (rectangular nr. 1)
    swept_solid_plate_1 = splinepy.helpme.create.swept(
        trajectory=trajectory,
        cross_section=cross_section_plate,
        cross_section_normal=None,
        anchor="auto",
        set_on_trajectory=False,
        rotation_adaption=None,
    )

    # create swept solid (rectangular nr. 2)
    # rotation adaption with 45 degrees
    swept_solid_plate_2 = splinepy.helpme.create.swept(
        trajectory=trajectory,
        cross_section=cross_section_plate,
        cross_section_normal=None,
        anchor="auto",
        set_on_trajectory=False,
        rotation_adaption=45 * np.pi / 180,
    )

    # create swept surfaces to compare anchor placement strategies
    swept_surface_rectangle_parametric = splinepy.helpme.create.swept(
        trajectory=trajectory,
        cross_section=cross_section_rectangle_line,
        cross_section_normal=None,
        anchor="parametric",
        set_on_trajectory=False,
        rotation_adaption=None,
    )

    swept_surface_rectangle_control_box = splinepy.helpme.create.swept(
        trajectory=trajectory,
        cross_section=cross_section_rectangle_line,
        cross_section_normal=None,
        anchor="control_box",
        set_on_trajectory=False,
        rotation_adaption=None,
    )

    swept_surface_rectangle_geometry_box = splinepy.helpme.create.swept(
        trajectory=trajectory,
        cross_section=cross_section_rectangle_line,
        cross_section_normal=None,
        anchor="geometry_box",
        set_on_trajectory=False,
        rotation_adaption=None,
    )

    swept_surface_rectangle_auto = splinepy.helpme.create.swept(
        trajectory=trajectory,
        cross_section=cross_section_rectangle_line,
        cross_section_normal=None,
        anchor="auto",
        set_on_trajectory=False,
        rotation_adaption=None,
    )

    ### VISUALIZATION ###

    # first window: swept surface
    gus.show(
        ["Trajectory", trajectory],
        ["1D Cross Section", cross_section_circle],
        ["Swept Surface", swept_surface_circle],
        resolution=51,
        control_mesh=False,
        control_point_ids=False,
    )

    # adjust show options
    swept_solid_disk_2.show_options["alpha"] = 0.3
    swept_solid_disk_3.show_options["alpha"] = 0.3
    trajectory.show_options["control_points"] = False

    # second window: swept solids (circular)
    gus.show(
        ["2D Cross Section", cross_section_disk],
        ["Swept Solid - Crooked Sweeping", swept_solid_disk_1],
        [
            "Swept Solid - Set on Control Points",
            swept_solid_disk_2,
            trajectory,
        ],
        [
            "Swept Solid - Set on Evaluation Points",
            swept_solid_disk_3,
            trajectory,
        ],
        resolution=51,
        control_mesh=False,
        control_point_ids=False,
    )

    # third window: swept solids (rectangular)
    gus.show(
        ["New Cross Section", cross_section_plate],
        ["Swept Solid without Rotation", swept_solid_plate_1],
        ["Swept Solid with 45° Rotation", swept_solid_plate_2],
        resolution=51,
        control_mesh=False,
        control_point_ids=False,
    )

    # fourth window: anchor placement for line cross-sections
    gus.show(
        ["Rectangle Line Cross Section", cross_section_rectangle_line],
        [
            "Anchor: parametric",
            swept_surface_rectangle_parametric,
            trajectory,
        ],
        [
            "Anchor: control_box",
            swept_surface_rectangle_control_box,
            trajectory,
        ],
        [
            "Anchor: geometry_box",
            swept_surface_rectangle_geometry_box,
            trajectory,
        ],
        [
            "Anchor: auto",
            swept_surface_rectangle_auto,
            trajectory,
        ],
        resolution=51,
        control_mesh=False,
        control_point_ids=False,
    )
