import gustaf as gus

import splinepy as sp

if __name__ == "__main__":
    # create mesh 2d - extracted from bspline
    mesh_2d = sp.BSpline(
        [2, 2],
        [[0, 0, 0, 0.33, 0.66, 1, 1, 1], [0, 0, 0, 1, 1, 1]],
        [
            [0, 0],
            [0.25, 0],
            [0.5, 0.05],
            [0.75, 0.1],
            [1, 0.1],
            [0, 0.2],
            [0.25, 0.2],
            [0.5, 0.2],
            [0.75, 0.2],
            [1, 0.2],
            [0, 0.4],
            [0.25, 0.4],
            [0.5, 0.35],
            [0.75, 0.3],
            [1, 0.3],
        ],
    ).extract.faces(resolutions=[50, 25])

    spline_2d = sp.BSpline(
        [2, 2],
        [[0, 0, 0, 4, 4, 4], [2, 2, 2, 3, 3, 3]],
        [
            [0, 0],
            [0.5, -0.2],
            [1, 0],
            [0, 0.2],
            [0.5, 0.2],
            [1, 0.2],
            [0, 0.4],
            [0.5, 0.4],
            [1, 0.4],
        ],
    )

    ffd_2d = sp.FFD(mesh_2d, spline_2d)
    ffd_2d.show(title="2D FFD - BSpline")

    # Now 3D - again, create a mesh
    volume_3d = sp.helpme.create.torus(3, 2.5).extract.faces([30, 30, 30])

    # create spline
    spline_3d = sp.helpme.create.box(4.5, 2.5, 5.5)
    spline_3d.elevate_degrees([0, 1, 2])

    # manipulate cps
    # center cps
    cp_bounds = spline_3d.control_point_bounds
    spline_3d.control_points -= (cp_bounds[1] - cp_bounds[0]) / 2
    # use multi_index to get mid slice ids
    z_slice_ids = spline_3d.multi_index[:, :, 1]
    # rotate
    spline_3d.control_points[z_slice_ids] = gus.utils.arr.rotate(
        spline_3d.control_points[z_slice_ids], [0, 0, 79]
    )
    # translate
    spline_3d.control_points[z_slice_ids] += [0, 0, 0]

    ffd_3d = sp.FFD(volume_3d, spline_3d)

    ffd_3d.show(title="3D FFD - Bezier")

    # if you only provide mesh, we will generate a bounding spline for you.
    ffd_without_spline = sp.FFD(volume_3d)

    ffd_without_spline.show(title="Default bounding spline.")
    # then, you can just access spline and apply deformation
    ffd_without_spline.spline.control_points[0] -= 3
    ffd_without_spline.show(title="Manipulated default spline")

    # Only provide spline and then mesh after initialization
    # this will put spline
    ffd_with_out_mesh = sp.FFD()
    ffd_with_out_mesh.spline = spline_3d
    ffd_with_out_mesh.mesh = volume_3d
    ffd_with_out_mesh.show(title="Setting spline first, then mesh")
