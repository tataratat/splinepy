import gmsh

import splinepy

try:
    import gustaf as gus

    gus_exists = True
except ModuleNotFoundError:
    gus_exists = False


if __name__ == "__main__":
    """Splines are a great parametrization for easy and complex shapes.
    Isogeometric analysis allows to calculate direct on spline parametrization.
    However, not every simulation is equipped with the options to use
    splines for analysis. Thus, it can be necessary to create computational
    meshes from splines. In this example, we demonstrate how to use gmsh
    for meshing splinepy splines.
    """

    ## For simple shapes

    # Init gmsh
    gmsh.initialize()

    # Create two splines
    disk = splinepy.helpme.create.disk(1, angle=180)
    box = splinepy.helpme.create.box(1, 1.2)

    # Add both splines to gmsh, it is necessary to use distinct point IDs
    dtag, point_tags = splinepy.io.gmsh.spline_to_gmsh(disk)
    ltag, lpoint_tags = splinepy.io.gmsh.spline_to_gmsh(
        box, startid=max(point_tags) + 5
    )

    # Generate mesh
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(2)
    gmsh.write("splines.brep")
    gmsh.write("splines.msh")

    gmsh.finalize()

    if gus_exists:
        gus.show(
            ["Splines", disk, box],
            ["Mesh", gus.io.load("splines.msh").to_edges()],
        )

    ## For complex geometries

    # In more complex geometry scripts, it might be necessary
    # to have more control about the associated point numbers and entities.
    # For this usecase we provide the data and recommend the user to
    # add it manually to the geometry

    gmsh.initialize()
    new_model = gmsh.model

    # Create line
    line = splinepy.helpme.create.line(
        [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0.1, 0.8, 0], [0, 0.0, 0.0]]
    )
    line.elevate_degrees(0)

    # Add line to gmsh
    gmsh_points, gmsh_dict = splinepy.io.gmsh.gmsh_data([line])[0]

    tags = list(range(100, 100 + len(gmsh_points)))
    for tag, point in zip(tags, gmsh_points):
        new_model.occ.addPoint(*point, tag=tag)
    linetag = new_model.occ.addBSpline(tags, **gmsh_dict)

    # Create curve line
    looptag = new_model.occ.addCurveLoop([linetag])
    # Create surface filling
    new_model.occ.addSurfaceFilling(looptag)

    new_model.occ.synchronize()
    new_model.mesh.generate(2)

    # Write files
    gmsh.write("line.brep")
    gmsh.write("line.msh")

    if gus_exists:
        gus.show(
            ["Spline", line], ["Mesh", gus.io.load("line.msh").to_edges()]
        )
