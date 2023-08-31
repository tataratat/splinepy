import gmsh
import numpy as np


def spline_to_gmsh(spline, model=None, startid=1):
    if model is None:
        model = gmsh.model

    points, gmsh_dict = gmsh_data([spline])[0]

    for i, point in enumerate(points):
        model.occ.addPoint(*point, tag=startid + i)

    if spline.para_dim == 1:
        tag = gmsh.model.occ.addBSpline(
            pointTags=np.arange(startid, startid + points.shape[0]),
            **gmsh_dict,
        )

    elif spline.para_dim == 2:
        tag = gmsh.model.occ.addBSplineSurface(
            pointTags=np.arange(startid, startid + points.shape[0]),
            **gmsh_dict,
        )
    else:
        raise ValueError(
            f"Spline has parametric dimension {spline.para_dim} and is not supported."
        )
    return tag, np.arange(startid, startid + points.shape[0])


def gmsh_data(splines):
    gmsh_data = []
    for spline in splines:
        # Line
        if spline.para_dim == 1:
            gmsh_data.append(gmsh_line(spline))
        # Surface
        elif spline.para_dim == 2:
            gmsh_data.append(gmsh_surface(spline))
        else:
            raise ValueError(
                f"""Gmsh meshing only supports line and surface splines,
                your spline is of degree {spline.para_dim}"""
            )
    return gmsh_data


def gmsh_line(spline):
    if spline.para_dim != 1:
        raise ValueError(
            f"Spline has parametric dimension {spline.para_dim} and is not supported."
        )
    if spline.control_points.shape[1] > 3:
        raise ValueError(
            "Gmsh export supports only three-dimensional geometries."
        )

    # Embed control points in 3D
    cps = np.hstack(
        [
            spline.control_points,
            np.zeros(
                [
                    spline.control_points.shape[0],
                    3 - spline.control_points.shape[1],
                ]
            ),
        ]
    )

    # Get knov values and multiplities
    kv0, km0 = np.unique(spline.nurbs.knot_vectors[0], return_counts=True)

    properties = {
        "degree": spline.degrees[0],
        "knots": kv0,
        "multiplicities": km0,
        "weights": spline.nurbs.weights.flatten(),
    }

    return cps, properties


def gmsh_surface(spline):
    if spline.para_dim != 2:
        raise ValueError(
            f"Spline has parametric dimension {spline.para_dim} and is not supported."
        )
    if spline.control_points.shape[1] > 3:
        raise ValueError(
            "Gmsh export supports only three-dimensional geometries."
        )

    # Embed control points in 3D
    cps = np.hstack(
        [
            spline.control_points,
            np.zeros(
                [
                    spline.control_points.shape[0],
                    3 - spline.control_points.shape[1],
                ]
            ),
        ]
    )

    # Number of control points in first direction
    len(spline.nurbs.knot_vectors[0]) - spline.degrees[0] - 1

    # Get knov values and multiplities
    kv0, km0 = np.unique(spline.nurbs.knot_vectors[0], return_counts=True)
    kv1, km1 = np.unique(spline.nurbs.knot_vectors[1], return_counts=True)

    properties = {
        "numPointsU": len(spline.nurbs.knot_vectors[0])
        - spline.degrees[0]
        - 1,
        "degreeU": spline.degrees[0],
        "degreeV": spline.degrees[1],
        "knotsU": kv0,
        "knotsV": kv1,
        "multiplicitiesU": km0,
        "multiplicitiesV": km1,
        "weights": spline.nurbs.weights.flatten(),
    }

    return cps, properties
