"""
Isophotes are a tool to visualize the continuity of Surfaces. They represent
areas of equal brightness, supposing a parallel light source. This can be used
to check interface connections between two patches. If isophotes are at least
G0 continuous, the connection is G1 continuous, further, if they are G1 smooth,
the connection is G2 continuous. This visual aid can help find kinks in
surfaces.

Further, this example shows how to compute the Mean and Gaussian Curvature of
a given surface based on the description in
    "Rogers, David F. An introduction to NURBS: with historical perspective.
     Morgan Kaufmann, 2001"
"""

import numpy as np

import splinepy as sp

# Create a random spline surface in 3D
a = sp.helpme.create.box(1, 1).bspline
a.elevate_degrees([0, 1, 0, 1])
a.insert_knots(0, [0.33333, 0.66667])
a.insert_knots(1, [0.33333, 0.66667])
a.control_points = np.hstack(
    [a.cps, np.random.random((a.cps.shape[0], 1)) * 0.5]
)


# Function to create a Callable for the SplineDataAdaptor
def isophotes_to_direction(normal=None, interval=0.02):
    if normal is None:
        normal = [1, 0, 0]
    normal = normal / np.linalg.norm(normal)

    def isophotes_scalars(data, on):
        v1 = data.derivative(on, [1, 0])
        v2 = data.derivative(on, [0, 1])
        normals = np.cross(v1, v2)
        test = np.arccos(
            np.matmul(normals, normal) / np.linalg.norm(normals, axis=1)
        )
        return np.floor(test / interval) % 2

    return isophotes_scalars


# Gaussian Curvature
def gaussian_curvature(data, on):
    v1 = data.derivative(on, [1, 0])
    v2 = data.derivative(on, [0, 1])
    v3 = data.derivative(on, [2, 0])
    v4 = data.derivative(on, [0, 2])
    v5 = data.derivative(on, [1, 1])
    # Using notation from "An introduction to NURBS historical perspective"
    normals = np.cross(v1, v2)
    A = np.einsum("ij,ij->i", normals, v3)
    B = np.einsum("ij,ij->i", normals, v4)
    C = np.einsum("ij,ij->i", normals, v5)

    k_g = (A * C - B**2) / np.linalg.norm(normals, axis=1) ** 4
    return k_g


# Mean Curvature
def mean_curvature(data, on=None):
    v1 = data.derivative(on, [1, 0])
    v2 = data.derivative(on, [0, 1])
    v3 = data.derivative(on, [2, 0])
    v4 = data.derivative(on, [0, 2])
    v5 = data.derivative(on, [1, 1])
    # Using notation from "An introduction to NURBS historical perspective"
    normals = np.cross(v1, v2)
    A = np.einsum("ij,ij->i", normals, v3)
    B = np.einsum("ij,ij->i", normals, v4)
    C = np.einsum("ij,ij->i", normals, v5)
    v1v1 = np.einsum("ij,ij->i", v1, v1)
    v1v2 = np.einsum("ij,ij->i", v1, v2)
    v2v2 = np.einsum("ij,ij->i", v2, v2)

    k_g = (A * v1v1 - 2 * v1v2 * B + C * v2v2) / (
        2 * np.linalg.norm(normals, axis=1) ** 3
    )
    return k_g


# Spline Adaptor to plot function
isophotes = sp.SplineDataAdaptor(
    a, function=isophotes_to_direction(normal=[0, 0, 1])
)
mean_curvature = sp.SplineDataAdaptor(a, function=mean_curvature)
gaussian_curvature = sp.SplineDataAdaptor(a, function=gaussian_curvature)

# Set Plot options
a.spline_data["isophotes"] = isophotes
a.spline_data["gaussian_curvature"] = gaussian_curvature
a.spline_data["mean_curvature"] = mean_curvature
a.show_options["lighting"] = "off"
a.show_options["control_points"] = False
a.show_options["knots"] = False

# Show Isophotes
a.show_options["data"] = "isophotes"
a.show(resolutions=1000, title="Isophotes")

# Show Gaussian Curvature
a.show_options["data"] = "gaussian_curvature"
a.show_options["scalarbar"] = True
a.show(resolutions=100, title="Gaussian Curvature")

# Show Curvature
a.show_options["data"] = "mean_curvature"
a.show_options["scalarbar"] = True
a.show(resolutions=100, title="Mean Curvature")
