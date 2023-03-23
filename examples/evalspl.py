import numpy as np
import splinepy

# abbreviation
# z: bezier
# r: rational bezier
# b: bspline
# n: nurbs

# initializing a spline should be a test itself, so provide `dict_spline`
# this is "iga-book"'s fig 2.15.
b2P2D = dict(
    degrees=[2, 2],
    knot_vectors=[
        [0, 0, 0, 0.5, 1, 1, 1],
        [0, 0, 0, 1, 1, 1],
    ],
    control_points=[
        [0, 0],
        [0, 1],
        [1, 1.5],
        [3, 1.5],
        [-1, 0],
        [-1, 2],
        [1, 4],
        [3, 4],
        [-2, 0],
        [-2, 2],
        [1, 5],
        [3, 5],
    ],
)

# half-half circle.
n2P2D = dict(
    degrees=[2, 1],
    knot_vectors=[
        [0, 0, 0, 1, 1, 1],
        [0, 0, 1, 1],
    ],
    control_points=[
        [-1.0, 0.0],
        [-1.0, 1.0],
        [0.0, 1.0],
        [-2.0, 0.0],
        [-2.0, 2.0],
        [0.0, 2.0],
    ],
    weights=[
        [1.0],
        [2**-0.5],
        [1.0],
        [1.0],
        [2**-0.5],
        [1.0],
    ],
)

#
z2P2D = dict(degrees=n2P2D["degrees"], control_points=n2P2D["control_points"])

#
r2P2D = dict(
    degrees=n2P2D["degrees"],
    control_points=n2P2D["control_points"],
    weights=n2P2D["weights"],
)

# 3D
z3P3D = dict(
    degrees=[1, 1, 1],
    control_points=[
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, -1.0, 1.0],
        [1.0, 0.0, 1.0],
        [-1.0, 1.0, 2.0],
        [2.0, 2.0, 2.0],
    ],
)

r3P3D = dict(
    **z3P3D,
    weights=[1.0] * len(z3P3D["control_points"]),
)

b3P3D = dict(
    **z3P3D,
    knot_vectors=[
        [0.0, 0.0, 1.0, 1.0],
        [0.0, 0.0, 1.0, 1.0],
        [0.0, 0.0, 1.0, 1.0],
    ],
)

n3P3D = dict(
    **r3P3D,
    knot_vectors=b3P3D["knot_vectors"],
)

# query points
q2D = [
    [0.01, 0.01],
    [0.01, 0.5],
    [0.9, 0.1],
    [0.8, 0.7],
    [0.4, 0.99],
]

q3D = [
    [0.1, 0.1, 0.1],
    [0.734, 0.525, 0.143],
    [0.9666, 0.991, 0.003],
    [0.5623, 0.0089, 0.99],
    [0.0431, 0.2, 0.523],
]


b = splinepy.BSpline(**b2P2D)
n = splinepy.NURBS(**n2P2D)