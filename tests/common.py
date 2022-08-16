from math import degrees
import unittest

import numpy as np

## abbreviation
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
        [-1.,  0.],
        [-1.,  1.],
        [ 0.,  1.],
        [-2.,  0.],
        [-2.,  2.],
        [ 0.,  2.],
    ],
    weights=[
        [1.],
        [(2 ** -.5)],
        [1.],
        [1.],
        [(2 ** -.5)],
        [1.],
    ],
)

# 
z2P2D = dict(
    degrees=n2P2D["degrees"], 
    control_points=n2P2D["control_points"]
)

# 
r2P2D = dict(
    degrees=n2P2D["degrees"], 
    control_points=n2P2D["control_points"], 
    weights=n2P2D["weights"]
)

# query points
q2D = [
    [.01, .01], 
    [.01, .5], 
    [.9, .1], 
    [.8, .7], 
    [.4, .99],
]
