"""Test Orientation between adjacent spline patches
         2--3
         |C |
         0--1
    0--1 2--3 1--3
    |D | |A | | B|
    2--3 0--1 0--2
         1--0
    |    |E |
    O-   3--2
"""

import numpy as np

import splinepy


def test_orientation():
    # Init Splines to be tested
    a_s = splinepy.Bezier(
        degrees=[1, 1], control_points=[[1, 1], [2, 1], [1, 2], [2, 2]]
    )
    b_s = splinepy.Bezier(
        degrees=[1, 1], control_points=[[2, 1], [2, 2], [3, 1], [3, 2]]
    )
    c_s = splinepy.Bezier(
        degrees=[1, 1], control_points=[[1, 2], [2, 2], [1, 3], [2, 3]]
    )
    d_s = splinepy.Bezier(
        degrees=[1, 1], control_points=[[0, 2], [1, 2], [0, 1], [1, 1]]
    )
    e_s = splinepy.Bezier(
        degrees=[1, 1], control_points=[[2, 1], [1, 1], [2, 0], [1, 0]]
    )
    # Expected Values
    expected_interfaces = [
        [13, 6, 18, 10],
        [-1, -1, 1, -1],
        [-1, -1, 3, -1],
        [-1, 0, -1, -1],
        [-1, -1, 2, -1],
    ]

    expected_orientations = [
        [0, 0, 3, 1, 0, 1, 1, 0],
        [0, 1, 1, 2, 1, 0, 1, 1],
        [0, 2, 4, 2, 0, 1, 0, 0],
        [0, 3, 2, 2, 0, 1, 1, 1],
    ]

    # Provide connectivity data
    mp = splinepy.Multipatch([a_s, b_s, c_s, d_s, e_s])

    # First compare interfaces to what is expected
    interfaces = mp.interfaces
    assert np.all(interfaces == expected_interfaces)

    # Then check orientations
    orientations = mp.interface_orientations()
    assert np.all(orientations == expected_orientations)
