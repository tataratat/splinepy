import numpy as np
import pytest

import splinepy


@pytest.fixture
def self_referencing_spline():
    # Spline that points to itself
    return splinepy.Bezier(
        degrees=[3, 1],
        control_points=[
            [0, 0],
            [2, 3],
            [-2, 3],
            [0, 0],
            [0, 1],
            [1, 2],
            [-1, 2],
            [0, 1],
        ],
    )


@pytest.fixture
def rotated_l_splines():
    # List of splinewith with this topology
    #
    #  o ---- o ---- o
    #  |   0  |   1  |
    #  o ---- o ---- o
    #  |   2  |
    #  o ---- o
    #
    rect_arc_1 = splinepy.Bezier(
        [2, 1], [[0, 0], [1, 0], [2, 0], [0, 1], [1, 2], [2, 1]]
    )
    rect_arc_2 = splinepy.Bezier([1, 1], [[2, 0], [3, 0], [2, 1], [3, 1]])
    rect_arc_3 = splinepy.Bezier(
        [2, 1], [[0, -1], [1, -2], [2, -1], [0, 0], [1, 0], [2, 0]]
    )
    return [
        rect_arc_1,
        rect_arc_2,
        rect_arc_3,
    ]


def test_interfaces(self_referencing_spline, rotated_l_splines):
    """ """
    # init multipatch with multiple splines
    multipatch = splinepy.Multipatch()
    multipatch.patches = rotated_l_splines

    # init multipatch with single spline
    single_p_multipatch = splinepy.Multipatch()
    single_p_multipatch.patches = [self_referencing_spline]

    # Determine connectivities
    single_p_multipatch.determine_interfaces()
    multipatch.determine_interfaces()
    assert np.all(
        single_p_multipatch.interfaces == np.array([[1, 0, -1, -1]], dtype=int)
    )

    assert np.all(
        multipatch.interfaces
        == np.array(
            # Global Face IDs and their connectivity
            # [0  1   2   3]  [4   5   6   7]   [8   9  10  11]
            [[-1, 4, 11, -1], [1, -1, -1, -1], [-1, -1, -1, 2]],
            dtype=int,
        )
    )


def test_boundaries(rotated_l_splines):
    """ """
    # init multipatch with multiple splines
    multipatch = splinepy.Multipatch()
    multipatch.patches = rotated_l_splines
    multipatch.determine_interfaces()

    # Using a function
    def west_side(points):
        return points[:, 0] < 0.1

    multipatch.boundary_from_function(west_side)
    assert np.all(
        multipatch.interfaces
        == np.array(
            # Global Face IDs and their connectivity
            # [0  1   2   3]  [4   5   6   7]   [8   9  10  11]
            [[-2, 4, 11, -1], [1, -1, -1, -1], [-2, -1, -1, 2]],
            dtype=int,
        )
    )

    multipatch.set_boundary([0, 1], [3, 3], 3)
    multipatch.set_boundary([1], [1])
    assert np.all(
        multipatch.interfaces
        == np.array(
            # Global Face IDs and their connectivity
            # [0  1   2   3]  [4   5   6   7]   [8   9  10  11]
            [[-2, 4, 11, -3], [1, -4, -1, -3], [-2, -1, -1, 2]],
            dtype=int,
        )
    )

    # Delete all boundaries and determine new ones based on continuity
    multipatch.boundaries_from_continuity()
    assert np.all(
        multipatch.interfaces
        == np.array(
            # Global Face IDs and their connectivity
            # [0  1   2   3]  [4   5   6   7]   [8   9  10  11]
            [[-1, 4, 11, -2], [1, -3, -4, -5], [-1, -6, -7, 2]],
            dtype=int,
        )
    )


def test_interfaces_and_boundaries(are_splines_equal):
    # 2 --- 3 1 --- 0
    # |  1  | |  3  |
    # 0 --- 1 3 --- 2
    # 3 --- 2 0 --- 2
    # |  2  | |  4  |
    # 1 --- 0 1 --- 3
    #
    # with 1 and three flipped in the 3rd axis

    # Create some splines in a random parametric swaps
    b1 = splinepy.Bezier(
        degrees=[1, 1, 1],
        control_points=[
            [-2, 0, 3],
            [0, 0, 3],
            [-2, 1, 3],
            [0, 1, 3],
            [-2, 0, 0],
            [0, 0, 0],
            [-2, 1, 0],
            [0, 1, 0],
        ],
    )
    b2 = splinepy.Bezier(
        degrees=[1, 1, 1],
        control_points=[
            [0, -1, 0],
            [-2, -1, 0],
            [0, 0, 0],
            [-2, 0, 0],
            [0, -1, 3],
            [-2, -1, 3],
            [0, 0, 3],
            [-2, 0, 3],
        ],
    )
    b3 = splinepy.Bezier(
        degrees=[1, 1, 1],
        control_points=[
            [1, 1, 3],
            [0, 1, 3],
            [1, 0, 3],
            [0, 0, 3],
            [1, 1, 0],
            [0, 1, 0],
            [1, 0, 0],
            [0, 0, 0],
        ],
    )
    b4 = splinepy.Bezier(
        degrees=[1, 1, 1],
        control_points=[
            [0, 0, 0],
            [0, -1, 0],
            [1, 0, 0],
            [1, -1, 0],
            [0, 0, 3],
            [0, -1, 3],
            [1, 0, 3],
            [1, -1, 3],
        ],
    )

    # Multipatch
    multipatch = splinepy.Multipatch([b1, b2, b3, b4])
    multipatch.boundaries_from_continuity()
    assert np.all(
        multipatch.interfaces
        == np.array(
            [
                [-1, 13, 9, -2, -3, -4],
                [20, -1, -5, 2, -4, -3],
                [-6, 1, -2, 18, -3, -4],
                [15, -5, 6, -6, -4, -3],
            ],
            dtype=int,
        )
    )

    # Test boundary multipatch extraction
    default_bmp = multipatch.boundary_multipatch()
    assert len(default_bmp.patches) == 16
    bmp_1 = multipatch.boundary_multipatch(1)
    assert len(bmp_1.patches) == 2
    boundary_1 = []
    for i_patch, i_face in zip(*multipatch.boundaries[0]):
        boundary_1.append(
            *multipatch.patches[i_patch].extract.boundaries([i_face])
        )
    for patch_0, patch_1 in zip(boundary_1, bmp_1.patches):
        assert are_splines_equal(patch_0, patch_1)

    assert len(multipatch.boundary_patch_ids(8)) == 0


def test_add_fields():
    """Test the evaluation and addition to field with setup checks"""
    # Create splines to form multipatch object and "solution" field
    b1 = splinepy.helpme.create.box(1, 2).bspline
    b2 = splinepy.helpme.create.box(2, 2).bspline
    b2.insert_knots(0, [0.2, 0.5])
    b2.cps[:] += [1, 0]

    # 2D field
    two_dim_field = [b1, b2]
    geometry = splinepy.Multipatch([b.create.embedded(3) for b in [b1, b2]])
    # Field dimension is 2 so should raise
    with pytest.raises(RuntimeError) as _:
        geometry.add_fields([two_dim_field], field_dim=1)
    # Success
    geometry.add_fields([two_dim_field], field_dim=2)
    geometry.add_fields([geometry.patches], field_dim=3)
