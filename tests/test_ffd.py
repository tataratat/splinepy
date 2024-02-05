from typing import List

import numpy as np
import pytest
from gustaf import Faces

from splinepy import FFD, BSpline


@pytest.fixture
def control_points_2d():
    return np.array(
        [
            [0, 0],
            [1, 0],
            [2, 0],
            [3, 0],
            [4, 0],
            [0, 1],
            [1, 2],
            [2, 1],
            [3, 2],
            [4, 1],
        ],
        dtype=np.float64,
    )


@pytest.fixture
def knot_vector_2d():
    return [[0, 0, 0, 0.3, 0.7, 1, 1, 1], [0, 0, 1, 1]]


@pytest.fixture
def degrees_2d_nu():
    return [2, 1]


@pytest.fixture
def quad_connec():
    return np.array(
        [
            [1, 0, 2, 3],
            [0, 1, 5, 4],
            [1, 3, 7, 5],
            [3, 2, 6, 7],
            [2, 0, 4, 6],
            [4, 5, 7, 6],
        ],
        dtype=np.int32,
    )


@pytest.fixture
def tri_connec():
    TF = np.array(
        [
            [1, 0, 2],
            [0, 1, 5],
            [1, 3, 7],
            [3, 2, 6],
            [2, 0, 4],
            [4, 5, 7],
            [2, 3, 1],
            [5, 4, 0],
            [7, 5, 1],
            [6, 7, 3],
            [4, 6, 2],
            [7, 6, 4],
        ],
        dtype=np.int32,
    )
    return TF


@pytest.fixture
def vertices_3d():
    V = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
        ],
        dtype=np.float64,
    )
    return V


@pytest.fixture
def faces_tri(vertices_3d, tri_connec):
    return Faces(vertices_3d, tri_connec)


@pytest.fixture
def bspline_2d(control_points_2d, degrees_2d_nu, knot_vector_2d):
    return BSpline(
        control_points=control_points_2d,
        degrees=degrees_2d_nu,
        knot_vectors=knot_vector_2d,
    )


@pytest.fixture
def control_points_3d():
    control_points = []
    for i in range(3):
        for j in range(3):
            for k in range(3):
                control_points.append([i * 0.5, j * 0.5, k * 0.5])
    return control_points


@pytest.fixture
def control_points_3d_deformed(control_points_3d):
    # change a control point so that there is a deformation
    control_points_3d[22] = [1.5, 0.5, 0.5]
    control_points_3d[2] = [0.25, 0.25, 0.75]
    return control_points_3d


@pytest.fixture
def knot_vector_3d():
    kv_3d = [
        [0, 0, 0, 4, 4, 4],
        [0, 0, 0, 1, 1, 1],
        [5, 5, 5, 10, 10, 10],
    ]
    return kv_3d


@pytest.fixture
def degrees_3d():
    d_3d = [2, 2, 2]
    return d_3d


@pytest.fixture
def faces_quad(vertices_3d, quad_connec):
    return Faces(vertices_3d, quad_connec)


@pytest.fixture
def bspline_3d(control_points_3d, knot_vector_3d, degrees_3d):
    spline_3d = BSpline(
        degrees_3d,
        knot_vector_3d,
        control_points_3d,
    )
    return spline_3d


@pytest.fixture
def bspline_3d_deformed(
    control_points_3d_deformed, knot_vector_3d, degrees_3d
):
    spline_3d = BSpline(
        degrees_3d,
        knot_vector_3d,
        control_points_3d_deformed,
    )
    return spline_3d


@pytest.mark.parametrize(
    "init_values, throw_error",
    [
        ([], False),  # test empty input
        (["jallo"], True),
        ([1], True),
        (["fixture_bspline_2d"], True),
        (["fixture_faces_tri"], False),
        (["fixture_faces_quad"], False),
        (["fixture_faces_tri", "fixture_bspline_3d"], False),
        (["fixture_faces_tri", "fixture_bspline_3d_deformed"], False),
        (["fixture_faces_tri", 1], True),
        (["fixture_faces_tri", "testing"], True),
    ],
)
def test_init_error(init_values: List, throw_error, request):
    for idx, value in enumerate(init_values):
        if isinstance(value, str) and "fixture_" in value:
            init_values[idx] = request.getfixturevalue(value[8:])
    if throw_error:
        with pytest.raises(TypeError):
            FFD(*init_values)
    else:
        FFD(*init_values)


def test_mesh_with_empty_init():
    _ = FFD()
