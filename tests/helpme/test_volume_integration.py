import numpy as np
import pytest

import splinepy


def test_volume_integration_1D(np_rng):
    """
    Test volume integration for splines using numerical integration of the
    Jacobi-Determinant
    """
    bezier = splinepy.Bezier(degrees=[1], control_points=[[0], [2]])

    assert np.allclose(bezier.integrate.volume(), 2.0)

    # test for other types same spline
    assert np.allclose(bezier.bspline.integrate.volume(), 2.0)
    assert np.allclose(bezier.rationalbezier.integrate.volume(), 2.0)
    assert np.allclose(bezier.nurbs.integrate.volume(), 2.0)

    # Check if equal after refinement
    bezier.elevate_degrees([0, 0, 0])

    assert np.allclose(bezier.integrate.volume(), 2.0)

    # Test knot insertion
    bspline = bezier.bspline
    bspline.insert_knots(0, np_rng.random(10))
    assert np.allclose(bezier.integrate.volume(), 2.0)


def test_volume_integration_embedded(np_rng):
    """
    Test volume integration for splines using numerical integration of the
    Jacobi-Determinant
    """
    # Test 1D -> 2D
    expected_result = 2.0**1.5
    bezier = splinepy.Bezier(degrees=[1], control_points=[[0, 0], [2, 2]])

    assert np.allclose(bezier.integrate.volume(), expected_result)

    # test for other types same spline
    assert np.allclose(bezier.bspline.integrate.volume(), expected_result)
    assert np.allclose(
        bezier.rationalbezier.integrate.volume(), expected_result
    )
    assert np.allclose(bezier.nurbs.integrate.volume(), expected_result)

    # Check if equal after refinement
    bezier.elevate_degrees([0, 0, 0])

    assert np.allclose(bezier.integrate.volume(), expected_result)

    # Test knot insertion
    bspline = bezier.bspline
    bspline.insert_knots(0, np_rng.random(10))
    assert np.allclose(bezier.integrate.volume(), expected_result)


def test_volume_integration_2D(np_rng):
    """
    Test volume integration for splines using numerical integration of the
    Jacobi-Determinant
    """
    # Test 2D
    bezier = splinepy.Bezier(
        degrees=[1, 1], control_points=[[0, 0], [2, 0], [0, 1], [2, 1]]
    )

    assert np.allclose(bezier.integrate.volume(), 2.0)

    # test for other types same spline
    assert np.allclose(bezier.bspline.integrate.volume(), 2.0)
    assert np.allclose(bezier.rationalbezier.integrate.volume(), 2.0)
    assert np.allclose(bezier.nurbs.integrate.volume(), 2.0)

    # Check if equal after refinement
    bezier.elevate_degrees([0, 0, 1])

    assert np.allclose(bezier.integrate.volume(), 2.0)

    # Test knot insertion
    bspline = bezier.bspline
    bspline.insert_knots(0, np_rng.random(10))
    assert np.allclose(bezier.integrate.volume(), 2.0)


def test_volume_integration_3D(np_rng):
    """
    Test volume integration for splines using numerical integration of the
    Jacobi-Determinant
    """
    # Test 3D
    bezier = splinepy.Bezier(
        degrees=[1, 1, 1],
        control_points=[
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [1, 1, 0],
            [0, 0, 3],
            [1, 0, 3],
            [0, 1, 3],
            [1, 1, 3],
        ],
    )

    assert np.allclose(bezier.integrate.volume(), 3.0)

    # Check if equal after refinement
    bezier.elevate_degrees([0, 0, 1, 1, 2, 2])

    assert np.allclose(bezier.integrate.volume(), 3.0)

    # Test knot insertion
    bspline = bezier.bspline
    bspline.insert_knots(0, np_rng.random(10))
    assert np.allclose(bezier.integrate.volume(), 3.0)

    # Move control points along axis does not change volume
    mi = bspline.multi_index
    bspline.cps[mi[3, :, :]][:, 1] += 10
    assert np.allclose(bezier.integrate.volume(), 3.0)


def test_complex_geometry(np_rng):
    """Test on a more complex geometry"""
    bezier = splinepy.Bezier(
        degrees=[1, 1, 1],
        control_points=[
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [1, 1, 0],
            [0, 0, 3],
            [1, 0, 3],
            [0, 1, 3],
            [1, 1, 3],
        ],
    )
    bezier.elevate_degrees([0, 1, 2, 0, 1, 2])
    bezier.cps += np_rng.random(bezier.cps.shape) * 0.1

    bezier_c = bezier.copy()

    # Move around internal control points

    mi = bezier.multi_index
    bezier.cps[mi[1:3, 1:3, 1:3]][:] += 0.05 * np_rng.random(
        bezier.cps[mi[1:3, 1:3, 1:3]][:].shape
    )
    assert np.allclose(bezier.integrate.volume(), bezier_c.integrate.volume())


def test_assertions(np_rng):
    """Test the assertions in volume function"""
    bezier = splinepy.Bezier(
        degrees=[1, 2], control_points=np_rng.random((6, 1))
    )

    with pytest.raises(
        ValueError, match=r"`Volume` not supported if para_dim > dim"
    ):
        bezier.integrate.volume()

    bezier = splinepy.Bezier(
        degrees=[2, 2], control_points=np_rng.random((9, 2))
    )

    with pytest.raises(
        ValueError, match=r"Integration order must be array of size para_dim"
    ):
        bezier.integrate.volume(orders=[2, 4, 5])


def test_function_integration(np_rng):
    col1_factor = 2

    def volume_function(x):
        vf = np.ones((x.shape[0], 2))
        # scale it with a factor to get a different value
        vf[:, 1] = col1_factor
        return vf

    bezier = splinepy.Bezier(
        degrees=[1, 2], control_points=np_rng.random((6, 2))
    )
    assert np.allclose(
        [bezier.integrate.volume(), col1_factor * bezier.integrate.volume()],
        bezier.integrate.parametric_function(volume_function),
    )
