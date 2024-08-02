import numpy as np

import splinepy


def test_transformation_class():
    """Test element transformation of single patch"""
    # Create quadratic spline
    spline = splinepy.BSpline(
        degrees=[2, 2],
        knot_vectors=[[0, 0, 0, 1, 1, 1], [0, 0, 0, 1, 1, 1]],
        control_points=splinepy.utils.data.cartesian_product(
            [np.linspace(0, 1, 3), np.linspace(0, 1, 3)]
        ),
    )

    # Randomly insert knots along both parametric dimensions
    spline.insert_knots(0, np.random.random(2))
    spline.insert_knots(1, np.random.random(2))

    # Create transformation class for spline
    trafo = splinepy.helpme.integrate.Transformation(spline)

    # Check whether element quadrature points lie inside element
    ukvs = spline.unique_knots
    trafo.compute_all_element_quad_points()
    # Check quadrature points of each element
    for element_id, quad_points in enumerate(trafo.all_quad_points):
        grid_id = trafo.get_element_grid_id(element_id)
        # Extract the corners of the current element
        element_corners = [
            ukv[grid_dim_id : grid_dim_id + 2]
            for ukv, grid_dim_id in zip(ukvs, grid_id)
        ]
        # Check if quadrature points lie within element
        for dim, corners in enumerate(element_corners):
            assert np.all(
                (quad_points[:, dim] > corners[0])
                & (quad_points[:, dim] < corners[1])
            ), f"Quadrature points do not lie within element for dimension {dim}"

    # For given spline, all Jacobians should be identity matrix
    eye = np.eye(spline.para_dim)
    trafo.compute_all_element_jacobian_inverses()
    for element_jacobians in trafo.all_jacobians:
        for jacobian_at_quad_point in element_jacobians:
            assert np.allclose(
                eye, jacobian_at_quad_point
            ), "All Jacobians should be identity matrix"

    # For created spline, all determinants should equal one
    trafo.compute_all_element_jacobian_determinants()
    assert np.allclose(
        trafo.all_jacobian_determinants,
        np.ones_like(trafo.all_jacobian_determinants),
    ), "All Jacobians' determinants should be equal to one"
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
def test_transformation_class():
    """Test element transformation of single patch"""
    # Create quadratic spline
    spline = splinepy.BSpline(
        degrees=[2, 2],
        knot_vectors=[[0, 0, 0, 1, 1, 1], [0, 0, 0, 1, 1, 1]],
        control_points=splinepy.utils.data.cartesian_product(
            [np.linspace(0, 1, 3), np.linspace(0, 1, 3)]
        ),
    )

    # Randomly insert knots along both parametric dimensions
    spline.insert_knots(0, np.random.random(2))
    spline.insert_knots(1, np.random.random(2))

    # Create transformation class for spline
    trafo = splinepy.helpme.integrate.Transformation(spline)

    # Check whether element quadrature points lie inside element
    ukvs = spline.unique_knots
    trafo.compute_all_element_quad_points()
    for element_id, quad_points in enumerate(trafo.all_quad_points):
        grid_id = trafo.get_element_grid_id(element_id)
        element_corners = [
            ukv[grid_dim_id : grid_dim_id + 2]
            for ukv, grid_dim_id in zip(ukvs, grid_id)
        ]
        # Check if quadrature points lie within element corners
        for dim, corners in enumerate(element_corners):
            assert np.all(
                (quad_points[:, dim] > corners[0])
                & (quad_points[:, dim] < corners[1])
            )

    # For given spline, all Jacobians should be identity matrix
    eye = np.eye(spline.para_dim)
    trafo.compute_all_element_jacobian_inverses()
    for element_jacobians in trafo.all_jacobians:
        for jacobian_at_quad_point in element_jacobians:
            assert np.allclose(eye, jacobian_at_quad_point)

    # For created spline, all determinants should equal one
    trafo.compute_all_element_jacobian_determinants()
    assert np.allclose(
        trafo.all_jacobian_determinants,
        np.ones_like(trafo.all_jacobian_determinants),
    )
