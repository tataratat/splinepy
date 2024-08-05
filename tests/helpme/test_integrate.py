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
