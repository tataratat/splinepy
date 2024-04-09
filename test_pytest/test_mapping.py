"""Test mapping of basis function derivatives"""

import numpy as np
import pytest

import splinepy


class Geometry:
    def __init__(self):
        """Define geometries the function will be mapped to

        Create three simple geometries, two where results are known, and one
        with no zero entries in the jacobian.

        Functions will be 2D and 3D only
        """
        self.scaling_factors = 1 / np.array([2.0, 0.3, 1.5])
        self.scaling3D = splinepy.Bezier(
            degrees=[1, 1, 1],
            control_points=[
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [0.0, 0.3, 0.0],
                [2.0, 0.3, 0.0],
                [0.0, 0.0, 1.5],
                [2.0, 0.0, 1.5],
                [0.0, 0.3, 1.5],
                [2.0, 0.3, 1.5],
            ],
        )
        # Spline that rotates field by 0.17 radiants
        self.cc, self.ss = np.cos(0.17), np.sin(0.17)
        self.rotation_matrix = np.array(
            ((self.cc, -self.ss), (self.ss, self.cc))
        )
        self.rotating2D = splinepy.Bezier(
            degrees=[1, 1],
            control_points=[
                [0, 0],
                [1, 0],
                [0, 1],
                [1, 1],
            ],
        )
        self.rotating2D.control_points = np.einsum(
            "ij,qj->qi", self.rotation_matrix, self.rotating2D.control_points
        )

        # Complicated spline
        self.askew_spline2D = splinepy.BSpline(
            degrees=[2, 2],
            control_points=[
                [0.0, 0.0],
                [1.0, 0.5],
                [2.0, 0.2],
                [0.5, 1.5],
                [1.0, 1.5],
                [1.5, 1.5],
                [0.0, 3.0],
                [1.0, 2.5],
                [2.0, 3.0],
            ],
            knot_vectors=[[0, 0, 0, 1, 1, 1], [0, 0, 0, 1, 1, 1]],
        )

        # Two solution fields
        self.solution_field_mono3D = splinepy.Bezier(
            degrees=[2, 1, 1], control_points=np.ones((12, 1)) * 0.2
        )
        self.solution_field_rando = splinepy.Bezier(
            degrees=[2, 1],
            control_points=np.random.default_rng().random((6, 1)),
        )
        self.solution_field_rando2D = splinepy.Bezier(
            degrees=[2, 1],
            control_points=np.random.default_rng().random((6, 2)),
        )

        self.query_points2D = np.random.default_rng().random((13, 2))
        self.query_points3D = np.random.default_rng().random((17, 3))


# initialize geo variable
geo = Geometry()


def test_cross_evaluation_of_different_implementations():
    """Divergence and Laplacian are implemented differently when gradients
    are called at the same time
    """

    # Test Basis functions first
    mapper = geo.solution_field_rando2D.mapper(geo.askew_spline2D)
    bf_results = mapper.basis_function_derivatives(
        geo.query_points2D, gradient=True, hessian=True, laplacian=True
    )
    bf_gradient, support_gradient = mapper.basis_gradient_and_support(
        geo.query_points2D
    )
    bf_hessian, support_hessian = mapper.basis_hessian_and_support(
        geo.query_points2D
    )
    bf_laplace, support_laplacian = mapper.basis_laplacian_and_support(
        geo.query_points2D
    )

    # Laplacian is computed differently depending on function call
    assert np.allclose(bf_results["laplacian"], bf_laplace)
    assert np.allclose(bf_results["support"], support_laplacian)

    # Check function calls
    assert np.allclose(bf_results["gradient"], bf_gradient)
    assert np.allclose(bf_results["support"], support_gradient)
    assert np.allclose(bf_results["hessian"], bf_hessian)
    assert np.allclose(bf_results["support"], support_hessian)

    # Try field values
    results = mapper.field_derivatives(
        geo.query_points2D,
        gradient=True,
        divergence=True,
        hessian=True,
        laplacian=True,
        basis_function_values=True,
    )
    laplacian = mapper.laplacian(geo.query_points2D)
    divergence = mapper.divergence(geo.query_points2D)
    gradient = mapper.gradient(geo.query_points2D)
    hessian = mapper.hessian(geo.query_points2D)

    # Divergence and laplacian have different implementations
    assert np.allclose(results["laplacian"], laplacian)
    assert np.allclose(results["divergence"], divergence)

    # Check function calls
    assert np.allclose(results["hessian"], hessian)
    assert np.allclose(results["gradient"], gradient)

    # Check with previous test
    assert np.allclose(
        results["basis_function_values"]["gradient"],
        bf_results["gradient"],
    )


def test_check_assertions():
    mapper = geo.solution_field_rando.mapper(geo.askew_spline2D)
    with pytest.raises(
        ValueError,
        match=r"Divergence can only be performed "
        r"on vector fields with para_dim = dim",
    ):
        mapper.divergence(geo.query_points2D)


def test_first_order_derivatives_analytical():
    mapper2D = geo.solution_field_rando.mapper(geo.rotating2D)
    mapper3D = geo.solution_field_mono3D.mapper(geo.scaling3D)
    bf_gradient, support = mapper2D.basis_gradient_and_support(
        geo.query_points2D
    )
    (
        bf_reference0,
        supportb,
    ) = geo.solution_field_rando.basis_derivative_and_support(
        geo.query_points2D, [1, 0]
    )
    (
        bf_reference1,
        supportb,
    ) = geo.solution_field_rando.basis_derivative_and_support(
        geo.query_points2D, [0, 1]
    )
    bf_reference = np.dstack((bf_reference0, bf_reference1))

    assert np.allclose(support, supportb)

    # Rotate bf_reference
    bf_reference = np.einsum("ij,qsj->qsi", geo.rotation_matrix, bf_reference)

    assert np.allclose(bf_gradient, bf_reference)

    bf_gradient, support = mapper3D.basis_gradient_and_support(
        geo.query_points3D
    )
    (
        bf_reference0,
        supportb,
    ) = geo.solution_field_mono3D.basis_derivative_and_support(
        geo.query_points3D, [1, 0, 0]
    )
    (
        bf_reference1,
        supportb,
    ) = geo.solution_field_mono3D.basis_derivative_and_support(
        geo.query_points3D, [0, 1, 0]
    )
    (
        bf_reference2,
        supportb,
    ) = geo.solution_field_mono3D.basis_derivative_and_support(
        geo.query_points3D, [0, 0, 1]
    )
    bf_reference = np.dstack((bf_reference0, bf_reference1, bf_reference2))
    bf_reference = np.einsum("qsi,i->qsi", bf_reference, geo.scaling_factors)
    assert np.allclose(bf_gradient, bf_reference)


def test_second_order_analytical():
    mapper2D = geo.solution_field_rando.mapper(geo.rotating2D)
    bf_hessian, support = mapper2D.basis_hessian_and_support(
        geo.query_points2D
    )
    bf_reference = np.zeros(
        (
            geo.query_points2D.shape[0],
            np.prod(geo.solution_field_rando.degrees + 1),
            2,
            2,
        )
    )

    (
        bf_reference[:, :, 0, 0],
        supportb,
    ) = geo.solution_field_rando.basis_derivative_and_support(
        geo.query_points2D, [2, 0]
    )
    (
        bf_reference[:, :, 1, 0],
        supportb,
    ) = geo.solution_field_rando.basis_derivative_and_support(
        geo.query_points2D, [1, 1]
    )
    (
        bf_reference[:, :, 1, 1],
        supportc,
    ) = geo.solution_field_rando.basis_derivative_and_support(
        geo.query_points2D, [0, 2]
    )
    bf_reference[:, :, 0, 1] = bf_reference[:, :, 1, 0]

    assert np.allclose(support, supportb) and np.allclose(support, supportc)

    # Rotate bf_reference
    bf_reference = np.einsum(
        "ij,qsjk,kl->qsil",
        geo.rotation_matrix,
        bf_reference,
        geo.rotation_matrix.T,
    )

    assert np.allclose(bf_hessian, bf_reference)


def test_second_order_fd(np_rng):
    "Use proximity to get points on askew geometry and approcimate hessian"
    mapper = geo.solution_field_rando2D.mapper(geo.askew_spline2D)
    center_point_reference = np_rng.random((1, 2)) * 0.8 + 0.1
    dx = 1e-4

    # Analytical solution
    bf_references = mapper.basis_function_derivatives(
        queries=center_point_reference,
        gradient=True,
        hessian=True,
        laplacian=True,
    )
    references = mapper.field_derivatives(
        queries=center_point_reference,
        gradient=True,
        hessian=True,
        laplacian=True,
        divergence=True,
    )

    # Compute aux values for FD
    center_point = geo.askew_spline2D.evaluate(center_point_reference)
    center_point = np.repeat(center_point, 9, axis=0)
    center_point += np.array(
        [
            [-dx, -dx],
            [0, -2 * dx],
            [dx, -dx],
            [-2 * dx, 0],
            [0, 0],
            [2 * dx, 0],
            [-dx, dx],
            [0, 2 * dx],
            [dx, dx],
        ]
    )
    # Approximate points in the physical domain
    center_point_parametric = geo.askew_spline2D.proximities(
        center_point,
        initial_guess_sample_resolutions=[10, 10],
        tolerance=1e-12,
    )
    bfv, support = geo.solution_field_rando2D.basis_and_support(
        center_point_parametric
    )
    values = geo.solution_field_rando2D.evaluate(center_point_parametric)

    assert np.allclose(
        center_point_parametric[4, :],
        center_point_reference,
        atol=1e-8,
    )

    # Use FD to approximate
    # Gradient of basis function
    bf_gradient = np.zeros(bf_references["gradient"].shape)
    bf_gradient[:, :, 0] = (bfv[5, :] - bfv[3, :]) / (4 * dx)
    bf_gradient[:, :, 1] = (bfv[7, :] - bfv[1, :]) / (4 * dx)

    assert np.allclose(bf_gradient, bf_references["gradient"], atol=1e-4)

    # Hessians of basis function
    bf_hessian = np.zeros(bf_references["hessian"].shape)
    bf_hessian[:, :, 0, 0] = (bfv[5, :] + bfv[3, :] - 2 * bfv[4, :]) / (
        4 * dx * dx
    )
    bf_hessian[:, :, 1, 1] = (bfv[1, :] + bfv[7, :] - 2 * bfv[4, :]) / (
        4 * dx * dx
    )
    bf_hessian[:, :, 1, 0] = (
        bfv[8, :] - bfv[6, :] - bfv[2, :] + bfv[0, :]
    ) / (4 * dx * dx)
    bf_hessian[:, :, 0, 1] = bf_hessian[:, :, 1, 0]

    assert np.allclose(bf_hessian, bf_references["hessian"], atol=1e-4)

    # Gradient
    gradient = np.zeros(references["gradient"].shape)
    gradient[:, :, 0] = (values[5, :] - values[3, :]) / (4 * dx)
    gradient[:, :, 1] = (values[7, :] - values[1, :]) / (4 * dx)

    assert np.allclose(gradient, references["gradient"], atol=dx)

    # Hessians
    hessian = np.zeros(references["hessian"].shape)
    hessian[:, :, 0, 0] = (values[5, :] + values[3, :] - 2 * values[4, :]) / (
        4 * dx * dx
    )
    hessian[:, :, 1, 1] = (values[1, :] + values[7, :] - 2 * values[4, :]) / (
        4 * dx * dx
    )
    hessian[:, :, 1, 0] = (
        values[8, :] - values[6, :] - values[2, :] + values[0, :]
    ) / (4 * dx * dx)
    hessian[:, :, 0, 1] = hessian[:, :, 1, 0]

    assert np.allclose(hessian, references["hessian"], atol=dx)

    # Check reduced values
    laplacian = (
        references["hessian"][:, :, 0, 0] + references["hessian"][:, :, 1, 1]
    )
    divergence = (
        references["gradient"][:, 0, 0] + references["gradient"][:, 1, 1]
    )

    assert np.allclose(laplacian, references["laplacian"])
    assert np.allclose(divergence, references["divergence"])
