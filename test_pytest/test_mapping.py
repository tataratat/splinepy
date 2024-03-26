"""Test mapping of basis function derivatives"""

import numpy as np

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


def check_assertions():
    mapper = geo.solution_field_rando.mapper(geo.askew_spline2D)
    assert mapper.divergence(geo.query_points2D)
