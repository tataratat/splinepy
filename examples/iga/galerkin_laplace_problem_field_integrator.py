r"""
This is a showcase of splinepy's FieldIntegrator class. It solves the Poisson
equation in the form
.. math::
        -\Delta u = f
with source term f=1.

First, homogeneous Dirichlet boundary conditions are applied. Then, inhomogeneous ones
are applied.
"""

import numpy as np

import splinepy as sp

# Test Case
n_refine = 15


def prepare_geometry_and_solution_field():
    # Define the Geometry
    geometry = sp.BSpline(
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

    # Define the solution field
    solution_field = geometry.copy()
    solution_field.control_points = np.ones(
        (geometry.control_points.shape[0], 1)
    )

    # Use refinement
    solution_field.elevate_degrees([0, 1, 0, 1])
    new_knots = np.linspace(1 / n_refine, 1, n_refine, endpoint=False)
    solution_field.insert_knots(0, new_knots)
    solution_field.insert_knots(1, new_knots)

    return geometry, solution_field


def poisson_lhs(mapper, quad_points, quad_weights, jacobian_det):
    """
    Assemble the system matrix of the Poisson equation.

    Parameters
    -----------
    mapper: splinepy.helpme.mapper.Mapper
        Mapper of solution field
    quad_points: np.ndarray
        Quadrature points
    quad_weights: np.ndarray
        Quadrature weights
    jacobian_det: np.ndarray
        Determinant of Jacobian, evaluated at quadrature points

    Returns
    -----------
    element_matrix: np.ndarray
        Element matrix
    """
    bf_gradient, _ = mapper.basis_gradient_and_support(quad_points)
    element_matrix = np.einsum(
        "qid,qjd,q,q->ij",
        bf_gradient,
        bf_gradient,
        quad_weights,
        jacobian_det,
        optimize=True,
    )
    return element_matrix.ravel()


def poisson_rhs(mapper, quad_points, quad_weights, jacobian_det):
    """
    Assemble the rhs of the Poisson equation with f=1

    Parameters
    -----------
    mapper: splinepy.helpme.mapper.Mapper
        Mapper of solution field
    quad_points: np.ndarray
        Quadrature points
    quad_weights: np.ndarray
        Quadrature weights
    jacobian_det: np.ndarray
        Determinant of Jacobian, evaluated at quadrature points

    Returns
    -----------
    element_vector: np.ndarray
        Element vector
    """
    bf = mapper._field_reference.basis(quad_points)
    element_vector = np.einsum(
        "qj,q,q->j", bf, quad_weights, jacobian_det, optimize=True
    )
    return element_vector


def show_solution(geometry, solution_field):
    """Visualize the solution"""
    geometry.spline_data["field"] = solution_field
    geometry.show_options["data"] = "field"
    geometry.show_options["cmap"] = "jet"
    geometry.show_options["lighting"] = "off"
    geometry.show_options["scalarbar"] = True
    geometry.show(knots=True, control_points=False)


if __name__ == "__main__":
    geometry, solution_field = prepare_geometry_and_solution_field()

    fi = sp.helpme.integrate.FieldIntegrator(
        geometry=geometry, solution_field=solution_field
    )

    fi.assemble_matrix(poisson_lhs)
    fi.assemble_vector(poisson_rhs)

    # Homogeneous Dirichlet boundary conditions
    fi.assign_homogeneous_dirichlet_boundary_conditions()
    fi.solve_linear_system()
    # Plot geometry and field
    show_solution(geometry, solution_field)

    # Inhomogeneous Dirichlet boundary conditions
    def dirichlet_function(points):
        """
        On the boundary apply: g(x,y) = x/10
        """
        return points[:, 0] / 10

    fi.apply_dirichlet_boundary_conditions(dirichlet_function)
    fi.solve_linear_system()
    show_solution(geometry, solution_field)
