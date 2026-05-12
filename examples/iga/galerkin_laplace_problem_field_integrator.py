r"""
This is a showcase of splinepy's FieldIntegrator class. It solves the Poisson
equation in the form
.. math::
        -\Delta u = f
with source term f=1.

First, homogeneous Dirichlet boundary conditions are applied. Then, a mix of homogeneous
and inhomogeneous Dirichlet boundary conditions are applied on 3 of 4 boundaries. A
zero Neumann boundary condition is implicitly applied on the fourth boundary.
"""

import numpy as np

import splinepy as sp

# Number of refinements for the solution field
n_refine = 15


def prepare_geometry_and_solution_field():
    """
    Creates single patch geometry and solution field with h- and p-refinement

    Returns
    -----------------
    geometry: splinepy.BSpline
        The geometry object
    solution_field: splinepy.BSpline
        The solution field with refinements. Its control points serve as the degree of
        freedoms (DoFs) of the solution. It is initialized to a vector of ones.
    """
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
    # Initialize solution vector
    solution_field.control_points = np.ones(
        (geometry.control_points.shape[0], 1)
    )

    # Apply p-refinement
    solution_field.elevate_degrees([0, 1, 0, 1])
    # Apply uniform h-refinement
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
    # Create geometry and solution field with h- and p-refinement
    geometry, solution_field = prepare_geometry_and_solution_field()

    fi = sp.helpme.integrate.FieldIntegrator(
        geometry=geometry, solution_field=solution_field
    )

    # Assemble the linear system for the Poisson equation
    fi.assemble_matrix(poisson_lhs)
    fi.assemble_vector(poisson_rhs)

    # Apply homogeneous Dirichlet boundary conditions
    fi.apply_homogeneous_dirichlet_boundary_conditions()
    # Solve linear system to obtain solution vector, which is stored as the
    # solution_field's control points
    fi.solve_linear_system()
    # Plot geometry and field
    show_solution(geometry, solution_field)

    # Function for inhomogeneous Dirichlet boundary conditions
    def dirichlet_function(points):
        """
        On the boundary apply: g(x,y) = x/4
        """
        return points[:, 0] / 4

    # Assemble again to override previous boundary conditions
    fi.assemble_matrix(poisson_lhs)
    fi.assemble_vector(poisson_rhs)

    # Apply boundary conditions on 3 boundaries (west, south and north boundary)
    fi.apply_homogeneous_dirichlet_boundary_conditions(west=True)
    fi.apply_dirichlet_boundary_conditions(
        dirichlet_function, south=True, north=True
    )
    # For zero Neumann boundary conditions we can keep the matrix and rhs as they are

    # Solve linear system
    fi.solve_linear_system()
    # Plot resulting solution
    show_solution(geometry, solution_field)
