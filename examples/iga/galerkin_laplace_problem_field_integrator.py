r"""
This is a showcase of splinepy's FieldIntegrator class and is doing the same as
'examples/galerkin_laplace_problem.py'.
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
    solution_field = sp.BSpline(
        degrees=geometry.degrees,
        control_points=np.ones((geometry.control_points.shape[0], 1)),
        knot_vectors=geometry.knot_vectors,
    )
    
    # Use refinement
    solution_field.elevate_degrees([0, 1, 0, 1])
    new_knots = np.linspace(1 / n_refine, 1, n_refine, endpoint=False)
    solution_field.insert_knots(0, new_knots)
    solution_field.insert_knots(1, new_knots)
    
    return geometry, solution_field

def poisson_lhs(mapper, quad_points, quad_weights, jacobian_det):
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

def poisson_rhs(solution_field, quad_points, quad_weights, jacobian_det, current_sol=None):
    bf, _ = solution_field.basis_and_support(quad_points)
    element_vector = np.einsum(
        "qj,q,q->j",
        bf,
        quad_weights,
        jacobian_det,
        optimize=True
    )
    return element_vector

if __name__ == "__main__":
    geometry, solution_field = prepare_geometry_and_solution_field()
    
    fi = sp.helpme.integrate.FieldIntegrator(geometry=geometry, solution_field=solution_field)
    
    fi.assemble_matrix(poisson_lhs)
    fi.assemble_vector(poisson_rhs)
    fi.assign_homogeneous_dirichlet_boundary_conditions()
    
    fi.solve_linear_system()
    
    # Plot geometry and field
    geometry.spline_data["field"] = solution_field
    geometry.show_options["data"] = "field"
    geometry.show_options["cmap"] = "jet"
    geometry.show_options["lighting"] = "off"
    geometry.show_options["scalarbar"] = True
    geometry.show(knots=True, control_points=False)
    