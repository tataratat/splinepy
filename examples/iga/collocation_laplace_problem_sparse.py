r"""
In this file we will prototype a collocation example (in the simplest way),
that shows how to implement a Poisson problem in the form
.. math::
        -\Delta u = f
with source term f=1

The refinement will only be applied to the solution field, to make the
calculations more efficient
"""

import numpy as np
from scipy.sparse import linalg

import splinepy as sp

# Test Case
n_refine = 15


# Source Function
def source_function(x):
    return np.ones(x.shape[0])


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


# Get greville points and geometric values
evaluation_points = solution_field.greville_abscissae()
mapper = solution_field.mapper(reference=geometry)
sp_laplacian = sp.utils.data.make_matrix(
    *mapper.basis_laplacian_and_support(evaluation_points),
    solution_field.cps.shape[0],
    as_array=False,
)


# Evaluate RHS function
physical_points = geometry.evaluate(evaluation_points)
rhs = source_function(physical_points)

# Set dirichlet values (identify boundary nodes and set matrix to identity rhs
# to 0)
boundary_dof_ids = np.concatenate(
    (
        solution_field.multi_index[0, :],
        solution_field.multi_index[-1, :],
        solution_field.multi_index[1:-1, 0],
        solution_field.multi_index[1:-1, -1],
    )
)
sp_laplacian[boundary_dof_ids] *= 0  # this does not change sparsity pattern
sp_laplacian[boundary_dof_ids, boundary_dof_ids] = 1

rhs[boundary_dof_ids] = 0.0

# Solve linear system
solution_field.control_points[:] = linalg.spsolve(-sp_laplacian, rhs).reshape(
    -1, 1
)


# Plot geometry and field
geometry.spline_data["field"] = solution_field
geometry.show_options["data"] = "field"
geometry.show_options["cmap"] = "jet"
geometry.show_options["lighting"] = "off"
geometry.show_options["scalarbar"] = True
geometry.show(knots=True, control_points=False)
