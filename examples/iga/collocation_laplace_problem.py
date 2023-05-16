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
import scipy
from scipy.sparse import csc_matrix

import splinepy as sp

try:
    import gustaf as gus

    has_gus = True
except ImportError:
    has_gus = False

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
evaluation_points = solution_field.greville_abscissae
mapper = solution_field.mapper(reference=geometry)
bf_laplacian, support = mapper.basis_laplacian_and_support(evaluation_points)
# Use the values and support to store it in global matrix (sparse)
laplacian = np.zeros(
    (evaluation_points.shape[0], solution_field.control_points.shape[0])
)

n_dofs = evaluation_points.shape[0]
row_ids = (
    np.arange(support.shape[0])
    .reshape(-1, 1)
    .repeat(support.shape[1], axis=1)
    .flatten()
)
col_ids = support.flatten()
laplacian = scipy.sparse.csc_array(
    (bf_laplacian.flatten(), (row_ids, col_ids)), shape=(n_dofs, n_dofs)
)

# Evaluate RHS function
physical_points = geometry.evaluate(evaluation_points)
rhs = source_function(physical_points)

# Set dirichlet values (identify boundary nodes and set matrix to identity rhs
# to 0)
indices = np.hstack(
    (
        solution_field.multi_index[0, :],
        solution_field.multi_index[-1, :],
        solution_field.multi_index[:, 0],
        solution_field.multi_index[:, -1],
    )
)
laplacian[indices, :] = 0
laplacian[indices, indices] = 1
rhs[indices] = 0

# Solve linear system
solution_field.control_points = scipy.sparse.linalg.spsolve(
    -laplacian, rhs
).reshape(-1, 1)

if has_gus:
    geometry = gus.spline.BSpline(spline=geometry)
    solution_field = gus.spline.BSpline(spline=solution_field)

    # Plot geometry and field
    geometry.spline_data["field"] = solution_field
    geometry.show_options["data_name"] = "field"
    geometry.show_options["cmap"] = "jet"
    geometry.show_options["lighting"] = "off"
    geometry.show_options["scalarbar"] = True
    gus.show(geometry, knots=True, control_points=False)
