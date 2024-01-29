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

import splinepy as sp

try:
    import scipy

    has_scipy = True
except ImportError:
    has_scipy = False

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

# Compute Laplacian matrix
laplacian = sp.utils.data.make_matrix(
    *mapper.basis_laplacian_and_support(evaluation_points),
    solution_field.control_points.shape[0],
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
if has_scipy:
    solution_field.control_points = scipy.sparse.linalg.spsolve(
        -laplacian, rhs
    ).reshape(-1, 1)
else:
    solution_field.control_points = np.linalg.solve(-laplacian, rhs).reshape(
        -1, 1
    )


# Plot geometry and field
geometry.spline_data["field"] = solution_field
geometry.show_options["data"] = "field"
geometry.show_options["cmap"] = "jet"
geometry.show_options["lighting"] = "off"
geometry.show_options["scalarbar"] = True
geometry.show(knots=True, control_points=False)
