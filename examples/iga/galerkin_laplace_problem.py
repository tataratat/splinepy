r"""
In this file we will prototype a galerkin IGA example (in the simplest way),
that shows how to implement a Poisson problem in the form
.. math::
        -\Delta u = f
with source term f=1

The refinement will only be applied to the solution field, to make the
calculations more efficient
"""

import numpy as np

import splinepy as sp

# Test Case
n_refine = 15


# Source Function
def source_function(x):
    return np.ones(x.shape[0])


# Auxiliary function to map positions into the element
def map_positions(positions, x_min, x_max, y_min, y_max):
    mapped_positions = np.empty(positions.shape, dtype=np.float64)
    mapped_positions[:, 0] = (x_max - x_min) * 0.5 * positions[:, 0] + (
        x_max + x_min
    ) * 0.5
    mapped_positions[:, 1] = (y_max - y_min) * 0.5 * positions[:, 0] + (
        y_max + y_min
    ) * 0.5
    return mapped_positions


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

# Retrieve integration points and weights
max_order = int(np.max(solution_field.degrees))
positions, weights = np.polynomial.legendre.leggauss(deg=max_order)
positions = sp.utils.data.cartesian_product(
    [positions for _ in range(geometry.para_dim)]
)
weights = sp.utils.data.cartesian_product(
    [weights for _ in range(geometry.para_dim)]
)
weights = np.prod(weights, axis=1)

# Assemble individual elements
ukv = solution_field.unique_knots
n_dofs = solution_field.control_points.shape[0]
system_matrix = np.zeros((n_dofs, n_dofs))
system_rhs = np.zeros(n_dofs)
mapper = solution_field.mapper(reference=geometry)

# Element Loop
for i in range(len(ukv[0]) - 1):
    for j in range(len(ukv[1]) - 1):
        # Auxiliary values
        mapped_positions = map_positions(
            positions,
            ukv[0][i],
            ukv[0][i + 1],
            ukv[1][j],
            ukv[1][j + 1],
        )
        det_jacs = np.linalg.det(geometry.jacobian(mapped_positions))

        # RHS
        # q : quadrature point | d: dim
        bf, support = solution_field.basis_and_support(mapped_positions)
        assert np.all(support[0, :] == support)
        system_rhs[support[0, :]] += np.einsum(
            "qj,q,q->j", bf, weights, det_jacs, optimize=True
        )

        # LHS
        bf_jacobian, support = mapper.basis_gradient_and_support(
            mapped_positions
        )
        assert np.all(support[0, :] == support)
        local_matrix = np.einsum(
            "qid,qjd,q,q->ij",
            bf_jacobian,
            bf_jacobian,
            weights,
            det_jacs,
            optimize=True,
        )
        support = sp.utils.data.cartesian_product(
            (support[0, :], support[0, :])
        )
        system_matrix[support[:, 0], support[:, 1]] += local_matrix.flatten()

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
system_matrix[indices, :] = 0
system_matrix[indices, indices] = 1
system_rhs[indices] = 0

# Solve linear system
solution_field.control_points = np.linalg.solve(
    system_matrix, system_rhs
).reshape(-1, 1)

# Plot geometry and field
geometry.spline_data["field"] = solution_field
geometry.show_options["data"] = "field"
geometry.show_options["cmap"] = "jet"
geometry.show_options["lighting"] = "off"
geometry.show_options["scalarbar"] = True
geometry.show(knots=True, control_points=False)
