r"""
This file implements a stokes test case based on the example stated in
"Finite Element Methods for Flow Problems", Donea and Huerta p.306
chapter 6.8.1 example with analytical solution.
"""

import numpy as np
from scipy.sparse import bmat, csr_matrix, linalg

import splinepy as sp

###
# Inputs
###

# Define the dynamic viscosity
viscosity = 1

# Define the Geometry (simple first order unit square)
geometry = sp.BSpline(
    degrees=[1, 1],
    control_points=[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
    knot_vectors=[[0, 0, 1, 1], [0, 0, 1, 1]],
)


solution_field_pressure = sp.BSpline(
    degrees=geometry.degrees,
    control_points=np.ones((geometry.control_points.shape[0], 1)),
    knot_vectors=geometry.knot_vectors,
)
solution_field_pressure.elevate_degrees([0, 1, 0, 1])
solution_field_pressure.uniform_refine([1], 7)
solution_field_pressure.uniform_refine([0], 7)
# Refinement leads to quadratic spline with 10x10cps

# Create Raviart-Thomas mixed style splines
solution_field_velocity_u = solution_field_pressure.copy()
solution_field_velocity_u.elevate_degrees([0])
solution_field_velocity_v = solution_field_pressure.copy()
solution_field_velocity_v.elevate_degrees([1])

###
# Source Functions
###


def source_function_u(x_vec):
    x = x_vec[:, 0]
    y = x_vec[:, 1]
    return (
        (12 - 24 * y) * x**4
        + (-24 + 48 * y) * x**3
        + (-48 * y + 72 * y**2 - 48 * y**3 + 12) * x**2
        + (-2 + 24 * y - 72 * y**2 + 48 * y**3) * x
        + (1 - 4 * y + 12 * y**2 - 8 * y**3)
    )


def source_function_v(x_vec):
    x = x_vec[:, 0]
    y = x_vec[:, 1]
    return (
        (8 - 48 * y + 48 * y**2) * x**3
        + (-12 + 72 * y - 72 * y**2) * x**2
        + (4 - 24 * y + 48 * y**2 - 48 * y**3 + 24 * y**4) * x
        + (-12 * y**2 + 24 * y**3 - 12 * y**4)
    )


###
# Computation of greville abscissae
###
greville_points_u = solution_field_velocity_u.greville_abscissae()
greville_points_v = solution_field_velocity_v.greville_abscissae()
greville_points_p = solution_field_pressure.greville_abscissae()

###
# Computation of the matrices
###
mapper_u = solution_field_velocity_u.mapper(reference=geometry)
mapper_v = solution_field_velocity_v.mapper(reference=geometry)
mapper_p = solution_field_pressure.mapper(reference=geometry)

# FIRST ROW ----------
# A_u_tu = mu * delta(u)
A_u_tu = -viscosity * sp.utils.data.make_matrix(
    *mapper_u.basis_laplacian_and_support(greville_points_u),
    solution_field_velocity_u.cps.shape[0],
    as_array=False,
)
# A_v_tu = 0
A_v_tu = None
#  A_p_tu = dpdx
gradient_p, support_p = mapper_p.basis_gradient_and_support(greville_points_u)
A_p_tu = sp.utils.data.make_matrix(
    gradient_p[:, :, 0],
    support_p,
    solution_field_pressure.cps.shape[0],
    as_array=False,
)

# SECOND ROW ----------
# A_u_tv = 0
A_u_tv = None
# A_v_tv = mu * delta(v)
A_v_tv = -viscosity * sp.utils.data.make_matrix(
    *mapper_v.basis_laplacian_and_support(greville_points_v),
    solution_field_velocity_v.cps.shape[0],
    as_array=False,
)
#  A_p_tv
gradient_p, support_p = mapper_p.basis_gradient_and_support(greville_points_v)
A_p_tv = sp.utils.data.make_matrix(
    gradient_p[:, :, 1],
    support_p,
    solution_field_pressure.cps.shape[0],
    as_array=False,
)
# THIRD ROW ----------
# A_u_tp dudx
gradient_u, support_u = mapper_u.basis_gradient_and_support(greville_points_p)
A_u_tp = sp.utils.data.make_matrix(
    gradient_u[:, :, 0],
    support_u,
    solution_field_velocity_u.cps.shape[0],
    as_array=False,
)
# A_v_tp dvdx
gradient_v, support_v = mapper_v.basis_gradient_and_support(greville_points_p)
A_v_tp = sp.utils.data.make_matrix(
    gradient_v[:, :, 1],
    support_v,
    solution_field_velocity_v.cps.shape[0],
    as_array=False,
)
# A_p_tp = 0
A_p_tp = None


###
# Evaluation of rhs
###
physical_points_u = geometry.evaluate(greville_points_u)
physical_points_v = geometry.evaluate(greville_points_v)
rhs_u = source_function_u(physical_points_u)
rhs_v = source_function_v(physical_points_v)
rhs_p = np.zeros(greville_points_p.shape[0])


###
# Imposition of strong bcs (There must be a better solution here...)
###
# identify boundary nodes u
boundary_dof_ids = np.concatenate(
    (
        solution_field_velocity_u.multi_index[0, :],
        solution_field_velocity_u.multi_index[-1, :],
        solution_field_velocity_u.multi_index[1:-1, 0],
        solution_field_velocity_u.multi_index[1:-1, -1],
    )
)
A_u_tu[boundary_dof_ids] *= 0  # this does not change sparsity pattern
A_p_tu[boundary_dof_ids] *= 0
A_u_tu[boundary_dof_ids, boundary_dof_ids] = 1
rhs_u[boundary_dof_ids] = 0.0

# identify boundary nodes v
boundary_dof_ids = None
boundary_dof_ids = np.concatenate(
    (
        solution_field_velocity_v.multi_index[0, :],
        solution_field_velocity_v.multi_index[-1, :],
        solution_field_velocity_v.multi_index[1:-1, 0],
        solution_field_velocity_v.multi_index[1:-1, -1],
    )
)
A_v_tv[boundary_dof_ids] *= 0  # this does not change sparsity pattern
A_v_tv[boundary_dof_ids, boundary_dof_ids] = 1
A_p_tv[boundary_dof_ids] *= 0
rhs_v[boundary_dof_ids] = 0.0

# Add condition for pressure (optional)
pressure_treatment = ""
if pressure_treatment.lower() == "normalize":
    n = greville_points_p.shape[0]
    row_idx = np.full(n, n - 1)  # all entries in last row
    col_idx = np.arange(n)  # columns 0..m-1
    data = np.ones(n)
    A_p_tp = csr_matrix((data, (row_idx, col_idx)), shape=(n, n))
    A_u_tp[-1, :] *= 0
    A_v_tp[-1, :] *= 0
elif pressure_treatment.lower() == "singlepoint":
    n = greville_points_p.shape[0]
    A_u_tp[-1, :] *= 0
    A_v_tp[-1, :] *= 0
    A_p_tp = csr_matrix(([1], ([n - 1], [n - 1])), shape=(n, n))
    rhs_p[n - 1] = 0.0


###
# Solve linear system
###
rhs_all = np.concatenate([rhs_u, rhs_v, rhs_p])
matrix_all = bmat(
    [
        [A_u_tu, A_v_tu, A_p_tu],
        [A_u_tv, A_v_tv, A_p_tv],
        [A_u_tp, A_v_tp, A_p_tp],
    ],
    format=A_u_tu.format,
)
solution_vector = linalg.spsolve(matrix_all, rhs_all)

###
# Assign solution to original splines
###
solution_field_velocity_u.control_points = np.reshape(
    solution_vector[: solution_field_velocity_u.cps.shape[0]], (-1, 1)
)
solution_field_velocity_v.control_points = np.reshape(
    solution_vector[
        solution_field_velocity_u.cps.shape[0] : -(
            solution_field_pressure.cps.shape[0]
        )
    ],
    (-1, 1),
)
solution_field_pressure.control_points = np.reshape(
    solution_vector[-(solution_field_pressure.cps.shape[0]) :], (-1, 1)
)


###
# Plot and analysis
###

# Analytical solutions


def analytical_solution_u(geometry, on):
    x_vec = geometry.evaluate(on)
    x = x_vec[:, 0]
    y = x_vec[:, 1]
    return (x**2) * (1.0 - x) * (1.0 - x) * (2 * y - 6 * y**2 + 4 * y**3)


def error_u(solution_u):
    def error_u(data, on):
        sol = solution_u.evaluate(on)
        ansol = analytical_solution_u(data, on)
        return np.abs(sol.flat - ansol)

    return error_u


def analytical_solution_v(geometry, on):
    x_vec = geometry.evaluate(on)
    x = x_vec[:, 0]
    y = x_vec[:, 1]
    return -(y**2) * (1.0 - y) * (1.0 - y) * (2 * x - 6 * x**2 + 4 * x**3)


def error_v(solution_v):
    def error_v(data, on):
        sol = solution_v.evaluate(on)
        ansol = analytical_solution_v(data, on)
        return np.abs(sol.flat - ansol)

    return error_v


def analytical_solution_p(geometry, on):
    x_vec = geometry.evaluate(on)
    x = x_vec[:, 0]
    return x * (1 - x)


def error_p(solution_p):
    def error_p(data, on):
        sol = solution_p.evaluate(on)
        ansol = analytical_solution_p(data, on)
        return np.abs(sol.flat - ansol)

    return error_p


# Plot geometry and fields
geometry.spline_data["field"] = solution_field_pressure
geometry.show_options["data"] = "field"
geometry.show_options["cmap"] = "jet"
geometry.show_options["lighting"] = "off"
geometry.show_options["scalarbar"] = True
solution_p = geometry.copy()
geometry.spline_data["field"] = solution_field_velocity_v
solution_v = geometry.copy()
geometry.spline_data["field"] = solution_field_velocity_u
solution_u = geometry.copy()
geometry.spline_data["field"] = sp.SplineDataAdaptor(
    geometry, function=error_u(solution_field_velocity_u)
)
error_field_u = geometry.copy()
geometry.spline_data["field"] = sp.SplineDataAdaptor(
    geometry, function=error_v(solution_field_velocity_v)
)
error_field_v = geometry.copy()
geometry.spline_data["field"] = sp.SplineDataAdaptor(
    geometry, function=error_p(solution_field_pressure)
)
error_field_p = geometry.copy()

sp.show(
    ["U-Velocity", solution_u],
    ["V-velocity", solution_v],
    ["Pressure", solution_p],
    ["|U - U_exp|", error_field_u],
    ["|V - V_exp|", error_field_v],
    ["|P - P_exp|", error_field_p],
    knots=True,
    control_points=False,
)


###
# Loss function u
###

# Check Loss function
sample_points = greville_points_p
mapper_u = solution_field_velocity_u.mapper(reference=geometry)
mapper_v = solution_field_velocity_v.mapper(reference=geometry)
mapper_p = solution_field_pressure.mapper(reference=geometry)


def loss_function_u(data, on):
    return (
        viscosity * mapper_v.laplacian(on)
        - mapper_p.gradient(on)[:, 0, 1]
        - source_function_v(data.evaluate(on))
    )


def loss_function_v(data, on):
    return (
        viscosity * mapper_u.laplacian(on)
        - mapper_p.gradient(on)[:, 0, 0]
        - source_function_u(data.evaluate(on))
    )


geometry.spline_data["field"] = sp.SplineDataAdaptor(
    geometry, function=loss_function_u
)
loss_u = geometry.copy()
loss_u.show()
