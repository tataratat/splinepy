r"""
This file implements a stokes test case based on the example stated in
"Finite Element Methods for Flow Problems", Donea and Huerta p.306
chapter 6.8.1 example with analytical solution.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import bmat, linalg

import splinepy as sp

###
# Inputs
###

# Define the dynamic viscosity
viscosity = 1
mass_factor = 1e-1  # Weight the incompressibility equations in the LSQ solve
impose_pressure_bcs = False  # Impose BCS for pressure

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
solution_field_pressure.elevate_degrees([0, 1])
solution_field_pressure.uniform_refine([1], 7)
solution_field_pressure.uniform_refine([0], 7)
# Refinement leads to quadratic spline with 10x10cps

# Create Raviart-Thomas mixed style splines
solution_field_velocity_u = solution_field_pressure.copy()
# For div-conforming tensor-product spaces in 2D, the x-velocity has
# one higher degree in x, and the y-velocity has one higher degree in y.
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
        + (12 - 48 * y + 72 * y**2 - 48 * y**3) * x**2
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


def boundary_conditions_u(x_vec):
    return np.zeros(x_vec.shape[0])


def boundary_conditions_v(x_vec):
    return np.zeros(x_vec.shape[0])


def boundary_conditions_p(x_vec):
    return np.zeros(x_vec.shape[0])


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
# RHS
rhs_u = source_function_u(geometry.evaluate(greville_points_u))

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
# RHS
rhs_v = source_function_v(geometry.evaluate(greville_points_v))

# THIRD ROW ----------
# A_u_tp dudx
gradient_u, support_u = mapper_u.basis_gradient_and_support(greville_points_p)
A_u_tp = sp.utils.data.make_matrix(
    mass_factor * gradient_u[:, :, 0],
    support_u,
    solution_field_velocity_u.cps.shape[0],
    as_array=False,
)
# A_v_tp dvdx
gradient_v, support_v = mapper_v.basis_gradient_and_support(greville_points_p)
A_v_tp = sp.utils.data.make_matrix(
    mass_factor * gradient_v[:, :, 1],
    support_v,
    solution_field_velocity_v.cps.shape[0],
    as_array=False,
)
# A_p_tp = 0
A_p_tp = None
# RHS
rhs_p = mass_factor * np.zeros(greville_points_p.shape[0])


###
# Imposition of BCs
###
# For U-Velocity
boundary_evaluation_point_ids_u = np.concatenate(
    (
        solution_field_velocity_u.multi_index[0, :],
        solution_field_velocity_u.multi_index[-1, :],
        solution_field_velocity_u.multi_index[1:-1, 0],
        solution_field_velocity_u.multi_index[1:-1, -1],
    )
)
boundary_evaluation_points_u = greville_points_u[
    boundary_evaluation_point_ids_u
]

# Create matrices and RHS
A_u_tpbcu = sp.utils.data.make_matrix(
    *solution_field_velocity_u.basis_and_support(boundary_evaluation_points_u),
    solution_field_velocity_u.cps.shape[0],
    as_array=False,
)
A_v_tpbcu = None
A_p_tpbcu = None
rhs_tpbcu = boundary_conditions_u(boundary_evaluation_points_u)

# For V-velocity
boundary_evaluation_point_ids_v = np.concatenate(
    (
        solution_field_velocity_v.multi_index[0, :],
        solution_field_velocity_v.multi_index[-1, :],
        solution_field_velocity_v.multi_index[1:-1, 0],
        solution_field_velocity_v.multi_index[1:-1, -1],
    )
)
boundary_evaluation_points_v = greville_points_v[
    boundary_evaluation_point_ids_v
]

# Create matrices and RHS
A_u_tpbcv = None
A_v_tpbcv = sp.utils.data.make_matrix(
    *solution_field_velocity_v.basis_and_support(boundary_evaluation_points_v),
    solution_field_velocity_v.cps.shape[0],
    as_array=False,
)
A_p_tpbcv = None
rhs_tpbcv = boundary_conditions_v(boundary_evaluation_points_v)

# For Pressure Field (redundant?)
if impose_pressure_bcs:
    boundary_evaluation_point_ids_p = np.concatenate(
        (
            solution_field_pressure.multi_index[0, :],
            solution_field_pressure.multi_index[-1, :],
        )
    )
    boundary_evaluation_points_p = greville_points_p[
        boundary_evaluation_point_ids_p
    ]

    # Create matrices and RHS
    A_u_tpbcp = None
    A_v_tpbcp = None
    A_p_tpbcp = sp.utils.data.make_matrix(
        *solution_field_pressure.basis_and_support(
            boundary_evaluation_points_p
        ),
        solution_field_pressure.cps.shape[0],
        as_array=False,
    )
    rhs_tpbcp = boundary_conditions_p(boundary_evaluation_points_p)


###
# Solve linear system
###
rhs_all = np.concatenate([rhs_u, rhs_v, rhs_p, rhs_tpbcu, rhs_tpbcv])
block_matrices = [
    [A_u_tu, A_v_tu, A_p_tu],
    [A_u_tv, A_v_tv, A_p_tv],
    [A_u_tp, A_v_tp, A_p_tp],
    [A_u_tpbcu, A_v_tpbcu, A_p_tpbcu],
    [A_u_tpbcv, A_v_tpbcv, A_p_tpbcv],
]
if impose_pressure_bcs:
    rhs_all = np.concatenate([rhs_all, rhs_tpbcp])
    block_matrices.append([A_u_tpbcp, A_v_tpbcp, A_p_tpbcp])
matrix_all = bmat(block_matrices, format=A_u_tu.format)
solution_vector, istop, itn, r1norm = linalg.lsqr(
    matrix_all,
    rhs_all,
    atol=1e-12,
    btol=1e-12,
    iter_lim=100 * matrix_all.shape[1],
)[:4]
print(f"LSQR istop={istop}, iterations={itn}, residual={r1norm:.3e}")


# Plot Matrix
plt.spy(matrix_all)
plt.show()


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


# if we do not impose pressure BCs, we need to align the computed pressure
# field with the analytical solution because the pressure field is only
# determined up to a constant
pressure_error_offset = 0.0
if not impose_pressure_bcs:
    # we choose the offset as the minimum value of the computed pressure field
    # such that the computed pressure field is (theoretically) non-negative
    pressure_error_offset = solution_field_pressure.evaluate(
        greville_points_p
    ).min()


def error_p(solution_p):
    def error_p(data, on):
        sol = solution_p.evaluate(on)
        corr_sol = sol.flat - pressure_error_offset
        ansol = analytical_solution_p(data, on)
        return np.abs(corr_sol - ansol)

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


###
# Loss function u
###

# Check Loss function
sample_points = greville_points_p
mapper_u = solution_field_velocity_u.mapper(reference=geometry)
mapper_v = solution_field_velocity_v.mapper(reference=geometry)
mapper_p = solution_field_pressure.mapper(reference=geometry)


def loss_function_u(data, on):
    return np.abs(
        -viscosity * mapper_u.laplacian(on).ravel()
        + mapper_p.gradient(on)[:, 0, 0]
        - source_function_u(data.evaluate(on))
    )


def loss_function_v(data, on):
    return np.abs(
        -viscosity * mapper_v.laplacian(on).ravel()
        + mapper_p.gradient(on)[:, 0, 1]
        - source_function_v(data.evaluate(on))
    )


def loss_function_p(data, on):
    return mass_factor * (
        mapper_u.gradient(on)[:, 0, 0] + mapper_v.gradient(on)[:, 0, 1]
    )


geometry.spline_data["field"] = sp.SplineDataAdaptor(
    geometry, function=loss_function_u
)
loss_u = geometry.copy()
geometry.spline_data["field"] = sp.SplineDataAdaptor(
    geometry, function=loss_function_v
)
loss_v = geometry.copy()
geometry.spline_data["field"] = sp.SplineDataAdaptor(
    geometry, function=loss_function_p
)
loss_p = geometry.copy()
sp.show(
    ["U-Velocity", solution_u],
    ["V-velocity", solution_v],
    ["Pressure", solution_p],
    ["|U - U_exp|", error_field_u],
    ["|V - V_exp|", error_field_v],
    ["|P - P_exp|", error_field_p],
    ["|R(u)|", loss_u],
    ["|R(v)|", loss_v],
    ["|R(p)|", loss_p],
    knots=True,
    control_points=False,
)
