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
n_refine = 50

# Parameters (steel)
Youngs_modulus = 210e9
poisson_ratio = 0.3
density = 7850.
lame_mu = 0.5 / (1 + poisson_ratio) * Youngs_modulus
lame_lambda = 0.5 / ((1 + poisson_ratio) *
                     (1 - 2 * poisson_ratio)) * Youngs_modulus


# Source Function
def source_function(x):
    return np.array([[0, -9.81]]).repeat(x.shape[0], axis=0)*density


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
geometry = sp.helpme.create.box(2, 0.02).bspline
geometry.elevate_degrees([0, 1, 0, 0, 1, 0])
new_knots = np.linspace(1 / n_refine, 1, n_refine, endpoint=False)
geometry.insert_knots(0, new_knots)
geometry.insert_knots(1, new_knots[1::20])

# Define the solution field
solution_field = sp.BSpline(
    degrees=geometry.degrees,
    control_points=np.ones((geometry.control_points.shape[0], 2)),
    knot_vectors=geometry.knot_vectors,
)


# # Use refinement
# solution_field.elevate_degrees([])
# new_knots = np.linspace(1 / n_refine, 1, n_refine, endpoint=False)
# solution_field.insert_knots(0, new_knots)
# solution_field.insert_knots(1, new_knots)


# Get greville points and geometric values
evaluation_points = solution_field.greville_abscissae()
mapper = solution_field.mapper(reference=geometry)
sp_hessian, supports = mapper.basis_hessian_and_support(evaluation_points)

# Calculate matrix coefficients (eye as kronecker is pretty inefficient, ...)
sigma = np.einsum("aqjj,ik->aqik", sp_hessian, np.eye((2))) * lame_mu
sigma += np.einsum("aqki->aqik", sp_hessian) * (lame_mu + lame_lambda)

# Calculate indices
indices = np.reshape(supports.repeat(
    4, axis=1) * 2, (supports.shape[0], supports.shape[1], 2, 2))
indices[:, :, :, 1] += 1

# This is not yet final
sigma_list = np.einsum("abcd->acbd", sigma).reshape(supports.shape[0] * 2, -1)
indices_list = np.einsum(
    "abcd->acbd", indices).reshape(supports.shape[0] * 2, -1)


sp_hessian = sp.utils.data.make_matrix(
    sigma_list,
    indices_list,
    solution_field.cps.size,
    as_array=False,
)

# Evaluate RHS function
physical_points = geometry.evaluate(evaluation_points)
rhs = source_function(physical_points).reshape(-1)

# Set dirichlet values (identify boundary nodes and set matrix to identity rhs
# to 0)
boundary_dof_ids = solution_field.multi_index[0, :].repeat(2) * 2
boundary_dof_ids[1::2] += 1

sp_hessian[boundary_dof_ids] *= 0  # this does not change sparsity pattern
sp_hessian[boundary_dof_ids, boundary_dof_ids] = 1

rhs[boundary_dof_ids] = 0.0

# Solve linear system
solution_field.control_points[:] = linalg.spsolve(-sp_hessian, rhs).reshape(
    -1, 2
)

# Errors


def error_strong_form(data, on):
    field_hessian = mapper.hessian(on)
    grad_sig = np.einsum("qijj->qi", field_hessian) * lame_mu
    grad_sig += np.einsum("qjij->qi", field_hessian) * (lame_mu + lame_lambda)
    return np.linalg.norm(grad_sig + source_function(on), axis=1)


errors = error_strong_form(1, evaluation_points)
print(f"Mean error : {np.mean(errors)}")
print(f"Max error  : {np.max(errors)}")


# Plot geometry and field
geometry.spline_data["field"] = solution_field
geometry.show_options["arrow_data"] = "field"
geometry.show_options["arrow_data_on"] = sp.utils.data.cartesian_product(
    [np.linspace(0, 1, 20)] * 2
)
geometry.show_options["arrow_data_scale"] = 1e3
geometry.show_options["scalarbar"] = True
geometry.show_options["cmap"] = "jet"
geometry.show_options["lighting"] = "off"
geometry.show_options["scalarbar"] = True
# geometry.show(knots=True, control_points=False)

# Error plot
geo_para = geometry.create.parametric_view()
geo_para.spline_data["error"] = sp.SplineDataAdaptor(
    geometry, function=error_strong_form)
geo_para.show_options["data_name"] = "error"
geo_para.show_options["cmap"] = "jet"
geo_para.show_options["lighting"] = "off"
geo_para.show_options["control_points"] = False
geo_para.show_options["knots"] = False
# geo_para.show()


deformation_plot = geometry.copy()
deformation_plot.cps += solution_field.cps * 1e6
deformation_plot.show(control_points=False, knots=False)
