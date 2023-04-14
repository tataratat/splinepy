r"""
In this file we will prototype a collocation example (in the simplest way),
that shows how to implement a Poisson problem in the form
.. math::
        \Delta u = f
with source term f=1
"""
import gustaf as gus
import numpy as np


# Source Function
def source_function(x):
    return np.ones(x.shape[0])


# Define the Geometry
geometry = gus.BSpline(
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
geometry = gus.Bezier(
    degrees=[1, 1], control_points=[[0, 0], [1, 0], [0, 1], [1, 1]]
).bspline

# Use refinement
geometry.elevate_degrees([0, 1, 0, 1])
geometry.insert_knots(0, [0.25, 0.5, 0.75])
geometry.insert_knots(1, [0.25, 0.5, 0.75])

# Define the solution field
solution_field = gus.BSpline(
    degrees=geometry.degrees,
    control_points=np.ones((geometry.control_points.shape[0], 1)),
    knot_vectors=geometry.knot_vectors,
)

# Get greville points and geometric values
evaluation_points = solution_field.greville_abscissae
geometric_jacobians = geometry.jacobian(evaluation_points)
inverse_geometric_jacobians = np.linalg.inv(geometric_jacobians)

# Basis function values
solution_jacobians = np.zeros(
    (evaluation_points.shape[0], solution_field.control_points.shape[0], 2)
)
solution_hessians = np.zeros(
    (evaluation_points.shape[0], solution_field.control_points.shape[0], 2, 2)
)

# Fill Jacobian & Hessian
# Jacobian
basis_f_der, supports = solution_field.basis_derivative_and_support(
    evaluation_points, orders=[1, 0]
)
np.put_along_axis(solution_jacobians[:, :, 0], supports, basis_f_der, axis=1)
basis_f_der, supports = solution_field.basis_derivative_and_support(
    evaluation_points, orders=[0, 1]
)
np.put_along_axis(solution_jacobians[:, :, 1], supports, basis_f_der, axis=1)
# Hessian
basis_f_der, supports = solution_field.basis_derivative_and_support(
    evaluation_points, orders=[2, 0]
)
np.put_along_axis(solution_hessians[:, :, 0, 0], supports, basis_f_der, axis=1)
basis_f_der, supports = solution_field.basis_derivative_and_support(
    evaluation_points, orders=[1, 1]
)
np.put_along_axis(solution_hessians[:, :, 1, 0], supports, basis_f_der, axis=1)
np.put_along_axis(solution_hessians[:, :, 0, 1], supports, basis_f_der, axis=1)
basis_f_der, supports = solution_field.basis_derivative_and_support(
    evaluation_points, orders=[0, 2]
)
np.put_along_axis(solution_hessians[:, :, 1, 1], supports, basis_f_der, axis=1)

# Evaluate RHS function
physical_points = geometry.evaluate(evaluation_points)
rhs = source_function(physical_points)

# a is control_point, n is evaluation point
physical_hessians = np.einsum(
    "eli,ealk,ekj->eaij",
    inverse_geometric_jacobians,
    solution_hessians,
    inverse_geometric_jacobians,
)
physical_hessians -= np.einsum(
    "eal,elm,bm,eni,ebnk,ekj->eaij",
    solution_jacobians,
    inverse_geometric_jacobians,
    geometry.control_points,
    inverse_geometric_jacobians,
    solution_hessians,
    inverse_geometric_jacobians,
)
# calculate trace for laplacian
laplacian = np.einsum("eajj->ea", physical_hessians)

# set dirichlet values
cmr = solution_field.control_mesh_resolutions
indices = np.arange(evaluation_points.shape[0])
index_list = np.zeros((evaluation_points.shape[0]), dtype=bool)
index_list[indices % cmr[0] == 0] = True
index_list[indices % cmr[0] == cmr[0] - 1] = True
index_list[np.floor_divide(indices, cmr[0]) == 0] = True
index_list[np.floor_divide(indices, cmr[0]) == cmr[1] - 1] = True
laplacian[index_list, :] = 0
index_list = np.where(index_list)[0]
laplacian[index_list, index_list] = 1
rhs[index_list] = 0

solution_field.control_points = np.linalg.solve(laplacian, rhs).reshape(-1, 1)

geometry.splinedata["field"] = solution_field
geometry.show_options["dataname"] = "field"
geometry.show_options["cmap"] = "jet"
geometry.show_options["lighting"] = "off"
gus.show(geometry, knots=False, control_points=False)
