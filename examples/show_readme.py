import gustaf as gus
import numpy as np
import vedo

import splinepy

# for consistent viewing for both dark and light modes
bg = "grey"

### 1.

# Initialize nurbs with any array-like input
nurbs = splinepy.NURBS(
    degrees=[2, 1],
    knot_vectors=[
        [0, 0, 0, 1, 1, 1],
        [0, 0, 1, 1],
    ],
    control_points=[
        [-1.0, 0.0],
        [-1.0, 1.0],
        [0.0, 1.0],
        [-2.0, 0.0],
        [-2.0, 2.0],
        [0.0, 2.0],
    ],
    weights=[
        [1.0],
        [2**-0.5],
        [1.0],
        [1.0],
        [2**-0.5],
        [1.0],
    ],
)

# vizusalize
nurbs.show(background=bg)


### 2.

# start with a copy of the original spline
modified = nurbs.copy()

# manipulate control points
# 1. all at once
modified.control_points /= 2.0
# 2. indexwise (flat indexing)
modified.control_points[[3, 4, 5]] *= [1.3, 2.0]
# 3. with grid-like indexing using multi_index helper
multi_index = modified.multi_index
modified.control_points[multi_index[0, 1]] = [-1.5, -0.3]
modified.control_points[multi_index[2, :]] += [2.0, 0.1]

# modified.show()  # visualize Nr. 1
view1 = modified.copy()

# elevate degrees and insert knots
modified.elevate_degrees([0, 1])
# modified.show()  # visualize Nr. 2
view2 = modified.copy()

modified.insert_knots(1, [0.5])
# modified.show()  # visualize Nr. 3
view3 = modified.copy()

splinepy.show(
    ["1", view1],
    ["2", view2],
    ["3", view3],
    background=bg,
)


### 3.

# first, create parametric coordinate queries
queries = [
    [0.1, 0.2],  # first query
    [0.4, 0.5],  # second query
    [0.1156, 0.9091],  # third query
]

# evaluate basis, spline and derivatives.
# for derivatives, specify order per parametric dimension.
basis = nurbs.basis(queries)
basis_derivative = nurbs.basis_derivative(queries, [1, 1])
physical_coordinates = nurbs.evaluate(queries)
physical_derivatives = nurbs.derivative(queries, [2, 0])

p_basis0 = nurbs.basis(queries, nthreads=2)
# or
splinepy.settings.NTHREADS = 3
p_basis1 = nurbs.basis(queries)

para_queries = []
phys_queries = []
colors = ["red", "yellow", "blue"]
for q, pc, c in zip(queries, physical_coordinates, colors):
    para_q, phys_q = gus.Vertices([q]), gus.Vertices([pc])
    para_q.show_options["c"] = c
    para_q.show_options["r"] = 10
    phys_q.show_options["c"] = c
    phys_q.show_options["r"] = 10
    para_queries.append(para_q)
    phys_queries.append(phys_q)

splinepy.show(
    ["Parametric", *para_queries, nurbs.create.parametric_view()],
    ["Physical", *phys_queries, nurbs],
    background=bg,
)

# see docs for options
para_coordinates = nurbs.proximities(physical_coordinates)
assert np.allclose(queries, para_coordinates)

(
    parametric_coordinates,
    physical_coordindates,
    physical_difference,
    distance,
    convergence_norm,
    first_derivatives,
    second_derivatives,
) = nurbs.proximities(physical_coordinates, return_verbose=True)


### 4.1

# basic splines
box = splinepy.helpme.create.box(1, 2, 3)  # length per dim
disk = splinepy.helpme.create.disk(outer_radius=3, inner_radius=2, angle=256)
torus = splinepy.helpme.create.torus(torus_radius=3, section_outer_radius=1.5)

splinepy.show(
    ["box", box],
    ["disk", disk],
    ["torus", torus],
    background=bg,
    resolutions=201,
)


# based on existing splines
extruded = nurbs.create.extruded(extrusion_vector=[1, 2, 3])
revolved = nurbs.create.revolved(axis=[1, 0, 0], center=[-1, -1, 0], angle=50)

splinepy.show(["extruded", extruded], ["revolved", revolved], background=bg)

### 4.2

# extract meshes as gustaf objects
control_mesh = nurbs.extract.control_mesh()
control_points = nurbs.extract.control_points()
mesh = nurbs.extract.faces([201, 33])  # specify sample resolutions

control_mesh_as_edges = control_mesh.to_edges()
control_mesh_as_edges.show_options.update(c="black", lw=8)
splinepy.show(
    ["control mesh", control_mesh_as_edges],
    ["control points", control_points],
    ["spline", mesh],
    background=bg,
)

# extract splines
boundaries = nurbs.extract.boundaries()
partial = nurbs.extract.spline(0, [0.5, 0.78])
partial_partial = nurbs.extract.spline(0, [0.1, 0.3]).extract.spline(
    1, [0.65, 0.9]
)
bases = nurbs.extract.bases()  # basis functions as splines
# insert knots to increase number of bezier patches
inserted = nurbs.copy()
inserted.insert_knots(0, [0.13, 0.87])
beziers_patches = inserted.extract.beziers()

splinepy.show(
    ["boundaries and part of splines", boundaries, partial, partial_partial],
    ["beziers", beziers_patches],
    ["bases", bases],
    background=bg,
)

### 4.3
# create gustaf mesh using extract.spline()
# or use gustaf's io functions (gustaf.io)
mesh = splinepy.helpme.create.torus(2, 1).extract.faces([100, 100, 100])

# initialize ffd and move control points
ffd = splinepy.FFD(mesh=mesh)
multi_index = ffd.spline.multi_index
ffd.spline.control_points[multi_index[-1, :, -1]] += [3, 0.5, 0.1]

cam = {
    "position": (-9.96434, -7.98279, 7.95655),
    "focal_point": (0.361918, -0.129347, 0.348385),
    "viewup": (0.440934, 0.254149, 0.860805),
    "roll": 72.6622,
    "distance": 15.0397,
    "clipping_range": (7.14564, 27.5698),
}

ffd.show(background=bg, cam=cam)

# get deformed mesh - FFD.mesh attribute deforms mesh before returning
deformed = ffd.mesh


### 4.4
data = [
    [-0.955, 0.293],
    [-0.707, 0.707],
    [-0.293, 0.955],
    [-1.911, 0.587],
    [-1.414, 1.414],
    [-0.587, 1.911],
]

curve, residual_curve = splinepy.helpme.fit.curve(data, degree=2)
# you can also use any existing spline's basis
surface, residual_surface = splinepy.helpme.fit.surface(
    data, size=[3, 2], fitting_spline=nurbs
)

# set visuals for data
d = gus.Vertices(data)
d.show_options.update(c="blue", r=15)

splinepy.show(
    ["curve fit", d, curve],
    ["surface fit", d, surface],
    background=bg,
)


### 4.5
# create solution spline
solution_field = nurbs.create.embedded(1)

# refine
solution_field.elevate_degrees([0, 1])
solution_field.uniform_refine(n_knots=[4, 4])

# create matrix using mapper
# collocation points at greville abcissae
mapper = solution_field.mapper(reference=nurbs)
laplacian, support = mapper.basis_laplacian_and_support(
    solution_field.greville_abscissae()
)
laplacian_matrix = splinepy.utils.data.make_matrix(
    laplacian,
    support,
    n_cols=solution_field.control_points.shape[0],
    as_array=True,
)

# for this one, we want zero to the mittle value to be zero to highlight
# sparse matrix
vmin = laplacian_matrix.min()
vmax = laplacian_matrix.max()
v_abs = max(vmin, vmax)
vmin = v_abs * np.sign(vmin)
vmax = v_abs * np.sign(vmax)

vedo_mat = vedo.pyplot.matrix(
    laplacian_matrix, cmap="seismic", vmin=vmin, vmax=vmax
)
splinepy.show(vedo_mat, background=bg)


### 5.

splinepy.microstructure.tiles.show(alpha=0.7, background=bg)

# get specific dimensions as dict
para_2_dim_2 = splinepy.microstructure.tiles.by_dim(para_dim=2, dim=2)
dim_2 = splinepy.microstructure.tiles.by_dim(dim=2)


# create microstructure generator
microstructure = splinepy.Microstructure()
# set outer spline and a (micro) tile
microstructure.deformation_function = nurbs
microstructure.microtile = splinepy.microstructure.tiles.get("Cross2D")
# tiling determines tile resolutions within each bezier patch
microstructure.tiling = [5, 3]

microstructure.show(background=bg)

# extract only generated parts as multipatch
generated = microstructure.create()


### 6.
# use previously generated microtiles
interface_info_array = generated.interfaces


# Mark boundaries to set boundary conditions
# In case of micro structure, you can use outer spline's boundary
def is_left_bdr(x):
    left = nurbs.extract.boundaries()[3]
    return (left.proximities(x, return_verbose=True)[3] < 1e-8).ravel()


generated.boundary_from_function(is_left_bdr, boundary_id=5)

splinepy.show(
    ["All", generated],
    ["Boundaries", generated.boundary_multipatch()],
    ["Boundary 5", generated.boundary_multipatch(5)],
    control_points=False,
    background=bg,
)

# export for the solver
splinepy.io.gismo.export("microstructure.xml", generated)


### 7.
# export
splinepy.io.mfem.export("quarter_circle.mesh", nurbs)

# load
quarter_circle = splinepy.io.mfem.load("quarter_circle.mesh")

# svg export
splinepy.io.svg.export("nurbs.svg", nurbs, background=bg)
