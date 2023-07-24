import gustaf as gus
import numpy as np

import splinepy

# First Test
generator = splinepy.microstructure.Microstructure()
generator.deformation_function = splinepy.Bezier(
    degrees=[2, 1],
    control_points=[[0, 0], [1, 0], [2, -1], [-1, 1], [1, 1], [3, 2]],
)
generator.microtile = [
    splinepy.Bezier(
        degrees=[3], control_points=[[0, 0.5], [0.5, 1], [0.5, 0], [1, 0.5]]
    ),
    splinepy.Bezier(
        degrees=[4],
        control_points=[
            [0.5, 0],
            [0.75, 0.5],
            [0.8, 0.8],
            [0.25, 0.5],
            [0.5, 1],
        ],
    ),
]
generator.tiling = [8, 8]
generator.show(
    knots=False, control_points=False, title="2D Lattice Microstructure"
)


# Second test
def parametrization_function(x):
    return (
        0.3 - 0.4 * np.maximum(abs(0.5 - x[:, 0]), abs(0.5 - x[:, 1]))
    ).reshape(-1, 1)


para_s = splinepy.BSpline(
    knot_vectors=[
        [0, 0, 1 / 5, 2 / 5, 3 / 5, 4 / 5, 1, 1],
        [0, 0, 1 / 3, 2 / 3, 1, 1],
    ],
    degrees=[1, 1],
    control_points=[
        [0.15],
        [0.15],
        [0.15],
        [0.15],
        [0.05],
        [0.05],
        [0.05],
        [0.05],
        [0.15],
        [0.15],
        [0.15],
        [0.05],
        [0.05],
        [0.1],
        [0.15],
        [0.05],
        [0.15],
        [0.15],
        [0.2],
        [0.2],
        [0.2],
        [0.15],
        [0.15],
        [0.24],
    ],
)

# Plot all available microtiles
micro_tiles = []
module_list = map(
    splinepy.microstructure.tiles.__dict__.get,
    splinepy.microstructure.tiles.__all__,
)
for mt in module_list:
    if hasattr(mt, "create_tile"):
        micro_tiles.append(
            [
                mt.__qualname__,
                mt().create_tile(),
            ]
        )

gus.show(
    *micro_tiles,
    resolutions=7,
    control_points=False,
    knots=True,
    lighting="off"
)

# Ellipsoid
a = splinepy.microstructure.tiles.Ellipsvoid().create_tile()
for ii in range(4):
    a, b = splinepy.microstructure.tiles.Ellipsvoid().create_tile(
        parameters=np.array([[0.5, 0.3, np.deg2rad(20), np.deg2rad(10)]]),
        parameter_sensitivities=np.eye(N=1, M=4, k=ii).reshape(1, 4, 1),
    )
    surfaces = []
    s_2_extract = [5, 4, 3, 2, 1, 0]
    for i, (aa, bb) in enumerate(zip(a, b[0])):
        aa.spline_data["field"] = bb
        aa.show_options["arrow_data"] = "field"
        aa.show_options[
            "arrow_data_on"
        ] = splinepy.utils.data.cartesian_product(
            [np.linspace(0, 1, 4) for _ in range(3)]
        )
        aa.show_options["alpha"] = 0
        aa.show_options["resolutions"] = 2
        aa.show_options["control_point_ids"] = False
        surface = aa.extract.boundaries()[s_2_extract[i]]
        surface.show_options["control_points"] = False
        surface.show_options["c"] = "grey"
        surface.show_options["resolutions"] = 20
        surfaces.append(surface)
    camera = dict(
        position=(1.9, 1.3, 3),
        focal_point=(0.5, 0.5, 0.5),
        viewup=(-0.1, 0.95, -0.3),
        distance=3.33943,
        clipping_range=(1.37, 5.0),
    )
    gus.show(surfaces + a, cam=camera, alpha=0.2)


# Cubevoid
a = splinepy.microstructure.tiles.Cubevoid().create_tile()
for ii in range(4):
    a, b = splinepy.microstructure.tiles.Cubevoid().create_tile(
        parameters=np.array([[0.5, 0.3, np.deg2rad(20), np.deg2rad(10)]]),
        parameter_sensitivities=np.eye(N=1, M=4, k=ii).reshape(1, 4, 1),
    )
    surfaces = []
    s_2_extract = [5, 4, 3, 2, 1, 0]
    for i, (aa, bb) in enumerate(zip(a, b[0])):
        aa.spline_data["field"] = bb
        aa.show_options["arrow_data"] = "field"
        aa.show_options[
            "arrow_data_on"
        ] = splinepy.utils.data.cartesian_product(
            [np.linspace(0, 1, 4) for _ in range(3)]
        )
        aa.show_options["alpha"] = 0
        aa.show_options["resolutions"] = 2
        aa.show_options["control_point_ids"] = False
        surface = aa.extract.boundaries()[s_2_extract[i]]
        surface.show_options["control_points"] = False
        surface.show_options["c"] = "grey"
        surface.show_options["resolutions"] = 20
        surfaces.append(surface)
    camera = dict(
        position=(1.9, 1.3, 3),
        focal_point=(0.5, 0.5, 0.5),
        viewup=(-0.1, 0.95, -0.3),
        distance=3.33943,
        clipping_range=(1.37, 5.0),
    )
    gus.show(surfaces + a, cam=camera)


def foo(x):
    return 0.1 * np.sum(x, axis=1).reshape(-1, 1)


# Cube 3D without closing face
generator = splinepy.microstructure.microstructure.Microstructure()
generator.microtile = splinepy.microstructure.tiles.Cube3D()

generator.parametrization_function = foo
generator.deformation_function = splinepy.Bezier(
    degrees=[1, 1], control_points=[[0, 0], [1, 0], [0, 1], [1, 1]]
).create.extruded(extrusion_vector=[0, 0, 1])
generator.tiling = [3, 3, 2]
generator.show(
    knots=True,
    control_points=False,
    title="3D Cube Microstructure",
    resolutions=3,
    alpha=0.05,
)


# Test new microstructure
def parameter_function_double_lattice(x):
    """
    Parametrization Function (determines thickness)
    """
    return para_s.evaluate(x)


generator = splinepy.microstructure.Microstructure()
# outer geometry
generator.deformation_function = splinepy.Bezier(
    degrees=[1, 1], control_points=[[0, 0], [2, 0], [0, 1], [2, 1]]
)
generator.microtile = splinepy.microstructure.tiles.DoubleLatticeTile()
# how many structures should be inside the cube
generator.tiling = [12, 6]
generator.parametrization_function = parameter_function_double_lattice
my_ms = generator.create(contact_length=0.4)
generator.show(
    use_saved=True,
    knots=True,
    control_points=False,
    title="2D Nuttile Parametrized Microstructure",
    contact_length=0.4,
    resolutions=2,
)
gus.show(my_ms, knots=True, control_points=False, resolution=2)


def parametrization_function_nut(x):
    return np.array([0.3]).reshape(-1, 1)


# Test new microstructure

generator = splinepy.microstructure.Microstructure()
# outer geometry
generator.deformation_function = splinepy.Bezier(
    degrees=[1, 1], control_points=[[0, 0], [1, 0], [0, 1], [1, 1]]
)
generator.microtile = splinepy.microstructure.tiles.NutTile2D()
# how many structures should be inside the cube
generator.tiling = [10, 10]
generator.parametrization_function = parametrization_function_nut
my_ms = generator.create(closing_face="x", contact_length=0.4)
generator.show(
    use_saved=True,
    knots=True,
    control_points=False,
    title="2D Nuttile Parametrized Microstructure",
    contact_length=0.4,
    resolutions=2,
)


# Second test
generator = splinepy.microstructure.Microstructure()
generator.deformation_function = splinepy.Bezier(
    degrees=[1, 1], control_points=[[0, 0], [1, 0], [0, 1], [1, 1]]
)
generator.microtile = splinepy.microstructure.tiles.CrossTile2D()
generator.tiling = [5, 5]
generator.parametrization_function = parametrization_function
ms = generator.create(closing_face="x", center_expansion=1.3)
generator.show(
    use_saved=True,
    knots=True,
    control_points=False,
    title="2D Crosstile Parametrized Microstructure",
)

# Third test
generator = splinepy.microstructure.Microstructure()
generator.deformation_function = splinepy.Bezier(
    degrees=[1, 1], control_points=[[0, 0], [1, 0], [0, 1], [1, 1]]
).create.extruded(extrusion_vector=[0, 0, 1])
generator.microtile = [
    splinepy.Bezier(
        degrees=[3],
        control_points=[
            [0, 0.5, 0.5],
            [0.5, 1, 0.5],
            [0.5, 0, 0.5],
            [1, 0.5, 0.5],
        ],
    ),
    splinepy.Bezier(
        degrees=[3],
        control_points=[
            [0.5, 0.5, 0.0],
            [0.5, 0, 0.5],
            [0.5, 1.0, 0.5],
            [0.5, 0.5, 1.0],
        ],
    ),
    splinepy.Bezier(
        degrees=[3],
        control_points=[
            [0.5, 0, 0.5],
            [1, 0.5, 0.5],
            [0, 0.5, 0.5],
            [0.5, 1, 0.5],
        ],
    ),
]
generator.tiling = [1, 2, 3]
generator.show(
    knots=False, control_points=False, title="3D Lattice Microstructure"
)

# Fourth test
generator = splinepy.microstructure.microstructure.Microstructure()
generator.deformation_function = splinepy.Bezier(
    degrees=[1, 1], control_points=[[0, 0], [1, 0], [0, 1], [1, 1]]
).create.extruded(extrusion_vector=[0, 0, 1])
generator.microtile = splinepy.microstructure.tiles.CrossTile3D()
generator.tiling = [2, 2, 3]
generator.show(
    control_points=False, resolutions=2, title="3D Crosstile Microstructure"
)

# Fifth test
# Non-uniform tiling
generator = splinepy.microstructure.microstructure.Microstructure()
generator.deformation_function = splinepy.BSpline(
    degrees=[2, 1],
    control_points=[
        [0, 0],
        [0.1, 0],
        [0.2, 0],
        [1, 0],
        [0, 1],
        [0.1, 1],
        [0.2, 2],
        [1, 3],
    ],
    knot_vectors=[[0, 0, 0, 0.15, 1, 1, 1], [0, 0, 1, 1]],
)
generator.microtile = [
    splinepy.Bezier(
        degrees=[3], control_points=[[0, 0.5], [0.5, 1], [0.5, 0], [1, 0.5]]
    ),
    splinepy.Bezier(
        degrees=[3], control_points=[[0.5, 0], [1, 0.5], [0, 0.5], [0.5, 1]]
    ),
]
generator.tiling = [5, 1]
generator.show(
    knot_span_wise=False,
    control_points=False,
    resolutions=20,
    title="2D Lattice with global tiling",
)


# Sixth test
# A Parametrized microstructure and its respective inverse structure
def foo(x):
    """
    Parametrization Function (determines thickness)
    """
    return (x[:, 0] * 0.05 + x[:, 1] * 0.05 + x[:, 2] * 0.1 + 0.1).reshape(
        -1, 1
    )


generator = splinepy.microstructure.Microstructure()
generator.deformation_function = splinepy.Bezier(
    degrees=[1, 1], control_points=[[0, 0], [1, 0], [0, 1], [1, 1]]
).create.extruded(extrusion_vector=[0, 0, 1])
generator.microtile = splinepy.microstructure.tiles.InverseCrossTile3D()
generator.tiling = [3, 3, 5]
generator.parametrization_function = foo

inverse_microstructure = generator.create(
    closing_face="z", seperator_distance=0.4, center_expansion=1.3
)

# Plot the results
_, showables_inverse = generator.show(
    closing_face="z",
    seperator_distance=0.4,
    center_expansion=1.3,
    title="Parametrized Inverse Microstructure",
    control_points=False,
    knots=True,
    return_showable_list=True,
    resolutions=5,
    c="hotpink",
)

# Corresponding Structure
generator.microtile = splinepy.microstructure.tiles.CrossTile3D()
microstructure = generator.create(
    closing_face="z", seperator_distance=0.4, center_expansion=1.3
)
_, showables = generator.show(
    closing_face="z",
    center_expansion=1.3,
    title="Parametrized Microstructure",
    control_points=False,
    knots=True,
    return_showable_list=True,
    resolutions=5,
)

gus.show(
    [*showables_inverse, *showables],
    title="Parametrized Microstructure and its inverse",
)
