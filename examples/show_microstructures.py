import gustaf as gus
import numpy as np

import splinepy

# Simple Microstructure consisting Bezier-splines only
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
ms = generator.create(macro_sensitivities=True)


# Plot the ms with some of its fields
for i in range(0, len(ms.fields), 3):
    ms.spline_data["field"] = i
    ms.show_options["arrow_data"] = "field"
    ms.show_options["arrow_data_on"] = np.linspace(0, 1, 5).reshape(-1, 1)
    ms.show_options["arrow_data_scale"] = 0.1
    ctps, dim = divmod(i, ms.dim)
    gus.show(
        [f"Derivative wrt C{ctps},{dim}", ms],
        lighting="off",
        control_points=False,
        resolutions=100,
    )

# Visualize all available microtiles
micro_tiles = []
module_list = map(
    splinepy.microstructure.tiles.__dict__.get,
    splinepy.microstructure.tiles.__all__,
)
for mt in module_list:
    if hasattr(mt, "create_tile"):
        if not isinstance(mt().create_tile(), tuple):
            mt().create_tile()
            raise ValueError("Must be tuple in updated version")
        micro_tiles.append(
            [
                mt.__qualname__,
                mt().create_tile()[0],
            ]
        )
gus.show(
    *micro_tiles,
    resolutions=7,
    control_points=False,
    knots=True,
    lighting="off",
)

# Parameter spline
para_spline = splinepy.BSpline(
    degrees=[1, 1],
    knot_vectors=[[0, 0, 2, 2], [0, 0, 1, 1]],
    control_points=[[0.25], [0.3], [0.1], [0.25]],
)


def para_function(x):
    return para_spline.evaluate(x)


def para_sens_function(x):
    return splinepy.utils.data.make_matrix(
        *para_spline.basis_and_support(x),
        para_spline.cps.shape[0],
        as_array=True,
    ).reshape(x.shape[0], 1, para_spline.cps.shape[0])


# Parametrized microstructure inner and outer derivatives
generator = splinepy.microstructure.Microstructure()
generator.deformation_function = splinepy.BSpline(
    degrees=[1, 1],
    knot_vectors=[[0, 0, 1, 2, 2], [0, 0, 1, 1]],
    control_points=[[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]],
)
generator.tiling = [2, 2]
generator.parameter_sensitivity_function = para_sens_function
generator.parametrization_function = para_function
generator.microtile = splinepy.microstructure.tiles.Cross2D()
ms = generator.create(macro_sensitivities=True)
len(ms.fields)

# Plot the ms with its fields
for i in range(0, len(ms.fields), 3):
    ms.spline_data["field"] = i
    ms.show_options["arrow_data"] = "field"
    ctps, dim = divmod(i - 4, ms.dim)
    title = (
        f"Derivative with respect to Parameter-Spline Coefficient {i}"
        if i < 4
        else f"Derivative wrt C{ctps},{dim}"
    )
    gus.show([title, ms], lighting="off", control_points=False, knots=False)


def foo(x):
    return 0.1 * np.sum(x, axis=1).reshape(-1, 1)


# Cube 3D without closing face
generator = splinepy.microstructure.microstructure.Microstructure()
generator.microtile = splinepy.microstructure.tiles.HollowCube()

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


# Test new microstructure
def parameter_function_double_lattice(x):
    """
    Parametrization Function (determines thickness)
    """
    return para_s.evaluate(x)


generator = splinepy.microstructure.Microstructure()
generator.deformation_function = splinepy.Bezier(
    degrees=[1, 1], control_points=[[0, 0], [2, 0], [0, 1], [2, 1]]
)
generator.microtile = splinepy.microstructure.tiles.DoubleLattice()
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
gus.show(my_ms, knots=True, control_points=False, resolutions=2)


def parametrization_function_nut(x):
    return np.ones((len(x), 1)) * 0.3


# Structure with Nut-tile
generator = splinepy.microstructure.Microstructure()
generator.deformation_function = splinepy.Bezier(
    degrees=[1, 1], control_points=[[0, 0], [1, 0], [0, 1], [1, 1]]
)
generator.microtile = splinepy.microstructure.tiles.HollowOctagon()
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


# Classical cross tile 2D
generator = splinepy.microstructure.Microstructure()
generator.deformation_function = splinepy.Bezier(
    degrees=[1, 1], control_points=[[0, 0], [1, 0], [0, 1], [1, 1]]
)
generator.microtile = splinepy.microstructure.tiles.Cross2D()
generator.tiling = [5, 5]
generator.parametrization_function = parameter_function_double_lattice
ms = generator.create(closing_face="x", center_expansion=1.3)
generator.show(
    use_saved=True,
    knots=True,
    control_points=False,
    title="2D Crosstile Parametrized Microstructure",
)

# Micro tile with three bezier lines
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

# Cross tile 3D in a microstructure
generator = splinepy.microstructure.microstructure.Microstructure()
generator.deformation_function = splinepy.Bezier(
    degrees=[1, 1], control_points=[[0, 0], [1, 0], [0, 1], [1, 1]]
).create.extruded(extrusion_vector=[0, 0, 1])
generator.microtile = splinepy.microstructure.tiles.Cross3D()
generator.tiling = [2, 2, 3]
generator.show(
    control_points=False, resolutions=2, title="3D Crosstile Microstructure"
)

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


def foo(x):
    """
    Parametrization Function (determines thickness)
    """
    return (x[:, 0] * 0.05 + x[:, 1] * 0.05 + x[:, 2] * 0.1 + 0.1).reshape(
        -1, 1
    )


# A Parametrized microstructure and its respective inverse structure
generator = splinepy.microstructure.Microstructure()
generator.deformation_function = splinepy.Bezier(
    degrees=[1, 1], control_points=[[0, 0], [1, 0], [0, 1], [1, 1]]
).create.extruded(extrusion_vector=[0, 0, 1])
generator.microtile = splinepy.microstructure.tiles.InverseCross3D()
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
generator.microtile = splinepy.microstructure.tiles.Cross3D()
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
