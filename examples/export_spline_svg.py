"""
Export spline plots as svg files
"""

import numpy as np

import splinepy as spp

np_rng = np.random.default_rng()

box = spp.helpme.create.box(10, 5).bspline
box.elevate_degrees([0, 1])
box.insert_knots(0, [0.5])

# Create a simple curve
spline_curve = spp.Bezier(degrees=[2], control_points=[[0, 0], [2, 0], [2, 2]])
spp.io.svg.export("spline_curve_1.svg", spline_curve)

# Make it more complex
rational_spline_curve = spline_curve.nurbs
rational_spline_curve.elevate_degrees([0, 0])
rational_spline_curve.insert_knots(0, [0.4, 0.6])
rational_spline_curve.weights = np_rng.random(
    rational_spline_curve.cps.shape[0]
)
spp.io.svg.export("spline_curve_3.svg", rational_spline_curve)

# Export some surfaces
circle = spp.helpme.create.surface_circle(1)
spp.io.svg.export("spline_circle_1.svg", circle)

circle = circle.nurbs
circle.insert_knots(0, [0.333, 0.666])
circle.insert_knots(1, [0.333, 0.666])
circle.spline_data["field"] = circle
# svg export also respects many of the show_option configurations
circle.show_options["arrow_data"] = "field"
circle.show_options["arrow_data_scale"] = 0.4
circle.show_options["arrow_data_on"] = spp.utils.data.cartesian_product(
    [np.linspace(0, 1, 10) for _ in range(2)]
)
spp.io.svg.export("spline_circle_2.svg", circle, box_margins=1.0)

# Create a spline with a scalar field
distorted_rectangle = spp.Bezier(
    degrees=[2, 2],
    control_points=[
        [0.0, 0.0],
        [1.0, 0.5],
        [2.0, 0.0],
        [0.5, 0.8],
        [1.2, 1.0],
        [2.3, 1.0],
        [0.0, 2.0],
        [1.0, 1.5],
        [2.0, 1.8],
    ],
).bspline
distorted_rectangle.insert_knots(0, [0.2, 0.4, 0.6, 0.8])
distorted_rectangle.insert_knots(1, [0.25, 0.5, 0.75])

# Extract a field (jacobian determinant spline projection)
jacobian_determinant_spline = distorted_rectangle.create.determinant_spline()
distorted_rectangle.spline_data["jac_dets"] = jacobian_determinant_spline
distorted_rectangle.show_options["data"] = "jac_dets"
spp.io.svg.export(
    "spline_rectangle_with_field.svg", distorted_rectangle, box_margins=1.0
)
