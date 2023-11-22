import numpy as np
import vedo.pyplot as plot

import splinepy as sp

# Composition
inner_function_polynomomial = sp.Bezier(
    degrees=[3],
    control_points=[[0.2, 0.2], [0.4, 0.2], [0.4, 0.8], [0.8, 0.8]],
)
outer_function_polynomial = sp.Bezier(
    degrees=[2, 2],
    control_points=[
        [0.0, 0.0],
        [0.5, 0.2],
        [1.0, 0.0],
        [0.0, 0.6],
        [0.7, 0.7],
        [1.4, 0.9],
        [-0.3, 1.0],
        [0.4, 1.2],
        [1.0, 2.0],
    ],
)

# sp.show(inner_function_polynomomial,[outer_function_polynomial, outer_function_polynomial.compose(inner_function_polynomomial)])


# Rational Surface
inner_function_rational = sp.RationalBezier(
    degrees=[2],
    control_points=[
        [0.2, 0.0],
        [0.5, 1.0],
        [0.8, 0.0]
    ],
    weights=[0.333, 1.0, 0.3333]
)

outer_function_rational =  sp.Bezier(
    degrees=[1],
    control_points=[[0.5,0],[1.5,0]]
    ).create.revolved(
        center=[0,0], angle=90, n_knot_spans=None, degree=True
    )

# sp.show(
#     inner_function_rational,
#     outer_function_rational,
#     [
#         outer_function_rational,
#         outer_function_rational.compose(inner_function_rational)
#     ]
# )

spline = sp.Bezier(
    degrees=[2], control_points=[[0.2, 0.9], [0.8, 0.7], [0.8, 0.2]]
).bspline
inner_function_polynomomial.show_options["c"] = "k"

sp.io.svg.export("test.svg", inner_function_polynomomial)

spline = sp.helpme.create.arc(10,90).nurbs
spline.insert_knots(0,[0.5])
spline.elevate_degrees([0,0])
sp.io.svg.export("test_rational.svg", spline)

quit()

# Create spline

nurbs = sp.NURBS(
    degrees=[2],
    control_points=[[0, 0], [1, 1], [3, 1], [4, 0], [5, 1]],
    knot_vectors=[[0, 0, 0, 0.25, 0.5, 1, 1, 1]],
    weights=[1, 1, 1, 0.8, 1],
)
nurbs.show()


def plot_basis_functions(spline, sample_rate=100, return_fig=False):
    """
    Plot basis functions using vedo plot.
    """

    if spline.para_dim != 1:
        raise ValueError("Only 1D basis functions supported")

    n_basis_functions = spline.control_points.shape[0]
    queries = np.linspace(*spline.parametric_bounds, sample_rate)
    basis, supports = spline.basis_and_support(queries)
    basis_function_matrix = np.zeros((sample_rate, n_basis_functions))
    np.put_along_axis(basis_function_matrix, supports, basis, axis=1)

    fig = plot.plot(
        queries,
        basis_function_matrix[:, 0],
        ylim=(0.0, 1.01),
        label="B0," + str(spline.degrees[0]),
        title="Basis functions",
        c=0,
    )
    for i in range(n_basis_functions):
        fig += plot.plot(
            queries,
            basis_function_matrix[:, i],
            ylim=(0.0, 1.01),
            label="B" + str(i) + "," + str(spline.degrees[0]),
            c=i,
        )

    fig.add_legend("top-right", s=0.7)
    if return_fig:
        return fig
    else:
        fig.show()


plot_basis_functions(nurbs, 1001, False)
