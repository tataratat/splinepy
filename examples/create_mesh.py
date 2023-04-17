import subprocess

import gustaf as gus
import numpy as np
import scipy

gismo_options = [
    {
        "tag": "Function",
        "attributes": {"type": "FunctionExpr", "id": "1", "dim": "2"},
        "text": "\n    ",
        "children": [
            {
                "tag": "c",
                "attributes": {"index": "0"},
                "text": "0",
            },
            {
                "tag": "c",
                "attributes": {"index": "1"},
                "text": "0",
            },
        ],
    },
    {
        "tag": "boundaryConditions",
        "attributes": {"multipatch": "0", "id": "2"},
        "children": [
            {
                "tag": "Function",
                "attributes": {
                    "type": "FunctionExpr",
                    "dim": "2",
                    "index": "0",
                },
                "text": "0",
            },
            {
                "tag": "Function",
                "attributes": {
                    "type": "FunctionExpr",
                    "dim": "2",
                    "index": "1",
                },
                "text": "\n      -1\n    ",
            },
            {
                "tag": "Function",
                "attributes": {
                    "type": "FunctionExpr",
                    "dim": "2",
                    "index": "2",
                    "c": "2",
                },
                "children": [
                    {
                        "tag": "c",
                        "attributes": {"index": "0"},
                        "text": "0",
                    },
                    {
                        "tag": "c",
                        "attributes": {"index": "1"},
                        "text": "-10000",
                    },
                ],
            },
            {
                "tag": "bc",
                "attributes": {
                    "type": "Dirichlet",
                    "function": "0",
                    "unknown": "0",
                    "name": "BID2",
                },
            },
            {
                "tag": "bc",
                "attributes": {
                    "type": "Dirichlet",
                    "function": "0",
                    "unknown": "1",
                    "name": "BID2",
                },
            },
            {
                "tag": "bc",
                "attributes": {
                    "type": "Neumann",
                    "function": "2",
                    "unknown": "0",
                    "name": "BID6",
                },
            },
        ],
    },
    {
        "tag": "OptionList",
        "attributes": {"id": "3"},
        "text": "\n    ",
        "children": [
            {
                "tag": "int",
                "attributes": {
                    "label": "DirichletStrategy",
                    "desc": "Method for enforcement of Dirichlet BCs [11..14]",
                    "value": "11",
                },
            },
            {
                "tag": "int",
                "attributes": {
                    "label": "DirichletValues",
                    "desc": "Method for computation of Dirichlet DoF values [1"
                    "00..103]",
                    "value": "101",
                },
            },
            {
                "tag": "int",
                "attributes": {
                    "label": "InterfaceStrategy",
                    "desc": "Method of treatment of patch interfaces [0..3]",
                    "value": "1",
                },
            },
            {
                "tag": "real",
                "attributes": {
                    "label": "bdA",
                    "desc": "Estimated nonzeros per column of the matrix: bdA*"
                    "deg + bdB",
                    "value": "2",
                },
            },
            {
                "tag": "int",
                "attributes": {
                    "label": "bdB",
                    "desc": "Estimated nonzeros per column of the matrix: bdA*"
                    "deg + bdB",
                    "value": "1",
                },
            },
            {
                "tag": "real",
                "attributes": {
                    "label": "bdO",
                    "desc": "Overhead of sparse mem. allocation: (1+bdO)(bdA*d"
                    "eg + bdB) [0..1]",
                    "value": "0.333",
                },
            },
            {
                "tag": "real",
                "attributes": {
                    "label": "quA",
                    "desc": "Number of quadrature points: quA*deg + quB",
                    "value": "1",
                },
            },
            {
                "tag": "int",
                "attributes": {
                    "label": "quB",
                    "desc": "Number of quadrature points: quA*deg + quB",
                    "value": "1",
                },
            },
            {
                "tag": "int",
                "attributes": {
                    "label": "quRule",
                    "desc": "Quadrature rule [1:GaussLegendre,2:GaussLobatto]",
                    "value": "1",
                },
            },
        ],
    },
]
gus.settings.NTHREADS = 8
length = 2
height = 1
tiling_x = 24
tiling_y = 12
nthreads = 24
load_b_tiles = 2
para_deg = 2

filename = "lattice_structure_" + str(tiling_x) + "x" + str(tiling_y) + ".xml"


def prepare_microstructure(parameters):
    generator = gus.spline.microstructure.Microstructure()
    deformation_function = gus.Bezier(
        degrees=[1, 1],
        control_points=[[0, 0], [length, 0], [0, height], [length, height]],
    )
    deformation_function = deformation_function.bspline
    deformation_function.insert_knots(
        0, [i * (1 / tiling_x) for i in range(1, tiling_x)]
    )
    deformation_function.insert_knots(
        1, [i * (1 / tiling_y) for i in range(1, tiling_y)]
    )
    generator.deformation_function = deformation_function
    generator.tiling = [1, 1]
    generator.microtile = gus.spline.microstructure.tiles.DoubleLatticeTile()

    iknots = [
        i * (1 / (tiling_x - para_deg)) for i in range(1, tiling_x - para_deg)
    ]
    jknots = [
        i * (1 / (tiling_y - para_deg)) for i in range(1, tiling_y - para_deg)
    ]
    pdegrees = [2, 2]
    pknot_vectors = [
        [0] * (para_deg + 1) + k + [1] * (para_deg + 1)
        for k in [iknots, jknots]
    ]
    pcontrol_points = parameters.reshape(tiling_x * tiling_y, 1)
    parameter_spline = gus.BSpline(
        degrees=pdegrees,
        knot_vectors=pknot_vectors,
        control_points=pcontrol_points,
    )

    def foo(x):
        """
        Parametrization Function (determines thickness)
        """
        return parameter_spline.evaluate(x)

    def foo_deriv(x):
        basis_function_matrix = np.zeros(
            (x.shape[0], parameter_spline.control_points.shape[0])
        )
        basis_functions, support = parameter_spline.basis_and_support(x)
        np.put_along_axis(
            basis_function_matrix, support, basis_functions, axis=1
        )
        return basis_function_matrix.reshape(1, 1, -1)

    generator.parametrization_function = foo
    generator.parameter_sensitivity_function = foo_deriv
    my_ms, my_ms_der = generator.create(contact_length=0.5)

    def identifier_function(deformation_function, face_id):
        boundary_spline = deformation_function.extract_boundaries(face_id)[0]

        def identifier_function(x):
            distance_2_boundary = boundary_spline.proximities(
                queries=x, initial_guess_sample_resolutions=4, tolerance=1e-9
            )[3]
            return distance_2_boundary.flatten() < 1e-8

        return identifier_function

    def identifier_function_neumann(x):
        return x[:, 0] >= (tiling_x - load_b_tiles) / tiling_x * length - 1e-12

    multipatch = gus.spline.splinepy.Multipatch(my_ms)
    multipatch.add_fields(
        *my_ms_der, check_compliance=True, check_conformity=True
    )
    multipatch.determine_interfaces()
    multipatch.boundary_from_function(
        identifier_function(generator.deformation_function, 0)
    )
    multipatch.boundary_from_function(
        identifier_function(generator.deformation_function, 1)
    )
    multipatch.boundary_from_function(
        identifier_function(generator.deformation_function, 2)
    )
    multipatch.boundary_from_function(
        identifier_function(generator.deformation_function, 3)
    )

    multipatch.boundary_from_function(identifier_function_neumann, mask=[5])
    gus.spline.io.gismo.export(
        filename,
        multipatch=multipatch,
        options=gismo_options,
        export_fields=True,
    )


def read_jacobians():
    jacs = np.genfromtxt(fname="sensitivities.out")
    return jacs


def read_objective_function():
    obj_val = float(np.genfromtxt(fname="objective_function.out"))
    return obj_val


def run_gismo(sensitivities=False, plot=False, refinement=False):
    process_call = [
        "./linear_elasticity_expressions",
        "-f",
        filename,
        "-p",
        str(nthreads),
        "--compute-objective-function",
        "--output-to-file",
    ]
    if sensitivities:
        process_call += [
            "--compute-sensitivities",
            "-x",
            filename + ".fields.xml",
        ]
    if plot:
        process_call += ["--plot"]
    if refinement:
        process_call += [
            "-r",
            "1",
        ]
    text = subprocess.run(process_call, capture_output=True, encoding="ascii")
    return text.returncode


def evaluate_iteration(x):
    prepare_microstructure(x)
    run_gismo(sensitivities=False, refinement=True)
    return read_objective_function()


def evaluate_jacobian(x):
    prepare_microstructure(x)
    run_gismo(sensitivities=True, refinement=True)
    return read_jacobians()


def call_back_optimization(x):
    print(
        f"{read_objective_function()}:"
        + "  ".join(str(xx) for xx in x.tolist())
    )


def main():
    initial_guess = np.ones((tiling_x * tiling_y, 1)) * 0.10

    # Mass constraint
    A = np.ones((1, tiling_x * tiling_y))
    c = 0
    d = tiling_x * tiling_y * 0.10
    C2 = scipy.optimize.LinearConstraint(A, c, d)

    optim = scipy.optimize.minimize(
        evaluate_iteration,
        initial_guess,
        method="SLSQP",
        jac=evaluate_jacobian,
        bounds=[(0.0111, 0.249) for _ in range(tiling_x * tiling_y)],
        constraints=C2,
        options={"disp": True},
        callback=call_back_optimization,
    )
    # Finalize
    prepare_microstructure(optim.x)
    run_gismo(sensitivities=False, plot=True, refinement=True)


if __name__ == "__main__":
    main()
