"""
Export to g+smo-compatible xml-file. Recreates poisson2d_bvp.xml from gismo repo
"""

import splinepy as spp
import splinepy.io.gismo as gismo

EPS = 1e-8

# Create multipatch geometry
arc_inner = spp.NURBS(
    degrees=[1,2],
    knot_vectors=[
        [0.0, 0.0, 1.0, 1.0],
        [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
    ],
    control_points=[
        [1.0, 0.0],
        [2.0, 0.0],
        [1.0, 1.0],
        [2.0, 2.0],
        [0.0, 1.0],
        [0.0, 2.0]
    ],
    weights=[1.0, 1.0, 0.707106781186548, 0.707106781186548, 1.0, 1.0]
)

arc_outer = arc_inner.copy()
arc_outer.control_points = [
    [2.0, 0.0],
    [3.0, 0.0],
    [2.0, 2.0],
    [3.0, 3.0],
    [0.0, 2.0],
    [0.0, 3.0]
]

multipatch = spp.Multipatch(
    splines=[arc_inner, arc_outer]
)
multipatch.determine_interfaces()

# Identifies the patches' boundaries which correspond to the Dirichlet BCs
def dirichlet_identifier(points):
    return (points[:,1] < EPS) | ((points[:,0] < EPS) & (points[:,1] > 2-EPS))

# Set Dirichlet boundaries to "BID2", Neumann BCs correspond to "BID1"
multipatch.boundary_from_function(dirichlet_identifier, boundary_id=2)

# Create function block for rhs
rhs = gismo.create_function_block(
    dim=2,
    function_id=1,
    function_string="2*pi^2*sin(pi*x)*sin(pi*y)",
    comment="Right-hand side function"
)

# Create function block for manufactured solution
manufactured_solution = gismo.create_function_block(
    dim=2,
    function_id=3,
    function_string="sin(pi*x) * sin(pi*y)",
    comment="The manufactured exact solution (for reference)"
)

# Create block for boundary conditions
boundary_conditions = gismo.create_boundary_conditions_block(
    bc_id=2,
    dim=2,
    function_list=[
        "sin(pi*x) * sin(pi*y)",
        ("pi*cos(pi*x) * sin(pi*y)", "pi*sin(pi*x) * cos(pi*y)"),
        "0"
    ],
    bc_list=[
        ("BID2", "Dirichlet", 0),
        ("BID1", "Neumann", 1)
    ],
    unknown_id=0,
    multipatch_id=0,
    comment="The boundary conditions (multipatch=number of patches)"
)

# Create dictionary for assembly options
assembly_options = gismo.create_assembly_options_block(
    options_id=4,
    comment="Assembler options"
)

# Visualize geometry and BCs
boundary_names = ["Neumann boundary", "Dirichlet boundary"]
spp.show(
    ["Multipatch", multipatch],
    *[[f"BID{i+1}: {boundary_names[i]}", multipatch.boundary_multipatch(i+1)] for i in range(len(multipatch.boundaries))],
    control_points=False
)

# Export to xml-file
gismo.export(
    fname="poisson2d_bvp_recreation.xml",
    multipatch=multipatch,
    indent=True,
    options=[rhs, boundary_conditions, manufactured_solution, assembly_options]
)