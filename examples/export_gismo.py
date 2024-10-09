"""
Export to g+smo-compatible xml-file.

Recreates poisson2d_bvp.xml from gismo repo
(see https://github.com/gismo/gismo/blob/stable/filedata/pde/poisson2d_bvp.xml)
"""

import splinepy as spp
from splinepy.io import gismo

EPS = 1e-8


def create_geometry():
    """Recreate multipatch geometry from g+smo repository

    Returns
    ---------
    multipatch: splinepy Multipatch
        Arc geometry consisting of two patches
    """
    # Create example's inner and outer arcs as NURBS patches with predefined values
    arc_inner = spp.NURBS(
        degrees=[1, 2],
        knot_vectors=[[0.0, 0.0, 1.0, 1.0], [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]],
        control_points=[
            [1.0, 0.0],
            [2.0, 0.0],
            [1.0, 1.0],
            [2.0, 2.0],
            [0.0, 1.0],
            [0.0, 2.0],
        ],
        weights=[1.0, 1.0, 0.707106781186548, 0.707106781186548, 1.0, 1.0],
    )

    arc_outer = arc_inner.copy()
    arc_outer.control_points = [
        [2.0, 0.0],
        [3.0, 0.0],
        [2.0, 2.0],
        [3.0, 3.0],
        [0.0, 2.0],
        [0.0, 3.0],
    ]

    # Initialize Multipatch object
    multipatch = spp.Multipatch(splines=[arc_inner, arc_outer])
    multipatch.determine_interfaces()

    # Identify the patches' boundaries which correspond to the Dirichlet BCs
    # Points with y=0 or (x=0 and y=2) are identified as Dirichlet boundaries
    def dirichlet_identifier(points):
        return (points[:, 1] < EPS) | (
            (points[:, 0] < EPS) & (points[:, 1] > 2 - EPS)
        )

    # Set Dirichlet boundaries to "BID2", Neumann BCs correspond to "BID1"
    multipatch.boundary_from_function(dirichlet_identifier, boundary_id=2)

    return multipatch


def generate_additional_xml_blocks():
    """
    Create blocks for rhs, manufactured solution, boundary conditions and assembly
    options.

    Returns
    -----------
    additional_blocks: list<dict>
        List of dictionaries used for splinepy's gismo export function
    """

    additional_blocks = gismo.AdditionalBlocks()

    # Create function block for rhs
    additional_blocks.add_function(
        dim=2,
        block_id=1,
        function_string="2*pi^2*sin(pi*x)*sin(pi*y)",
        comment="Right-hand side function",
    )

    # Create function block for manufactured solution
    additional_blocks.add_function(
        dim=2,
        block_id=3,
        function_string="sin(pi*x) * sin(pi*y)",
        comment="The manufactured exact solution (for reference)",
    )

    # Create block for boundary conditions
    additional_blocks.add_boundary_conditions(
        block_id=2,
        dim=2,
        function_list=[
            "sin(pi*x) * sin(pi*y)",
            ("pi*cos(pi*x) * sin(pi*y)", "pi*sin(pi*x) * cos(pi*y)"),
            "0",
        ],
        bc_list=[("BID2", "Dirichlet", 0), ("BID1", "Neumann", 1)],
        cv_list=[
            ("0", "0", "1", "0"),
            ("0", "0", "2", "sin(x)"),
        ],  # unknown, patch, corner, function
        unknown_id=0,
        multipatch_id=0,
        comment="The boundary conditions (multipatch=number of patches)",
    )

    # Create dictionary for assembly options
    additional_blocks.add_assembly_options(
        block_id=4, comment="Assembler options"
    )

    return additional_blocks.to_list()


if __name__ == "__main__":
    multipatch = create_geometry()

    additional_blocks = generate_additional_xml_blocks()

    # Visualize geometry and BCs
    boundary_names = ["Neumann boundary", "Dirichlet boundary"]
    # spp.show(
    #     ["Multipatch", multipatch],
    #     *[
    #         [
    #             f"BID{i+1}: {boundary_names[i]}",
    #             multipatch.boundary_multipatch(i + 1),
    #         ]
    #         for i in range(len(multipatch.boundaries))
    #     ],
    #     control_points=False,
    # )

    # Export multipatch geometry and additional options/functions to an XML file
    gismo.export(
        fname="poisson2d_bvp_recreation.xml",
        multipatch=multipatch,
        indent=True,
        additional_blocks=additional_blocks,
    )
