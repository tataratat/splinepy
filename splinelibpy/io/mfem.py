"""
Supports single nurbs export for mfem.

Currently hardcoded for 2D-single-patch-splines.
"""

import numpy as np

def form_lines(*args):
    """
    Formulate a string, taking each *args as a line.

    Parameters
    -----------
    *args: *str

    Returns
    --------
    line_separated_str: str
    """
    line_separated_str = ""
    for a in args:
        line_separated_str += a + "\n"

    return line_separated_str

def mfem_index_mapping(
        para_dim,
        degrees,
        knot_vectors,
        flat_list=True,
):
    """
    Returns index to properly reorganize control points/weights.


    2D:
            <---------------------- 1.2
    0.3 *-------------------------------* 0.2
        |                               |
    2.1 |   --------------------> 3.... | 2.2
     ^  |   --------------------> 3.8   |  ^
     |  |   --------------------> 3.7   |  |
     |  |   --------------------> 3.6   |  |
     |  |   --------------------> 3.5   |  |
     |  |   --------------------> 3.4   |  |
     |  |   --------------------> 3.3   |  |
     |  |   --------------------> 3.2   |  |
     |  |   --------------------> 3.1   |  |
        |                               |
    0.0 *-------------------------------* 0.1
           -----------------------> 1.1

    Parameters
    -----------
    para_dim: int
    degrees: (1D) array-like
    knot_vectors: list
    flat_list: bool
      If false, each sections are returned in a separate list

    Returns
    --------
    to_mfem: (n,) np.ndarray
    inverse: (n,) np.ndarray
    """
    def flatten(l):
        """unrolls any nested list"""
        if len(l) == 0:
            return l
        if isinstance(l[0], list):
           return flatten(l[0]) + flatten(l[1:])
        return l[:1] + flatten(l[1:])


    if int(para_dim) == 2:
        cps_per_dim = [
            len(knot_vectors[i]) - degrees[i] - 1 for i in range(para_dim)
        ] # degrees = order - 1
        cp_ids = np.arange(np.prod(cps_per_dim)).tolist()

        # group0 - list of int
        group0 = [
            cp_ids[0],
            cp_ids[cps_per_dim[0] - 1],
            cp_ids[-1],
            cp_ids[-cps_per_dim[0]],
        ]

        # group1 - list of list
        group1 = [
            cp_ids[group0[0]+1:group0[1]],
            cp_ids[group0[3]+1:group0[2]][::-1],
        ]

        # group2 - list of list
        group2 = [
            cp_ids[cps_per_dim[0]::cps_per_dim[0]][:-1],
            cp_ids[group0[1] + cps_per_dim[0]::cps_per_dim[0]][:-1],
        ]

        # group3 - list of list
        group3 = [
            cp_ids[i + 1:j] for i, j in zip(group2[0], group2[1])
        ]

        groups = [group0, group1, group2, group3]
        flat_groups = flatten(groups)
        print(groups)
        print(flat_groups)

        # Quick check
        assert len(flat_groups) == len(np.unique(flat_groups)),\
            "Something went wrong during reorganizing indices for MFEM."

        return (
            flat_groups if flat_list else groups, # to_mfem
            np.argsort(flat_groups) # inverse
        )
        

    else:
        raise NotImplementedError
    


def mfem(nurbs, fname):
    """
    Exports current nurbs in `mfem` format.

    IDs of MFEM Geometry types are:
    - Segment = 1
    - Square = 3
    - CUBS = 5

    Sections required are:
    - dimension  (space dim)
    - elements
    - boundary
    - edges
    - vertices
    - knotvectors
    - weights
    - FiniteElementSpace
    - FiniteElementCollection: <name>
    - VDim
    - Ordering (1?)

    Parameters
    -----------
    nurbs: NURBS
    fname: str

    Returns
    --------
    None
    """
    intro_sec = (
        "MFEM NURBS mesh v1.0 - Generated with splinelibpy\n"
        + "\n"
        + "#\n"
        + "# MFEM Geometry Types (see mesh/geom.hpp)\n"
        + "#\n"
        + "# SEGMENT = 1\n"
        + "# SQUARE  = 3\n"
        + "# CUBE    = 3\n"
        + "#\n"
        + "\n"
    )

    dimension_sec = (
        "dimension\n"
        + str(s)
        + "\n"
        + "\n"
    )

    elements_sec = (
        "elements\n"
        + "1"
    )























"""
general curvlinear form
MFEM mesh v1.0

# Space dimension: 2 or 3
dimension
<dimension>

# Mesh elements, e.g. tetrahedrons (4)
elements
<number of elements>
<element attribute> <geometry type> <vertex index 1> ... <vertex index m>
...

# Mesh faces/edges on the boundary, e.g. triangles (2)
boundary
<number of boundary elements>
<boundary element attribute> <geometry type> <vertex index 1> ... <vertex index m>
...

# Number of vertices (no coordinates)
vertices
<number of vertices>

# Mesh nodes as degrees of freedom of a finite element grid function
nodes
FiniteElementSpace
FiniteElementCollection: <finite element collection>
VDim: <dimension>
Ordering: 0
<x-coordinate degrees of freedom>
...
<y-coordinate degrees of freedom>
...
<z-coordinate degrees of freedom>
...
"""
