"""
Supports single nurbs export for mfem.

Currently hardcoded for 2D-single-patch-splines.
"""

import numpy as np

# single function imports
from splinepy.utils import make_c_contiguous
from splinepy.io.utils import (form_lines,
                                  next_line,
                                  make_meaningful)

# keywords : possible assert value
_mfem_meaningful_keywords = {
    "MFEM NURBS mesh v1.0": "intro",
    "dimension": "para_dim",
    "elements": 1,
    "boundary": 4,
    "edges": 4,
    "vertices": 4,
    "knotvectors": 2,
    "weights": None,
    "FiniteElementSpace": None,
    "FiniteElementCollection": None,
    "VDim": "dim",
    "Ordering": 1,
}


def read_mfem(fname,):
    """
    Reads mfem spline and returns a spline.
    Again, only supports 2D single patch.

    Parameters
    -----------
    fname: str

    Returns
    --------
    nurbs: dict
      dict ready to be used for init. ex) NURBS(**nurbs)
    """
    from copy import deepcopy

    mk = deepcopy(_mfem_meaningful_keywords)
    # they follow a strict order or keywords, so just gather those in order
    # Ordering is a hotkey because control points comes right after
    hotkeys = ["knotvectors", "weights", "Ordering",]
    nurbs_dict = dict(knotvectors=[], weights=[], Ordering=[])
    collect = False
    with open(fname, "r") as f:
        for l in f:
            line = make_meaningful(l)
            if not line:
                continue

            keyword_hit = [k for k in mk.keys() if line.startswith(k)]
            if len(keyword_hit) > 1:
                raise ValueError("double keyword hit!")

            elif (
                len(keyword_hit) == 1
                and keyword_hit[0] in hotkeys
            ):
                collect = True
                current_key = deepcopy(keyword_hit[0])
                continue

            elif (
                len(keyword_hit) == 1
                and keyword_hit[0] not in hotkeys
            ):
                collect = False
                current_key = None

            elif (
                len(keyword_hit) == 0
                and collect
            ):
                pass

            elif(
                len(keyword_hit) == 0
                and not collect
            ):
                continue

            if collect:
                nurbs_dict[current_key].append(line)
            else:
                continue

    # parse nurbs_dict
    # kvs
    knot_vectors = nurbs_dict["knotvectors"][1:]
    knot_vectors[0] = eval("[" + knot_vectors[0].replace(" ", ",") + "]")
    knot_vectors[1] = eval("[" + knot_vectors[1].replace(" ", ",") + "]")
    # pop some values
    # ds
    degrees = [knot_vectors[0].pop(0), knot_vectors[1].pop(0)]
    ncps = knot_vectors[0].pop(0) * knot_vectors[1].pop(0)
    # ws
    weights = nurbs_dict["weights"]
    weights = [eval(w) for w in weights] # hopefully not too slow
    # cps
    control_points = nurbs_dict["Ordering"]
    control_points = [
        eval(f"[{cp.replace(' ', ',')}]") for cp in control_points
    ]

    # double check
    # maybe separate them
    if (
        ncps != len(control_points)
        or ncps != len(weights)
        or int(nurbs_dict["knotvectors"][0]) != 2
    ):
        raise ValueError("Inconsistent spline info in " + fname)

    _, reorder = mfem_index_mapping(
        len(knot_vectors),
        degrees,
        knot_vectors,
        flat_list=True,
    )

    return dict(
        degrees=degrees,
        knot_vectors=knot_vectors,
        control_points=make_c_contiguous(control_points)[reorder],
        weights=make_c_contiguous(weights)[reorder],
    )


def read_solution(fname, reference_nurbs,):
    """
    Given solution and reference_nurbs, returns solution-nurbs-dict.

    Parameters
    -----------
    fname: str
    reference_nurbs: NURBS

    Returns
    --------
    solution_nurbs: dict
    """
    from copy import deepcopy

    with open(fname, "r") as f:
        hotkey = "Ordering"
        line = next_line(f)
        solution = []
        collect = False
        while line is not None:
            if collect: # check this first. should be true most time
                solution.append(eval(f"[{line.replace(' ', ',')}]"))
            elif hotkey in line:
                collect = True

            line = next_line(f)

    if len(reference_nurbs.control_points) != len(solution):
        raise ValueError("Solution length does not match reference nurbs")

    _, reorder = mfem_index_mapping(
        reference_nurbs.para_dim,
        reference_nurbs.degrees,
        reference_nurbs.knot_vectors,
        flat_list=True,
    )

    return dict(
        degrees=deepcopy(reference_nurbs.degrees),
        knot_vectors=deepcopy(reference_nurbs.knot_vectors),
        control_points=make_c_contiguous(solution)[reorder],
        weights=deepcopy(reference_nurbs.weights),
    )


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

        # Quick check
        assert len(flat_groups) == len(np.unique(flat_groups)),\
            "Something went wrong during reorganizing indices for MFEM."

        return (
            flat_groups if flat_list else groups, # to_mfem
            np.argsort(flat_groups) # inverse
        )
        

    else:
        raise NotImplementedError
    


def write_mfem(fname, nurbs, precision=10):
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
    if nurbs.whatami.startswith("NURBS"):
        pass
    elif nurbs.whatami.startswith("BSpline"):
        nurbs = nurbs.nurbs # if bspline, turn it into nurbs
    else:
        raise TypeError("Sorry, invalid spline object.")

    intro_sec = form_lines(
        "MFEM NURBS mesh v1.0",
        "",
        "#",
        "# Generated with splinepy",
        "#",
        "# MFEM Geometry Types (see mesh/geom.hpp)",
        "#",
        "# SEGMENT = 1",
        "# SQUARE  = 3",
        "# CUBE    = 5",
        "#",
        "",
    )

    dimension_sec = form_lines(
        "dimension",
        str(nurbs.para_dim),
        "",
    )

    if nurbs.para_dim == 2:
        elements_sec = form_lines(
            "elements",
            "1",
            "1 3 0 1 2 3",
            "",
        )

        boundary_sec = form_lines(
            "boundary",
            "4",
            "1 1 0 1",
            "2 1 2 3",
            "3 1 3 0",
            "4 1 1 2",
            "",
        )

        edges_sec = form_lines(
            "edges",
            "4",
            "0 0 1",
            "0 3 2",
            "1 0 3",
            "1 1 2",
            "",
        )

        vertices_sec = form_lines(
            "vertices",
            "4",
            "",
        )

        # I am not sure if mixed order is allowed, but incase not, let's
        # match orders
        max_degree = max(nurbs.degrees)
        for i, d in enumerate(nurbs.degrees):
            d_diff = int(max_degree - d)
            if d_diff > 0:
                for _ in range(d_diff):
                    nurbs.elevate_degree(i)

        cnr = nurbs.control_mesh_resolutions

        # double-check
        if not (nurbs.degrees == nurbs.degrees[0]).all():
            raise RuntimeError(
                "Something went wrong trying to match degrees of nurbs "
                + "before export."
            )

        # This is reusable
        def kv_sec(spline):
            kvs = form_lines(
                "knotvectors",
                str(len(spline.knot_vectors)),
            )
            kvs2 = ""
            for i in range(spline.para_dim):
                kvs2 += form_lines(
                    str(spline.degrees[i])
                    + " "
                    + str(cnr[i])
                    + " "
                    + str(spline.knot_vectors[i])[1:-1].replace(",", "")
                )

            kvs = kvs + kvs2

            kvs += "\n" # missing empty line

            return kvs

        knotvectors_sec = kv_sec(nurbs)

    else:
        raise NotImplementedError

    # disregard inverse
    reorder_ids, _ = mfem_index_mapping(
        para_dim=nurbs.para_dim,
        degrees=nurbs.degrees,
        knot_vectors=nurbs.knot_vectors,
        flat_list=True,
    )

    with np.printoptions(
        formatter=dict(float_kind=lambda x: f"{x:.{precision}f}")
    ):
        # weights - string operation
        weights_sec = str(nurbs.weights.flatten()[reorder_ids].tolist())
        weights_sec = weights_sec[1:-1] # remove []
        weights_sec = weights_sec.replace("\n", "") # remove \n
        weights_sec = weights_sec.replace(",", "") # remove ,
        weights_sec = weights_sec.replace(" ", "\n") # add \n
        weights_sec = "weights\n" + weights_sec # add title
        weights_sec += "\n\n" # empty line 

        # cps
        cps_sec = str(nurbs.control_points[reorder_ids].tolist())
        cps_sec = cps_sec.replace("[", "") # remove [
        cps_sec = cps_sec.replace(",", "") # remove ,
        cps_sec = cps_sec.replace("\n", "") # remove \n
        cps_sec = cps_sec.replace("]", "\n") # replace ] with \n

    fe_space_sec = form_lines(
        "FiniteElementSpace",
        "FiniteElementCollection: NURBS" + str(nurbs.degrees[0]),
        "VDim: " + str(nurbs.dim),
        "Ordering: 1",
        "",
    )

    with open(fname, "w") as f:
        f.write(intro_sec)
        f.write(dimension_sec)
        f.write(elements_sec)
        f.write(boundary_sec)
        f.write(edges_sec)
        f.write(vertices_sec)
        f.write(knotvectors_sec)
        f.write(weights_sec)
        f.write(fe_space_sec)
        f.write(cps_sec)
