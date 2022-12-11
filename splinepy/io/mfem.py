"""Supports single nurbs export for mfem.

Currently hardcoded for 2D-single-patch-splines.
"""

import numpy as np

# single function imports
from splinepy.io.ioutils import (form_lines, next_line, make_meaningful)
from splinepy.splinepy_core import retrieve_mfem_information

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


def load(fname):
    """Reads mfem spline and returns a spline.

    Again, only supports 2D single patch.

    Parameters
    -----------
    fname: str

    Returns
    --------
    nurbs: dict
      dict ready to be used for init.
    """
    from copy import deepcopy

    mk = deepcopy(_mfem_meaningful_keywords)
    # they follow a strict order or keywords, so just gather those in order
    # Ordering is a hotkey because control points comes right after
    hotkeys = [
            "knotvectors",
            "weights",
            "Ordering",
    ]
    nurbs_dict = dict(knotvectors=[], weights=[], Ordering=[])
    collect = False
    with open(fname, "r") as f:
        for single_line in f:
            line = make_meaningful(single_line)
            if not line:
                continue

            keyword_hit = [k for k in mk.keys() if line.startswith(k)]
            if len(keyword_hit) > 1:
                raise ValueError("double keyword hit!")

            elif (len(keyword_hit) == 1 and keyword_hit[0] in hotkeys):
                collect = True
                current_key = deepcopy(keyword_hit[0])
                continue

            elif (len(keyword_hit) == 1 and keyword_hit[0] not in hotkeys):
                collect = False
                current_key = None

            elif (len(keyword_hit) == 0 and collect):
                pass

            elif (len(keyword_hit) == 0 and not collect):
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
    weights = [eval(w) for w in weights]  # hopefully not too slow
    # cps
    control_points = nurbs_dict["Ordering"]
    control_points = [
            eval(f"[{cp.replace(' ', ',')}]") for cp in control_points
    ]

    # double check
    # maybe separate them
    if (
            ncps != len(control_points) or ncps != len(weights)
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
            control_points=np.ascontiguousarray(control_points)[reorder],
            weights=np.ascontiguousarray(weights)[reorder],
    )


def read_solution(
        fname,
        reference_nurbs,
):
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
            if collect:  # check this first. should be true most time
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
            control_points=np.ascontiguousarray(solution)[reorder],
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

    .. code-block::

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

    def flatten(list_):
        """unrolls any nested list"""
        if len(list_) == 0:
            return list_
        if isinstance(list_[0], list):
            return flatten(list_[0]) + flatten(list_[1:])
        return list_[:1] + flatten(list_[1:])

    if int(para_dim) == 2:
        cps_per_dim = [
                len(knot_vectors[i]) - degrees[i] - 1 for i in range(para_dim)
        ]  # degrees = order - 1
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
                cp_ids[group0[0] + 1:group0[1]],
                cp_ids[group0[3] + 1:group0[2]][::-1],
        ]

        # group2 - list of list
        group2 = [
                cp_ids[cps_per_dim[0]::cps_per_dim[0]][:-1],
                cp_ids[group0[1] + cps_per_dim[0]::cps_per_dim[0]][:-1],
        ]

        # group3 - list of list
        group3 = [cp_ids[i + 1:j] for i, j in zip(group2[0], group2[1])]

        groups = [group0, group1, group2, group3]
        flat_groups = flatten(groups)

        # Quick check
        assert len(flat_groups) == len(np.unique(flat_groups)),\
            "Something went wrong during reorganizing indices for MFEM."

        return (
                flat_groups if flat_list else groups,  # to_mfem
                np.argsort(flat_groups)  # inverse
        )

    else:
        raise NotImplementedError


def export_cartesian(
        fname,
        spline_list,
        tolerance=None,
        boundary_functions=None,
        connectivity_only=False
):
    """
    Export list of bezier splines in mfem export.

    Parameters
    ----------
    fname: string
      file name
    bezier_list: list
      list of bezier spline objects
    tolerance : float
      tolerance to collapse two neighboring points
    boundary_functions : list(Callable)
      Functions that define boundaries
    connectivity_only : bool
      For use in other function, only retrieve connectivity

    Returns
    -------
    None
    """
    from splinepy import Spline

    # Check first spline
    if not type(spline_list) == list:
        raise ValueError("export_cartesian expects list for export.")
    if not issubclass(type(spline_list[0]), Spline):
        raise ValueError("Unsupported type in list")

    # Set Tolerance
    if tolerance is None:
        tolerance = 1e-5
    # Set boundary list
    if boundary_functions is None:
        boundary_functions = []

    # Set Problem dimensions
    para_dim = spline_list[0].para_dim
    dim = spline_list[0].dim
    if not ((para_dim == dim) and (dim == 3 or dim == 2)):
        raise ValueError("Only 2D2D or 3D3D splines are supported")

    # For use in loop
    def _sanity_check(spline):
        if not issubclass(type(spline), Spline):
            raise ValueError("Unsupported type in list")
        if not (spline.dim == dim and spline.para_dim == para_dim):
            raise ValueError("Dimension Mismatch between splines")

    # Auxiliary function to identify corner vertices of the underlying spline
    #  representation. As this function requires the numeration of its
    # underlying hypercube Element, the following part needs to be different
    # for the 2D and the 3D case
    if para_dim == 2:
        geometry_type = 3
        boundary_type = 1
        n_vertex_per_element = 4
        n_vertex_per_boundary = 2

        def _corner_vertex_ids(spline):
            cmr = spline.control_mesh_resolutions
            return np.array(
                    [
                            0,
                            cmr[0] - 1,
                            cmr[0] * cmr[1] - 1,
                            cmr[0] * (cmr[1] - 1),
                    ],
                    dtype="int32"
            )
    elif para_dim == 3:
        geometry_type = 5
        boundary_type = 3
        n_vertex_per_element = 8
        n_vertex_per_boundary = 4

        def _corner_vertex_ids(spline):
            cmr = spline.control_mesh_resolutions
            return np.array(
                    [
                            0,
                            cmr[0] - 1,
                            cmr[0] * cmr[1] - 1,
                            cmr[0] * (cmr[1] - 1),
                            (cmr[2] - 1) * cmr[1] * cmr[0],
                            (cmr[2] - 1) * cmr[1] * cmr[0] + cmr[0] - 1,
                            (cmr[2] - 1) * cmr[1] * cmr[0] + cmr[0] * cmr[1]
                            - 1,
                            (cmr[2] - 1) * cmr[1] * cmr[0] + cmr[0]
                            * (cmr[1] - 1),
                    ],
                    dtype="int32"
            )

    # Create a list of all corner vertices ordered by spline patch number
    corner_vertices = np.empty((0, dim))
    for spline in spline_list:
        _sanity_check(spline)
        corner_vertex_ids = _corner_vertex_ids(spline)
        corner_vertices = np.vstack(
                (corner_vertices, spline.control_points[corner_vertex_ids, :])
        )

    # Retrieve information using bezman
    (connectivity, vertex_ids, edges, boundaries,
     success) = retrieve_mfem_information(corner_vertices, tolerance)

    if connectivity_only:
        return connectivity

    boundary_ids = np.ones(boundaries.shape[0], dtype=int)
    if not len(boundary_functions) == 0:
        # Retrieve all boundary points (to minimize computations)
        boundary_vertex_ids = np.unique(boundaries)
        # Boundary_vertex_ids refer to the new enumeration
        _, vertex_ids_unique = np.unique(vertex_ids, return_index=True)
        boundary_corner_verts = corner_vertices[
                vertex_ids_unique[boundary_vertex_ids], :]
        # Loop over all boundary functions to identify boundaries
        for i_bound, b_func in enumerate(boundary_functions):
            try:
                is_on_boundary = b_func(boundary_corner_verts)
            except BaseException:
                raise ValueError("Boundary Function Layout incompatible")
            set_boundary = np.sum(
                    np.isin(boundaries, boundary_vertex_ids[is_on_boundary]),
                    axis=1
            ) == n_vertex_per_boundary
            boundary_ids[set_boundary] = i_bound + 2

    # Check if algorithm was successfull
    if not success:
        # todo : only connectivity stores data, see issue #33
        # Use <SplineBase>.permute to potentially save mesh export
        raise NotImplementedError(
                "MFEM export not implemented for unstructured meshes"
        )

    # Write all gathered information into a file
    with open(fname, "w") as f:
        # Header
        f.write(f"MFEM NURBS mesh v1.0\n\ndimension\n{dim}\n\n")

        # Elements
        n_elements = len(spline_list)
        f.write(f"elements\n{n_elements}\n")
        f.write(
                '\n'.join(
                        f"1 {geometry_type} "
                        + ' '.join(str(id) for id in row) for row in
                        vertex_ids.reshape(-1, n_vertex_per_element).tolist()
                )
        )

        # Boundaries
        n_boundaries = boundaries.shape[0]
        f.write(f"\n\nboundary\n{n_boundaries}\n")
        # Here currently all boudaries are set to 1
        f.write(
                '\n'.join(
                        f"{boundary_id} {boundary_type} "
                        + ' '.join(str(id)
                                   for id in row)
                        for row, boundary_id in zip(
                                boundaries.reshape(-1, n_vertex_per_boundary).
                                tolist(), boundary_ids.tolist()
                        )
                )
        )

        # Edges
        n_edges = edges.shape[0]
        f.write(f"\n\nedges\n{n_edges}\n")
        # Here currently all boudaries are set to 1
        f.write(
                '\n'.join(
                        ' '.join(str(id)
                                 for id in row)
                        for row in edges.reshape(-1, 3).tolist()
                )
        )

        # Write Number Of vertices
        f.write(f"\n\nvertices\n{np.max(vertex_ids)+1}\n\n")

        # Export Splines
        f.write("patches\n\n")
        for spline in spline_list:
            f.write(f"knotvectors\n{para_dim}\n")
            cmr = spline.control_mesh_resolutions
            for i_para_dim in range(para_dim):
                f.write(f"{spline.degrees[i_para_dim]} ")
                f.write(f"{cmr[i_para_dim]} ")
                if "knot_vectors" not in spline.required_properties:
                    f.write("0.0 " * int(spline.degrees[i_para_dim] + 1))
                    f.write(
                            "1.0 " * int(spline.degrees[i_para_dim] + 1) + "\n"
                    )
                else:
                    f.write(
                            ' '.
                            join(str(kvi)
                                 for kvi in spline.kvs[i_para_dim]) + "\n"
                    )
            f.write(f"dimension\n{dim}\n")
            f.write("controlpoints_cartesian\n")
            if "weights" not in spline.required_properties:
                f.write(
                        '\n'.join(
                                (' '.join(str(x_i)
                                          for x_i in row) + " 1.0")
                                for row in spline.control_points.tolist()
                        ) + "\n\n"
                )
            else:
                f.write(
                        '\n'.join(
                                (
                                        ' '.join(str(x_i) for x_i in coords)
                                        + " " + str(weight[0])
                                ) for (coords, weight) in zip(
                                        spline.control_points.tolist(),
                                        spline.weights.tolist()
                                )
                        ) + "\n\n"
                )

        f.close()


def export(fname, nurbs, precision=10):
    """Exports current nurbs in `mfem` format.

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
        nurbs = nurbs.nurbs  # if bspline, turn it into nurbs
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
                        str(spline.degrees[i]) + " " + str(cnr[i]) + " "
                        + str(spline.knot_vectors[i])[1:-1].replace(",", "")
                )

            kvs = kvs + kvs2

            kvs += "\n"  # missing empty line

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
        weights_sec = weights_sec[1:-1]  # remove []
        weights_sec = weights_sec.replace("\n", "")  # remove \n
        weights_sec = weights_sec.replace(",", "")  # remove ,
        weights_sec = weights_sec.replace(" ", "\n")  # add \n
        weights_sec = "weights\n" + weights_sec  # add title
        weights_sec += "\n\n"  # empty line

        # cps
        cps_sec = str(nurbs.control_points[reorder_ids].tolist())
        cps_sec = cps_sec.replace("[", "")  # remove [
        cps_sec = cps_sec.replace(",", "")  # remove ,
        cps_sec = cps_sec.replace("\n", "")  # remove \n
        cps_sec = cps_sec.replace("]", "\n")  # replace ] with \n

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
