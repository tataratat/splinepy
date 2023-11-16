"""Supports single nurbs export for mfem.

Currently hardcoded for 2D-single-patch-splines.
"""

import gustaf as _gus
import numpy as _np

from splinepy.io.ioutils import dict_to_spline as _dict_to_spline
from splinepy.io.ioutils import form_lines as _form_lines
from splinepy.io.ioutils import make_meaningful as _make_meaningful
from splinepy.io.ioutils import next_line as _next_line

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
    from copy import deepcopy as _deepcopy

    mk = _deepcopy(_mfem_meaningful_keywords)
    # they follow a strict order or keywords, so just gather those in order
    # Ordering is a hotkey because control points comes right after
    hotkeys = [
        "knotvectors",
        "weights",
        "Ordering",
    ]
    nurbs_dict = {"knotvectors": [], "weights": [], "Ordering": []}
    collect = False
    with open(fname) as f:
        for single_line in f:
            line = _make_meaningful(single_line)
            if not line:
                continue

            keyword_hit = [k for k in mk if line.startswith(k)]
            if len(keyword_hit) > 1:
                raise ValueError("double keyword hit!")

            elif len(keyword_hit) == 1 and keyword_hit[0] in hotkeys:
                collect = True
                current_key = _deepcopy(keyword_hit[0])
                continue

            elif len(keyword_hit) == 1 and keyword_hit[0] not in hotkeys:
                collect = False
                current_key = None

            elif len(keyword_hit) == 0 and collect:
                pass

            elif len(keyword_hit) == 0 and not collect:
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
        ncps != len(control_points)
        or ncps != len(weights)
        or int(nurbs_dict["knotvectors"][0]) != 2
    ):
        raise ValueError("Inconsistent spline info in " + fname)

    mfem_dict_spline = {
        "degrees": degrees,
        "knot_vectors": knot_vectors,
        "control_points": _np.ascontiguousarray(control_points),
        "weights": _np.ascontiguousarray(weights),
    }

    spline = _dict_to_spline([mfem_dict_spline])[0]

    _, reorder = dof_mapping(spline)

    spline.cps[:] = spline.cps[reorder]
    spline.ws[:] = spline.ws[reorder]

    return spline


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

    with open(fname) as f:
        hotkey = "Ordering"
        line = _next_line(f)
        solution = []
        collect = False
        while line is not None:
            if collect:  # check this first. should be true most time
                solution.append(eval(f"[{line.replace(' ', ',')}]"))
            elif hotkey in line:
                collect = True

            line = _next_line(f)

    if len(reference_nurbs.control_points) != len(solution):
        raise ValueError("Solution length does not match reference nurbs")

    mfem_dict_spline = reference_nurbs.todict()

    spline = _dict_to_spline([mfem_dict_spline])[0]

    _, reorder = dof_mapping(spline)

    spline.cps = _np.ascontiguousarray(solution)[reorder]
    spline.ws[:] = spline.ws[reorder]

    return spline


def _mfem_dof_map_2d(spline):
    """
    Dof mapping for 2D splines
    """
    mi = spline.multi_index

    to_m = []

    # corner vertices
    x_ids = [0, -1, -1, 0]
    y_ids = [0, 0, -1, -1]
    vertices = mi[x_ids, y_ids].tolist()

    # edge vertices
    edges = []
    edges.extend(mi[1:-1:, 0].tolist())  # bottom
    edges.extend(mi[1:-1, -1].tolist()[::-1])  # top
    edges.extend(mi[0, 1:-1].tolist())  # left
    edges.extend(mi[-1, 1:-1].tolist())  # right

    # internal
    internal = mi[1:-1, 1:-1].tolist()

    to_m.extend(vertices)
    to_m.extend(edges)
    to_m.extend(internal)

    return to_m, _np.argsort(to_m).tolist()


def _mfem_dof_map_3d(spline):
    """
    Dof mapping for 3D splines
    """
    mi = spline.multi_index

    to_m = []

    # corner vertices
    x_ids = [0, -1, -1, 0, 0, -1, -1, 0]
    y_ids = [0, 0, -1, -1, 0, 0, -1, -1]
    z_ids = [0, 0, 0, 0, -1, -1, -1, -1]

    vertices = mi[x_ids, y_ids, z_ids].tolist()

    # edges
    edges = []
    edges.extend(mi[1:-1, 0, 0].tolist())  # 1
    edges.extend(mi[1:-1, -1, 0].tolist()[::-1])  # 2
    edges.extend(mi[1:-1, 0, -1].tolist())  # 3
    edges.extend(mi[1:-1, -1, -1].tolist()[::-1])  # 4

    edges.extend(mi[0, 1:-1, 0].tolist())  # 5
    edges.extend(mi[-1, 1:-1, 0].tolist())  # 6
    edges.extend(mi[0, 1:-1, -1].tolist())  # 7
    edges.extend(mi[-1, 1:-1, -1].tolist())  # 8

    edges.extend(mi[0, 0, 1:-1].tolist())  # 9
    edges.extend(mi[-1, 0, 1:-1].tolist())  # 10
    edges.extend(mi[-1, -1, 1:-1].tolist())  # 11
    edges.extend(mi[0, -1, 1:-1].tolist())  # 12

    # faces - this one is extra funky
    faces = []
    faces.extend(mi[1:-1, -2:0:-1, 0].tolist())  # 1
    faces.extend(mi[1:-1, 0, 1:-1].tolist())  # 2
    faces.extend(mi[-1, 1:-1, 1:-1].tolist())  # 3
    faces.extend(mi[-2:0:-1, -1, 1:-1].tolist())  # 4
    faces.extend(mi[0, -2:0:-1, 1:-1].tolist())  # 5
    faces.extend(mi[1:-1, 1:-1, -1].tolist())  # 6

    # internal - thank you for being so friendly
    internal = mi[1:-1, 1:-1, 1:-1].tolist()

    # merge them
    to_m.extend(vertices)
    to_m.extend(edges)
    to_m.extend(faces)
    to_m.extend(internal)

    return to_m, _np.argsort(to_m).tolist()


def dof_mapping(spline):
    """
    MFEM stores vertices in a special order. This function gives dof index
    mapping between splinepy and mfem

    Parameters
    ----------
    spline: Spline

    Returns
    -------
    to_mfem: list
    from_mfem: array
    """
    if spline.para_dim == 2:
        return _mfem_dof_map_2d(spline)
    elif spline.para_dim == 3:
        return _mfem_dof_map_3d(spline)
    else:
        raise ValueError("supports para_dim of 2 and 3")


def export_cartesian(
    fname,
    spline_list,
    tolerance=None,
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

    Returns
    -------
    None
    """
    from splinepy import Multipatch as _Multipatch

    # Check first spline
    if not isinstance(spline_list, (list, _Multipatch)):
        raise ValueError("export_cartesian expects list for export.")
    if isinstance(spline_list, list):
        spline_list = _Multipatch(splines=spline_list)

    # Set Tolerance
    if tolerance is None:
        tolerance = 1e-5

    # Set Problem dimensions
    para_dim = spline_list.para_dim
    dim = spline_list.dim
    if not ((para_dim == dim) and (dim in {2, 3})):
        raise ValueError("Only 2D2D or 3D3D splines are supported")

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
            return _np.array(
                [
                    0,
                    cmr[0] - 1,
                    cmr[0] * cmr[1] - 1,
                    cmr[0] * (cmr[1] - 1),
                ],
                dtype="int32",
            )

        # Face enumeration is different from splinepy's
        sub_element_vertices = _np.array([[0, 1], [1, 2], [3, 2], [0, 3]])
        # splinepy's face 3 corresponds to mfem's face 0
        splinepy_face_id_2_mfem_face_id = _np.array([3, 1, 0, 2])

    elif para_dim == 3:
        geometry_type = 5
        boundary_type = 3
        n_vertex_per_element = 8
        n_vertex_per_boundary = 4

        def _corner_vertex_ids(spline):
            cmr = spline.control_mesh_resolutions
            return _np.array(
                [
                    0,
                    cmr[0] - 1,
                    cmr[0] * cmr[1] - 1,
                    cmr[0] * (cmr[1] - 1),
                    (cmr[2] - 1) * cmr[1] * cmr[0],
                    (cmr[2] - 1) * cmr[1] * cmr[0] + cmr[0] - 1,
                    (cmr[2] - 1) * cmr[1] * cmr[0] + cmr[0] * cmr[1] - 1,
                    (cmr[2] - 1) * cmr[1] * cmr[0] + cmr[0] * (cmr[1] - 1),
                ],
                dtype="int32",
            )

        sub_element_vertices = _np.array(
            [
                [0, 1, 2, 3],
                [1, 2, 6, 5],
                [3, 2, 6, 7],
                [0, 3, 7, 4],
                [0, 1, 5, 4],
                [4, 5, 6, 7],
            ]
        )
        splinepy_face_id_2_mfem_face_id = _np.array([3, 1, 4, 2, 0, 5])

    # Create a list of all corner vertices ordered by spline patch number
    corner_vertices = _np.vstack(
        [
            spline.cps[_corner_vertex_ids(spline), :]
            for spline in spline_list.patches
        ]
    )

    # Retrieve information using bezman
    connectivity = spline_list.interfaces

    (_, _, inverse_numeration, _) = _gus.utils.arr.close_rows(
        corner_vertices, tolerance, return_intersection=False
    )
    vertex_ids = inverse_numeration.reshape(-1, n_vertex_per_element)
    # Get boundaries from interfaces
    boundary_elements, boundary_faces = _np.where(connectivity < 0)
    boundary_ids = -connectivity[boundary_elements, boundary_faces]
    # Convert from splinepy enumeration to mfem enumeration
    boundaries = _np.take_along_axis(
        vertex_ids[boundary_elements, :],
        sub_element_vertices[splinepy_face_id_2_mfem_face_id[boundary_faces]],
        axis=1,
    )

    # Write all gathered information into a file
    with open(fname, "w") as f:
        # Header
        f.write(f"MFEM NURBS mesh v1.0\n\ndimension\n{dim}\n\n")

        # Elements
        n_elements = len(spline_list.patches)
        f.write(f"elements\n{n_elements}\n")
        f.write(
            "\n".join(
                f"1 {geometry_type} " + " ".join(str(id) for id in row)
                for row in vertex_ids.reshape(
                    -1, n_vertex_per_element
                ).tolist()
            )
        )

        # Boundaries
        n_boundaries = boundaries.shape[0]
        f.write(f"\n\nboundary\n{n_boundaries}\n")
        # Here currently all boundaries are set to 1
        f.write(
            "\n".join(
                f"{boundary_id} {boundary_type} "
                + " ".join(str(id) for id in row)
                for row, boundary_id in zip(
                    boundaries.reshape(-1, n_vertex_per_boundary).tolist(),
                    boundary_ids.tolist(),
                )
            )
        )

        # Write number of edges 0, for auto knot2edge
        f.write(f"\n\nedges\n{0}\n")

        # Write Number Of vertices
        f.write(f"\nvertices\n{int(_np.max(vertex_ids)+1)}\n\n")

        # Export Splines
        f.write("patches\n\n")
        for spline in spline_list.patches:
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
                        " ".join(str(kvi) for kvi in spline.kvs[i_para_dim])
                        + "\n"
                    )
            f.write(f"dimension\n{dim}\n")
            f.write("controlpoints_cartesian\n")
            if "weights" not in spline.required_properties:
                f.write(
                    "\n".join(
                        (" ".join(str(x_i) for x_i in row) + " 1.0")
                        for row in spline.control_points.tolist()
                    )
                    + "\n\n"
                )
            else:
                f.write(
                    "\n".join(
                        (
                            " ".join(str(x_i) for x_i in coords)
                            + " "
                            + str(weight[0])
                        )
                        for (coords, weight) in zip(
                            spline.control_points.tolist(),
                            spline.weights.tolist(),
                        )
                    )
                    + "\n\n"
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
    else:
        nurbs = nurbs.nurbs  # turns it into nurbs

    intro_sec = _form_lines(
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

    dimension_sec = _form_lines(
        "dimension",
        str(nurbs.para_dim),
        "",
    )

    if nurbs.para_dim == 2:
        elements_sec = _form_lines(
            "elements",
            "1",
            "1 3 0 1 2 3",
            "",
        )

        boundary_sec = _form_lines(
            "boundary",
            "4",
            "1 1 0 1",
            "2 1 2 3",
            "3 1 3 0",
            "4 1 1 2",
            "",
        )

        edges_sec = _form_lines(
            "edges",
            "4",
            "0 0 1",
            "0 3 2",
            "1 0 3",
            "1 1 2",
            "",
        )

        vertices_sec = _form_lines(
            "vertices",
            "4",
            "",
        )

        # I am not sure if mixed order is allowed, but in case not, let's
        # match orders
        max_degree = max(nurbs.degrees)
        for i, d in enumerate(nurbs.degrees):
            d_diff = int(max_degree - d)
            if d_diff > 0:
                for _ in range(d_diff):
                    nurbs.elevate_degrees(i)

        cnr = nurbs.control_mesh_resolutions

        # double-check
        if not (nurbs.degrees == nurbs.degrees[0]).all():
            raise RuntimeError(
                "Something went wrong trying to match degrees of nurbs "
                + "before export."
            )

        # This is reusable
        def kv_sec(spline):
            kvs = _form_lines(
                "knotvectors",
                str(len(spline.knot_vectors)),
            )
            kvs2 = ""
            for i in range(spline.para_dim):
                kvs2 += _form_lines(
                    str(spline.degrees[i])
                    + " "
                    + str(cnr[i])
                    + " "
                    + str(list(spline.knot_vectors[i]))[1:-1].replace(",", "")
                )

            kvs = kvs + kvs2

            kvs += "\n"  # missing empty line

            return kvs

        knotvectors_sec = kv_sec(nurbs)

    else:
        raise NotImplementedError

    # disregard inverse
    reorder_ids, _ = dof_mapping(nurbs)

    with _np.printoptions(
        formatter={"float_kind": lambda x: f"{x:.{precision}f}"}
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

    fe_space_sec = _form_lines(
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
