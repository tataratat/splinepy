from functools import wraps as _wraps

import numpy as _np
from gustaf import Edges as _Edges
from gustaf import Faces as _Faces
from gustaf import Vertices as _Vertices
from gustaf import Volumes as _Volumes
from gustaf.create.edges import from_data as _from_data
from gustaf.utils import connec as _connec
from gustaf.utils.arr import enforce_len as _enforce_len

from splinepy import settings as _settings
from splinepy.helpme import visualize as _visualize
from splinepy.utils import log as _log
from splinepy.utils.data import cartesian_product as _cartesian_product


def edges(
    spline,
    resolution=100,
    all_knots=True,
):
    """Extract edges (lines) from a given spline, or multipatch object. Only
    entity you can extract without dimension limit.

    Parameters
    -----------
    spline: Spline/ Multipatch
      (`self`-argument if called via extract member of a spline)
    resolution: int or list
      samples per parametric dimension
    all_knots: bool
      Switch to allow all knot-line extraction or just contour (para_dim>0)

    Returns
    --------
    edges: Edges
    """
    from splinepy import Multipatch as _Multipatch

    if spline.para_dim == 1:
        # get a single value
        if isinstance(resolution, (_np.ndarray, tuple, list)):
            resolution = max(resolution)
        # ensure int
        resolution = int(resolution)

        vertices = spline.sample(resolution)
        e = _connec.range_to_edges(
            (0, resolution),
            closed=False,
        )
        if isinstance(spline, _Multipatch):
            e = _np.vstack(
                [e + i * resolution for i in range(len(spline.patches))]
            )
        return _Edges(
            vertices=vertices,
            edges=e,
        )

    else:
        # for para_dim > 1, we support array-like resolutions
        resolution = _enforce_len(resolution, spline.para_dim)

        # recursion for Multipatch splines
        if isinstance(spline, _Multipatch):
            return _Edges.concat(
                [
                    edges(spline=s, resolution=resolution, all_knots=all_knots)
                    for s in spline.patches
                ]
            )

        # Single patch case
        if all_knots:
            relevant_knots = spline.unique_knots
        else:
            relevant_knots = list(spline.parametric_bounds.T)

        temp_edges = []  # edges' is not a valid syntax
        for i in range(spline.para_dim):
            split_knots = relevant_knots.copy()
            split_knots.pop(i)
            n_knot_lines = _np.prod([len(s) for s in split_knots], dtype=int)
            # Create query points
            # Sampling along a line can be done using cartesian product,
            # however, to preserve the correct order we need to permute the
            # columns, example:
            #
            # [y0, x, z]    [x, y0, z]
            # [y1, x, z] -> [x, y1, z]
            # [y2, x, z]    [x, y2, z]
            reorder_mask = (
                [*range(1, i + 1)] + [0] + [*range(1 + i, spline.para_dim)]
            )
            extract_knot_queries = _cartesian_product(
                [_np.linspace(*spline.parametric_bounds[:, i], resolution[i])]
                + split_knots,
                reverse=True,
            )[:, reorder_mask]

            # Retrieve connectivity to be repeated
            single_line_connectivity = _connec.range_to_edges(
                (0, resolution[i]),
                closed=False,
            )

            temp_edges.append(
                _Edges(
                    vertices=spline.evaluate(extract_knot_queries),
                    edges=_np.vstack(
                        [
                            single_line_connectivity + resolution[i] * j
                            for j in range(n_knot_lines)
                        ]
                    ),
                )
            )

        return _Edges.concat(temp_edges)


def _uniform_3d_faces(spline, res):
    # Extract boundaries and sample from boundaries
    if isinstance(spline, _settings.NAME_TO_TYPE["Multipatch"]):
        boundaries = spline.boundary_multipatch()
    else:
        boundaries = _settings.NAME_TO_TYPE["Multipatch"](
            splines=spline.extract.boundaries()
        )

    n_faces = (res - 1) ** 2
    vertices = res**2
    f_loc = _connec.make_quad_faces(_enforce_len(res, 2))
    face_connectivity = _np.empty((n_faces * len(boundaries.patches), 4))

    # Create Connectivity for Multipatches
    for i in range(len(boundaries.patches)):
        face_connectivity[i * n_faces : (i + 1) * n_faces] = f_loc + (
            i * vertices
        )

    # make faces and merge vertices before returning
    return _Faces(vertices=boundaries.sample(res), faces=face_connectivity)


def faces(
    spline,
    resolution,
    watertight=True,
):
    """Extract faces from spline. Valid iff para_dim is one of the following:
    {2, 3}. In case of {3}, it will return only surfaces. If internal faces are
    desired, used `spline.extract.volumes().faces()`. Note that dimension
    higher than 3 is not showable.

    Parameters
    -----------
    spline: Spline / Multipatch
      (`self`-argument if called via extract member of a spline)
    resolution: int or list
      samples per parametric dimension.
    watertight: bool
      Default is True. Only related to para_dim = 3 splines. If False,
      overlapping vertices at boundary edges won't be merged.

    Returns
    --------
    faces: faces
    """
    from splinepy import Multipatch as _Multipatch

    if spline.para_dim == 2:
        if isinstance(spline, _Multipatch):
            n_faces = (resolution - 1) ** spline.para_dim
            vertices = resolution**spline.para_dim
            f_loc = _connec.make_quad_faces(
                _enforce_len(resolution, spline.para_dim)
            )
            face_connectivity = _np.empty((n_faces * len(spline.patches), 4))
            # Create Connectivity for Multipatches
            for i in range(len(spline.patches)):
                face_connectivity[i * n_faces : (i + 1) * n_faces] = f_loc + (
                    i * vertices
                )
        else:
            face_connectivity = _connec.make_quad_faces(
                _enforce_len(resolution, spline.para_dim)
            )
        faces = _Faces(
            vertices=spline.sample(resolution),
            faces=face_connectivity,
        )

    elif spline.para_dim == 3:
        full_res = _enforce_len(resolution, spline.para_dim)
        if (full_res[0] == full_res[1] == full_res[2]) or isinstance(
            spline, _Multipatch
        ):
            faces = _uniform_3d_faces(spline, int(max(full_res)))
        else:
            # single patch should support varying resolution sample
            b_faces = []
            for i, bs in enumerate(spline.extract.boundaries()):
                res = full_res.tolist()
                pop_dim = i // 2
                res.pop(pop_dim)

                b_faces.append(
                    _Faces(bs.sample(res), _connec.make_quad_faces(res))
                )

            faces = _Faces.concat(b_faces)

    else:
        raise ValueError("Invalid spline to make faces.")

    if watertight:
        faces.merge_vertices()

    return faces


def volumes(spline, resolution, watertight=False):
    """Extract volumes from spline. Valid iff spline.para_dim == 3.

    Parameters
    -----------
    spline: Spline / Multipatch
      (`self`-argument if called via extract member of a spline)
    resolution: int or list
      Sample resolution
    watertight : bool
      Return watertight mesh

    Returns
    --------
    volumes: Volumes
    """
    from splinepy import Multipatch as _Multipatch

    if spline.para_dim != 3 or spline.dim != 3:
        raise ValueError(
            "Volume extraction from a spline is only valid for "
            "para_dim: 3 dim: 3 splines."
        )

    if isinstance(spline, _Multipatch):
        if not isinstance(resolution, int):
            raise ValueError(
                "Resolution for sampling must be integer type for Multipatch "
                "objects"
            )
        p_connect = _connec.make_hexa_volumes(_enforce_len(resolution, 3))
        n_elements_per_patch = p_connect.shape[0]
        n_vertices_per_patch = resolution**spline.para_dim
        connectivity = _np.empty(
            (n_elements_per_patch * len(spline.patches), p_connect.shape[1])
        )

        for i in range(len(spline.patches)):
            ids = slice(
                i * n_elements_per_patch,
                (i + 1) * n_elements_per_patch,
                None,
            )
            connectivity[ids] = p_connect + i * n_vertices_per_patch
    else:
        connectivity = _connec.make_hexa_volumes(_enforce_len(resolution, 3))

    # Create volumes
    volumes = _Volumes(
        vertices=spline.sample(resolution),
        volumes=connectivity,
    )

    if watertight:
        volumes.merge_vertices()

    return volumes


def control_points(spline):
    """Extracts control points and return as vertices. Same can be achieved by
    doing `gustaf.Vertices(spline.control_points)`

    Parameters
    -----------
    spline: Spline / Multipatch
      (`self`-argument if called via extract member of a spline)

    Returns
    --------
    cps_as_Vertices: Vertices
    """
    return _Vertices(spline.control_points)


def control_edges(spline):
    """Extract control edges (mesh). Valid iff para_dim is 1.

    Parameters
    -----------
    spline: Spline / Multipatch
      (`self`-argument if called via extract member of a spline)

    Returns
    --------
    edges: Edges
    """
    from splinepy import Multipatch as _Multipatch

    if spline.para_dim != 1:
        raise ValueError("Invalid spline type!")

    if isinstance(spline, _Multipatch):
        # @todo avoid loop and transfer range_to_edges to cpp
        return _Edges.concat([control_edges(s) for s in spline.patches])
    else:
        return _Edges(
            vertices=spline.control_points,
            edges=_connec.range_to_edges(
                len(spline.control_points), closed=False
            ),
        )


def control_faces(spline):
    """Extract control face (mesh). Valid iff para_dim is 2.

    Parameters
    -----------
    spline: Spline / Multipatch
      (`self`-argument if called via extract member of a spline)

    Returns
    --------
    faces: Faces
    """
    from splinepy import Multipatch as _Multipatch

    if spline.para_dim != 2:
        raise ValueError("Invalid spline type!")

    if isinstance(spline, _Multipatch):
        # @todo avoid loop and transfer range_to_edges to cpp
        return _Faces.concat([control_faces(s) for s in spline.patches])
    else:
        return _Faces(
            vertices=spline.control_points,
            faces=_connec.make_quad_faces(spline.control_mesh_resolutions),
        )


def control_volumes(spline):
    """Extract control volumes (mesh). Valid iff para_dim is 3.

    Parameters
    -----------
    spline: Spline / Multipatch
      (`self`-argument if called via extract member of a spline)
    Returns
    --------
    volumes: Volumes
    """
    from splinepy import Multipatch as _Multipatch

    if spline.para_dim != 3:
        raise ValueError("Invalid spline type!")

    if isinstance(spline, _Multipatch):
        # @todo avoid loop and transfer range_to_edges to cpp
        return _Volumes.concat([control_volumes(s) for s in spline.patches])
    else:
        return _Volumes(
            vertices=spline.control_points,
            volumes=_connec.make_hexa_volumes(spline.control_mesh_resolutions),
        )


def control_mesh(spline):
    """Calls control_edges, control_faces, control_volumes based on current
    spline.

    Parameters
    -----------
    spline: Spline / Multipatch
      (`self`-argument if called via extract member of a spline)

    Returns
    --------
    control_mesh: Edges or Faces or Volumes
    """
    if spline.para_dim == 1:
        return control_edges(spline)
    elif spline.para_dim == 2:
        return control_faces(spline)
    elif spline.para_dim == 3:
        return control_volumes(spline)
    else:
        raise ValueError(
            "Invalid para_dim to extract control_mesh. Supports 1 to 3."
        )


def spline(spline, para_dim, split_plane):
    """Extract a sub spline from a given representation.

    Parameters
    ----------
    spline: Spline / Multipatch
      (`self`-argument if called via extract member of a spline)
    para_dim : int
      parametric dimension to be extracted
    split_plane : float / tuple<float, float>
      interval or value in parametric space to be extracted from the spline
      representation

    Returns
    -------
    spline
    """
    from splinepy.spline import Spline as _Spline

    # Check type
    if not isinstance(spline, _Spline):
        raise TypeError("Unknown spline representation passed to sub spline")

    # Check arguments for sanity
    if para_dim > spline.para_dim:
        raise ValueError(
            "Requested parametric dimension exceeds spline's parametric"
            " dimensionality."
        )

    # create 1d array for easy type / size check
    split_plane = _np.asanyarray(split_plane).ravel()
    if not _np.issubdtype(split_plane.dtype, _np.floating):
        raise ValueError("split_plane must be floating type.")

    if split_plane.size > 2:
        raise ValueError("split_plane must have size of 1 or 2.")

    is_interval = len(split_plane) == 2

    # Check if is bezier-type
    if not spline.has_knot_vectors:
        spline_copy = spline.nurbs if spline.is_rational else spline.bspline
    else:
        spline_copy = spline.copy()

    for _ in range(spline_copy.degrees[para_dim] + 1):
        # Will do nothing if spline already has sufficient number of knots
        # at given position go to c^-1 extraction
        spline_copy.insert_knots(para_dim, split_plane)

    # start and end id. indices correspond to [first dim][first appearance]
    start_id = (
        _np.where(
            abs(spline_copy.knot_vectors[para_dim].numpy() - split_plane[0])
            < _settings.TOLERANCE
        )[0][-1]
        - spline_copy.degrees[para_dim]
    )

    end_id = (
        _np.where(
            abs(spline_copy.knot_vectors[para_dim].numpy() - split_plane[-1])
            < _settings.TOLERANCE
        )[0][-1]
        - spline_copy.degrees[para_dim]
    )

    # Use MultiIndex for extraction
    mi = spline_copy.multi_index
    indices = [slice(None, None, None) for _ in range(spline_copy.para_dim)]
    if not is_interval:
        # Check if boundary
        if (
            abs(spline_copy.knot_vectors[para_dim][-1] - split_plane[-1])
            < _settings.TOLERANCE
        ):
            start_id -= 1
        indices[para_dim] = start_id
    else:
        indices[para_dim] = slice(start_id, end_id, None)
    indices = mi[tuple(indices)]

    # Return new_spline
    spline_info = {}
    spline_info["control_points"] = spline_copy.cps[indices]
    spline_info["degrees"] = spline_copy.degrees.tolist()
    if not is_interval:
        spline_info["degrees"].pop(para_dim)
    if spline.has_knot_vectors:
        spline_info["knot_vectors"] = spline_copy.knot_vectors.copy()
        if not is_interval:
            spline_info["knot_vectors"].pop(para_dim)
        else:
            spline_info["knot_vectors"][para_dim] = spline_copy.knot_vectors[
                para_dim
            ][start_id : (end_id + spline_copy.degrees[para_dim] + 1)]

    if spline.is_rational:
        spline_info["weights"] = spline_copy.weights[indices]

    return type(spline)(**spline_info)


def boundaries(spline, boundary_ids=None):
    r"""
    Extracts boundary spline.

    The boundaries deducted from the parametric axis which is normal to the
    boundary (:math:`j`), if the boundary is at parametric axis position
    :math:`x_{j}=x_{j_{min}}` the corresponding boundary is :math:`2j`, else at
    parametric axis position :math:`x_j=x_{j_{min}}` the boundary is
    :math:`2j+1`


    Parameters
    -----------
    spline: Spline / Multipatch
      (`self`-argument if called via extract member of a spline)
    boundary_ids: list
      Only considered for Spline. Default is None and returns all boundaries.

    Returns
    -------
    boundary_spline: type(self)
        boundary spline, which has one less para_dim
    """
    from splinepy import Multipatch as _Multipatch
    from splinepy import splinepy_core as _core

    # Pass to respective c++ implementation
    if isinstance(spline, _Multipatch):
        return spline.boundary_multipatch()
    else:
        bids = [] if boundary_ids is None else list(boundary_ids)
        return [
            type(spline)(spline=c)
            for c in _core.extract_boundaries(spline, bids)
        ]


def arrow_data(spline, adata_name):
    """
    Creates edges that represent arrow_data.
    This function respects entries in show_options.

    Parameters
    ----------
    spline: Spline
    adata_name: str

    Returns
    -------
    adata_edges: gus.Edges
    """
    if adata_name not in spline.spline_data:
        raise KeyError(
            f"given splines does not have arrow data - {adata_name}."
        )

    default_res = 10 if spline.name.startswith("Multi") else 100
    res = spline.show_options.get("resolutions", default_res)
    sampled_spline = _visualize._sample_spline(spline, res)

    adata = _visualize._sample_arrow_data(
        spline, adata_name, sampled_spline, _enforce_len(res, spline.para_dim)
    )

    # if `on` is specified, arrow_data will be a vertices, else it's processed
    # within sampled_spline
    gus_mesh = sampled_spline if adata is None else adata

    arrow_data_value = gus_mesh.vertex_data.as_arrow(adata_name, None, True)

    # if data and origin does not have save dimension, they can't be
    # represented as edges.
    # We can either raise error or return gus_mesh. We will do the latter
    # One can extract origin and value using gus_mesh.vertices and
    # gus_mesh.vertex_data[adata_name]
    if arrow_data_value.shape[1] != spline.dim:
        _log.warning(
            "dimension mismatch between arrow_data and spline."
            "returning sampled mesh as is (not as gus.Edges)."
        )
        return gus_mesh

    as_edges = _from_data(
        gus_mesh,
        arrow_data_value,
        gus_mesh.show_options.get("arrow_data_scale", None),
        data_norm=gus_mesh.vertex_data.as_scalar(adata),
    )

    return as_edges


def bases(spline, ids=None, parametric_view=True):
    """
    Creates spline that represents basis functions of given ids.
    Parametric bounds equals to the given spline.

    Parameters
    ----------
    spline: Spline
    ids: list
      Supports ids of basis functions. Default is None and returns all supports

    Returns
    -------
    basis_splines: list
    """
    # if none, extract all
    if ids is None:
        ids = _np.arange(len(spline.control_points))

    # type check
    if not isinstance(ids, (list, tuple, _np.ndarray)):
        raise TypeError(f"Bases ids ({type(ids)}) should be a list.")

    # first, create conforming parametric view
    if parametric_view:
        conform_para_view = spline.create.parametric_view(
            axes=False, conform=True
        )
        base_basis = conform_para_view.create.embedded(spline.para_dim + 1)
    else:
        base_basis = spline.create.embedded(spline.dim + 1)

    # for template, the last axis should be all zero
    # also, set control point show_option to False
    base_basis.control_points[:, -1] = 0.0
    base_basis.show_options["control_points"] = False

    # loop and e
    basis_splines = []
    for i in ids:
        basis_i = base_basis.copy()
        basis_i.control_points[i, -1] = 1.0
        basis_i.show_options["c"] = int(i)
        basis_splines.append(basis_i)

    return basis_splines


class Extractor:
    """Helper class to allow direct extraction from spline obj (BSpline or
    NURBS). Internal use only.

    Examples
    ---------
    >>> spline_faces = my_spline.extract.faces()
    """

    __slots__ = ("_helpee",)

    def __init__(self, spl):
        from splinepy import Multipatch as _Multipatch
        from splinepy import Spline as _Spline

        if not isinstance(spl, (_Spline, _Multipatch)):
            raise ValueError("Extractor expects a Spline or Multipatch type")
        self._helpee = spl

    @_wraps(edges)
    def edges(self, *args, **kwargs):
        return edges(self._helpee, *args, **kwargs)

    @_wraps(faces)
    def faces(self, *args, **kwargs):
        return faces(self._helpee, *args, **kwargs)

    @_wraps(volumes)
    def volumes(self, *args, **kwargs):
        return volumes(self._helpee, *args, **kwargs)

    @_wraps(control_points)
    def control_points(self):
        return control_points(self._helpee)

    @_wraps(control_edges)
    def control_edges(self):
        return control_edges(self._helpee)

    @_wraps(control_faces)
    def control_faces(self):
        return control_faces(self._helpee)

    @_wraps(control_points)
    def control_volumes(self):
        return control_volumes(self._helpee)

    @_wraps(control_mesh)
    def control_mesh(self):
        return control_mesh(self._helpee)

    def beziers(self):
        """
        Extract all (rational) Bezier patches within a given spline

        Parameters
        ----------
        None

        Returns
        -------
        beziers : list
          List of all individual bezier patches representing the non-zero
          knot-spans
        """
        if not self._helpee.has_knot_vectors:
            return [self._helpee]
        return self._helpee.extract_bezier_patches()

    @_wraps(boundaries)
    def boundaries(self, *args, **kwargs):
        return boundaries(self._helpee, *args, **kwargs)

    def spline(self, splitting_plane=None, interval=None):
        """Extract a spline from a spline.

        Use a (number of) splitting planes to extract a subsection from the
        parametric domain of it.

        Parameters
        ----------
        splitting_plane : int or dictionary {int : (float, float)}
          if int : parametric dimension to be extracted
          if dict : list of splitting planes and ranges to be passed
        interval : float or tuple<float,float>
          interval or value in parametric space to be extracted from the
          spline representation

        Returns
        -------
        spline
        """
        if isinstance(splitting_plane, dict):
            if interval is not None:
                raise ValueError("Arguments incompatible expect dictionary")
            splitting_plane = dict(
                sorted(splitting_plane.items(), key=lambda x: x[0])[::-1]
            )
            spline_copy = self._helpee.copy()
            for key, item in splitting_plane.items():
                spline_copy = spline(spline_copy, key, item)
            return spline_copy
        else:
            return spline(self._helpee, splitting_plane, interval)

    @_wraps(arrow_data)
    def arrow_data(self, *args, **kwargs):
        return arrow_data(self._helpee, *args, **kwargs)

    @_wraps(bases)
    def bases(self, *args, **kwargs):
        return bases(self._helpee, *args, **kwargs)
