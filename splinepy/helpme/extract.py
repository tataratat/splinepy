import numpy as _np
from gustaf import Edges as _Edges
from gustaf import Faces as _Faces
from gustaf import Vertices as _Vertices
from gustaf import Volumes as _Volumes
from gustaf.utils import connec as _connec
from gustaf.utils.arr import enforce_len as _enforce_len

from splinepy import settings as _settings
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
    resolution: int
      samples per parametric dimension
    all_knots: bool
      Switch to allow all knot-line extraction or just contour (para_dim>0)

    Returns
    --------
    edges: Edges
    """
    from splinepy import Multipatch as _Multipatch

    # Check resolution input
    if not isinstance(resolution, int):
        raise ValueError("Resolution must be integer-type")

    if spline.para_dim == 1:
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
                [_np.linspace(*spline.parametric_bounds[:, i], resolution)]
                + split_knots,
                reverse=True,
            )[:, reorder_mask]

            # Retrieve connectivity to be repeated
            single_line_connectivity = _connec.range_to_edges(
                (0, resolution),
                closed=False,
            )

            temp_edges.append(
                _Edges(
                    vertices=spline.evaluate(extract_knot_queries),
                    edges=_np.vstack(
                        [
                            single_line_connectivity + resolution * j
                            for j in range(n_knot_lines)
                        ]
                    ),
                )
            )

        return _Edges.concat(temp_edges)


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
    spline: BSpline or NURBS
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
        # Extract boundaries and sample from boundaries
        if isinstance(spline, _Multipatch):
            boundaries = spline.boundary_multipatch()
        else:
            boundaries = _Multipatch(splines=spline.extract.boundaries())

        n_faces = (resolution - 1) ** 2
        vertices = resolution**2
        f_loc = _connec.make_quad_faces(_enforce_len(resolution, 2))
        face_connectivity = _np.empty((n_faces * len(boundaries.patches), 4))

        # Create Connectivity for Multipatches
        for i in range(len(boundaries.patches)):
            face_connectivity[i * n_faces : (i + 1) * n_faces] = f_loc + (
                i * vertices
            )

        # make faces and merge vertices before returning
        faces = _Faces(
            vertices=boundaries.sample(resolution), faces=face_connectivity
        )

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
      spline to be meshed into volumes
    resolution: int
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
    None

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
    para_dim : int
      parametric dimension to be extract ted
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
    if isinstance(split_plane, list):
        if not (
            (len(split_plane) == 2)
            and (isinstance(split_plane[0], float))
            and (isinstance(split_plane[0], float))
        ):
            raise ValueError(
                "Range must be float or tuple of floats with length 2"
            )
    elif not isinstance(split_plane, float):
        raise ValueError(
            "Range must be float or tuple of floats with length 2"
        )
    else:
        # Convert float to tuple to facilitate
        split_plane = [split_plane]

    # Check if is bezier-type
    is_bezier = "Bezier" in spline.whatami
    is_rational = "weights" in spline.required_properties
    if is_bezier:
        spline_copy = spline.nurbs if is_rational else spline.bspline
    else:
        spline_copy = spline.copy()

    for _ in range(spline_copy.degrees[para_dim]):
        # Will do nothing if spline already has sufficient number of knots
        # at given position
        spline_copy.insert_knots(para_dim, split_plane)

    # Start extraction
    cps_res = spline_copy.control_mesh_resolutions
    # start and end id. indices correspond to [first dim][first appearance]
    start_id = _np.where(
        abs(spline_copy.knot_vectors[para_dim].numpy() - split_plane[0])
        < _settings.TOLERANCE
    )[0][0]
    end_id = _np.where(
        abs(spline_copy.knot_vectors[para_dim].numpy() - split_plane[-1])
        < _settings.TOLERANCE
    )[0][0]
    para_dim_ids = _np.arange(_np.prod(cps_res))
    for i_pd in range(para_dim):
        para_dim_ids -= para_dim_ids % cps_res[i_pd]
        para_dim_ids = para_dim_ids // cps_res[i_pd]
    # indices are shifted by one
    para_dim_ids = para_dim_ids % cps_res[para_dim] + 1

    # Return new_spline
    spline_info = {}
    spline_info["control_points"] = spline_copy.cps[
        (para_dim_ids >= start_id) & (para_dim_ids <= end_id)
    ]
    spline_info["degrees"] = spline_copy.degrees.tolist()
    if start_id == end_id:
        spline_info["degrees"].pop(para_dim)
    if not is_bezier:
        spline_info["knot_vectors"] = spline_copy.knot_vectors.copy()
        if start_id == end_id:
            spline_info["knot_vectors"].pop(para_dim)
        else:
            start_knot = spline_copy.knot_vectors[para_dim][start_id]
            knots_in_between = spline_copy.knot_vectors[para_dim][
                start_id : (end_id + spline_copy.degrees[para_dim])
            ]
            end_knot = spline_copy.knot_vectors[para_dim][
                (end_id + spline_copy.degrees[para_dim] - 1)
            ]

            spline_info["knot_vectors"][para_dim] = _np.concatenate(
                ([start_knot], knots_in_between, [end_knot])
            )

    if is_rational:
        spline_info["weights"] = spline_copy.weights[
            (para_dim_ids >= start_id) & (para_dim_ids <= end_id)
        ]

    return type(spline)(**spline_info)


def boundaries(spline, boundary_ids=None):
    """
    Extracts boundary spline.

    The boundaries deducted from the parametric axis which is normal to the
    boundary (j), if the boundary is at parametric axis position x_j=x_jmin
    the corresponding boundary is 2*j, else at parametric axis position
    x_j=x_jmin the boundary is 2*j+1


    Parameters
    -----------
    spline: Spline / Multipatch
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


class Extractor:
    """Helper class to allow direct extraction from spline obj (BSpline or
    NURBS). Internal use only.

    Examples
    ---------
    >>> spline_faces = my_spline.extract.faces()
    """

    def __init__(self, spl):
        from splinepy import Multipatch as _Multipatch
        from splinepy import Spline as _Spline

        if not isinstance(spl, (_Spline, _Multipatch)):
            raise ValueError("Extractor expects a Spline or Multipatch type")
        self._spline = spl

    def edges(self, *args, **kwargs):
        return edges(self._spline, *args, **kwargs)

    def faces(self, *args, **kwargs):
        return faces(self._spline, *args, **kwargs)

    def volumes(self, *args, **kwargs):
        return volumes(self._spline, *args, **kwargs)

    def control_points(self):
        return control_points(self._spline)

    def control_edges(self):
        return control_edges(self._spline)

    def control_faces(self):
        return control_faces(self._spline)

    def control_volumes(self):
        return control_volumes(self._spline)

    def control_mesh(self):
        return control_mesh(self._spline)

    def beziers(self):
        if not self._spline.has_knot_vectors:
            return [self._spline]
        return self._spline.extract_bezier_patches()

    def boundaries(self, *args, **kwargs):
        return boundaries(self._spline, *args, **kwargs)

    def spline(self, splitting_plane=None, interval=None):
        """Extract a spline from a spline.

        Use a (number of) splitting planes to extract a subsection from the
        parametric domain of it.

        Parameters
        ----------
        splitting_plane : int / dictionary (int : (float))
          if integer : parametric dimension to be extracted
          if dictionary : list of splitting planes and ranges to be passed
        interval : float / tuple<float,float>
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
            spline_copy = self._spline.copy()
            for key, item in splitting_plane.items():
                spline_copy = spline(spline_copy, key, item)
            return spline_copy
        else:
            return spline(self._spline, splitting_plane, interval)
