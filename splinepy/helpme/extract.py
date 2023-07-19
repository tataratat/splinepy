import numpy as np
from gustaf import Edges, Faces, Vertices, Volumes
from gustaf.utils import arr, connec

from splinepy import settings
from splinepy.utils.data import cartesian_product


def edges(
    spline,
    resolution=100,
    extract_dim=None,
    extract_knot=None,
    all_knots=False,
):
    """Extract edges (lines) from a given spline. Only entity you can extract
    without dimension limit.

    Parameters
    -----------
    spline: Spline
    resolution: int
    extract_dim: int
      Parametric dimension to extract.
    extract_knot: list
      (spline.para_dim - 1,) shaped knot location along extract_dim
    all_knots: bool
      Switch to allow all knot-line extraction.

    Returns
    --------
    edges: Edges
    """
    if not all_knots:
        resolution = int(resolution)

    if spline.para_dim == 1:
        return Edges(
            vertices=spline.sample(resolution),
            edges=connec.range_to_edges(
                (0, resolution),
                closed=False,
            ),
        )

    else:
        # This should be possible for spline of any dimension.
        # As long as it satisfies the following condition
        if extract_knot is not None:
            if len(extract_knot) != spline.para_dim - 1:
                raise ValueError(
                    "Must satisfy len(extract_knot) == spline.para_dim -1."
                )

        # This may take awhile.
        if all_knots:
            temp_edges = []  # edges' is not a valid syntax
            unique_knots = np.array(spline.unique_knots, dtype=object)
            for i in range(spline.para_dim):
                mask = np.ones(spline.para_dim, dtype=bool)
                mask[i] = False
                # gather knots along current knot
                extract_knot_queries = cartesian_product(
                    unique_knots[mask], reverse=False
                )

                for ekq in extract_knot_queries:
                    temp_edges.append(
                        edges(spline, resolution[i], i, ekq, False)
                    )

            return Edges.concat(temp_edges)

        # Get parametric points to extract
        queries = np.empty(
            (resolution, spline.para_dim),
            dtype="float64",  # hardcoded for splinelibpy
            order="C",  # hardcoded for splinelibpy
        )
        # get ~extract_dim
        not_ed = np.arange(spline.para_dim).tolist()
        not_ed.pop(extract_dim)
        queries[:, not_ed] = extract_knot

        # get knot extrema
        uniq_knots = spline.unique_knots[extract_dim]
        min_knot_position = min(uniq_knots)
        max_knot_position = max(uniq_knots)

        queries[:, extract_dim] = np.linspace(
            min_knot_position,
            max_knot_position,
            resolution,
        )

        return Edges(
            vertices=spline.evaluate(queries),
            edges=connec.range_to_edges(
                (0, resolution),
                closed=False,
            ),
        )


def faces(
    spline,
    resolutions,
    watertight=True,
):
    """Extract faces from spline. Valid iff para_dim is one of the followings:
    {2, 3}. In case of {3}, it will return only surfaces. If internal faces are
    desired, used `spline.extract.volumes().faces()`. Note that dimension
    higher than 3 is not showable.

    Parameters
    -----------
    spline: BSpline or NURBS
    resolutions: int or list
    watertight: bool
      Default is True. Only related to para_dim = 3 splines. If False,
      overlapping vertices at boundary edges won't be merged.

    Returns
    --------
    faces: faces
    """
    resolutions = arr.enforce_len(resolutions, spline.para_dim)

    if spline.para_dim == 2:
        return Faces(
            vertices=spline.sample(resolutions),
            faces=connec.make_quad_faces(resolutions),
        )

    elif spline.para_dim == 3:
        # extract boundaries and sample from them
        # alternatively, as each groups share basis and will have same
        # resolution query, we could reuse the basis.
        boundaries = spline.extract_boundaries()
        grouped_boundaries = [
            boundaries[i * 2 : (i + 1) * 2] for i in range(spline.para_dim)
        ]

        list_res = list(resolutions)
        vertices = []
        faces = []
        offset = 0
        for i, gb in enumerate(grouped_boundaries):
            this_res = list_res.copy()  # shallow
            this_res.pop(i)

            tmp_faces = connec.make_quad_faces(this_res)
            offset_size = np.prod(this_res)

            for g in gb:  # each spline
                vertices.append(g.sample(this_res))
                faces.append(tmp_faces + int(offset))
                offset += offset_size

        # make faces and merge vertices before returning
        f = Faces(vertices=np.vstack(vertices), faces=np.vstack(faces))

        if watertight:
            f.merge_vertices()

        return f

    else:
        raise ValueError("Invalid spline to make faces.")


def volumes(spline, resolutions):
    """Extract volumes from spline. Valid iff spline.para_dim == 3.

    Parameters
    -----------
    spline: BSpline or NURBS
    resolutions:

    Returns
    --------
    volumes: Volumes
    """
    if spline.para_dim != 3:
        raise ValueError(
            "Volume extraction from a spline is only valid for "
            "para_dim: 3 dim: 3 splines."
        )

    return Volumes(
        vertices=spline.sample(resolutions),
        volumes=connec.make_hexa_volumes(resolutions),
    )


def control_points(spline):
    """Extracts control points and return as vertices. Same can be achieved by
    doing `gustaf.Vertices(spline.control_points)`

    Parameters
    -----------
    spline: BSpline or NURBS

    Returns
    --------
    cps_as_Vertices: Vertices
    """
    return Vertices(spline.control_points)


def control_edges(spline):
    """Extract control edges (mesh). Valid iff para_dim is 1.

    Parameters
    -----------
    edges: BSpline or NURBS

    Returns
    --------
    edges: Edges
    """
    if spline.para_dim != 1:
        raise ValueError("Invalid spline type!")

    return Edges(
        vertices=spline.control_points,
        edges=connec.range_to_edges(len(spline.control_points), closed=False),
    )


def control_faces(spline):
    """Extract control face (mesh). Valid iff para_dim is 2.

    Parameters
    -----------
    spline: BSpline or NURBS

    Returns
    --------
    faces: Faces
    """
    if spline.para_dim != 2:
        raise ValueError("Invalid spline type!")

    return Faces(
        vertices=spline.control_points,
        faces=connec.make_quad_faces(spline.control_mesh_resolutions),
    )


def control_volumes(spline):
    """Extract control volumes (mesh). Valid iff para_dim is 3.

    Parameters
    -----------
    spline: BSpline or NURBS

    Returns
    --------
    volumes: Volumes
    """
    if spline.para_dim != 3:
        raise ValueError("Invalid spline type!")

    return Volumes(
        vertices=spline.control_points,
        volumes=connec.make_hexa_volumes(spline.control_mesh_resolutions),
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
            "Invalid para_dim to extract control_mesh. " "Supports 1 to 3."
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
    from splinepy.spline import Spline

    # Check type
    if not isinstance(spline, Spline):
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
        split_plane = list([split_plane])

    # Check if is bezier-type
    is_bezier = "Bezier" in spline.whatami
    is_rational = "weights" in spline.required_properties
    if is_bezier:
        if is_rational:
            spline_copy = spline.nurbs
        else:
            spline_copy = spline.bspline
    else:
        spline_copy = spline.copy()

    for _ in range(spline_copy.degrees[para_dim]):
        # Will do nothing if spline already has sufficient number of knots
        # at given position
        spline_copy.insert_knots(para_dim, split_plane)

    # Start extraction
    cps_res = spline_copy.control_mesh_resolutions
    # start and end id. indices correspond to [first dim][first appearance]
    start_id = np.where(
        abs(spline_copy.knot_vectors[para_dim].numpy() - split_plane[0])
        < settings.TOLERANCE
    )[0][0]
    end_id = np.where(
        abs(spline_copy.knot_vectors[para_dim].numpy() - split_plane[-1])
        < settings.TOLERANCE
    )[0][0]
    para_dim_ids = np.arange(np.prod(cps_res))
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

            spline_info["knot_vectors"][para_dim] = np.concatenate(
                ([start_knot], knots_in_between, [end_knot])
            )

    if is_rational:
        spline_info["weights"] = spline_copy.weights[
            (para_dim_ids >= start_id) & (para_dim_ids <= end_id)
        ]

    return type(spline)(**spline_info)


class Extractor:
    """Helper class to allow direct extraction from spline obj (BSpline or
    NURBS). Internal use only.

    Examples
    ---------
    >>> my_spline = <your-spline>
    >>> spline_faces = my_spline.extract.faces()
    """

    def __init__(self, spl):
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
        return self._spline.extract_boundaries(*args, **kwargs)

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
