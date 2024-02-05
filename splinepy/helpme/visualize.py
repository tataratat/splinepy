import gustaf as _gus
import numpy as _np
from gustaf import Vertices as _Vertices
from gustaf.helpers import options as _options
from gustaf.utils.arr import enforce_len as _enforce_len

from splinepy import settings as _settings
from splinepy.utils import log as _log

_vedo_spline_common_options = (
    _options.Option("vedo", "knots", "Show spline's knots.", (bool,)),
    _options.Option("vedo", "knot_c", "Knot color.", (str, tuple, list, int)),
    _options.Option(
        "vedo",
        "knot_lw",
        "Line width of knots. Number of pixels. "
        "Applicable to para_dim > 1.",
        (float, float),
    ),
    _options.Option(
        "vedo",
        "knot_alpha",
        "Transparency of knots in range [0, 1]. Applicable to para_dim > 1",
        (float, int),
    ),
    _options.Option(
        "vedo",
        "control_points",
        "Show spline's control points."
        "Options propagates to control mesh, unless it is specified.",
        (bool,),
    ),
    _options.Option(
        "vedo",
        "control_point_r",
        "Control point radius",
        (float, int),
    ),
    _options.Option(
        "vedo",
        "control_point_c",
        "Color of control_points in {rgb, RGB, str of (hex, name), int}",
        (str, tuple, list, int),
    ),
    _options.Option(
        "vedo",
        "control_point_alpha",
        "Transparency of control points in range [0, 1].",
        (float, int),
    ),
    _options.Option(
        "vedo",
        "control_point_ids",
        "Show ids of control_points.",
        (bool,),
    ),
    _options.Option(
        "vedo",
        "control_mesh",
        "Show spline's control mesh.",
        (bool,),
    ),
    _options.Option(
        "vedo",
        "control_mesh_c",
        "Color of control_mesh in {rgb, RGB, str of (hex, name), int}",
        (str, tuple, list, int),
    ),
    _options.Option(
        "vedo",
        "control_mesh_lw",
        "Line width of control mesh. Number of pixels",
        (float, int),
    ),
    _options.Option(
        "vedo",
        "resolutions",
        "Sampling resolution for spline.",
        (int, list, tuple, _np.ndarray),
    ),
)


class SplineShowOption(_options.ShowOption):
    """
    Show options for splines.
    """

    __slots__ = ()

    # if we start to support more backends, most of this options should become
    # some sort of spline common.
    _valid_options = _options.make_valid_options(
        *_options.vedo_common_options,
        *_vedo_spline_common_options,
        _options.Option(
            "vedo",
            "fitting_queries",
            "Shows fitting queries if they are locally saved in splines.",
            (bool,),
        ),
        _options.Option(
            "vedo",
            "arrow_data_on",
            "Specify parametric coordinates to place arrow_data.",
            (list, tuple, _np.ndarray),
        ),
    )

    _helps = "Spline"

    def __init__(self, helpee):
        """
        Parameters
        ----------
        helpee: GustafSpline
        """
        self._helpee = helpee
        # checks if helpee inherits from GustafSpline
        if self._helps not in str(type(helpee).__mro__):
            raise TypeError(
                f"This show option if for {self._helps}.",
                f"Given helpee is {type(helpee)}.",
            )
        self._options = {}


class MultipatchShowOption(_options.ShowOption):
    """
    Show options for Multipatches.
    """

    __slots__ = ()

    # if we start to support more backends, most of this options should become
    # some sort of spline common.
    _valid_options = _options.make_valid_options(
        *_options.vedo_common_options,
        *_vedo_spline_common_options,
    )

    _helps = "Multipatch"

    def __init__(self, helpee):
        """
        Parameters
        ----------
        helpee: GustafSpline
        """
        self._helpee = helpee
        # checks if helpee inherits from GustafSpline
        if self._helps not in str(type(helpee).__mro__):
            raise TypeError(
                f"This show option if for {self._helps}.",
                f"Given helpee is {type(helpee)}.",
            )
        self._options = {}


def _sample_knots_1d(spline):
    """
    Samples knots for splines of para_dim=1.
    Assumes this is one of the showable splines.
    """

    def sample_uks(s):
        """samples unique knots from single patch spline"""
        return s.evaluate(_np.asanyarray(s.unique_knots[0]).reshape(-1, 1))

    # sample unique knots in physical space and create vertices
    # for multi-patch, we will call this function for its patches
    if spline.name.startswith("Multi"):
        evaluated_knots = _np.vstack([sample_uks(p) for p in spline.patches])
    else:
        evaluated_knots = sample_uks(spline)

    # create vertices - for multipatch, we could remove duplicating vertices
    # it gets annoying
    knots = _Vertices(evaluated_knots)

    # apply show options to show this vertices as 'x'
    knots.show_options["labels"] = ["x"] * len(knots.vertices)
    knots.show_options["label_options"] = {
        "justify": "center",
        "c": spline.show_options.get("knot_c", "green"),
    }

    return knots


def _sample_knots_2d(spline, res):
    """
    Samples knots for splines of para_dim=2.
    Assumes this is one of the showable splines.
    """
    # extract knot lines as edges
    knot_lines = spline.extract.edges(res, all_knots=True)

    # apply show options
    knot_lines.show_options["c"] = spline.show_options.get("knot_c", "black")
    knot_lines.show_options["lw"] = spline.show_options.get("knot_lw", 3)
    knot_lines.show_options["alpha"] = spline.show_options.get("knot_alpha", 1)

    return knot_lines


def _sample_knots(spline, res=None):
    """
    Samples knots for splines.
    Assumes this is one of the showable splines
    """
    # call sample functions for corresponding para_dim
    para_dim = spline.para_dim
    if para_dim > 1:
        if res is None:
            raise ValueError("sample resolution not given.")
        # 2d and 3d splines has same knot sampling routine
        return _sample_knots_2d(spline, res)
    else:
        return _sample_knots_1d(spline)


def _sample_spline(spline, res):
    """
    Sample splines.
    Assumes this is one of the showable splines

    Parameters
    ----------
    None

    Returns
    -------
    gus_dict: dict
      dict of sampled spline and knots.
    """
    # sampling spline is quite simple, so we don't split functions here
    para_dim = spline.para_dim
    if para_dim == 1:
        sp = spline.extract.edges(res)
        sp.show_options["lw"] = 8
    elif para_dim > 1:
        # we don't return water tight meshes here
        # see (gustaf #120)
        sp = spline.extract.faces(res, watertight=False)

    return sp


def _process_scalar_field(spline, data, sampled_spline, res):
    """
    Assuming field data is defined, extracts field value and applies
    to sampled spline.
    """
    # get data name and try to sample scalar field
    sampled_spline_data = spline.spline_data.as_scalar(data, res, default=None)

    # apply sampled data as vertex data to sampled_spline
    if sampled_spline_data is not None:
        # transfer data
        sampled_spline.vertex_data[data] = sampled_spline_data

        # transfer options - maybe vectorized query?
        keys = (
            "vmin",
            "vmax",
            "scalarbar",
            "scalarbar3d",
            "cmap",
            "cmap_alpha",
        )
        spline.show_options.copy_valid_options(
            sampled_spline.show_options, keys
        )

        # mark that we want to see this data
        sampled_spline.show_options["data"] = data

    else:
        raise ValueError(f"No spline_data named ({data}) for {spline}.")


def _process_spline_color(spline, sampled_spline, res):
    """
    Processes color and scalar field.
    If both are defined, you will see the color of the scalar field.
    """
    # first, get defaults
    default_lighting = "glossy" if spline.dim > 2 else "off"
    default_color = "green" if spline.para_dim > 1 else "black"

    # apply basic color
    sampled_spline.show_options["c"] = spline.show_options.get(
        "c", default_color
    )

    # apply alpha and lighting
    sampled_spline.show_options["alpha"] = spline.show_options.get("alpha", 1)
    sampled_spline.show_options["lighting"] = spline.show_options.get(
        "lighting", default_lighting
    )

    # process scalar field if data is defined
    data = spline.show_options.get("data", None)
    if data is not None:
        _process_scalar_field(spline, data, sampled_spline, res)


def _sample_arrow_data(spline, adata, sampled_spline, res):
    """
    Process arrow_data. Also known as vector_data.
    """
    # early exit?
    if adata is None:
        return None

    # use getitem to retrieve adapted_adata. this will raise key error if
    # it doesn't exist
    adapted_adata = spline.spline_data[adata]

    # if location is specified, this will be a separate Vertices obj with
    # configured arrow_data. Multipatch might just return a multipatch,
    # so a careful step here.
    has_locations = getattr(adapted_adata, "has_locations", False)
    adata_on = "arrow_data_on" in spline.show_options.keys()  # noqa SIM118
    create_vertices = has_locations or adata_on

    # this case causes conflict of interest. raise
    if has_locations and adata_on:
        raise ValueError(
            f"arrow_data-({adata}) has fixed location, "
            "but and `arrow_data_on` is set.",
        )

    # two cases of arrow data sampling
    # 1 - same as sampled spline
    #      -> add arrow as vertex_data
    # 2 - at specified locations
    #      -> create vertices at specified locations
    #         then, add arrow as vertex_data to them
    if not create_vertices:
        # sample arrows and append to sampled spline.
        sampled_spline.vertex_data[adata] = spline.spline_data.as_arrow(
            adata, resolutions=res
        )

        # transfer options to sampled_spline
        keys = ("arrow_data_scale", "arrow_data_color", "arrow_data")
        spline.show_options.copy_valid_options(
            sampled_spline.show_options, keys
        )

    else:
        # prepare corresponding queries
        if getattr(adapted_adata, "has_locations", False):
            queries = adapted_adata.locations
            on = None
        else:
            queries = spline.show_options["arrow_data_on"]
            on = queries

        # bound /  dim check - only for splines
        if not spline.name.startswith("Multi"):
            bounds = spline.parametric_bounds
            if queries.shape[1] != len(bounds[0]):
                raise ValueError(
                    "Dimension mismatch: arrow_data locations-"
                    f"{queries.shape[1]} / para_dim-{spline.para_dim}."
                )
            # tolerance padding. may still cause issues in splinepy.
            # in that case, we will have to scale queries.
            lb_diff = queries.min(axis=0) - bounds[0] + _settings.TOLERANCE
            ub_diff = queries.max(axis=0) - bounds[1] - _settings.TOLERANCE
            if any(lb_diff < 0) or any(ub_diff > 0):
                raise ValueError(
                    f"Given locations for ({adata}) are outside the "
                    f"parametric bounds ({bounds})."
                )

        # get arrow
        adata_values = spline.spline_data.as_arrow(adata, on=on)

        # create vertices that can be shown as arrows
        loc_vertices = _Vertices(spline.evaluate(queries))
        loc_vertices.vertex_data[adata] = adata_values

        # transfer options
        keys = ("arrow_data_scale", "arrow_data_color", "arrow_data")
        spline.show_options.copy_valid_options(loc_vertices.show_options, keys)

        return loc_vertices


def _create_control_points(spline):
    """
    Creates control points.
    """
    # control points (vertices)
    cps = spline.extract.control_points()
    cps.show_options["c"] = spline.show_options.get("control_point_c", "red")
    cps.show_options["r"] = spline.show_options.get("control_point_r", 10)
    cps.show_options["alpha"] = spline.show_options.get(
        "control_point_alpha", 1.0
    )
    return cps


def _create_control_point_ids(spline):
    """
    Creates control point ids.
    """
    # a bit redundant, but nicely separable
    cp_ids = spline.extract.control_points()
    cp_ids.show_options["labels"] = _np.arange(len(cp_ids.vertices))
    cp_ids.show_options["label_options"] = {"font": "VTK"}

    return cp_ids


def _create_control_mesh(spline):
    """
    Creates control mesh.
    """
    # pure control mesh
    c_mesh = spline.extract.control_mesh()  # either edges or faces
    if spline.para_dim != 1:
        c_mesh = c_mesh.to_edges(unique=True)

    c_mesh.show_options["c"] = spline.show_options.get("control_mesh_c", "red")
    c_mesh.show_options["lw"] = spline.show_options.get("control_mesh_lw", 4)
    c_mesh.show_options["alpha"] = spline.show_options.get(
        "control_points_alpha", 0.8
    )

    return c_mesh


def _create_fitting_queries(spline):
    """
    Create fitting queries. This will raise if spline is not BSpline or
    fitting queries does not exist.
    """
    fqs = _Vertices(spline._fitting_queries)
    fqs.show_options["c"] = "blue"
    fqs.show_options["r"] = 10
    return fqs


def make_showable(spline):
    """
    Goes through common procedures for preparing showable splines.

    Parameters
    ----------
    None

    Returns
    -------
    gus_dict: dict
      dict of sampled spline as gustaf objects
    """
    # get sampling resolution
    default_res = 10 if spline.name.startswith("Multi") else 100
    res = spline.show_options.get(
        "resolutions", default_res
    )  # , spline.para_dim)

    # get spline and knots
    sampled_gus = {}
    sampled_spline = _sample_spline(spline, res)
    sampled_gus["spline"] = sampled_spline  # save to dict
    if spline.show_options.get("knots", True):
        sampled_gus["knots"] = _sample_knots(spline, res)

    # apply spline color - this includes scalar field
    _process_spline_color(spline, sampled_spline, res)

    # axis?
    sampled_spline.show_options["axes"] = spline.show_options.get(
        "axes", False
    )

    # process arrow data, this may return something if arrow_data_on is defined
    arrow_data = _sample_arrow_data(
        spline,
        spline.show_options.get("arrow_data", None),
        sampled_spline,
        _enforce_len(res, spline.para_dim),
    )
    if arrow_data is not None:
        sampled_gus["arrow_data"] = arrow_data

    # process control_points
    show_cps = spline.show_options.get("control_points", True)
    if show_cps:
        sampled_gus["control_points"] = _create_control_points(spline)

    # control_point_ids
    if spline.show_options.get("control_point_ids", show_cps):
        sampled_gus["control_point_ids"] = _create_control_point_ids(spline)

    # control_mesh
    if spline.show_options.get("control_mesh", show_cps):
        sampled_gus["control_mesh"] = _create_control_mesh(spline)

    # fitting queries
    if hasattr(spline, "_fitting_queries") and spline.show_options.get(
        "fitting_queries", True
    ):
        sampled_gus["fitting_queries"] = _create_fitting_queries(spline)

    return sampled_gus


def show(spline, **kwargs):
    """Shows splines with various options.

    They are excessively listed, so that it can be adjustable.

    Parameters
    -----------
    spline: Spline / Multipatch
    resolutions: int or (spline.para_dim,) array-like
    control_points: bool
    knots: bool
    fitting_queries: bool
    return_gustaf: bool
      Return dict of gustaf discrete objects, for example,
      {Vertices, Edges, Faces}, instead of opening a window
    return_showable: bool
      Return dict of showable objects.
    parametric_space: bool
      Only relevant for `vedo` backend.
    c: str
      Default is None. Black for curves, else green.
    alpha: float
    lighting: str
    control_point_ids: bool
    color_spline: str
    cmap: str

    Returns
    --------
    things_to_show: dict
      iff return_discrete==True, dict of gustaf objects that are showable.
      iff return_showable==True, dict of backend objects that are showable.
    """
    # Showing is only possible for following splines
    allowed_dim_combo = (
        (1, 2),
        (1, 3),
        (2, 2),
        (2, 3),
        (3, 3),
    )
    if (spline.para_dim, spline.dim) not in allowed_dim_combo:
        raise ValueError("Sorry, can't show given spline.")

    orig_show_options = None
    if kwargs:
        orig_show_options = spline.show_options
        spline._show_options = spline.__show_option__(spline)
        orig_show_options.copy_valid_options(spline.show_options)
        for key, value in kwargs.items():
            try:
                spline.show_options[key] = value
            except BaseException:
                _log.debug(
                    f"Skipping invalid option {key} for "
                    f"{spline.show_options._helps}."
                )
                continue

    # Prepare things to show dict.
    things_to_show = make_showable(spline)

    # set original options back
    if orig_show_options is not None:
        spline._show_options = orig_show_options

    para_space = kwargs.pop("parametric_space", False)
    if para_space:
        p_spline = spline.create.parametric_view()
        things_to_show["parametric_spline"] = p_spline

    if kwargs.get("return_gustaf", False):
        return things_to_show

    if kwargs.get("return_showable", False):
        return {key: value.showable() for key, value in things_to_show.items()}

    if para_space:
        things_to_show.pop("parametric_spline")

        return _gus.show(
            ["Parametric Space", p_spline],
            ["Physical Space", *things_to_show.values()],
            **kwargs,
        )

    return _gus.show(things_to_show, **kwargs)
