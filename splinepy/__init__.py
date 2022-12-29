from splinepy import (
    bezier,
    bspline,
    io,
    load,
    multipatch,
    nurbs,
    rational_bezier,
    settings,
    spline,
    splinepy_core,
    utils,
)
from splinepy._version import __version__
from splinepy.bezier import Bezier
from splinepy.bspline import BSpline
from splinepy.load import load_solution, load_splines
from splinepy.multipatch import Multipatch
from splinepy.nurbs import NURBS
from splinepy.rational_bezier import RationalBezier
from splinepy.spline import Spline

# set NAME_TO_TYPE
settings.NAME_TO_TYPE = settings.__splinepy_name_to_type__()

# configure logging
utils.log.configure()

__all__ = [
    "__version__",
    "bezier",
    "bspline",
    "io",
    "load",
    "multipatch",
    "nurbs",
    "rational_bezier",
    "spline",
    "splinepy_core",
    "utils",
    "settings",
    "Bezier",
    "BSpline",
    "Multipatch",
    "NURBS",
    "RationalBezier",
    "Spline",
    "load_splines",
    "load_solution",
]
