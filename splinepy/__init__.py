from splinepy._version import __version__
from splinepy import splinepy_core

from splinepy import spline
from splinepy import bezier
from splinepy import rational_bezier
from splinepy import bspline
from splinepy import nurbs
from splinepy import utils
from splinepy import load
from splinepy import io
from splinepy import settings

from splinepy.spline import Spline
from splinepy.bezier import Bezier
from splinepy.rational_bezier import RationalBezier
from splinepy.bspline import BSpline
from splinepy.nurbs import NURBS
from splinepy.load import (load_splines, load_solution)

# set NAME_TO_TYPE
settings.NAME_TO_TYPE = settings.__splinepy_name_to_type__()

# configure logging
utils.log.configure()

__all__ = [
        "__version__",
        "splinepy_core",
        "spline",
        "bezier",
        "rational_bezier",
        "bspline",
        "nurbs",
        "load",
        "utils",
        "io",
        "settings",
        "Spline",
        "Bezier",
        "RationalBezier",
        "BSpline",
        "NURBS",
        "load_splines",
        "load_solution",
]
