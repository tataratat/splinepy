from splinepy._version import __version__
from splinepy import splinepy_core

from splinepy import bezier
from splinepy import bspline
from splinepy import io
from splinepy import load
from splinepy import multipatch
from splinepy import nurbs
from splinepy import rational_bezier
from splinepy import spline
from splinepy import settings
from splinepy import utils

from splinepy.bezier import Bezier
from splinepy.bspline import BSpline
from splinepy.multipatch import Multipatch
from splinepy.nurbs import NURBS
from splinepy.rational_bezier import RationalBezier
from splinepy.spline import Spline
from splinepy.load import (load_splines, load_solution)

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
