from gustaf import show  # shortcut for frequently used module / func

from splinepy import (
    bezier,
    bspline,
    helpme,
    io,
    microstructure,
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
from splinepy.helpme.ffd import FFD
from splinepy.microstructure.microstructure import Microstructure
from splinepy.multipatch import Multipatch
from splinepy.nurbs import NURBS
from splinepy.rational_bezier import RationalBezier
from splinepy.spline import Spline
from splinepy.utils.data import (
    SplineDataAdaptor,
    cartesian_product,
    make_matrix,
    uniform_query,
)

# set NAME_TO_TYPE
settings.NAME_TO_TYPE = settings.__splinepy_name_to_type__()

# configure logging
utils.log.configure()


def to_derived(spline):
    """
    Returns derived spline type based on NAME_TO_TYPE conversion.

    Parameters
    ----------
    spline: CoreSpline

    Returns
    -------
    derived_spline: DerivedSpline
    """
    return settings.NAME_TO_TYPE[spline.name](spline=spline)


__all__ = [
    "__version__",
    "bezier",
    "bspline",
    "io",
    "multipatch",
    "nurbs",
    "rational_bezier",
    "spline",
    "splinepy_core",
    "utils",
    "microstructure",
    "settings",
    "helpme",
    "Bezier",
    "BSpline",
    "Multipatch",
    "NURBS",
    "RationalBezier",
    "Spline",
    "SplineDataAdaptor",
    "FFD",
    "to_derived",
    "show",
    "Microstructure",
    "TileLib",
    "cartesian_product",
    "make_matrix",
    "uniform_query",
]
