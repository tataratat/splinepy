"""splinepy/splinepy/settings.py
Convenient global variables.
"""


def __splinepy_name_to_type__():
    """Workaround to provide flexible string to type conversion without
    causing circular import
    """
    from splinepy import NURBS, Bezier, BSpline, Multipatch, RationalBezier

    return dict(
        Bezier=Bezier,
        RationalBezier=RationalBezier,
        BSpline=BSpline,
        NURBS=NURBS,
        MultiPatch=Multipatch,
    )


TOLERANCE = 1e-11
"""
Default tolerance for any operation that involves floats
"""

NTHREADS = 1
"""
Default number of threads.
"""

NAME_TO_TYPE = None
"""
String to Type conversion. Any function that may return different types should
use this. If you have derived classes, replace this with your own conversion
dict. In splinepy, this is set in __init__.
"""

CHECK_BOUNDS = True
"""
Bool to check bounds of queries if requested. Can be set to false to
accelerate process
"""
