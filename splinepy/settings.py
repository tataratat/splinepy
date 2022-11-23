"""splinepy/splinepy/settings.py
Convenient global variables.
"""


def __splinepy_name_to_type__():
    """Workaround to provide flexible string to type conversion without
    causing circular import
    """
    from splinepy import Bezier, RationalBezier, BSpline, NURBS
    return dict(
            Bezier=Bezier,
            RationalBezier=RationalBezier,
            BSpline=BSpline,
            NURBS=NURBS,
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
