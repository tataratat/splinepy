"""splinepy/splinepy/settings.py
Convenient global variables.
"""

from splinepy import Bezier, RationalBezier, BSpline, NURBS

TOLERANCE = 1e-11
"""
Default tolerance for any operation that involves floats
"""

NTHREADS = 1
"""
Default number of threads.
"""

NAME_TO_TYPE = dict(
        Bezier=Bezier,
        RationalBezier=RationalBezier,
        BSpline=BSpline,
        NURBS=NURBS,
)
"""
String to Type conversion. Any function that may return different types should
use this. If you have derived classes, replace this with your own conversion
dict.
"""
