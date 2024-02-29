from _typeshed import Incomplete

from splinepy.bezier import BezierBase as _BezierBase

class RationalBezier(_BezierBase):
    def __init__(
        self,
        degrees: Incomplete | None = None,
        control_points: Incomplete | None = None,
        weights: Incomplete | None = None,
        spline: Incomplete | None = None,
    ) -> None: ...
    @property
    def rationalbezier(self): ...
    @property
    def nurbs(self): ...
