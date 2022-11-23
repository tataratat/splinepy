from splinepy.bezier import BezierBase


class RationalBezier(BezierBase):

    __slots__ = ()

    def __init__(
            self,
            degrees=None,
            control_points=None,
            weights=None,
            spline=None
    ):
        """
        RationalBezier (Spline).

        Parameters
        -----------
        degrees: (para_dim,) array-like
        control_points: (m, dim) array-like
        weights: (m) array-like
        spline: Spline

        Returns
        --------
        None
        """
        super().__init__(
                spline=spline,
                degrees=degrees,
                control_points=control_points,
                weights=weights,
        )
