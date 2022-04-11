import logging
import copy

import numpy as np

from splinepy import utils
from splinepy._splinepy import *
from splinepy._spline import Spline

class Bezier(Spline):

    def __init__(self, degrees=None, control_points=None):
        """
        Bezier (Spline).

        Parameters
        -----------
        degrees: (para_dim,) list-like
        control_points: (m, dim) list-like

        Returns
        --------
        None
        """
        super().__init__(
            degrees=degrees,
            control_points=control_points,
        )

    def _update_c(self,):
        """
        """
        pass
