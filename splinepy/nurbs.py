import logging
import copy

import numpy as np

from splinepy import utils
from splinepy._splinepy import *
from splinepy._spline import Spline

class NURBS(Spline):

    # Required Properties
    _required_properties = [
        "degrees",
        "knot_vectors",
        "control_points",
        "weights",
    ]

    def __init__(self,
            degrees=None,
            knot_vectors=None,
            control_points=None,
            weights=None,
    ):
        """
        NURBS.

        Parameters
        -----------
        degrees: (para_dim,) list-like
        knot_vectors: (para_dim, n) list
        control_points: (m, dim) list-like
        weights: (m,) list-like

        Returns
        --------
        None
        """
        super().__init__(
            degrees=degrees,
            knot_vectors=knot_vectors,
            control_points=control_points,
            weights=weights,
        )


    def _update_c(self,):
        """
        Updates/Init cpp spline, if it is ready to be updated.
        Checks if all the entries are filled before updating.

        Parameters
        -----------
        None

        Returns
        --------
        None
        """
        if (
            self.dim is None
            or self.para_dim is None
            or self.degrees is None
            or self.knot_vectors is None
            or self.control_points is None
            or self.weights is None
        ):
            logging.debug(
                "Spline - Not enough information to update cpp spline. "
                "Skipping update."
            )
            if hasattr(self, "_c_spline"):
                delattr(self, "_c_spline")
            return None

        c_spline_class = f"NURBS{self.para_dim}P{self.dim}D"
        self._c_spline = eval(c_spline_class)(
            degrees=self.degrees,
            knot_vectors=self.knot_vectors,
            control_points=self.control_points,
            weights=self.weights,
        )

        logging.debug(f"Spline - Your spline is {self.whatami}.")

    def _update_p(self,):
        """
        Reads cpp spline and writes it here.
        Probably get an error if cpp isn't ready for this.

        Parameters
        -----------
        None

        Returns
        --------
        None
        """
        self._degrees = self._c_spline.degrees
        self._knot_vectors = self._c_spline.knot_vectors
        self._control_points = self._c_spline.control_points
        self._weights = self._c_spline.weights
        # do the setters work.
        self._para_dim = self._c_spline.para_dim
        self._dim = self._c_spline.dim
        logging.debug(
            "Spline - Updated python spline. CPP spline and python spline are "
            "now identical."
        )

    def copy(self,):
        """
        Returns freshly initialized Nurbs of self.

        Parameters
        -----------
        None

        Returns
        --------
        new_nurbs: `NURBS`
        """
        new_nurbs = NURBS(**self.todict())

        return new_nurbs

    def extract_bezier_patches(self):
        return self._c_spline.extract_bezier_patches()


    def todict(self, tolist=False):
        """
        Returns copy of degrees, knot_vectors, control_points, weights in dict.

        Parameters
        -----------
        tolist : bool
          Default is False. Convert numpy properties into lists

        Returns
        --------
        dict_spline: dict
          Keys are {degrees, knot_vectors, control_points, weights}.
        """
        if tolist:
            return dict(
                degrees=self.degrees.tolist(),
                knot_vectors=self.knot_vectors.tolist(),
                control_points=self.control_points.tolist(),
                weights=self.weights.tolist()
            )
        else:
            return dict(
                degrees=copy.deepcopy(self.degrees),
                knot_vectors=copy.deepcopy(self.knot_vectors),
                control_points=copy.deepcopy(self.control_points),
                weights=copy.deepcopy(self.weights)
            )
