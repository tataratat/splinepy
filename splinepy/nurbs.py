import logging

import numpy as np

from splinepy import utils
from splinepy._splinepy import *
from splinepy._spline import Spline
from splinepy.rational_bezier import RationalBezier


class NURBS(Spline):

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

    def extract_bezier_patches(self):
        """
        Extract all knot spans as Bezier patches to perform further operations
        such as compositions and multiplications

        Parameters
        ----------
        None

        Returns 
        -------
        extracted Beziers : list
        """
        # Extract bezier patches and create PyRationalBezier objects
        list_of_c_object_beziers = self._c_spline.extract_bezier_patches()

        # Transform into Rational Bezier Splinepy objects
        extracted_bezier_patches = []
        for c_object_spline in list_of_c_object_beziers:
            i_rational_bezier = RationalBezier()
            i_rational_bezier._c_spline = c_object_spline
            i_rational_bezier._update_p()
            extracted_bezier_patches.append(i_rational_bezier)

        return extracted_bezier_patches
