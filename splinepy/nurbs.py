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
