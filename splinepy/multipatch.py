"""
Multipatch Spline Configuration
"""

import numpy as np

from utils.log import logging

class Multipatch():
    """
    System of patches to store information such as boundaries and
    connectivities
    """

    def __init__(self,
      splines=None,
      connectivity=None
    ):
        """
        Multipatch

        Parameters
        ----------
        splines : list-like
          List of splines to store as multipatch

        Returns
        -------
        None
        """
        # Init values
        self._connectivity = None
        self._splines = None
        self._boundaries = None

        logging.debug("Instantiated Multipatch object")

        # Set properties
        if splines is not None:
            self.splines = splines

        if connectivity is not None:
            self.connectivity = connectivity



    @property
    def splines(self):
      """
      List of splines in splinepy format

      Returns a list of splines that are stored in the multipatch system
      """
      return self._spline_list

    @splines.setter
    def splines(self, list_of_splines):
        from splinepy.spline import Spline

        if isinstance(list_of_splines, list):
            if issubclass(type(spline), Spline):
               list_of_splines= list(list_of_splines)
            else:
                raise ValueError("Wrong input format")

        # Check if all entries in the list are splines
        spline_para_dim = list_of_splines[0].paradim
        spline_dim = list_of_splines[0].dim
        for spline in list_of_splines:
            if not issubclass(type(spline), Spline):
               raise ValueError(
                   "Only Splinepy-Spline types are allowed as list entries")
            if spline.dim == spline_dim:
               raise ValueError(
                   "Dimension mismatch between splines in list of splines")
            if spline.paradim == spline_para_dim:
               raise ValueError(
                   "Parametric dimension mismatch between splines in list of splines")

        self._spline_list = list_of_splines

    @property
    def connectivity(self):
        """
        Determine the connectivity of a given multipatch system

        Returns its connectivity in the form as an array of size n_patches x n_sides_per_patch
        """
        return self._connectivity

    @connectivity.setter
    def connectivity(self, con):
        """
        No Actual implementation only some checks

        TBD
        """
        if self.splines is None:
          raise ValueError("Connectivity set for unknown list of splines.")

        # Check types
        if not isinstance(con, np.ndarray) or (not len(con.shape)):
            raise ValueError(
                "Connectivity must be stored in a numpy 2D array.")
        # One boundary for min and max for each parametric dimension
        n_boundaries_per_spline = self.splines[0].paradim * 2
        if (con.shape[1] == n_boundaries_per_spline) or (
            con.shape[0] == len(self.splines)):
            raise ValueError(
                "Connectivity size mismatch. Expected size is n_patch x n_boundaries")

        # Assignment
        self._connectivity = con

    @property
    def boundaries(self):
        """
        IMO, boundaries should be stored as negative values in the connectivity
        """
        pass

    @boundaries.setter
    def boundaries(self, bids):
        pass

    @property
    def para_dim(self):
        """
        Parametric dimension of the splineset
        """
        if self.splines is None:
            raise ValueError("No Splines provided")
        return self.splines[0].para_dim

    @property
    def dim(self):
        """
        Physical dimension of the splineset
        """
        if self.splines is None:
            raise ValueError("No Splines provided")
        return self.splines[0].dim


    def determine_connectivity(self):
        """
        Uses MFEM export routine to retrieve connectivity info

        Stores new informaton as connectivity

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        from splinepy.splinepy_core import get_connectivity_from_face_centers
        from splinepy import settings

        if self.splines is None:
            raise ValueError("No Splines provided")

        face_centers = np.vstack(
            [s.evaluate_face_centers for s in self.splines])

        self.connectivity = get_connectivity_from_face_centers(
            face_centers, settings.TOLERANCE, self.para_dim)

    def add_boundary_with_function(function):
        """
        Uses all faces, that were identified as boundaries in the connectivity array and checks if they fit a boundary function

        See boundary
        """
        pass








