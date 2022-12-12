"""
Multipatch Spline Configuration
"""

import numpy as np

from splinepy.utils.log import logging


class Multipatch():
    """
    System of patches to store information such as boundaries and
    connectivities
    """

    def __init__(
            self,
            splines=None,
            connectivity=None,
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

        if not isinstance(list_of_splines, list):
            if issubclass(type(list_of_splines), Spline):
                list_of_splines = list(list_of_splines)
            else:
                raise ValueError("Wrong input format")

        # Check if all entries in the list are splines
        spline_para_dim = list_of_splines[0].para_dim
        spline_dim = list_of_splines[0].dim
        for spline in list_of_splines:
            if not issubclass(type(spline), Spline):
                raise ValueError(
                        "Only Splinepy-Spline types are allowed as list "
                        "entries"
                )
            if not spline.dim == spline_dim:
                raise ValueError(
                        "Dimension mismatch between splines in list of splines"
                )
            if not spline.para_dim == spline_para_dim:
                raise ValueError(
                        "Parametric dimension mismatch between splines in list"
                        " of splines"
                )

        self._spline_list = list_of_splines

    @property
    def connectivity(self):
        """
        Determine the connectivity of a given multipatch system

        Returns its connectivity in the form as an array of size
        n_patches x n_sides_per_patch
        """
        if self._connectivity is None:
            raise ValueError(
                    "No connectivity available, "
                    "try calling determine_connectivity"
            )
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
                    "Connectivity must be stored in a numpy 2D array."
            )

        # One boundary for min and max for each parametric dimension
        n_boundaries_per_spline = self.splines[0].para_dim * 2
        if not (
                (con.shape[1] == n_boundaries_per_spline)
                and (con.shape[0] == len(self.splines))
        ):
            raise ValueError(
                    "Connectivity size mismatch. Expected size is "
                    "n_patch x n_boundaries"
            )

        # Assignment
        self._connectivity = con

    @property
    def boundaries(self):
        """
        Boundaries are stored within the connectivity array as negativ entries

        Negativ entries mean, that there is a boundary, the absolute value
        holds the boundary ID
        """
        # Get minimum boundary id
        max_BID = self.connectivity.min()

        boundary_list = []
        for i_bid in range(-1, max_BID, -1):
            logging.debug(f"Extracting boundary with ID {abs(i_bid)}")
            boundary_list.append(self.connectivity)
            logging.debug(
                    f"Found {boundary_list[-1][1].size} boundary "
                    f"elements on boundary {abs(i_bid)}"
            )

        return boundary_list

    def set_boundary(self, spline_ids, boundary_faces, boundary_id=None):
        """
        Adds a boundary to a specific set of spline and spline-

        Parameters
        ----------
        spline_ids : np.ndarray
          array of spline-ids that are assigned a boundary
        boundary_faces : np.ndarray
          array of corresponding faces to be assigned a boundary
        boundary_id : int
          boundary_id to be assigned. If not chosen set to new lowest value

        Returns
        -------
        None
        """
        # Get minimum boundary id
        max_BID = self.connectivity.min()

        if boundary_id is None:
            new_BID = max_BID - 1
            logging.debug(f"Creating new boundary with ID {abs(new_BID)}")
        else:
            # Make sure its negative
            new_BID = -abs(int(boundary_id))
            if new_BID < max_BID:
                logging.debug(f"Creating new boundary with ID {abs(new_BID)}")
            else:
                logging.debug(
                        "Adding new boundary elements to existing "
                        f"boundary new boundary with ID {abs(new_BID)}"
                )

        try:
            old_indices = self.connectivity[spline_ids, boundary_faces]

            # Check if all old indices are negativ
            if (old_indices < 0).all():
                self.connectivity[spline_ids, boundary_faces] = new_BID
            else:
                raise ValueError(
                        "One or more of the assigned boundary elements do not"
                        " ly on the patch surface, please check topology"
                )
        except BaseException:
            raise ValueError(
                    "spline_ids and boundary_faces need to be one-dimensional."
                    "\nIf this error proceeds please check if connectivity "
                    "exists by calling, \nprint(<>.connectivity)"
            )

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

    @property
    def spline_face_centers(self):
        """
        This property function is a placeholder for once i figure out how to
        use tracked-arrays

        talk to at @j042
        """
        # If spline list is empty will throw excetion
        return np.vstack([s.evaluate_face_centers for s in self.splines])

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

        # Using the property instead of the the member, all necessery
        # checks will be performed
        self.connectivity = get_connectivity_from_face_centers(
                self.spline_face_centers, settings.TOLERANCE, self.para_dim
        )

        logging.debug(
                "Successfully provided new connectivity using uff algorithm"
        )

    def add_boundary_with_function(self, function, only_unassigned=False):
        """
        Uses all faces, that were identified as boundaries in the
        connectivity array and checks if they fit a boundary function

        See boundary
        """
        # Get minimum boundary id
        max_BID = self.connectivity.min()

        # retrieve all boundary elements
        if only_unassigned:
            boundary_ids = (self.connectivity < 0)
        else:
            boundary_ids = self.connectivity == -1

        # Check if there is a boundary
        if not boundary_ids.any():
            raise ValueError("No boundary elements could be identified")

        # Check all face centers
        relevant_face_centers = self.spline_face_centers[
                boundary_ids.flatten(), :]

        # Cols and Rows
        row_ids, col_ids = np.where(boundary_ids)

        # Function is applied to all face centers
        try:
            new_boundary_bools = function(relevant_face_centers)
        except BaseException:
            ValueError(
                    "Function is not applicable to array. Function layout must"
                    " f(array<n_face_points,dim>) -> bools<n_face_points>"
            )

        # Assign new boundary ID to connectivity array
        self.connectivity[row_ids[new_boundary_bools],
                          col_ids[new_boundary_bools]] = max_BID - 1
