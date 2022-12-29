"""
Multipatch Spline Configuration
"""

import numpy as np

from splinepy._base import SplinepyBase
from splinepy.settings import TOLERANCE
from splinepy.utils.data import enforce_contiguous


class Multipatch(SplinepyBase):
    """
    System of patches to store information such as boundaries and
    interfaces
    """

    def __init__(
        self,
        splines=None,
        interfaces=None,
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
        self._interfaces = None
        self._spline_list = None
        self._boundaries = None

        self._logd("Instantiated Multipatch object")

        # Set properties
        if splines is not None:
            self.splines = splines

        if interfaces is not None:
            self.interfaces = interfaces

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
                    "Only Splinepy-Spline types are allowed as list " "entries"
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
    def interfaces(self):
        """
        Determine the interfaces of a given multipatch system

        Returns its interfaces in the form as an array of size
        n_patches x n_sides_per_patch
        """
        if self.splines is None:
            raise ValueError("Connectivity set for unknown list of splines.")

        if self._interfaces is None:
            self._logd("No interfaces available, calculating on the fly")
            self.determine_interfaces()
        return self._interfaces

    @interfaces.setter
    def interfaces(self, con):
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
        con = enforce_contiguous(con, np.int32)

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
        self._interfaces = con

    @property
    def boundaries(self):
        """
        Boundaries are stored within the interfaces array as negativ entries

        Negativ entries mean, that there is a boundary, the absolute value
        holds the boundary ID
        """
        # Get minimum boundary id
        max_BID = self.interfaces.min()

        boundary_list = []
        for i_bid in range(-1, max_BID, -1):
            self._logd(f"Extracting boundary with ID {abs(i_bid)}")
            boundary_list.append(np.where(self.interfaces == i_bid))
            self._logd(
                f"Found {boundary_list[-1][1].size} boundary "
                f"elements on boundary {abs(i_bid)}"
            )

        return boundary_list

    def set_boundary(
        self,
        spline_ids,
        boundary_faces,
        boundary_id=None,
    ):
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
        max_BID = self.interfaces.min()

        if boundary_id is None:
            new_BID = max_BID - 1
            self._logd(f"Creating new boundary with ID {abs(new_BID)}")
        else:
            # Make sure its negative
            new_BID = -abs(int(boundary_id))
            if new_BID < max_BID:
                self._logd(f"Creating new boundary with ID {abs(new_BID)}")
            else:
                self._logd(
                    "Adding new boundary elements to existing "
                    f"boundary new boundary with ID {abs(new_BID)}"
                )

        try:
            old_indices = self.interfaces[spline_ids, boundary_faces]

            # Check if all old indices are negativ
            if (old_indices < 0).all():
                self.interfaces[spline_ids, boundary_faces] = new_BID
            else:
                raise ValueError(
                    "One or more of the assigned boundary elements do not"
                    " ly on the patch surface, please check topology"
                )
        except BaseException:
            raise ValueError(
                "spline_ids and boundary_faces need to be one-dimensional."
                "\nIf this error proceeds please check if interfaces "
                "exists by calling, \nprint(<>.interfaces)"
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
    def spline_boundary_centers(self):
        """
        This property function is a placeholder for once i figure out how to
        use tracked-arrays

        talk to at @j042
        """
        # If spline list is empty will throw excetion
        return np.vstack([s.boundary_centers for s in self.splines])

    def determine_interfaces(self):
        """
        Retrieve interfaces info

        Stores new informaton as interfaces

        Parameters
        ----------
        None

        Returns
        -------
        interfaces : array-like (n_patch x n_boundaries)
        """
        from splinepy import settings
        from splinepy.splinepy_core import interfaces_from_boundary_centers

        # Using the setter instead of the the member, all necessery
        # checks will be performed
        self.interfaces = interfaces_from_boundary_centers(
            self.spline_boundary_centers, settings.TOLERANCE, self.para_dim
        )

        self._logd("Successfully provided new interfaces using uff algorithm")

        return self.interfaces

    def add_boundary_with_function(
        self,
        function,
        only_unassigned=False,
        boundary_id=None,
    ):
        """
        Uses all faces, that were identified as boundaries in the
        interfaces array and checks if they fit a boundary function

        Parameters
        ----------
        function : Callable
          Function called on every boundary center point to check if it is on
          the boundary, returns bool-type
        only_unassigned : bool
          Uses only previously unassigned boundaries
          (i.e. on boundary 1)
        boundary_id : int
          boundary_id to be assigned. If not chosen set to new lowest value

        Returns
        -------
        None
        """
        # Get minimum boundary id
        max_BID = self.interfaces.min()

        if boundary_id is None:
            new_BID = max_BID - 1
            self._logd(f"Creating new boundary with ID {abs(new_BID)}")
        else:
            # Make sure its negative
            new_BID = -abs(int(boundary_id))
            if new_BID < max_BID:
                self._logd(f"Creating new boundary with ID {abs(new_BID)}")
            else:
                self._logd(
                    "Adding new boundary elements to existing "
                    f"boundary new boundary with ID {abs(new_BID)}"
                )

        # retrieve all boundary elements
        if only_unassigned:
            boundary_ids = self.interfaces < 0
        else:
            boundary_ids = self.interfaces == -1

        # Check if there is a boundary
        if not boundary_ids.any():
            raise ValueError("No boundary elements could be identified")

        # Check all face centers
        relevant_boundary_centers = self.spline_boundary_centers[
            boundary_ids.flatten(), :
        ]

        # Cols and Rows
        row_ids, col_ids = np.where(boundary_ids)

        # Function is applied to all face centers
        try:
            new_boundary_bools = function(relevant_boundary_centers)
        except BaseException:
            ValueError(
                "Function is not applicable to array. Function layout must"
                " f(array<n_face_points,dim>) -> bools<n_face_points>"
            )

        # Assign new boundary ID to interfaces array
        self.interfaces[
            row_ids[new_boundary_bools], col_ids[new_boundary_bools]
        ] = new_BID

    def add_boundary_from_seed(self, seed_position, tolerance=None):
        """WIP
        Starting from a seed position, the splines are propagated until they
        reach a kink (no g1 continuity). This uses the spline boundary
        information and determines the interface information.

        Parameters
        ----------
        seed_position : array-like (n_dim)
          Seed position from where the propagation starts
        tolerance : double
          Maximum angle (rad) between two normal vectors on boundary which is
          accepted as g1
        """
        if tolerance is None:
            tolerance = TOLERANCE

        pass
