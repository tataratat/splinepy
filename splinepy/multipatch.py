"""
Multipatch Spline Configuration
"""

import numpy as np

from splinepy import settings
from splinepy._base import SplinepyBase
from splinepy.spline import Spline
from splinepy.splinepy_core import (
    boundaries_from_continuity,
    boundary_centers,
    extract_all_boundary_splines,
    interfaces_from_boundary_centers,
)
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
        as_boundary=False,
    ):
        """
        Multipatch

        Parameters
        ----------
        splines : list-like
          List of splines to store as multipatch
        interfaces : array-like
          Defines the connectivity inbetween patches
        as_boundary : bool
          Multipatch is a boundary object of a higher dimensional geometry. If
          set to true, additional checks are performed on the interfaces,
          requiring strict interconnectivity between all patches

        Returns
        -------
        None
        """
        # Init values
        self._init_members()

        self._logd("Instantiated Multipatch object")

        # Set properties
        self._as_boundary = as_boundary

        if splines is not None:
            self.splines = splines

        if interfaces is not None:
            self.interfaces = interfaces

    def _init_members(self):
        """Defaults all relevant members to None

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        self._boundaries = None
        self._spline_list = None
        self._interfaces = None
        self._boundary_splines = None
        self._spline_boundary_centers = None

    @property
    def splines(self):
        """
        List of splines in splinepy format

        Returns
        -------
        splines : list
         list of splines that are stored in the multipatch system


        """
        return self._spline_list

    @splines.setter
    def splines(self, list_of_splines):
        # We need to start from scratch if new splines are added
        self._init_members()

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
        # Updating the spline set will erase all precalculated data
        self._boundaries = None

        self._spline_list = list_of_splines

    @property
    def interfaces(self):
        """
        Determine the interfaces of a given multipatch system

        Returns its interfaces in the form as an array of size
        n_patches x n_sides_per_patch

        Returns
        -------
        interfaces : np.ndarray
          inter-patch connectivitiy and boundaries
        """
        if self.splines is None:
            raise ValueError("Connectivity set for unknown list of splines.")

        if self._interfaces is None:
            self._logd("No interfaces available, calculating on the fly")
            self.determine_interfaces()
        return self._interfaces

    @interfaces.setter
    def interfaces(self, con):
        # Performes some checks prior to assigning the interface connectivity
        # to the multipatch object
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

        # If the multipatch is a boundary representation, than it must not
        # contain any negative entries as lower dimensional objects must always
        # be interconnected
        if self._as_boundary:
            if np.any(con < 0):
                raise ValueError(
                    "Interfaces are not interconnected, but interconnection is"
                    " required for boundary representations"
                )

        # Assignment
        self._interfaces = con

    @property
    def boundaries(self):
        """
        Boundaries are stored within the interfaces array as negativ entries

        Negativ entries mean, that there is a boundary, the absolute value
        holds the boundary ID

        Returns
        -------
        boundaries : list<np.ndarray>
          list of 2D arrays, representing the patch and face ID
        """
        # Get minimum boundary id - boundary ids are stored as negativ values
        max_BID = self.interfaces.min()

        boundary_list = []
        for i_bid in range(-1, max_BID - 1, -1):
            self._logd(f"Extracting boundary with ID {abs(i_bid)}")
            boundary_list.append(np.where(self.interfaces == i_bid))
            self._logd(
                f"Found {boundary_list[-1][1].size} boundary "
                f"elements on boundary {abs(i_bid)}"
            )

        return boundary_list

    def boundary_patches(self, use_saved=True, nthreads=None):
        """Extract all boundary patches of a given Multipatch system as splines

        Parameters
        ----------
        use_saved : bool
          Reuse an already computed patch system
        nthreads : int
          Number of threads to be used for extraction, defaults to
          settings.NTHREADS
        Returns
        -------
        boundary_patches : Multipatch
          Embedded splines representing boundaries of system
        """
        if nthreads is None:
            nthreads = settings.NTHREADS
        if (not use_saved) or (self._boundary_splines is None):
            self._logd("Determining boundary spline patches")
            patches = extract_all_boundary_splines(
                self.splines, self.interfaces, nthreads
            )
            self._boundary_splines = Multipatch(
                splines=[
                    settings.NAME_TO_TYPE[p.name](spline=p) for p in patches
                ],
                as_boundary=True,
            )

        return self._boundary_splines

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
        # Get minimum boundary id - boundary ids are stored as negativ values
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

        Returns
        -------
        para_dim : int
          Parametric dimensionality of the multipatch system
        """
        if self.splines is None:
            raise ValueError("No Splines provided")
        return self.splines[0].para_dim

    @property
    def dim(self):
        """
        Physical dimension of the splineset
        Returns
        -------
        dim : int
          Physical dimensionality of the multipatch system
        """
        if self.splines is None:
            raise ValueError("No Splines provided")

        return self.splines[0].dim

    @property
    def spline_boundary_centers(self):
        """
        Evaluated centers of the individual patches, used to determine
        connectivity and identify boundaries based on position

        Returns
        -------
        spline_centers : np.ndarray
          coordinates of the patch-boundary centers
        """
        # If spline list is empty will throw excetion
        if self._spline_boundary_centers is None:
            self._spline_boundary_centers = np.vstack(
                [boundary_centers(s) for s in self.splines]
            )

        return self._spline_boundary_centers

    def determine_interfaces(self, tolerance=None):
        """
        Retrieve interfaces info

        Stores new informaton as interfaces

        Parameters
        ----------
        tolerance : double
          normed distance between two neighboring points to be considered equal

        Returns
        -------
        interfaces : array-like (n_patch x n_boundaries)
        """
        if tolerance is None:
            tolerance = settings.TOLERANCE

        # Using the setter instead of the the member, all necessery
        # checks will be performed
        self.interfaces = interfaces_from_boundary_centers(
            self.spline_boundary_centers, tolerance, self.para_dim
        )

        self._logd("Successfully provided new interfaces using uff algorithm")

        return self.interfaces

    def boundary_from_function(
        self,
        function,
        mask=None,
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
        mask : array-like
          If assigned, takes only boundaries with matching ids
        boundary_id : int
          boundary_id to be assigned. If not chosen set to new lowest value

        Returns
        -------
        None
        """
        # Get minimum boundary id - boundary ids are stored as negativ values
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
        if mask is None:
            boundary_ids = self.interfaces < 0
        else:
            boundary_ids = np.isin(-self.interfaces, np.abs(mask))

        # Check if there is a boundary
        if not boundary_ids.any():
            self._logd(
                "No boundary elements could be identified that match "
                "requirements"
            )
            return None

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
            raise ValueError(
                "Function is not applicable to array. Function layout must"
                " f(array<n_face_points,dim>) -> bools<n_face_points>"
            )

        # Assign new boundary ID to interfaces array
        self.interfaces[
            row_ids[new_boundary_bools], col_ids[new_boundary_bools]
        ] = new_BID

    def boundaries_from_continuity(
        self,
        use_saved=True,
        nthreads=None,
        tolerance=None,
    ):
        """
        Starting from a seed position, the splines are propagated until they
        reach a kink (no g1 continuity). This uses the spline boundary
        information and determines the interface information.

        Parameters
        ----------
        use_saved : bool
          Allow function to reuse precomputed values
        nthreads : int
          Number of threads to be used to determine boundaries and aux values
        tolerance : double
          Tolerance between two normal vectors [cos(angle)] to be considered G1

        Returns
        -------
        None
        """
        if tolerance is None:
            tolerance = settings.TOLERANCE
        if nthreads is None:
            nthreads = settings.NTHREADS
        b_patches = self.boundary_patches(
            use_saved=use_saved, nthreads=nthreads
        )

        # Pass information to c++ backend
        self._logd("Start propagation of information...")
        n_new_boundaries = boundaries_from_continuity(
            b_patches.splines,
            b_patches.interfaces,
            self.interfaces,
            tolerance,
            nthreads,
        )
        self._logd(f"{n_new_boundaries} new boundaries were assigned")
        return None

    def combine_boundaries(self, mask=None):
        """
        Combines all boundaries that match an id from the mask to a single
        boundary

        Parameters
        ----------
        mask : array-like
          List of boundary ids to be reassigned a new combined ID

        Returns
        -------
        None
        """
        # retrieve all boundary elements
        boundary_ids = np.isin(-self.interfaces, np.abs(mask))

        if not boundary_ids.any():
            self._logd(
                "No boundary elements could be identified that match "
                "requirements"
            )
            return None

        self.interfaces[boundary_ids] = np.min(mask)
