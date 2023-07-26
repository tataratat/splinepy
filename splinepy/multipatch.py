"""
Multipatch Spline Configuration
"""

import numpy as np

from splinepy import settings
from splinepy._base import SplinepyBase
from splinepy.helpme import visualize
from splinepy.helpme.extract import Extractor
from splinepy.spline import _default_if_none, _get_helper
from splinepy.splinepy_core import PyMultiPatch, boundaries_from_continuity
from splinepy.utils.data import MultipatchData


class Multipatch(SplinepyBase, PyMultiPatch):
    """
    System of patches to store information such as boundaries and
    interfaces
    """

    __slots__ = (
        "_extractor",
        "_show_options",
        "_spline_data",
    )

    __show_option__ = visualize.SplineShowOption

    def __init__(self, splines=None, interfaces=None, *, spline=None):
        """
        Multipatch

        Parameters
        ----------
        splines : list-like
          List of splines to store as multipatch
        interfaces : array-like
          Defines the connectivity inbetween patches
        spline : PyMultiPatch
          keyword only argument, implemented to support to_derived() interface.
          calls move constructor.

        Returns
        -------
        None
        """
        # Init values
        if spline is not None:
            if not isinstance(spline, PyMultiPatch):
                raise TypeError("spline must be PyMultiPatch.")
            super().__init__(spline)
        elif splines is not None:
            super().__init__(splines, settings.NTHREADS, False)
        else:
            super().__init__()

        self._logd("Initialized Multipatch object")

        if interfaces is not None:
            self.interfaces = interfaces

    @property
    def splines(self):
        """
        List of splines in splinepy format

        Returns
        -------
        splines : list
         list of splines that are stored in the multipatch system


        """
        return self.patches

    @splines.setter
    def splines(self, list_of_splines):
        self.patches = list_of_splines

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
        # empty list as input will compute interfaces based on input
        return super().interfaces([])

    @interfaces.setter
    def interfaces(self, con):
        """super() checks validity of input"""
        # Assignment
        super().interfaces(con)

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

    def boundary_patches(self, nthreads=None):
        """Extract all boundary patches of a given Multipatch system as splines

        Parameters
        ----------
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

        # apply nthreads
        previous_nthreads = self.n_default_threads
        self.n_default_threads = nthreads

        b_patches = super().boundary_multi_patch()

        self.n_default_threads = previous_nthreads

        return b_patches

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
        return self.sub_patch_centers()

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

        # save previous default tol
        previous_tol = self.tolerance

        # set given tolerance and compute interfaces
        self.tolerance = tolerance
        interfaces = self.interfaces

        self._logd("Successfully provided new interfaces using uff algorithm")

        # set default tol back
        self.tolerance = previous_tol

        return interfaces

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

    @property
    def spline_data(self):
        """
        Spline data helper for splines. @todo (does not do anything at the
        moment)

        Parameters
        ----------
        None

        Returns
        -------
        spline_data: MultipatchData
        """
        return _get_helper(self, "_spline_data", MultipatchData)

    @property
    def show_options(self):
        """
        Show option manager for MultiPatches.

        Parameters
        ----------
        None

        Returns
        -------
        show_options: SplineShowOption
        """
        return _get_helper(self, "_show_options", self.__show_option__)

    def showable(self, **kwargs):
        """Equivalent to

        .. code-block:: python

           splinepy.helpme.visualize.show(
               spline, return_showable=True, **kwargs
            )

        Parameters
        ----------
        kwargs: kwargs
          see splinepy.helpme.visualize.show

        Returns
        -------
        spline_showable: dict
        """
        return visualize.show(self, return_showable=True, **kwargs)

    def show(self, **kwargs):
        return visualize.show(self, **kwargs)

    def boundaries_from_continuity(
        self,
        nthreads=None,
        tolerance=None,
    ):
        """
        Starting from a seed position, the splines are propagated until they
        reach a kink (no g1 continuity). This uses the spline boundary
        information and determines the interface information.

        Parameters
        ----------
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
        b_patches = self.boundary_patches(nthreads=nthreads)

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

    def add_fields(
        self,
        fields,
        check_compliance=False,
        check_conformity=False,
        nthreads=None,
    ):
        """
        Add fields using lists of splines

        Parameters
        ----------
        fields : list
          Any number of list of splines to represent n-dimensional field with
          equal parametric dimensionality
        check_compliance : bool (False)
          Check if field list is admissible, by comparing the parametric
          dimensionality of the field entries with the spline list, and compare
          patch sizes
        check_conformity : bool (False)
          Check for conformity between patches and fields by comparing degrees
          and control-mesh-resolutions
        nthreads : int

        Returns
        -------
        None
        """
        if nthreads is None:
            nthreads = settings.NTHREADS

        super().add_fields(
            fields,
            check_compliance,
            check_compliance,
            check_conformity,
            check_conformity,
            nthreads,
        )

    @property
    def fields(self):
        """Save fields as individul splines on patches

        Parameters
        ----------
        None

        Returns
        -------
        fields : list<list>
          List of all field representation in the form of list of splines
        """
        return super().fields()

    def sample(self, resolutions, nthreads=None):
        """
        Uniformly sample along each parametric dimensions from spline.

        Parameters
        -----------
        resolutions: int
        nthreads: int

        Returns
        --------
        results: (math.product(resolutions), dim) np.ndarray
        """

        if not isinstance(resolutions, int) and hasattr(
            resolutions, "__getitem__"
        ):
            self._logd(
                "sample() only supports uniform sample. Taking first entry."
            )
            # for now, just take the first elem
            resolutions = int(resolutions[0])

        self._logd(
            f"Sampling {resolutions ** self.para_dim} points from spline."
        )
        return super().sample(
            resolutions,
            nthreads=_default_if_none(nthreads, settings.NTHREADS),
            same_parametric_bounds=False,
        )

    def evaluate(self, queries, nthreads=None):
        """
        Evaluate each individual spline at specific parametric positions. To be
        used with caution, as there is no check if the queries are within the
        parametric bounds if settings.CHECK_BOUNDS is set to false.

        Parameters
        -----------
        queries: (n, para_dim) array-like
        nthreads: int

        Returns
        --------
        results: (math.product(resolutions), dim) np.ndarray
        """

        return super().evaluate(
            queries,
            nthreads=_default_if_none(nthreads, settings.NTHREADS),
        )

    @property
    def extract(self):
        """Return Extractor object to provide extract functionality for
        multiple splines

        Parameters
        -----------
        None

        Returns
        --------
        extractor: Extractor
        """
        return _get_helper(self, "_extractor", Extractor)
