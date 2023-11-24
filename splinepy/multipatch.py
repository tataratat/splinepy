"""
Multipatch Spline Configuration
"""

import numpy as _np

from splinepy import settings as _settings
from splinepy._base import SplinepyBase as _SplinepyBase
from splinepy.helpme import visualize as _visualize
from splinepy.helpme.extract import Extractor as _Extractor
from splinepy.spline import _default_if_none, _get_helper
from splinepy.splinepy_core import PyMultipatch as _PyMultipatch
from splinepy.splinepy_core import (
    boundaries_from_continuity as _boundaries_from_continuity,
)
from splinepy.utils.data import MultipatchData as _MultipatchData


class Multipatch(_SplinepyBase, _PyMultipatch):
    """
    System of patches to store information such as boundaries and
    interfaces
    """

    __slots__ = (
        "_extractor",
        "_show_options",
        "_spline_data",
    )

    __show_option__ = _visualize.SplineShowOption

    def __init__(self, splines=None, interfaces=None, *, spline=None):
        """
        Multipatch

        Parameters
        ----------
        splines : list-like
          List of splines to store as multipatch
        interfaces : array-like
          Defines the connectivity in between patches
        spline : PyMultipatch
          keyword only argument, implemented to support to_derived() interface.
          calls move constructor.

        Returns
        -------
        None
        """
        # Init values
        if spline is not None:
            if not isinstance(spline, _PyMultipatch):
                raise TypeError("spline must be PyMultipatch.")
            super().__init__(spline)
        elif splines is not None:
            super().__init__(splines, _settings.NTHREADS, False)
        else:
            super().__init__()

        self._logd("Initialized Multipatch object")

        if interfaces is not None:
            self.interfaces = interfaces

    @property
    def patches(self):
        """
        List of splines in splinepy format

        Returns
        -------
        patches : list
         list of splines that are stored in the multipatch system


        """
        return super().patches

    @patches.setter
    def patches(self, list_of_splines):
        _PyMultipatch.patches.fset(self, list_of_splines)

    @property
    def interfaces(self):
        """
        Determine the interfaces of a given multipatch system

        Returns its interfaces in the form as an array of size
        n_patches x n_sides_per_patch

        The interfaces are stored as indices. If the index is negative, the ID
        is referring to a boundary (e.g. for boundary conditions). If the
        number is positive, the ID refers to the global face-id, which is
        composed of the element ID and the local face ID :
        global_face = local_face + number_of_faces_per_element * element_id

        Example:

        .. code-block ::

          O -- 3 -- O O -- 3 -- O
          |         | |         |
          0    0    1 0    1    1
          |         | |         |
          O -- 2 -- O O -- 2 -- O

        would be stored as:

        .. code-block ::

          [
            [-1, 4, -1, -1],
            [1, -1, -1, -1]
          ]

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
        Boundaries are stored within the interfaces array as negative entries

        Negative entries mean, that there is a boundary, the absolute value
        holds the boundary ID

        Returns
        -------
        boundaries : list<np.ndarray>
          list of 2D arrays, representing the patch and face ID
        """
        # Get minimum boundary id - boundary ids are stored as negative values
        max_BID = self.interfaces.min()

        boundary_list = []
        for i_bid in range(-1, max_BID - 1, -1):
            self._logd(f"Extracting boundary with ID {abs(i_bid)}")
            boundary_list.append(_np.where(self.interfaces == i_bid))
            self._logd(
                f"Found {boundary_list[-1][1].size} boundary "
                f"elements on boundary {abs(i_bid)}"
            )

        return boundary_list

    def boundary_multipatch(self, nthreads=None):
        """Extract all boundary patches of a given Multipatch system as splines

        Parameters
        ----------
        nthreads : int
          Number of threads to be used for extraction, defaults to
          settings.NTHREADS

        Returns
        -------
        boundary_multipatch : Multipatch
          Embedded splines representing boundaries of system
        """
        if nthreads is None:
            nthreads = _settings.NTHREADS

        # apply nthreads
        previous_nthreads = self.n_default_threads
        self.n_default_threads = nthreads

        b_patches = super().boundary_multipatch()

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
        # Get minimum boundary id - boundary ids are stored as negative values
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

        except BaseException:
            raise ValueError(
                "spline_ids and boundary_faces need to be one-dimensional, "
                "array-like objects to access 2D interface numpy array"
            )

        # raise if any old indices include non-negative entries
        # as it indicates that it's not a boundary sub-patch
        if (old_indices > -1).any():
            raise ValueError(
                "One or more of the assigned boundary elements are not"
                " on the patch surface, please check topology"
            )

        # set boundary id
        self.interfaces[spline_ids, boundary_faces] = new_BID

    @property
    def para_dim(self):
        """
        Parametric dimension of the splineset

        Returns
        -------
        para_dim : int
          Parametric dimensionality of the multipatch system
        """
        return super().para_dim

    @property
    def dim(self):
        """
        Physical dimension of the splineset
        Returns
        -------
        dim : int
          Physical dimensionality of the multipatch system
        """
        return super().dim

    @property
    def sub_patch_centers(self):
        """
        Evaluated centers of the individual patches, used to determine
        connectivity and identify boundaries based on position

        Returns
        -------
        spline_centers : np.ndarray
          coordinates of the patch-boundary centers
        """
        return super().sub_patch_centers()

    def determine_interfaces(self, tolerance=None):
        """
        Retrieve interfaces info

        Stores new information as interfaces

        Parameters
        ----------
        tolerance : double
          normed distance between two neighboring points to be considered equal

        Returns
        -------
        interfaces : array-like (n_patch x n_boundaries)
        """
        if tolerance is None:
            tolerance = _settings.TOLERANCE

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
        # Get minimum boundary id - boundary ids are stored as negative values
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
            boundary_ids = _np.isin(-self.interfaces, _np.abs(mask))

        # Check if there is a boundary
        if not boundary_ids.any():
            self._logd(
                "No boundary elements could be identified that match "
                "requirements"
            )
            return None

        # Check all face centers
        relevant_boundary_centers = self.sub_patch_centers[
            boundary_ids.flatten(), :
        ]

        # Cols and Rows
        row_ids, col_ids = _np.where(boundary_ids)

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
        return _get_helper(self, "_spline_data", _MultipatchData)

    @property
    def show_options(self):
        """
        Show option manager for Multipatches.

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
        return _visualize.show(self, return_showable=True, **kwargs)

    def show(self, **kwargs):
        return _visualize.show(self, **kwargs)

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
            tolerance = _settings.TOLERANCE
        if nthreads is None:
            nthreads = _settings.NTHREADS
        b_patches = self.boundary_multipatch(nthreads=nthreads)

        # Pass information to c++ backend
        self._logd("Start propagation of information...")
        n_new_boundaries = _boundaries_from_continuity(
            b_patches.patches,
            b_patches.interfaces,
            self.interfaces,
            tolerance,
            nthreads,
        )
        self._logd(f"{n_new_boundaries} new boundaries were assigned")

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
        boundary_ids = _np.isin(-self.interfaces, _np.abs(mask))

        if not boundary_ids.any():
            self._logd(
                "No boundary elements could be identified that match "
                "requirements"
            )
            return None

        self.interfaces[boundary_ids] = _np.min(mask)

    def add_fields(
        self,
        fields,
        field_dim,
        check_compliance=True,
        check_conformity=True,
        nthreads=None,
    ):
        """
        Add fields using lists of splines

        Parameters
        ----------
        fields : list
          Any number of list of splines to represent n-dimensional field with
          equal parametric dimensionality
        field_dim: int
          Dimensionality of the individual fields
        check_compliance : bool (True)
          Check if field list is admissible, by comparing the parametric
          dimensionality of the field entries with the spline list, and compare
          patch sizes
        check_conformity : bool (True)
          Check for conformity between patches and fields by comparing degrees
          and control-mesh-resolutions
        nthreads : int

        Returns
        -------
        None
        """
        if nthreads is None:
            nthreads = _settings.NTHREADS

        super().add_fields(
            fields,
            field_dim,
            check_compliance,
            check_compliance,
            check_conformity,
            check_conformity,
            nthreads,
        )

    @property
    def fields(self):
        """Save fields as individual splines on patches

        Parameters
        ----------
        None

        Returns
        -------
        fields : list<Multipatch>
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
            nthreads=_default_if_none(nthreads, _settings.NTHREADS),
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
            nthreads=_default_if_none(nthreads, _settings.NTHREADS),
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
        return _get_helper(self, "_extractor", _Extractor)
