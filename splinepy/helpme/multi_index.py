"""
Helps you find control point ids.
"""
import numpy as np


class MultiIndex:
    """
    Performs same task as `np.ravel_multi_index`. However, you can use it with
    __getitem__ like queries.
    """

    __slots__ = ("_raveled_indices",)

    def __init__(self, grid_resolutions):
        """
        Initialize using grid's resolution per dimension. For splines, use
        `Spline.control_mesh_resolutions`.

        Parameters
        ----------
        grid_resolutions: array-like
        """
        if not isinstance(grid_resolutions, (tuple, list, np.ndarray)):
            raise TypeError(
                "grid_resolutions must be array-like type, not"
                f"{type(grid_resolutions)}"
            )

        # create raveled indices
        raveled = np.arange(np.prod(grid_resolutions), dtype="int32")

        # to allow general __getitem__ like query, turn them into
        # grid's shape
        self._raveled_indices = raveled.reshape(*grid_resolutions, order="F")

    def __getitem__(self, args):
        """
        Retuns raveled indices.

        Parameters
        ----------
        args: int, slice, ellipsis, np.newaxis, array-like
          Any input np.ndarray.__getitem__ would take is valid.

        Returnes
        --------
        raveled_indices: np.ndarray
          raveled indices, which are in a raveled array.
        """
        return self._raveled_indices.__getitem__(args).T.ravel()
