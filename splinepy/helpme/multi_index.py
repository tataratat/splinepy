"""
Helps you find control point ids.
"""
import numpy as np


class MultiIndex:
    """
    Performs same task as `np.ravel_multi_index`. However, you can use it with
    __getitem__ like queries.
    """
    __slots__ = ("_raveled_indices", "_ndim")
    _ellipsis_type = type(Ellipsis)

    def __init__(self, grid_resolutions):
        """
        Initialize using grid's resolution per dimension. For splines, use
        `Spline.control_mesh_resolutions`.

        Parameters
        ----------
        grid_resolutions: array-like
        """
        if not isinstance(grid_resolutions, (tuple, list, np.ndarray)):
            raise TypeError("grid_resolutions must be array-like type, not" f"{type(grid_resolutions)}")

        # create raveled indices
        raveled = np.arange(np.product(grid_resolutions), dtype="int32")

        # to allow general __getitem__ like query, turn them into
        # grid's shape
        self._raveled_indices = raveled.reshape(*grid_resolutions[::-1])
        self._ndim = self._raveled_indices.ndim


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
        # only interested in __iter__, not __getitem__
        if not hasattr(args, "__iter__"):
            args = [args]
        else:
            args = list(args)

        # Ellipsis? don't do anything
        # (Ellipsis not in args) raises ValueError during comparison
        # if one of them is an array. so, loop and check if any of the entry
        # are Ellipsis
        fill_slices = True
        for a in args:
            if isinstance(a, self._ellipsis_type):
                fill_slices = False
                break

        # no ellipsis? fill in with slices
        if fill_slices:
            s = slice(None)
            for _ in range(len(args), self._ndim):
                args.append(s)

        return self._raveled_indices.__getitem__(tuple(args[::-1])).ravel()


