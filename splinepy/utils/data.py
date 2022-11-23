"""splinepy/utils/data.py.

Helps helpee to manage data. Some useful data structures.
"""

import numpy as np


class TrackedArray(np.ndarray):
    """Taken from nice implementations of `trimesh` (see LICENSE).
    `https://github.com/mikedh/trimesh/blob/main/trimesh/caching.py`. Minor
    adaption, since we don't have hashing functionalities.

    All the inplace functions will set modified flag and if some operations
    has potential to cause un-trackable behavior, writeable flags will be set
    to False.

    Note, if you really really want, it is possible to change the tracked
    array without setting modified flag.
    """

    __slots__ = ("_modified", "_source")

    def __array_finalize__(self, obj):
        """Sets default flags for any arrays that maybe generated based on
        tracked array."""
        self._modified = True
        self._source = int(0)

        if isinstance(obj, type(self)):
            if isinstance(obj._source, int):
                self._source = obj
            else:
                self._source = obj._source

    @property
    def mutable(self):
        return self.flags['WRITEABLE']

    @mutable.setter
    def mutable(self, value):
        self.flags.writeable = value

    def _set_modified(self):
        """set modified flags to itself and to the source."""
        self._modified = True
        if isinstance(self._source, type(self)):
            self._source._modified = True

    def copy(self, *args, **kwargs):
        """copy gives np.ndarray.

        no more tracking.
        """
        return np.array(self, copy=True)

    def view(self, *args, **kwargs):
        """Set writeable flags to False for the view."""
        v = super(self.__class__, self).view(*args, **kwargs)
        v.flags.writeable = False
        return v

    def __iadd__(self, *args, **kwargs):
        self._set_modified()
        return super(self.__class__, self).__iadd__(*args, **kwargs)

    def __isub__(self, *args, **kwargs):
        self._set_modified()
        return super(self.__class__, self).__isub__(*args, **kwargs)

    def __imul__(self, *args, **kwargs):
        self._set_modified()
        return super(self.__class__, self).__imul__(*args, **kwargs)

    def __idiv__(self, *args, **kwargs):
        self._set_modified()
        return super(self.__class__, self).__idiv__(*args, **kwargs)

    def __itruediv__(self, *args, **kwargs):
        self._set_modified()
        return super(self.__class__, self).__itruediv__(*args, **kwargs)

    def __imatmul__(self, *args, **kwargs):
        self._set_modified()
        return super(self.__class__, self).__imatmul__(*args, **kwargs)

    def __ipow__(self, *args, **kwargs):
        self._set_modified()
        return super(self.__class__, self).__ipow__(*args, **kwargs)

    def __imod__(self, *args, **kwargs):
        self._set_modified()
        return super(self.__class__, self).__imod__(*args, **kwargs)

    def __ifloordiv__(self, *args, **kwargs):
        self._set_modified()
        return super(self.__class__, self).__ifloordiv__(*args, **kwargs)

    def __ilshift__(self, *args, **kwargs):
        self._set_modified()
        return super(self.__class__, self).__ilshift__(*args, **kwargs)

    def __irshift__(self, *args, **kwargs):
        self._set_modified()
        return super(self.__class__, self).__irshift__(*args, **kwargs)

    def __iand__(self, *args, **kwargs):
        self._set_modified()
        return super(self.__class__, self).__iand__(*args, **kwargs)

    def __ixor__(self, *args, **kwargs):
        self._set_modified()
        return super(self.__class__, self).__ixor__(*args, **kwargs)

    def __ior__(self, *args, **kwargs):
        self._set_modified()
        return super(self.__class__, self).__ior__(*args, **kwargs)

    def __setitem__(self, *args, **kwargs):
        self._set_modified()
        super(self.__class__, self).__setitem__(*args, **kwargs)

    def __setslice__(self, *args, **kwargs):
        self._set_modified()
        super(self.__class__, self).__setslice__(*args, **kwargs)

    def __getslice__(self, *args, **kwargs):
        """
        return slices I am pretty sure np.ndarray does not have __*slice__
        """
        self._set_modified()
        slices = super(self.__class__, self).__getitem__(*args, **kwargs)
        if isinstance(slices, np.ndarray):
            slices.flags.writeable = False
        return slices


def make_tracked_array(array, dtype=None, copy=True):
    """Taken from nice implementations of `trimesh` (see LICENSE).
    `https://github.com/mikedh/trimesh/blob/main/trimesh/caching.py`.

    ``Properly subclass a numpy ndarray to track changes.
    Avoids some pitfalls of subclassing by forcing contiguous
    arrays and does a view into a TrackedArray.``

    Factory-like wrapper function for TrackedArray.

    Parameters
    ------------
    array: array-like object
      To be turned into a TrackedArray
    dtype: np.dtype
      Which dtype to use for the array
    copy: bool
      Default is True. copy if True.

    Returns
    ------------
    tracked : TrackedArray
      Contains input array data
    """
    # if someone passed us None, just create an empty array
    if array is None:
        array = []
    # make sure it is contiguous then view it as our subclass
    tracked = np.ascontiguousarray(array, dtype=dtype)
    if copy:
        tracked = tracked.copy().view(TrackedArray)
    else:
        tracked = tracked.view(TrackedArray)

    # should always be contiguous here
    assert tracked.flags['C_CONTIGUOUS']

    return tracked


def is_modified(array):
    """
    Returns true if:
      1. array is list of tracked arrays and any of them is modified.
      2. array is modified.

    Parameters
    ----------
    array: spline or array-like

    Returns
    --------
    modified: bool
    """
    if isinstance(array, list):
        modified = False
        for a in array:
            # recusively check
            modified |= is_modified(a)

        return modified

    elif isinstance(array, TrackedArray):
        return array._modified

    else:
        raise TypeError(f"{array} is not trackable.")


def without_none_values(dict_):
    """
    Returns a new dict without None value items.

    Parameters
    ----------
    dict_: dict

    Returns
    -------
    dict_without_none: dict
    """
    dict_without_none = dict()
    for key, value in dict_.items():
        if value is not None:
            dict_without_none[key] = value

    return dict_without_none
