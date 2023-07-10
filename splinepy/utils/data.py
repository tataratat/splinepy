"""splinepy/utils/data.py.

Helps helpee to manage data. Some useful data structures.
"""
from itertools import accumulate, chain, repeat

import numpy as np
from gustaf.helpers.data import DataHolder, TrackedArray, make_tracked_array

from splinepy._base import SplinepyBase

__all__ = [
    "TrackedArray",
    "make_tracked_array",
    "DataHolder",
    "is_modified",
    "enforce_contiguous_values",
    "cartesian_product",
    "SplineDataAdaptor",
    "SplineData",
]


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


def enforce_contiguous(array, dtype=None, asarray=False):
    """
    If input is an instance / subclass of np.ndarray, this will check
    if they are configuous. If so, returns same object, else turns makes it
    contiguous and returns.

    Parameters
    ----------
    array: array-like
      data to be transformed
    dtype : type
      base data type (default float)
    asarray : bool
      also transform list or tuple to np.ndarray

    Returns
    -------
    contiguous_array: array-like
    """
    if isinstance(array, np.ndarray):
        if array.flags["C_CONTIGUOUS"] and (
            dtype is None or dtype is array.dtype
        ):
            return array
        return np.ascontiguousarray(array, dtype=dtype)

    if asarray and isinstance(array, (list, tuple)):
        return np.ascontiguousarray(array, dtype=dtype)

    return array


def enforce_contiguous_values(dict_):
    """
    Returns a new dict where values are contiguous. If the value is
    an array, enforces contiguous elements for one recursive level.
    This also removes None values.

    Parameters
    ----------
    dict_: dict

    Returns
    -------
    dict_with_contiguous: dict
    """
    dict_with_contiguous = dict()
    for key, value in dict_.items():
        # most likely will be asking about np.ndarray
        if isinstance(value, np.ndarray):
            value = enforce_contiguous(value)

        # or, list. only check if the first element is np.ndarray
        elif isinstance(value, list):
            if isinstance(value[0], np.ndarray):
                value = [enforce_contiguous(v) for v in value]

        # abandon Nones
        elif value is None:
            continue

        # set item
        dict_with_contiguous[key] = value

    return dict_with_contiguous


def cartesian_product(arrays, reverse=True):
    """
    Fast Cartesian product.
    Adapted from stackoverflow entry:
    https://stackoverflow.com/questions/11144513/
    cartesian-product-of-x-and-y-array-points-into-single-array-of-2d-points/
    49445693#49445693
    Answer from `Paul Panzer`.

    In splinepy, we use the ordering, where lowest dimension changes fastest.
    Thus, a reverse option.

    Parameters
    ----------
    arrays: tuple or list
      tuple or list of array-like
    reverse: bool
      Default is True. Reverses order of fastest iterating dimension.
      If False, highest indexed dimension iterates fastest.

    Returns
    -------
    cartesian: (n, len(arrays)) np.ndarray
    """
    n_arr = len(arrays)
    dtype = np.result_type(*arrays)

    # for reverse, use view of reverse ordered arrays
    if reverse:
        arrays = arrays[::-1]

    shape = *map(len, arrays), n_arr
    cartesian = np.empty(shape, dtype=dtype)

    # reverse order the view of output array, so that reverse
    # is again in original order, but just reversed fastest iterating dim
    if reverse:
        _cartesian = cartesian[..., ::-1]
    else:
        _cartesian = cartesian

    dim_reducing_views = (
        *accumulate(
            chain((_cartesian,), repeat(0, n_arr - 1)), np.ndarray.__getitem__
        ),
    )
    idx = slice(None), *repeat(None, n_arr - 1)
    for i in range(n_arr - 1, 0, -1):
        dim_reducing_views[i][..., i] = arrays[i][idx[: n_arr - i]]
        dim_reducing_views[i - 1][1:] = dim_reducing_views[i]
    _cartesian[..., 0] = arrays[0][idx]

    return cartesian.reshape(-1, n_arr)


class SplineDataAdaptor(SplinepyBase):
    """
    Prepares data to be presentable on spline. To support both
    scalar-data and vector-data, which are representable with colors and
    arrows respectively, this class will prepare data accordingly.
    """

    __slots__ = (
        "data",
        "function",
        "locations",
        "is_spline",
        "has_function",
        "has_locations",
        "has_evaluate",
        "arrow_data_only",
        "_user_created",
    )

    def __init__(self, data, locations=None, function=None):
        """ """
        # default
        self._user_created = True
        self.data = data
        self.is_spline = False
        self.has_function = False
        self.has_locations = False
        self.has_evaluate = False
        self.arrow_data_only = False

        # is spline we know?
        if "PySpline" in str(type(data).__mro__):
            self.is_spline = True

        # data has evaluate?
        if hasattr(data, "evaluate"):
            self.has_evaluate = callable(data.evaluate)

        # has function?
        if function is not None:
            self.has_function = True
            if not callable(function):
                raise ValueError("Given function isn't callable")
            self.function = function

        # locations? - keep this compatible with functions. maybe
        # we want to have some state dependent value at certain locations
        if locations is not None:
            # set what holds true
            self.has_locations = True
            self.arrow_data_only = True
            self.locations = np.asanyarray(locations)

            # if this is not a spline we know, it doesn't have a function,
            # it should:
            # -> `data.evaluate` is callable, or
            # -> len(data) == len(locations)
            if not self.is_spline and not self.has_function:
                len_matches = False
                if hasattr(data, "__len__"):
                    len_matches = len(locations) == len(data)
                if not (self.has_evaluate or len_matches):
                    raise ValueError(
                        "Data cannot be represented at specified locations."
                        "Requires one of the following requirements: "
                        "1) is a spline derived from splinepy's spline; "
                        "2) data has `data.evaluate()`; "
                        "3) length of the data and location should match."
                    )
            # location is specified, meaning we don't need sample()
            return None

        # can call sample or has a function?
        if not self.has_function and not self.is_spline:
            raise ValueError(
                "None spline data should at least have an accompanying "
                "function."
            )

    def as_vertex_data(self, resolutions=None, on=None):
        """
        Parameters
        ----------
        resolutions: list or tuple
        at: (n, d) array-like

        Returns
        -------
        vertex_data: (m, r) array-like
        """
        if resolutions is not None and on is not None:
            raise ValueError(
                "Please only specify either `resolutions` or `on`"
            )

        if self.has_locations and (resolutions is not None or on is not None):
            raise ValueError(
                "Location dependent data can't be evaluated with `resolutions`"
                " or `at`."
            )

        # if resolutions is specified, this is not a location query
        if resolutions is not None:
            if self.has_function:
                return self.function(self.data, resolutions=resolutions)
            elif self.is_spline and self.data.para_dim > 2:
                # TODO: replace this with generalized query helpers.
                return self.data.extract.faces(
                    resolutions, watertight=False
                ).vertices
            else:
                return self.data.sample(resolutions)

        # runtime location query
        if on is not None:
            if self.has_function:
                return self.function(self.data, on=on)
            elif self.has_evaluate:
                return self.data.evaluate(on)
            else:
                raise ValueError(
                    "Given data can't support data extraction on specified "
                    f"locations ({on})."
                )

        # location specified - either evaluate function at the
        if self.has_locations:
            if self.has_function:
                # function may want locations
                try:
                    return self.function(self.data, self.locations)
                except TypeError:  # maybe too many args
                    return self.function(self.data)
            elif self.has_evaluate:
                return self.data.evaluate(self.locations)
            else:
                return self.data

        # should be returned by now
        raise RuntimeError("Something went wrong while preparing spline data.")


class SplineData(DataHolder):
    """
    Data manager for splines.
    """

    def __init__(self, helpee):
        """ """
        if "Spline" not in str(type(helpee).__mro__):
            raise AttributeError("Helpee is not a Spline")

        super().__init__(helpee)

    def __setitem__(self, key, value):
        """
        Selectively accept spline data.

        Parameters
        ----------
        key: str
        value: object

        Returns
        -------
        None
        """
        if isinstance(value, SplineDataAdaptor):
            self._saved[key] = value
        else:
            adapted = SplineDataAdaptor(value)  # will test usability
            adapted._user_created = False  # mark for __getitem__
            self._saved[key] = adapted

    def __getitem__(self, key):
        """
        Returns value from __setitem__

        Parameters
        ----------
        key: str

        Returns
        -------
        value: object
        """
        saved = super().__getitem__(key)
        if saved._user_created:
            return saved
        else:
            return saved.data

    def as_scalar(self, key, resolutions, default=None):
        """
        Return scalar value at given resolutions

        Parameters
        ----------
        key: str
        resolutions: list or tuple
        default: object
          Default is None and will return is key doesn't exist

        Returns
        -------
        value: np.ndarray
        """
        if key not in self._saved:
            return default

        saved = super().__getitem__(key)
        # will raise
        return saved.as_vertex_data(resolutions=resolutions)

    def as_arrow(self, key, resolutions=None, on=None, default=None):
        """
        Returns as-arrow-representable data on certain places, with given
        resolution, or on predefined places.

        Parameters
        ----------
        key: str
        resolutions: list or tuple
        on: array-like
        """
        if key not in self._saved:
            return default

        saved = super().__getitem__(key)
        # will raise
        return saved.as_vertex_data(resolutions=resolutions, on=on)
