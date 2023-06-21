"""splinepy/utils/data.py.

Helps helpee to manage data. Some useful data structures.
"""
from itertools import accumulate, chain, repeat

import numpy as np
from gustaf.helpers.data import TrackedArray, make_tracked_array

__all__ = [
    "TrackedArray",
    "make_tracked_array",
    "is_modified",
    "enforce_contiguous_values",
    "cartesian_product",
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


def enforce_contiguous(array, dtype=None):
    """
    If input is an instance / subclass of np.ndarray, this will check
    if they are configuous. If so, returns same object, else turns makes it
    contiguous and returns.

    Parameters
    ----------
    array: array-like

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
