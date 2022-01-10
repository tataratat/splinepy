"""splinelibpy/splinelibpy/utils.py

Utility functions.
"""

import os
import logging

import numpy as np


def configure_logging(debug=False, logfile=None):
    """
    Logging configurator.

    Parameters
    -----------
    debug: bool
    logfile: str

    Returns
    --------
    None
    """
    logger = logging.getLogger()
    if debug:
        logger.setLevel(logging.DEBUG)

    else:
        logger.setLevel(logging.INFO)

    if logfile is not None:
        file_logger_handler = logging.FileHandler(logfile)


def is_property(property_dict, key, class_name):
    """
    Checks if property exist in given dict. If not add a debug log.

    Parameters
    -----------
    property_dict: dict
    key: str
    class_name: str

    Returns
    --------
    is_property: bool
    """
    if key in property_dict:
        return True

    else:
        logging.debug(
            class_name
            + " - `"
            + key
            + "` does not exist yet."
        )
        return False


def make_c_contiguous(array, dtype=None):
    """
    Make given array like object a c contiguous np.ndarray.
    dtype is optional.

    Parameters
    -----------
    array: array-like
    dtype: type or str
      (Optional) `numpy` interpretable type or str, describing type.
      Difference is that type will always return a copy and str will only copy
      if types doesn't match.

    Returns
    --------
    c_contiguous_array: np.ndarray
    """
    if isinstance(array, np.ndarray):
        if array.flags.c_contiguous:
            if dtype is not None:
                if isinstance(dtype, type):
                    return array.astype(dtype)

                elif isinstance(dtype, str):
                    if array.dtype.name != dtype:
                        return array.astype(dtype)

            return array

    if dtype:
        return np.ascontiguousarray(array, dtype=dtype)

    else:
        return np.ascontiguousarray(array)


def abs_fname(fname):
    """
    Checks if fname is absolute. If not, returns abspath. Tilde safe.

    Parameters
    ----------
    fname: str

    Returns
    --------
    abs_fname: str
      Maybe same to fname, maybe not.
    """
    if os.path.isabs(fname):
        pass

    elif "~" in fname:
        fname = os.path.expanduser(fname)

    else:
        fname = os.path.abspath(fname)

    return fname
