import numpy as np

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
    dtype: type
      (Optional) `numpy` interpretable type.

    Returns
    --------
    c_contiguous_array: np.ndarray
    """
    if isinstance(array, np.ndarray):
        if array.flags.c_contiguous:
            if dtype is not None:
                if array.dtype is not dtype:
                    return array.astype(dtype)

            return array

    if dtype:
        return np.ascontiguousarray(array, dtype=dtype)

    else:
        return np.ascontiguousarray(array)
