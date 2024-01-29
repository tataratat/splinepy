"""splinepy/utils/data.py.

Helps helpee to manage data. Some useful data structures.
"""

from itertools import accumulate as _accumulate
from itertools import chain as _chain
from itertools import repeat as _repeat

import numpy as _np
from gustaf.helpers.data import DataHolder as _DataHolder
from gustaf.utils.arr import enforce_len as _enforce_len

from splinepy import splinepy_core as _splinepy_core
from splinepy._base import SplinepyBase as _SplinepyBase

try:
    import scipy as _scipy

    has_scipy = True
except ImportError:
    has_scipy = False

__all__ = [
    "PhysicalSpaceArray",
    "_DataHolder",
    "enforce_contiguous_values",
    "cartesian_product",
    "SplineDataAdaptor",
    "SplineData",
]


class PhysicalSpaceArray(_np.ndarray):
    """numpy array object that keeps mirroring inplace changes to the source.
    Meant to help control_points.
    """

    __slots__ = (
        "_source_ptr",
        "_row_indices",
        "_full_row_indices",
        "_super_arr",
    )

    def __array_finalize__(self, obj):
        """Sets default flags for any arrays that maybe generated based on
        physical space array. For more information,
        see https://numpy.org/doc/stable/user/basics.subclassing.html"""
        self._source_ptr = None
        self._super_arr = None

        # for arrays created based on this subclass
        if isinstance(obj, type(self)):
            # this is copy. nothing to worry here
            if self.base is None:
                return None

            # first child array
            if self.base is obj:
                # make sure this is not a recursively born child
                # for example, `arr[[1,2]][:,2]`
                if obj._source_ptr is not None:
                    self._super_arr = obj
                return None

            # multi generation child array
            if obj._super_arr is not None and self.base is obj.base:
                self._super_arr = obj._super_arr
                return None

            return None

    def _sync_source_ptr(self):
        # this is super arr. super arr only has _source_ptr because we should've
        # set it manually
        if self._source_ptr is not None:  # and self._super_arr is None
            self._source_ptr.sync(self)
            return None

        # this is a sub arr viewing super arr
        # this will perform a full sync
        if self._super_arr is not None and self._source_ptr is None:
            self._super_arr._source_ptr.sync(self._super_arr)
            return None

    def copy(self, *args, **kwargs):
        """copy creates regular numpy array"""
        return _np.array(self, *args, copy=True, **kwargs)

    def view(self, *args, **kwargs):
        """Set writeable flags to False for the view."""
        v = super(self.__class__, self).view(*args, **kwargs)
        v.flags.writeable = False
        return v

    def __iadd__(self, *args, **kwargs):
        sr = super(self.__class__, self).__iadd__(*args, **kwargs)
        self._sync_source_ptr()
        return sr

    def __isub__(self, *args, **kwargs):
        sr = super(self.__class__, self).__isub__(*args, **kwargs)
        self._sync_source_ptr()
        return sr

    def __imul__(self, *args, **kwargs):
        sr = super(self.__class__, self).__imul__(*args, **kwargs)
        self._sync_source_ptr()
        return sr

    def __idiv__(self, *args, **kwargs):
        sr = super(self.__class__, self).__idiv__(*args, **kwargs)
        self._sync_source_ptr()
        return sr

    def __itruediv__(self, *args, **kwargs):
        sr = super(self.__class__, self).__itruediv__(*args, **kwargs)
        self._sync_source_ptr()
        return sr

    def __imatmul__(self, *args, **kwargs):
        sr = super(self.__class__, self).__imatmul__(*args, **kwargs)
        self._sync_source_ptr()
        return sr

    def __ipow__(self, *args, **kwargs):
        sr = super(self.__class__, self).__ipow__(*args, **kwargs)
        self._sync_source_ptr()
        return sr

    def __imod__(self, *args, **kwargs):
        sr = super(self.__class__, self).__imod__(*args, **kwargs)
        self._sync_source_ptr()
        return sr

    def __ifloordiv__(self, *args, **kwargs):
        sr = super(self.__class__, self).__ifloordiv__(*args, **kwargs)
        self._sync_source_ptr()
        return sr

    def __ilshift__(self, *args, **kwargs):
        sr = super(self.__class__, self).__ilshift__(*args, **kwargs)
        self._sync_source_ptr()
        return sr

    def __irshift__(self, *args, **kwargs):
        sr = super(self.__class__, self).__irshift__(*args, **kwargs)
        self._sync_source_ptr()
        return sr

    def __iand__(self, *args, **kwargs):
        sr = super(self.__class__, self).__iand__(*args, **kwargs)
        self._sync_source_ptr()
        return sr

    def __ixor__(self, *args, **kwargs):
        sr = super(self.__class__, self).__ixor__(*args, **kwargs)
        self._sync_source_ptr()
        return sr

    def __ior__(self, *args, **kwargs):
        sr = super(self.__class__, self).__ior__(*args, **kwargs)
        self._sync_source_ptr()
        return sr

    def __setitem__(self, key, value):
        # set first. invalid setting will cause error
        sr = super(self.__class__, self).__setitem__(key, value)

        # for child array, we will just sync and exit.
        if self._source_ptr is None:
            self._sync_source_ptr()
            return sr

        # if set item above worked.
        # then, mirror / copy rowwise
        # if this was called due to inplace operations of subset of entries
        # this would be redundant
        #
        # here, we will try to parse very simple, yet, quite often used
        # indexings
        contig = _np.ascontiguousarray  # frequently used
        if isinstance(key, int) or (
            isinstance(key, _np.ndarray) and key.ndim == 0
        ):
            # single get item is always contiguous
            self._source_ptr.set_row(key, self[key])

        elif isinstance(key, _np.ndarray) and key.ndim == 1:
            if key.dtype.kind.startswith("b"):
                key = contig(self.row_indices()[key])

            self._source_ptr.sync_rows(contig(key), self)

        elif isinstance(key, list):
            if isinstance(key[0], bool):
                key = contig(self.row_indices()[key])
            self._source_ptr.sync_rows(key, self)

        # tuple is a bit tricky
        elif isinstance(key, tuple) and (
            isinstance(key[0], int)
            or (isinstance(key[0], _np.ndarray) and key[0].ndim == 0)
        ):
            self._source_ptr.set_row(key[0], self[key[0]])

        elif isinstance(key, tuple) and (
            isinstance(key[0], (list, _np.ndarray))
        ):
            if isinstance(key[0][0], (_np.bool_, bool)):
                key = (contig(self.row_indices()[key[0]]),)

            self._source_ptr.sync_rows(key[0], self)

        elif isinstance(key, slice) or (
            isinstance(key, tuple) and isinstance(key[0], slice)
        ):
            key = key if isinstance(key, slice) else key[0]

            start = 0 if key.start is None else key.start
            # flip negatives
            start = start + len(self) if start < 0 else start
            stop = len(self) if key.stop is None else key.stop
            # flip negatives
            stop = stop + len(self) if stop < 0 else stop
            step = 1 if key.step is None else key.step

            ids = _np.arange(
                start, stop, step, dtype="int32"
            )  # should be contig
            self._source_ptr.sync_rows(ids, self)

        elif isinstance(key, type(Ellipsis)) or (
            isinstance(key, tuple) and isinstance(key[0], type(Ellipsis))
        ):
            self._source_ptr.sync(self)

        else:
            # even more complex? then
            # TODO - just sync if this seems to take too long
            ids = _np.unique(self.full_row_indices()[key])
            self._source_ptr.sync_rows(contig(ids), self)

        return sr

    def row_indices(self):
        """
        Returns row_indices. Same results as np.arange(len(arr))
        """
        if hasattr(self, "_row_indices"):
            return self._row_indices

        len_, d = self._source_ptr.len(), self._source_ptr.dim()
        self._row_indices = _np.arange(len_, dtype="int32")
        self._full_row_indices = self._row_indices.repeat(d).reshape(len_, d)
        return self._row_indices

    def full_row_indices(self):
        """
        Returns full_row_indices. Same results as np.indices(arr.shape)[0]
        """
        if hasattr(self, "_full_row_indices"):
            return self._full_row_indices

        len_, d = self._source_ptr.len(), self._source_ptr.dim()
        self._row_indices = _np.arange(len_, dtype="int32")
        self._full_row_indices = self._row_indices.repeat(d).reshape(len_, d)
        return self._full_row_indices


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
    if isinstance(array, _np.ndarray):
        if array.flags["C_CONTIGUOUS"] and (
            dtype is None or _np.dtype(dtype) is array.dtype
        ):
            return array
        return _np.ascontiguousarray(array, dtype=dtype)

    if asarray and isinstance(array, (list, tuple, _splinepy_core.KnotVector)):
        return _np.ascontiguousarray(array, dtype=dtype)

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
    dict_with_contiguous = {}
    for key, value in dict_.items():
        # most likely will be asking about np.ndarray
        if isinstance(value, _np.ndarray):
            contig_value = enforce_contiguous(value)

        # or, list. only check if the first element is np.ndarray
        elif isinstance(value, list) and isinstance(value[0], _np.ndarray):
            contig_value = [enforce_contiguous(v) for v in value]

        # abandon Nones
        elif value is None:
            continue

        else:
            contig_value = value

        # set item
        dict_with_contiguous[key] = contig_value

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
    dtype = _np.result_type(*arrays)

    # for reverse, use view of reverse ordered arrays
    if reverse:
        arrays = arrays[::-1]

    shape = *map(len, arrays), n_arr
    cartesian = _np.empty(shape, dtype=dtype)

    # reverse order the view of output array, so that reverse
    # is again in original order, but just reversed fastest iterating dim
    _cartesian = cartesian[..., ::-1] if reverse else cartesian

    dim_reducing_views = (
        *_accumulate(
            _chain((_cartesian,), _repeat(0, n_arr - 1)),
            _np.ndarray.__getitem__,
        ),
    )
    idx = slice(None), *_repeat(None, n_arr - 1)
    for i in range(n_arr - 1, 0, -1):
        dim_reducing_views[i][..., i] = arrays[i][idx[: n_arr - i]]
        dim_reducing_views[i - 1][1:] = dim_reducing_views[i]
    _cartesian[..., 0] = arrays[0][idx]

    return cartesian.reshape(-1, n_arr)


def make_matrix(values, supports, n_cols, as_array=False):
    r"""
    Create a matrix from values and supports.

    Used for basis functions their derivatives and everything mapped.
    Uses scipy if available. If `as_array` is true, dense matrix (numpy) is
    enforced

    This matrix can be used for approximations and IGA-applications. With given
    queries :math:`\pmb{\xi}`, control points :math:`\mathbf{C}`, associated to
    spline :math:`\mathcal{S}`, this function can be used to return matrix
    :math:`\mathbf{A}`, such that

    .. math::
        \mathcal{S}(\pmb{\xi}) = A(\pmb{\xi}) \cdot \mathbf{C}

    Parameters
    ----------
    values : np.ndarray
      Values to be inserted into the matrix
    supports : np.ndarray
      Column values corresponding to the given values
    n_cols : int
      Number of columns of the matrix (matches n_ctps for basis-functions)
    as_array : bool
      Return as numpy / dense type

    Returns
    -------
    matrix : np.ndarray / scipy.sparse.csr_array
      Matrix
    """
    # Return in Matrix shape (scipy if available)
    if has_scipy and not as_array:
        return _scipy.sparse.csr_array(
            (
                values.ravel(),
                (
                    _np.arange(values.shape[0]).repeat(values.shape[1]),
                    supports.ravel(),
                ),
            ),
            shape=(values.shape[0], n_cols),
        )
    else:
        matrix = _np.zeros((values.shape[0], n_cols))
        _np.put_along_axis(matrix, supports, values, axis=1)
        return matrix


def uniform_query(bounds, resolutions):
    """
    Creates uniform query within the given bounds and resolutions.
    Same as `gus.create.vertices.raster(bounds, resolution).vertices`.

    Parameters
    ----------
    bounds: 2D array-like
      [[lower_1, lower_2, ...], [upper_1, upper_2, ...]]
    resolutions: 1D array-like
      [resolution_1, resolution_2, ...]
    """
    # unpack bounds
    lower_b, upper_b = bounds

    if not (len(lower_b) == len(upper_b) == len(resolutions)):
        raise ValueError(
            "lower and upper bounds should have same len as resolutions"
        )

    # create per-dimension queries
    queries_per_dim = []
    for lb, ub, r in zip(lower_b, upper_b, resolutions):
        queries_per_dim.append(_np.linspace(lb, ub, r))

    return cartesian_product(queries_per_dim, reverse=True)


class SplineDataAdaptor(_SplinepyBase):
    """
    Prepares data to be presentable on spline. To support both
    scalar-data and vector-data, which are representable with colors and
    arrows respectively, this class will prepare data accordingly.

    Parameters
    ----------
    data: any
      Any data that you want to plot on to the spline. If this is not a
      spline, it requires a `function` that takes this data to create
      appropriate values.
    locations: 2D array-like
      Optional. If specified, used to evaluate data at specific locations.
      Applicable for arrow data
    function: callable
      Optional. A callable used to evaluate values using data at query points.
      function(data, queries) -> values
    parametric_bounds: 2D array-like
      Optional. Parametric bounds of supporting spline. Required for stand-alone
      use. Otherwise, set by helping spline. Used to compute query points for
      resolution based sampling.
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
        "parametric_bounds",
        "_user_created",
    )

    def __init__(
        self, data, locations=None, function=None, parametric_bounds=None
    ):
        """ """
        from splinepy import Multipatch as _Multipatch
        from splinepy import Spline as _Spline

        # default
        self._user_created = True
        self.data = data
        self.is_spline = False
        self.has_function = False
        self.has_locations = False
        self.has_evaluate = False
        self.arrow_data_only = False
        self.parametric_bounds = parametric_bounds

        # is spline we know?
        if isinstance(data, (_Spline, _Multipatch)):
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
            self.locations = _np.asanyarray(locations)

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
                # parametric bounds needs to be set by now.
                queries = uniform_query(
                    self.parametric_bounds,
                    _enforce_len(resolutions, len(self.parametric_bounds[0])),
                )
                # pass to evaluate.
                return self.function(
                    self.data,
                    on=queries,
                )
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


class SplineData(_DataHolder):
    """
    Data manager for splines.
    """

    __slots__ = ()

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

        # set parametric_bounds
        self._saved[key].parametric_bounds = self._helpee.parametric_bounds

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


class MultipatchData(SplineData):
    """
    Data manager for multipatch
    """

    __slots__ = ()

    def __init__(self, helpee):
        """
        Strictly, we don't need to inherit from SplineData and only implement
        interfaces, since _saved equivalent is managed in multipatch.
        However, it's an index based holder and we want to support the same
        syntax.
        Here, we will use _saved as key_to_id lookup
        """
        if "Multipatch" not in str(type(helpee).__mro__):
            raise AttributeError("Helpee is not a multipatch")

        _DataHolder.__init__(self, helpee)

    def __setitem__(self, key, value):
        """
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
            # set parametric_bounds of the first patch
            value.parametric_bounds = self._helpee.patches[0].parametric_bounds
        elif "PyMultipatch" in str(type(value).__mro__):
            # get id of this field that's about to be added
            self._saved[key] = len(self._helpee.fields)
            # add fields expects list of list
            # TODO @jzwar add type check to except both?
            self._helpee.add_fields([value.patches], value.dim)
        elif isinstance(value, int):
            # also accept integers as a reference to the field
            self._saved[key] = value
        else:
            raise TypeError(
                "MultipatchData supports SplineDataAdapter or Multipatch."
                f"given type: {type(value)}"
            )

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

        # integer input refers to saved field's number
        if isinstance(key, int):
            return self._helpee.fields[key]
        elif isinstance(key, str):
            saved = _DataHolder.__getitem__(self, key)  # will raise KeyError
            if isinstance(saved, int):
                return self._helpee.fields[saved]
            elif isinstance(saved, SplineDataAdaptor):
                return saved
        raise RuntimeError(
            "Invalid saved data / key type. Please help us by writing an "
            "issue at github.com/tataratat/splinepy. Thank you!"
        )

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

        saved = self.__getitem__(key)
        if isinstance(saved, SplineDataAdaptor):
            # explicitly checks that saved IS data adaptor
            pass
        else:
            # or will create one on the fly
            saved = SplineDataAdaptor(saved)

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

        saved = self.__getitem__(key)
        if isinstance(saved, SplineDataAdaptor):
            # explicitly checks that saved IS data adaptor
            pass
        else:
            # or will create one on the fly
            saved = SplineDataAdaptor(saved)

        return saved.as_vertex_data(resolutions=resolutions, on=on)
