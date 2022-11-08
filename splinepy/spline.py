"""
Abstract Spline
"""

import copy
import os

import numpy as np

from splinepy import utils
from splinepy import io
from splinepy import splinepy_core as core
from splinepy._base import SplinepyBase


class RequiredProperties:
    """
    Helper class to hold required properties of each spline.
    """

    # name mangling to make direct access very annoying
    __required_spline_properties = {
        "Bezier": ["degrees", "control_points"],
        "RationalBezier": ["degrees", "control_points", "weights"],
        "BSpline": ["degrees", "knot_vectors", "control_points"],
        "NURBS": ["degrees", "knot_vectors", "control_points", "weights"],
    }

    @classmethod
    def __call__(cls, spline):
        """
        Given splinepy.Spline, returns required properties.
        Wraps `__getitem__`.

        Parameters
        -----------
        spline: Spline or type

        Returns
        --------
        required_spline_properties: list
        """
        if isinstance(spline, str):
            return cls.__getitem__(spline)

        # extract pure `type`
        s_type = spline if isinstance(spline, type) else type(spline)

        # __getitem__ checks if it is supported
        return cls.__getitem__(s_type.__qualname__)

    @classmethod
    def __getitem__(cls, spline_type):
        """
        Given a name of spline type in str, returns required properties.
        Checks if given spline type is supported.

        Parameters
        -----------
        spline_type: str

        Returns
        --------
        required_spline_properties: list
        """
        # do we support this spline type?
        if spline_type not in cls.__required_spline_properties:
            raise ValueError(
                f"Sorry, we don't have support for ({spline_type})-types."
                "Supported spline types are:"
                f"{list(cls.__required_spline_properties.keys())}"
            )

        # copy to avoid inplace manipulation
        return cls.__required_spline_properties[spline_type].copy()

    @classmethod
    def of(cls, spline):
        """
        Wrapper of __call__ and __getitem__ as classmethod.

        Parameters
        -----------
        spline: Spline or str

        Returns
        --------
        required_spline_properties: list

        Examples
        ---------
        >>> myspline = NURBS(...)
        >>> requirements = RequiredSplineProperties.of(myspline)
        >>> requirements
        ['degrees', 'knot_vectors', 'control_points', 'weights']
        """
        if isinstance(spline, str):
            return cls.__getitem__(spline)
        else:
            return cls.__call__(spline)

    @classmethod
    def union(cls, *splines):
        """
        Checks union of required properties of given splines.
        If nothing is given, returns union of all supported splines.

        Parameters
        ----------
        *splines: Spline or str

        Returns
        -------
        rp_union: list
        """
        if len(splines) == 0:
            splines = list(cls.__required_spline_properties.keys())

        rp_union = set()
        for s in splines:
            rp_union.update(cls.of(s))

        return list(rp_union)

    @classmethod
    def intersection(cls, *splines):
        """
        Checks intersection of required properties of given splines.
        If nothing is given, returns intersection of all supported splines.

        Parameters
        ----------
        *splines: Spline or str

        Returns
        -------
        rp_intersection: list
        """
        if len(splines) == 0:
            splines = list(cls.__required_spline_properties.keys())

        rp_intersection = set(cls.of(splines[-1])
        for s in splines[:-1]:
            rp_intersection.intersection_update(cls.of(s))

        return list(rp_intersection)


# initialize one for direct use
required_properties = RequiredProperties()


def sync(spl, c2p=False):
    """
    Helper function to sync python spline with cpp core spline.
    Data syncing direction is python to cpp, unless c2p is set to True.

    Parameters
    -----------
    spl: Spline
    c2p: bool
      Default is False. Transfers core spline properties to this spline.

    Returns
    --------
    None
    """
    if c2p:
        # properly sync current properties as tracked array
        for key, prop in spl._core_spline.current_properties().items():
            # use setters
            setattr(spl, key, prop)

        spl._data["identities"] = spl._core_spline.current_identities()

        # at this point, c and p are synced. 
        return None

    # p2c
    required_props = spl.required_properties

    # gather properties into dict while we loop through each props.
    # will be used  to init core class later.
    dict_spline = dict()
    modified = False

    # Check if all data is available for spline type
    for i, rp in enumerate(required_props):
        tmp_rp = getattr(self, rp)
        if tmp_rp is None:
            spl._logd(
                f"Missing {rp} to initialize core spline.",
                "Skipping or removing existing core spline."
            )
            if hasattr(spl, "_core_spline"):
                spl._core_spline = None
                delattr(spl, "_core_spline")

            return None

        else:
            dict_spline[rp] = tmp_rp
            modified = utils.data.is_modified(tmp_rp)

    # update core spline
    if modified:
        spl._core_spline = core.CoreSpline(**dict_spline)

    else:
        # we are here because everthing is synced
        return None


def is_modified(spl):
    """
    Checks if there are inplace changes in spline properties.
    If properties are incomplete, it will return answers to existing
    properties. Default is False - for empty splines, it will return False.

    Parameters
    -----------
    spl: Spline

    Returns
    --------
    modified: bool
    """
    modified = False
    for rp in spl.required_properties:
        prop = getattr(spl, rp, None)
        if prop is not None:
            modified |= utils.data.is_modified(prop)

    return modified


def _set_modified_false(spl):
    """
    Sets all modified flags to False. Internal use only.

    Parameters
    ----------
    spl: Spline

    Returns
    -------
    None
    """
    for rp in spl.required_properties:
        prop = getattr(spl, rp, None)

        if prop is not None:
            if isinstance(prop, list):
                for p in prop:
                    p._modified = False

            prop._modified = False


def identities(spl):
    """
    Returns identity of given spline
    """

def properties(spl):
    """
    """


class Spline:
    """
    Spline base class. Wraps CoreSpline.
    """

    __slots__ = (
        "_logi",`
        "_logd",`
        "_logw",`
        "_core_spline",
        "_data",
    )

    def __init__(
            self,
            degrees=None,
            knot_vectors=None,
            control_points=None,
            weights=None,
            c_spline=None,
    ):
        """
        Spline.

        Parameters
        -----------
        degrees: (para_dim,) list-like
        knot_vectors: (para_dim, n) list
        control_points: (m, dim) list-like
        weights: (m,) list-like

        Attributes
        -----------
        whatami: str
        para_dim: int
        dim: int
        degrees: np.ndarray
        knot_vectors: list
        control_points: np.ndarray
        parametric_bounds: np.ndarray
        control_point_bounds: np.ndarray

        Returns
        --------
        None
        """
        # logger shortcut
        self._logi = log.prepend_log(type(self).__qualname__, utils.log.info)
        self._logd = log.prepend_log(type(self).__qualname__, utils.log.debug)
        self._logw = log.prepend_log(
            type(self).__qualname__, utils.log.warning
        )

        # initialize rest of the slots
        self._data = dict()
        self._core_spline = None

        self.degrees = degrees
        self.control_points = control_points

        # Bezier has no kvs
        if knot_vectors is not None:
            self.knot_vectors = knot_vectors

        # rational splines has weights
        if weights is not None:
            self.weights = weights

        # if there's c_spline, we will share the c_spline.
        if c_spline is not None:
            self._core_spline = c_spline

    def clear(self):
        """
        Clear all properties.

        Parameters
        -----------
        None

        Returns
        --------
        None
        """
        self._core_spline = None
        self._data = dict()

        self._logd("All attributes are cleared!")

    @property
    def _c_spline(self):
        """
        Returns core spline. Updates core spline before returning if
        any property is modified.

        Parameters
        -----------
        None

        Returns
        --------
        core_spline: CoreSpline
        """
        if self._core_spline is not None:
            return self._core_spline

        else:
            props_set = dict()
            # prepare error message to hint missing attr
            for rp in self.required_properties:
                props_set[rp] = "set"  if getattr(self, rp) else "missing"

            raise ValueError(
                    "Spline is not initialized yet."
                    f"Property status: {props_set}"
            ) 

    @_c_spline.setter
    def _c_spline(self, spline):
        """
        Sets spline using the core spline.

        Parameters
        ----------
        spline: dict or CoreSpline
          dict describing spline or CoreSpline. In case of CoreSpline, this
          object will share core spline.

        Returns
        -------
        None
        """
        spline_t = type(spline)
        # maybe this can be a helper
        #if isinstance(spline_t, dict):
        #    # let's at least check if those are what we support
        #    req_prop_all = RequiredProperties.union()
        #    for k in spline.keys():
        #        if k not in req_prop_all:
        #            raise KeyError(
        #                    f"{k} is not supported property of any spline"
        #            )
        #    self._core_spline = core.CoreSpline(**spline)

        #elif isinstance(spline_t, core.CoreSpline):
        if isinstance(spline_t, core.CoreSpline):
            self._core_spline = spline # take given ref
            sync(self, c2p=True)

        elif spline is None:
            self.clear()

        else:
            raise TypeError(f"{spline_t} is not supported _c_spline.")

    @property
    def required_properties(self):
        """
        Returns required property of current spline.

        Parameters
        -----------
        None

        Returns
        --------
        requied_properties: dict
        """
        return RequiredProperties.of(self)

    @property
    def whatami(self):
        """
        Answers a deep philosophical question of "what am i?"

        Parameters
        -----------
        None

        Returns
        --------
        whatami: str
        """
        return self._data.get("identities", dict()).get("whatami", "Nothing")

    @property
    def para_dim(self):
        """
        Returns parametric dimension.

        Parameters
        -----------
        None

        Returns
        --------
        para_dim: int
        """
        return self._data.get()
        if hasattr(self, "_c_spline"):
            return self._para_diml

        else:
            return None

    @property
    def dim(self):
        """
        Returns physical dimension.

        Parameters
        -----------
        None

        Returns
        --------
        dim: int
        """
        if hasattr(self, "_dim"):
            return self._dim 

        else:
            return None

    @propertyr
    def degrees(self):
        """
        Returns Degrees.

        Parameters
        -----------
        None

        Returns
        --------
        degrees: (para_dim,) np.ndarray
        """
        if hasattr(self, "_degrees"):
            return self._degrees

        else:
            return None

    @degrees.setter
    def degrees(self, degrees):
        """
        Sets Degrees.

        Parameters
        -----------
        degrees: list-like

        Returns
        --------
        None
        """
        if degrees is None:
            if hasattr(self, "_degrees"):
                delattr(self, "_degrees")
                if hasattr(self, "_para_dim") and self.knot_vectors is None:
                    delattr(self, "_para_dim")
            return None

        # flatten copies
        degrees = utils.make_c_contiguous(degrees, "int32").flatten()

        if self.para_dim is None:
            self._para_dim = degrees.size

        else:
            if self.para_dim != degrees.size:
                raise InputDimensionError(
                    "Input dimension does not match spline's para_dim "
                )

        self._degrees = degrees

        self._logd(f"Degrees set: {self.degrees}")

        self._check_and_update_c()

    @property
    def knot_vectors(self):
        """
        Returns knot vectors.

        Parameters
        -----------
        None

        Returns
        --------
        knot_vectors: list
        """
        if hasattr(self, "_knot_vectors"):
            return self._knot_vectors

        else:
            return None

    @knot_vectors.setter
    def knot_vectors(self, knot_vectors):
        """
        Set knot vectors.

        Parameters
        -----------
        knot_vectors: list

        Returns
        --------
        None
        """
        if knot_vectors is None:
            if hasattr(self, "_knot_vectors"):
                delattr(self, "_knot_vectors")
                if hasattr(self, "_para_dim") and self.degrees is None:
                    delattr(self, "_para_dim")
            return None

        if not isinstance(knot_vectors, list):
            raise TypeError("knot_vectors should be a `list`!")

        if self.para_dim is None:
            self._para_dim = len(knot_vectors)

        else:
            if self.para_dim != len(knot_vectors):
                raise InputDimensionError(
                    "Input dimension does not match spline's para_dim "
                )

        self._knot_vectors = copy.deepcopy(knot_vectors)

        self._logd("Knot vectors set:")
        for i, kv in enumerate(self.knot_vectors):
            self._logd(
                f"  {i}"
                ". knot vector length: "
                f"{len(kv)}"
            )

        self._check_and_update_c()

    @property
    def unique_knots(self,):
        """
        Returns unique knots.
        Computed using `np.unique`.
        Does not store results.

        Parameters
        -----------
        None

        Returns
        --------
        unique_knots: list
        """
        self._logd("Computing unique knots using `np.unique`.")
        unique_knots = []
        if "Bezier" in self.whatami:
            unique_knots = [[0, 1]] * self.para_dim
        else:
            for k in self.knot_vectors:
                unique_knots.append(np.unique(k).tolist())

        return unique_knots

    @property
    def parametric_bounds(self,):
        """
        Returns bounds of parametric_space.

        Parameters
        -----------
        None

        Returns
        --------
        parametric_bounds: (2, para_dim) np.ndarray
        """
        self._logd("Computing parametric_bounds")
        # beziers
        if "knot_vectors" not in self.required_properties:
            return [[0, 1] * self.para_dim]

        # bsplines
        kvs = self.knot_vectors
        if kvs is None:
            return None

        lower_bounds = []
        upper_bounds = []

        use_minmax = False
        if not hasattr(self, "_c_spline"):
            # in this case, `_check_and_update_c()` wasn't called.
            self._logd(
                "Entries of `knot_vectors` has not been checked. "
                "Values of `parametric_bounds` will be min and max."
            )
            use_minmax = True

        for kv in kvs:
            lower_bounds.append(min(kv) if use_minmax else kv[0])
            upper_bounds.append(max(kv) if use_minmax else kv[-1])

        return np.vstack((lower_bounds, upper_bounds))

    # keep this until gustaf stops using knot_vector_bounds
    knot_vector_bounds = parametric_bounds

    @property
    def control_points(self,):
        """
        Returns control points.

        Parameters
        -----------
        None

        Returns
        --------
        control_points_: (n, dim) np.ndarray
        """
        if hasattr(self, "_control_points"):
            return self._control_points

        else:
            return None

    @control_points.setter
    def control_points(self, control_points):
        """
        Set control points.

        Parameters
        -----------
        control_points: (n, dim) list-like

        Returns
        --------
        None
        """
        if control_points is None:
            if hasattr(self, "_control_points"):
                delattr(self, "_control_points")
                delattr(self, "_dim")
            return None

        control_points = utils.make_c_contiguous(
            control_points,
            "float64"
        ).copy()

        if self.dim is None:
            self._dim = control_points.shape[1]

        else:
            if self.dim != control_points.shape[1]:
                raise InputDimensionError(
                    "Input dimension does not match spline's dim "
                )

        self._control_points = control_points
        self._logd(
            f"{self.control_points.shape[0]} Control points set."
        )

        self._check_and_update_c()

    def _lexsort_control_points(self, order):
        """
        Lexsorts control points.
        Not always the solution.
        Leading underscore to mark development/internal use.

        Parameters
        -----------
        order: list

        Returns
        --------
        None
        """
        ind = np.lexsort([self.control_points[:, i] for i in order])
        self._logd(f"`lexsort` control points ({order})")
        self.control_points = self.control_points[ind]

    @property
    def control_point_bounds(self,):
        """
        Returns bounds of control points.

        Parameters
        -----------
        None

        Returns
        --------
        control_point_bounds: (2, dim) np.ndarray
        """
        self._logd("Computing control_point_bounds")
        cps = self.control_points

        return np.vstack(
            (
                cps.min(axis=0),
                cps.max(axis=0),
            )
        )

    @property
    def weights(self,):
        """
        Returns weights.

        Parameters
        -----------
        None

        Returns
        --------
        self._weights: (n, 1) list-like
        """
        if hasattr(self, "_weights"):
            return self._weights

        else:
            return None

    @weights.setter
    def weights(self, weights):
        """
        Set weights.

        Parameters
        -----------
        weights: (n,) list-like

        Returns
        --------
        None
        """
        if weights is None:
            if hasattr(self, "_weights"):
                delattr(self, "_weights")
            return None

        weights = utils.make_c_contiguous(
            weights,
            dtype=np.float64
        ).reshape(-1, 1)

        if self.control_points is not None:
            if self.control_points.shape[0] != weights.shape[0]:
                raise ValueError(
                    "Number of control points and number of weights does not "
                    "match."
                )

        self._weights = weights

        self._logd(f"{self.weights.shape[0]} Weights set.")

        self._check_and_update_c()

    @property
    def control_mesh_resolutions(self,):
        """
        Returns control mesh resolutions.

        Parameters
        -----------
        None

        Returns
        --------
        control_mesh_resolutions: list

        Raises
        -------
        TypeError: if one of the required properties to compute
          control_mesh_resolutions is missing.
        """
        cmr = []

        # Special case Bezier
        if "Bezier" in type(self).__qualname__:
            cmr = [d + 1 for d in self.degrees]
        else:
            for kv, d in zip(self.knot_vectors, self.degrees):
                cmr.append(len(kv) - d - 1)

        return cmr

    def _check_and_update_c(self,):
        """
        Check all available information before updating the backend

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        required_props = self.required_properties



        # Check if enough knot vectors were given
        if "knot_vectors" in required_props:
            if len(self.degrees) != len(self.knot_vectors):
                raise InputDimensionError(
                    "Dimension mis-match between `degrees` and `knot_vectors`"
                )

            # Check if knot vectors are large enough
            for i, (d, kv) in enumerate(zip(self.degrees, self.knot_vectors)):
                if len(kv) < int(2 * (d+1)):
                    raise InputDimensionError(
                        "Not enough knots in knot vector along parametric"
                        f" dimension {i}"
                    )

                # check if knots are in increasing order
                if (np.diff(kv) < 0).any():
                    raise ValueError(
                        f"Knots of {i}. "
                        "knot vector are not in increasing order."
                    )
                # we know that knot vector is sorted by now.
                # check if knots are >= 0.
                if kv[0] < 0:
                    raise ValueError(
                        f"{i}. knot vector includes negative knots."
                    )
            
        # Check if required number of control points is present
        n_required_cps = np.prod(self.control_mesh_resolutions)
        n_defined_cps = self.control_points.shape[0]
        if n_required_cps != n_defined_cps:
            raise InputDimensionError(
                "Number of control points invalid: "
                f"expected {n_required_cps}, but given {n_defined_cps}"
            )

        # template splines are instantiated with following naming rules.
        c_spline_cls = (
            f"core.{type(self).__qualname__}{self.para_dim}P{self.dim}D"
        )
        self._c_spline = eval(c_spline_cls)(**dict_spline)

    def _update_p(self,):
        """
        Reads cpp spline and writes it here.
        Probably get an error if cpp isn't ready for this.

        Parameters
        -----------
        None

        Returns
        --------
        None
        """
        required_props = self.required_properties

        # don't use setters here, as it will call `check_and_update_c` each
        # time and it causes us to lose array ref.
        for rp in required_props:
            setattr(self, f"_{rp}", getattr(self._c_spline, rp))

        # but, we still need to do setter's job.
        self._para_dim = self._c_spline.para_dim
        self._dim = self._c_spline.dim
        self._logd(
            "Updated python spline. CPP spline and python spline are"
            "now identical."
        )

    def evaluate(self, queries, n_threads=1):
        """
        Evaluates spline.

        Parameters
        -----------
        queries: (n, para_dim) list-like
        n_threads: int
          Default is 1. If n_thread > 1, it calls multithreading version of
          `evaluate`. Parallel loop implemented in cpp.

        Returns
        --------
        results: (n, dim) np.ndarray
        """
        if self.whatami == "Nothing":
            return None

        queries = utils.make_c_contiguous(queries, dtype="float64")

        if queries.shape[1] != self.para_dim:
            raise InputDimensionError(
                "`queries` does not match current pametric dimension."
            )

        self._logd("Evaluating spline...")

        if int(n_threads) > 1:
            return self._c_spline.p_evaluate(
                queries=queries,
                n_threads=n_threads
            )

        else:
            return self._c_spline.evaluate(queries=queries)

    def derivative(self, queries, orders, n_threads=1):
        """
        Evaluates derivatives of spline.

        Parameters
        -----------
        queries: (n, para_dim) list-like
        orders: (para_dim,) list-like
        n_threads: int
          Default is 1. If n_thread > 1, it calls multithreading version of
          `derivate`. Parallel loop implemented in cpp.

        Returns
        --------
        results: (n, dim) np.ndarray
        """
        if self.whatami == "Nothing":
            return None

        queries = utils.make_c_contiguous(queries, dtype="float64")
        orders = utils.make_c_contiguous(orders, dtype="float64")

        if queries.shape[1] != self.para_dim:
            raise InputDimensionError(
                "`queries` does not match current pametric dimension."
            )
        if orders.size != self.para_dim:
            raise InputDimensionError(
                "`orders` does not match current pametric dimension."
            )

        self._logd("Evaluating derivatives of the spline...")

        if int(n_threads) > 1:
            return self._c_spline.p_derivative(
                queries=queries,
                orders=orders,
                n_threads=n_threads
            )

        else:
            return self._c_spline.derivative(
                queries=queries,
                orders=orders
            )

    def basis_functions(self, queries, n_threads=1):
        """
        Returns basis function values and their support ids of given queries.

        Parameters
        -----------
        queries: (n, para_dim) list-like
        n_threads: int
          Default is 1. Higher number is currently not supported. #TODO: DOIT

        Returns
        --------
        results: tuple
          tuple of two elements.
          first: (prod(degrees .+ 1)) array of basis function values.
          second: support ids.
        """
        if self.whatami == "Nothing":
            return None

        queries = utils.make_c_contiguous(queries, dtype="float64")

        if queries.shape[1] != self.para_dim:
            raise InputDimensionError(
                "`queries` does not match current pametric dimension."
            )

        self._logd("Evaluating basis functions")

        return self._c_spline.basis_functions(queries)

    def insert_knots(self, parametric_dimension, knots):
        """
        Inserts knots.

        Parameters
        -----------
        parametric_dimension: int
        knots: list or float

        Returns
        --------
        None
        """
        if self.whatami == "Nothing":
            return None

        if parametric_dimension >= self.para_dim:
            raise ValueError(
                "Invalid parametric dimension to insert knots."
            )

        if isinstance(knots, float):
            knots = [knots]

        elif isinstance(knots, np.ndarray):
            knots = knots.tolist()

        if not isinstance(knots, list):
            raise TypeError(
                "We couldn't convert input to `list`. Please give us `list`."
            )

        if max(knots) > max(self.knot_vectors[parametric_dimension]):
            raise ValueError(
                "One of the query knots not in valid knot range. (Too big)"
            )

        if min(knots) < min(self.knot_vectors[parametric_dimension]):
            raise ValueError(
                "One of the query knots not in valid knot range. (Too small)"
            )

        self._c_spline.insert_knots(
            int(parametric_dimension),
            knots
        )

        self._logd(f"Inserted {len(knots)} knot(s).")

        self._update_p()

    def remove_knots(self, parametric_dimension, knots, tolerance=1e-8):
        """
        Tries to removes knots. If you've compiled `splinepy` in `Debug`
        and your removal request is not "accepted", you will get an error.
        See the comments for `Nurbs::remove_knots` @
        `splinepy/src/nurbs.hpp` for more info.

        Parameters
        -----------
        parametric_dimension: int
        knots: list or float
        tolerance: float

        Returns
        --------
        None
        """
        if self.whatami == "Nothing":
            return None

        if parametric_dimension >= self.para_dim:
            raise ValueError(
                "Invalid parametric dimension to remove knots."
            )

        if isinstance(knots, float):
            knots = [knots]

        elif isinstance(knots, np.ndarray):
            knots = knots.tolist()

        if not isinstance(knots, list):
            raise TypeError(
                "We couldn't convert input to `list`. Please give us `list`."
            )

        if max(knots) > max(self.knot_vectors[parametric_dimension]):
            raise ValueError(
                "One of the query knots not in valid knot range. (Too big)"
            )

        if min(knots) < min(self.knot_vectors[parametric_dimension]):
            raise ValueError(
                "One of the query knots not in valid knot range. (Too small)"
            )

        total_knots_before = len(self.knot_vectors[int(parametric_dimension)])
        self._c_spline.remove_knots(
            int(parametric_dimension),
            knots,
            tolerance,
        )

        self._update_p()

        self._logd(
            f"Tried to remove {len(knots)} knot(s)."
        )
        self._logd(
            "Actually removed {nk} knot(s).".format(
                nk=(
                    total_knots_before
                    - len(self.knot_vectors[int(parametric_dimension)])
                )
            )
        )

    def normalize_knot_vectors(self,):
        """
        Sets all knot vectors into a range of [0,1], if applicable

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        if "knot_vectors" in self.required_properties:
            for i, kv in enumerate(self.knot_vectors):
                offset = kv[0]
                scale = 1 / (kv[-1] - offset)
                n_kv = [(k - offset) * scale for k in kv]
                self._knot_vectors[i] = n_kv

            self._check_and_update_c()

    def permute_parametric_axes(self, permutation_list, inplace=True):
        """
        Permutates the parametric dimensions

        This function can be used, e.g., to  interchange the parametric
        dimensions xi and eta in order to have them in the right orientation
        (applications in boundary condition definition or mfem export)

        Parameters
        ----------
        permutation_list : list
            New order of parametric dimensions
        inplace: bool
            Default is True. If True, modifies spline inplace, else, returns
            a modified_spline.

        Returns
        -------
        modified_spline : type(self)
            spline with reordered parametric dimensions. iff `inplace=True`.
        """
        # Data collector for new spline object
        spline_data_dict = {}

        # Sanity checks
        if not isinstance(permutation_list, list):
            raise ValueError("Permutation list incomprehensive")
        if not set(range(self.para_dim)) == set(permutation_list):
            raise ValueError("Permutation list invalid")

        self._logd("Permuting parametric axes...")

        # Update knot_vectors where applicable
        if "knot_vectors" in self.required_properties:
            spline_data_dict["knot_vectors"] = [self.knot_vectors[permutation_list[i]]
                                                for i in range(self.para_dim)]
        # Update degrees
        spline_data_dict["degrees"] = [self.degrees[permutation_list[i]]
                                       for i in range(self.para_dim)]

        # Retrieve control mesh resolutions
        ctps_dims = self.control_mesh_resolutions
        new_ctps_dims = [ctps_dims[permutation_list[i]]
                         for i in range(self.para_dim)]
        n_ctps = self.control_points.shape[0]

        # Map global to local index
        # i_glob = i + n_i * j  + n_i * n_j * k ...
        local_indices = np.empty([n_ctps, self.para_dim], dtype=int)
        global_indices = np.arange(n_ctps, dtype=int)
        for i_p in range(self.para_dim):
            local_indices[:, i_p] = global_indices % ctps_dims[i_p]
            global_indices -= local_indices[:, i_p]
            global_indices = np.floor_divide(global_indices, ctps_dims[i_p])

        # Reorder indices
        local_indices[:] = local_indices[:, permutation_list]

        # Rearange global to local
        global_indices = np.matmul(
            local_indices, np.cumprod([1] + new_ctps_dims)[0:-1])
        # Get inverse mapping
        global_indices = np.argsort(global_indices)
        if "weights" in self.required_properties:
            spline_data_dict["weights"] = self.weights[global_indices]
        spline_data_dict["control_points"] = self.control_points[global_indices, :]

        if inplace:
            self._logd("  applying permutation inplace")
            self.clear()
            for rp in self.required_properties:
                setattr(self, rp, spline_data_dict[rp])

            return None

        else:
            self._logd("  returning permuted spline")
            return type(self)(**spline_data_dict)

    def elevate_degree(self, parametric_dimension):
        """
        Elevate degree.

        Parameters
        -----------
        parametric_dimension: int

        Returns
        --------
        None
        """
        if self.whatami == "Nothing":
            return None

        if parametric_dimension >= self.para_dim:
            raise ValueError(
                "Invalid parametric dimension to elevate_degree."
            )

        self._c_spline.elevate_degree(parametric_dimension)
        self._logd(
            f"Elevated {parametric_dimension}.-dim. "
            "degree of the spline."
        )

        self._update_p()

    def reduce_degree(self, parametric_dimension, tolerance=1e-8):
        """
        Tries to reduce degree.

        Parameters
        -----------
        parametric_dimension: int
        tolerance: float

        Returns
        --------
        reduced: bool
        """
        if self.whatami == "Nothing":
            return None

        if parametric_dimension >= self.para_dim:
            raise ValueError(
                "Invalid parametric dimension to reduce degree."
            )

        reduced = self._c_spline.reduce_degree(
            parametric_dimension,
            tolerance
        )

        self._logd(
            f"Tried to reduce {parametric_dimension}.-dim. "
            "degree of the spline."
        )

        if reduced:
            self._logd(
                f"Successfully reduced {parametric_dimension}.-dim. "
                "degree"
            )
            self._update_p()
        else:
            self._logd(
                f"Could not reduce {parametric_dimension}.-dim. "
                "degree"
            )

        return reduced

    def sample(self, query_resolutions):
        """
        Uniformly sample along each parametric dimensions from spline.

        Parameters
        -----------
        query_resolutions: (n,) list-like

        Returns
        --------
        results: (math.product(query_resolutions), dim) np.ndarray
        """
        if self.whatami == "Nothing":
            return None

        query_resolutions = utils.make_c_contiguous(
            query_resolutions,
            dtype="int32",
        ).flatten()

        if query_resolutions.size != self.para_dim:
            raise InputDimensionError(
                "Query resolutions don't match current parametric dimension."
            )

        is_one_or_less = [int(qr) <= 1 for qr in query_resolutions]
        if any(is_one_or_less):
            self._logd(
                "You cannot sample less than 2 points per each "
                "parametric dimension."
            )
            self._logd("Applying minimum sampling resolution 2.")

            query_resolutions[is_one_or_less] = int(2)

        self._logd(
            f"Sampling {np.product(query_resolutions)} "
            "points from spline."
        )

        return self._c_spline.sample(query_resolutions)

    def nearest_pcoord(
            self,
            queries,
            kdt_resolutions=None,
            n_threads=1
    ):
        """
        Given physical coordinate, finds a parametric coordinate that maps to
        the nearest physical coordinate. Also known as "point inversion".
        Makes initial guess by searching for the nearest points with kd-tree.
        With `kdt_resolutions` value, physical points will be sampled to plant
        a kd-tree. Probably useful if you have a lot of queries. The planted
        kd-tree is saved in c_spline internally and will re-plant only if
        `kdt_resolutions` differ from previous function call. Setting any
        entry in `kdt_resolutions` will make sure no tree is planted. Be aware,
        if there's no tree, this will raise runtime error from c++ side.
        Turns out, it is very important to have a good initial guess for
        newton method. Hence, this method is default.

        Parameters
        -----------
        queries: (n, dim) list-like
        kdt_resolutions: (para_dim) list-like
          Default is None and will build [10] * para_dim tree.
          If any entry is negative, it won't plant a tree.
        nthreads: int
          number of threads for parallel execution.

        Returns
        --------
        nearest_pcoord: (n, para_dim) np.ndarray
        """
        if self.whatami == "Nothing":
            return None

        queries = utils.make_c_contiguous(queries, dtype="float64")

        if kdt_resolutions is None:
            self._logd(
                "`kdt_resolutions` is None, "
                "setting default resolution ([10] * para_dim)."
            )
            kdt_resolutions = [10] * self.para_dim

        elif isinstance(kdt_resolutions, (list, np.ndarray)):
            if len(kdt_resolutions) != self.para_dim:
                raise InputDimensionError(
                    "`kdt_resolutions` does not match current para_dim"
                )

        else:
            raise TypeError("`kdt_resolutions` should be list or np.ndarray.")

        if (
            queries.shape[1] != self.dim
            and len(kdt_resolutions) != self.para_dim
        ):
            raise InputDimensionError(
                "`queries` does not match current dimension."
            )

        self._logd("Searching for nearest parametric coord...")

        return self._c_spline.nearest_pcoord_kdt(
            queries=queries,
            resolutions=kdt_resolutions,
            nthreads=n_threads
        )

    def nearest_pcoord_midpoint(self, queries, n_threads=1):
        """
        Given physical coordinate, finds a parametric coordinate that maps to
        the nearest physical coordinate. Also known as "point inversion".
        Initial guess is mid point of parametric space. This tends to fail.
        `nearest_pcoord` is recommended.


        Parameters
        -----------
        queries: (n, dim) list-like
        nthreads: int
          number of threads for parallel execution.

        Returns
        --------
        nearest_pcoord: (n, para_dim) np.ndarray
        """
        if self.whatami == "Nothing":
            return None

        queries = utils.make_c_contiguous(queries, dtype="float64")

        if queries.shape[1] != self.dim:
            raise InputDimensionError(
                "`queries` does not match current dimension."
            )

        self._logd("Searching for nearest parametric coord...")

        return self._c_spline.nearest_pcoord_midpoint(
            queries=queries,
            nthreads=n_threads
        )

    def export(self, fname):
        """
        Export spline. Please be aware of the limits of `.iges`

        Parameters
        -----------
        fname: str

        Returns
        --------
        None
        """
        self._check_and_update_c()

        if self.whatami == "Nothing":
            return None

        fname = str(fname)

        dirname = os.path.dirname(fname)
        if not os.path.isdir(dirname) and dirname != "":
            os.makedirs(dirname)

        ext = os.path.splitext(fname)[1]

        if ext == ".iges":
            self._c_spline.write_iges(fname)

        elif ext == ".xml":
            self._c_spline.write_xml(fname)

        elif ext == ".itd":
            self._c_spline.write_irit(fname)

        elif ext == ".npz":
            io.npz.export(fname, self)

        elif ext == ".json":
            io.json.export(fname, self)

        elif ext == ".mesh":
            # mfem nurbs mesh.
            io.mfem.export(fname, self, precision=12)

        else:
            raise ValueError(
                "We can only export "
                "< .iges | .xml | .itd | .npz | .mesh | .json> extentions"
            )

        self._logi(f"Exported current spline as {fname}.")

    def todict(self, tolist=False):
        """
        Return current spline as dict. Copies all the dict values.
        Does not check current status

        Parameters
        -----------
        tolist: bool
          Default is False. Converts `np.ndarray` into lists.

        Returns
        --------
        dict_spline: dict
        """
        self._logd("Preparing dict_spline...")
        dict_spline = dict()
        # loop and copy entries.
        for p in self.required_properties:
            # no prop? no prob. default is None
            tmp_prop = None

            # prepare copy if there's prop
            if hasattr(self, p):
                tmp_prop = getattr(self, p)
                should_copy = True
                # attr are either list or np.ndarray
                # prepare list if needed.
                if isinstance(tmp_prop, np.ndarray) and tolist:
                    tmp_prop = tmp_prop.tolist()  # copies
                    should_copy = False

                if should_copy:
                    tmp_prop = copy.deepcopy(tmp_prop)
            # update
            dict_spline[p] = tmp_prop

        return dict_spline

    def copy(self,):
        """
        Returns deepcopy of self.

        Parameters
        -----------
        None

        Returns
        --------
        new_spline: type(self)
        """
        # all the properties are deepcopyable
        return copy.deepcopy(self)

    # member function aliases for urgent situations
    # `weights` is defined in `NURBS`
    ds = degrees
    kvs = knot_vectors
    cps = control_points
    ws = weights
