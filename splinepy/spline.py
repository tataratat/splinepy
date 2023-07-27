"""
Abstract Spline
"""
import copy
import os

import numpy as np

from splinepy import helpme, io, settings
from splinepy import splinepy_core as core
from splinepy import utils
from splinepy._base import SplinepyBase
from splinepy.helpme import visualize
from splinepy.helpme.check import Checker
from splinepy.helpme.create import Creator
from splinepy.helpme.extract import Extractor
from splinepy.helpme.integrate import Integrator
from splinepy.utils.data import SplineData


class RequiredProperties(SplinepyBase):
    """
    Helper class to hold required properties of each spline.
    """

    # name mangling to make direct access very annoying
    __required_spline_properties = {
        "Bezier": ("degrees", "control_points"),
        "RationalBezier": ("degrees", "control_points", "weights"),
        "BSpline": ("degrees", "knot_vectors", "control_points"),
        "NURBS": ("degrees", "knot_vectors", "control_points", "weights"),
    }

    # union of all req props. nice for check support
    __union_all = {
        rp for rps in __required_spline_properties.values() for rp in rps
    }

    # intersection of all req props. nice for checking minimal support
    __intersection_all = __union_all.intersection(
        *[set(rps) for rps in __required_spline_properties.values()]
    )

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

        return cls.__required_spline_properties[spline_type]

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
            return cls.__union_all

        rp_union = set()
        for s in splines:
            rp_union.update(cls.of(s))

        return rp_union

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
            return cls.__intersection_all

        rp_intersection = set(cls.of(splines[-1]))
        for s in splines[:-1]:
            rp_intersection.intersection_update(cls.of(s))

        return rp_intersection


def _default_if_none(arg, default):
    """
    Returns default if arg is None

    Parameters
    ----------
    arg: object
    default: object

    Returns
    -------
    default_if_none: object
    """
    return default if arg is None else arg


def _prepare_coordinates(spl):
    """
    Prepares physical space array. Internally called when saved control points
    are  not physicalspacearrays.

    Parameters
    ----------
    spl: Spline

    Returns
    -------
    None
    """
    if not core.has_core(spl):
        return None

    ws = True  # dummy init
    cps = spl._data.get("control_points", None)
    if spl.is_rational:
        ws = spl._data.get("weights", None)

    # make sure we have all the properties we need
    if cps is None or ws is None:
        current_properties = spl.current_core_properties()
        cps = current_properties["control_points"]
        if spl.is_rational:
            ws = current_properties["weights"]

    # get coordinate pointers
    coordinate_pointers = spl._coordinate_pointers()

    cps = cps.view(utils.data.PhysicalSpaceArray)
    cps._source_ptr = coordinate_pointers[0]
    spl._data["control_points"] = cps

    if spl.is_rational:
        ws = ws.view(utils.data.PhysicalSpaceArray)
        ws._source_ptr = coordinate_pointers[1]
        spl._data["weights"] = ws


def _get_helper(spl, attr_name, helper_class):
    """
    Returns a helper object from given spline. If it doesn't exist, it will
    create one first.

    Parameters
    ----------
    spl: Spline
    attr_name: str
    helper_class: type

    Returns
    -------
    helper: object
    """
    attr = getattr(spl, attr_name, None)
    if attr is not None:  # hasattr
        return attr

    # don't have attr - create and set
    helper_obj = helper_class(spl)
    setattr(spl, attr_name, helper_obj)

    return helper_obj


def _call_required_properties(spl, exclude="?"):
    """
    Calls all getters of spline's required_properties.
    This guarantees that spline's saved data (spl._data) contains local
    copy/reference of core spline's properties, as getters properly processes
    and saves spline properties.

    This is useful to call within the setters before calling spl._new_core.
    For example, spline manipulation calls like spl.elevate_degrees() deletes
    all the saved data as most of the saved data is invalid afterwards.
    If you then call a setter, it will set a new value and naturally try to
    create a new spline using cached properties in spl._data, which is
    incomplete at this point.

    Parameters
    ----------
    spl: Spline
    exclude: str
      Default is "?" as no member name can start with a question mark.

    Returns
    -------
    None
    """
    for rp in spl.required_properties:
        if not rp.startswith(exclude):
            _ = getattr(spl, rp, None)


class Spline(SplinepyBase, core.PySpline):
    r"""
    Spline base class. Extends :class:`.PySpline` with documentation.

    Generally, all types of splines can be seen as mappings from a
    :math:`N_{param}`-dimensional parametric domain
    :math:`\Omega_{param} \subset \mathbb{R}^{N_{param}}` to a
    :math:`N_{phys}`-dimensional physical domain
    :math:`\Omega_{phys} \subset \mathbb{R}^{N_{phys}}`, i.e.,

    .. math::
            M : \Omega_{param} \to \Omega_{phys}

    :code:`splinepy` generally supports different combinations of embeddings,
    that is different combinations of dimensionalities :math:`N_{param}` and
    :math:`N_{phys}`, up to dimension 10.
    The supported spline types differ in the way how the spline function
    :math:`M` is defined. For an overview of the mathematical theories for the
    different types of supported splines, we refer to the documentation of the
    classes :class:`.BSplineBase` and :class:`.BezierBase` for a more in-depth
    discussion of the theory on B-Spline/Bezier families, as well as to their
    children :class:`.BSpline`, :class:`.NURBS`, :class:`.Bezier`, and
    :class:`.RationalBezier` for usage examples.

    Parameters
    -----------
    spline: Spline
        Initialize using another spline. Will SHARE core spline.
    degrees: (para_dim,) array-like
        Keyword only parameter.
    knot_vectors: (para_dim, n) list of array-like
        Keyword only parameter.
    control_points: (m, dim) array-like
        Keyword only parameter.
    weights: (m,) array-like
        Keyword only parameter.

    Returns
    -------
    None
    """

    __slots__ = (
        "_extractor",
        "_checker",
        "_creator",
        "_integrator",
        "_show_options",
        "_spline_data",
    )

    __show_option__ = visualize.SplineShowOption

    def __init__(self, spline=None, **kwargs):
        # return if this is an empty init
        if spline is None and len(kwargs) == 0:
            super().__init__()
            return None

        # current core spline based init if given spline has no local changes
        if spline is not None and isinstance(spline, core.PySpline):
            # spline based init - takes SplinepyBase, even the nullptr
            super().__init__(spline)
            # this should keep all the cps and kv references if there's any
            # make shallow copy
            self._data = spline._data.copy()
            return None

        else:
            super().__init__()  # alloc

            # use setters to check minimal requirements and enforce
            # contiguous array
            # in depth checks are done in core
            req_properties = self.required_properties
            if type(self).__qualname__.startswith("Spline"):
                # for Spline class, which can take anything,
                # we want to set it in this specific order that
                # it won't keep creating core
                set_order = (
                    "weights",
                    "control_points",
                    "knot_vectors",
                    "degrees",
                )
            else:
                set_order = req_properties

            for to_set in set_order:
                if to_set in kwargs:
                    setattr(self, to_set, kwargs[to_set])

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
        try:
            return RequiredProperties.of(self)
        except ValueError:
            # Spline
            self._logd(
                f"`required_properties` of {type(self)} is undefined."
                "Returning maximal set of properties."
            )
            return RequiredProperties.union()

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
        return super().para_dim

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
        return super().dim

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
        return super().whatami

    @property
    def name(self):
        """
        Name of the spline provided by cpp side as a str.
        Should be same as type(spline).__qualname__ if you are using
        a specified spline, i.e., not `Spline`.

        Parameters
        -----------
        None

        Returns
        -------
        name: str
        """
        return super().name

    @property
    def has_knot_vectors(self):
        """
        Returns True iff spline has knot vectors. Bezier splines don't.

        Parameters
        ----------
        None

        Returns
        -------
        has_knot_vectors: bool
        """
        return super().has_knot_vectors

    @property
    def is_rational(self):
        """
        Returns True iff spline is rational. NURBS is rational, for example.

        Parameters
        ----------
        None

        Returns
        -------
        is_rational: bool
        """
        return super().is_rational

    @property
    def extract(self):
        """Returns spline extractor. Can directly perform extractions available
        at `splinepy/helpme/extract.py`.

        Examples
        ---------
        >>> spline_faces = spline.extract.faces()

        Parameters
        -----------
        None

        Returns
        --------
        extractor: Extractor
        """
        return _get_helper(self, "_extractor", Extractor)

    @property
    def check(self):
        """Returns spline checker. Can directly perform type and value checks
        at `splinepy/helpme/check.py`.

        Examples
        ---------
        >>> spline.check.valid_queries(queries)

        Parameters
        -----------
        None

        Returns
        --------
        extractor: Extractor
        """
        return _get_helper(self, "_checker", Checker)

    @property
    def integrate(self):
        """Returns spline integrator. Can directly perform integrations
        available at `splinepy/helpme/integrate.py`.

        Examples
        ---------

        .. code-block :: python

          spline_faces = spline.integrate.volume()

        Parameters
        -----------
        None

        Returns
        --------
        integrator: Integrator
        """
        return _get_helper(self, "_integrator", Integrator)

    @property
    def create(self):
        """Returns spline creator Can be used to create new splines using
        geometric relations.

        Examples
        --------
        >>> prism = spline.create.extrude(axis=[1,4,1])

        Parameters
        ----------
        None

        Returns
        -------
        creator: spline.Creator
        """
        return _get_helper(self, "_creator", Creator)

    @property
    def show_options(self):
        """
        Show option manager for splines.

        Parameters
        ----------
        None

        Returns
        -------
        show_options: SplineShowOption
        """
        return _get_helper(self, "_show_options", self.__show_option__)

    @property
    def spline_data(self):
        """
        Spline data helper for splines.

        Parameters
        ----------
        None

        Returns
        -------
        spline_data: SplineData
        """
        return _get_helper(self, "_spline_data", SplineData)

    def clear(self):
        """
        Clears core spline and saved data

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        # annulment of core
        core.annul_core(self)
        self._data = dict()

    def _new_core(self, *, keep_properties=False, raise_=True, **kwargs):
        """
        Creates a new core spline based on given kwargs.
        Clears any saved data afterwards.

        Parameters
        -----------
        keep_properties: bool
          Default is False. If True, keeps current properties. Otherwise,
          all internal temporary data will be deleted, including current
          properties.
        raise_: bool
          Default is True. Raises if kwargs don't contain required
          properties.
        **kwargs: kwargs
          Keyword only argument. Takes properties.

        Returns
        -------
        core_created: bool
        """
        # let's remove None valued items and make values contiguous array
        kwargs = utils.data.enforce_contiguous_values(kwargs)

        # type check for knot vectors
        kvs = kwargs.get("knot_vectors", None)
        if kvs is not None:
            kwargs["knot_vectors"] = [
                kv.numpy() if isinstance(kv, core.KnotVector) else kv
                for kv in kvs
            ]

        # Spline, you do whatever
        if type(self).__qualname__ == "Spline":
            # maybe minimal set check?
            super()._new_core(**kwargs)

        # specified ones needs specific sets of kwargs
        # in case of an incomplete set of kwargs, nothing will happen
        elif set(kwargs.keys()).issuperset(set(self.required_properties)):
            rp = self.required_properties
            rp_dict = {}
            for k, v in kwargs.items():
                if k in rp:
                    rp_dict[k] = v

            super()._new_core(**rp_dict)

        elif raise_:
            raise RuntimeError(
                f"{kwargs.keys()} is not valid set of keywords to create "
                f"a new core for {type(self).__qualname__}."
                f"Valid ones should include: {self.required_properties}."
            )

        else:
            self._logd(
                "Couldn't initialize a new_core with given keywords."
                f"Keywords without None values are {kwargs.keys()}."
            )

        # clear saved data
        if not keep_properties:
            self._data = dict()

        return core.has_core(self)

    @property
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
        ds = self._data.get("degrees", None)

        if core.has_core(self):
            if ds is None:
                ds = self._current_core_degrees()
                ds.flags.writeable = False
                self._data["degrees"] = ds
            else:
                ds.flags.writeable = False

        return ds

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
        # clear core and return
        if degrees is None:
            core.annul_core(self)
            self._data = dict()
            return None

        # len should match knot_vector's len
        if self.knot_vectors is not None:
            if len(self.knot_vectors) != len(degrees):
                raise ValueError(
                    f"len(degrees) ({len(degrees)}) should match "
                    f"len(self.knot_vectors) ({len(self.knot_vectors)})."
                )

        # set - copies
        self._data["degrees"] = np.asarray(
            degrees, dtype="int32", order="C"
        ).copy()

        # make sure _new_core call works
        if core.has_core(self):
            _call_required_properties(self, exclude="degrees")

        # try to sync core with current status
        self._new_core(
            keep_properties=True,
            raise_=False,
            **self._data,
        )

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
        kvs = self._data.get("knot_vectors", None)

        if core.has_core(self):
            if not self.has_knot_vectors:
                return None

            if kvs is not None:
                # check type and early exit if it is core type
                if isinstance(kvs[0], core.KnotVector):
                    return kvs

            # get one from the core
            kvs = [self._knot_vector(i) for i in range(self.para_dim)]
            self._data["knot_vectors"] = kvs

        return kvs

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
        # clear core and return
        if knot_vectors is None:
            core.annul_core(self)
            self._data = dict()
            return None

        # len should match degrees's len
        if self.degrees is not None:
            if len(knot_vectors) != len(self.degrees):
                raise ValueError(
                    f"len(knot_vectors) ({len(knot_vectors)}) should "
                    f"match len(self.degrees) ({len(self.degrees)})."
                )

        # set - copies
        new_kvs = []
        for kv in knot_vectors:
            if isinstance(kv, core.KnotVector):
                new_kvs.append(kv.numpy())
            else:
                new_kvs.append(utils.data.enforce_contiguous(kv, "float64"))

        self._data["knot_vectors"] = new_kvs

        # make sure _new_core call works
        if core.has_core(self):
            _call_required_properties(self, exclude="knot_vectors")

        # try to sync core with current status
        self._new_core(
            keep_properties=True,
            raise_=False,
            **self._data,
        )

    @property
    def unique_knots(self):
        """
        Returns unique knots.
        Does not store results.

        Parameters
        -----------
        None

        Returns
        --------
        unique_knots: list
        """
        if "Bezier" in self.name:
            self._logd(
                "Returning parametric_bounds as "
                "Bezier spline's unique knots."
            )
            return list(self.parametric_bounds.T)

        else:
            self._logd("Computing unique knots")
            unique_knots = list()
            for kv in self.knot_vectors:
                unique_knots_per_dim = [kv[0]]
                for k, d in zip(kv[1:], np.diff(kv)):
                    if d > settings.TOLERANCE:
                        unique_knots_per_dim.append(k)
                unique_knots.append(np.array(unique_knots_per_dim))

            return unique_knots

    @property
    def parametric_bounds(self):
        """
        Returns bounds of parametric_space.

        Parameters
        -----------
        None

        Returns
        --------
        parametric_bounds: (2, para_dim) np.ndarray
        """
        return super().parametric_bounds

    @property
    def control_points(self):
        """
        Returns control points.

        Parameters
        -----------
        None

        Returns
        --------
        control_points_: (n, dim) np.ndarray
        """
        cps = self._data.get("control_points", None)

        if core.has_core(self):
            if not isinstance(cps, utils.data.PhysicalSpaceArray):
                _prepare_coordinates(self)

        return self._data.get("control_points", None)

    @control_points.setter
    def control_points(self, control_points):
        """
        Set control points.

        Parameters
        -----------
        control_points: (n, dim) array-like

        Returns
        --------
        None
        """
        # clear core and return
        if control_points is None:
            core.annul_core(self)
            self._data = dict()
            return None

        # len should match weights' len
        if self.weights is not None:
            if len(self.weights) != len(control_points):
                raise ValueError(
                    f"len(control_points) ({len(control_points)}) "
                    f"should match len(weights) ({len(self.weights)})."
                )
            # we need to remove exising weights so that pointers don't get
            # mixed
            if isinstance(self.weights, utils.data.PhysicalSpaceArray):
                self._data["weights"] = self._data["weights"].copy()

        # set - copies
        if isinstance(control_points, utils.data.PhysicalSpaceArray):
            # don't want to have two arrays updating
            # same cps
            self._data["control_points"] = control_points.copy(order="C")
        else:
            self._data["control_points"] = np.asarray(
                control_points, dtype="float64", order="C"
            ).copy()

        # make sure _new_core call works
        if core.has_core(self):
            _call_required_properties(self, exclude="control_points")

        # try to sync core with current status
        self._new_core(
            keep_properties=True,
            raise_=False,
            **self._data,
        )

    @property
    def control_point_bounds(self):
        """
        Returns bounds of control points.

        Parameters
        -----------
        None

        Returns
        --------
        control_point_bounds: (2, dim) np.ndarray
        """
        cps = self.control_points

        return np.vstack((cps.min(axis=0), cps.max(axis=0)))

    @property
    def control_mesh_resolutions(self):
        """
        Returns control mesh resolutions.

        Parameters
        -----------
        None

        Returns
        --------
        control_mesh_resolutions: (para_dim) np.ndarray
        """
        return super().control_mesh_resolutions

    def mapper(self, reference):
        """Retrieve a mapper that can be used to get physical derivatives such
        as a gradient or hessian in physical space

        Parameters
        ----------
        reference : spline
          Spline that represents the geometry of the field

        Returns
        -------
        mapper : Mapper
          Mapper to calculate physical gradients and hessians
        """
        return helpme.mapper.Mapper(self, reference=reference)

    @property
    def greville_abscissae(self):
        """
        Returns greville abscissae.

        Parameters
        -----------
        None

        Returns
        --------
        greville_abscissae: (para_dim) np.ndarray
        """
        return super().greville_abscissae

    @property
    def multi_index(self):
        """
        Easy control point / weights access using (unraveled) multi index.
        Useful for getting indices if you want to "slice" a control mesh,
        or find ids that relate to a certain boundary.

        Parameters
        ----------
        None

        Returns
        -------
        multi_index_helper: MultiIndex
        """
        mi = self._data.get("multi_index", None)
        if mi is not None:
            return mi

        # not stored, create one
        # this is okay to be copied together in copy()
        self._data["multi_index"] = helpme.multi_index.MultiIndex(
            self.control_mesh_resolutions
        )

        return self._data["multi_index"]

    @property
    def weights(self):
        """
        Returns weights.

        Parameters
        -----------
        None

        Returns
        --------
        self._weights: (n, 1) array-like
        """
        ws = self._data.get("weights", None)

        if core.has_core(self):
            if not isinstance(ws, utils.data.PhysicalSpaceArray):
                _prepare_coordinates(self)

        return self._data.get("weights", None)

    @weights.setter
    def weights(self, weights):
        """
        Set weights.

        Parameters
        -----------
        weights: (n,) array-like

        Returns
        --------
        None
        """
        # clear core and return
        if weights is None:
            core.annul_core(self)
            self._data = dict()
            return None

        # len should match control_points' len
        if self.control_points is not None:
            if len(self.control_points) != len(weights):
                raise ValueError(
                    f"len(weights) ({len(weights)}) should match "
                    f"len(control_points) ({len(self.control_points)})."
                )

            # we need to remove exising cps so that pointers don't get mixed
            if isinstance(self.control_points, utils.data.PhysicalSpaceArray):
                self._data["control_points"] = self._data[
                    "control_points"
                ].copy()

        # set - copies
        if isinstance(weights, utils.data.PhysicalSpaceArray):
            # don't want to have two arrays updating
            # same cps
            self._data["weights"] = weights.copy(order="C").reshape(-1, 1)
        else:
            self._data["weights"] = (
                np.asarray(weights, dtype="float64", order="C")
                .copy()
                .reshape(-1, 1)
            )

        # make sure _new_core call works
        if core.has_core(self):
            _call_required_properties(self, exclude="weights")

        # try to sync core with current status
        self._new_core(
            keep_properties=True,
            raise_=False,
            **self._data,
        )

    def evaluate(self, queries, nthreads=None):
        """
        Evaluates spline.

        Parameters
        -----------
        queries: (n, para_dim) array-like
        n_threads: int

        Returns
        --------
        results: (n, dim) np.ndarray
        """
        self._logd("Evaluating spline")

        queries = utils.data.enforce_contiguous(
            queries, dtype="float64", asarray=settings.CHECK_BOUNDS
        )

        if settings.CHECK_BOUNDS:
            self.check.valid_queries(queries)

        return super().evaluate(
            queries,
            nthreads=_default_if_none(nthreads, settings.NTHREADS),
        )

    def sample(self, resolutions, nthreads=None):
        """
        Uniformly sample along each parametric dimensions from spline.

        Parameters
        -----------
        resolutions: (n,) array-like
        nthreads: int

        Returns
        --------
        results: (math.product(resolutions), dim) np.ndarray
        """
        from gustaf.utils import arr

        resolutions = arr.enforce_len(resolutions, self.para_dim)

        self._logd(f"Sampling {np.prod(resolutions)} " "points from spline.")

        return super().sample(
            resolutions, nthreads=_default_if_none(nthreads, settings.NTHREADS)
        )

    def derivative(self, queries, orders, nthreads=None):
        """
        Evaluates derivatives of spline.

        Parameters
        -----------
        queries: (n, para_dim) array-like
        orders: (para_dim,) or (n, para_dim) array-like
        nthreads: int

        Returns
        --------
        results: (n, dim) np.ndarray
        """
        self._logd("Evaluating derivatives of the spline")

        queries = utils.data.enforce_contiguous(
            queries, dtype="float64", asarray=settings.CHECK_BOUNDS
        )

        if settings.CHECK_BOUNDS:
            self.check.valid_queries(queries)

        orders = utils.data.enforce_contiguous(orders, dtype="int32")

        return super().derivative(
            queries=queries,
            orders=orders,
            nthreads=_default_if_none(nthreads, settings.NTHREADS),
        )

    def jacobian(self, queries, nthreads=None):
        """
        Evaluates jacobians on spline.

        Parameters
        -----------
        queries: (n, para_dim) array-like
        n_threads: int

        Returns
        --------
        results: (n, para_dim, dim) np.ndarray
        """
        self._logd("Determining spline jacobians")

        queries = utils.data.enforce_contiguous(
            queries, dtype="float64", asarray=settings.CHECK_BOUNDS
        )

        if settings.CHECK_BOUNDS:
            self.check.valid_queries(queries)

        return super().jacobian(
            queries,
            nthreads=_default_if_none(nthreads, settings.NTHREADS),
        )

    def support(self, queries, nthreads=None):
        """
        Returns basis support ids of given queries.

        Parameters
        -----------
        queries: (n, para_dim) array-like
        n_threads: int

        Returns
        --------
        support: (n, prod(degrees + 1)) np.ndarray
        """
        self._logd("Evaluating support ids")
        queries = utils.data.enforce_contiguous(
            queries, dtype="float64", asarray=settings.CHECK_BOUNDS
        )

        if settings.CHECK_BOUNDS:
            self.check.valid_queries(queries)

        return super().support(
            queries=queries,
            nthreads=_default_if_none(nthreads, settings.NTHREADS),
        )

    def basis(self, queries, nthreads=None):
        """
        Returns basis function values on the supports of given queries.

        Parameters
        -----------
        queries: (n, para_dim) array-like
        n_threads: int

        Returns
        --------
        basis: (n, prod(degrees + 1)) np.ndarray
        """
        self._logd("Evaluating basis functions")
        queries = utils.data.enforce_contiguous(
            queries, dtype="float64", asarray=settings.CHECK_BOUNDS
        )

        if settings.CHECK_BOUNDS:
            self.check.valid_queries(queries)

        return super().basis(
            queries=queries,
            nthreads=_default_if_none(nthreads, settings.NTHREADS),
        )

    def basis_and_support(self, queries, nthreads=None):
        """
        Returns basis function values and their support ids of given queries.
        Same as calling `basis` and `support` at the same time.

        Parameters
        -----------
        queries: (n, para_dim) array-like
        n_threads: int

        Returns
        --------
        basis: (n, prod(degrees + 1)) np.ndarray
        support: (n, prod(degrees + 1)) np.ndarray
        """
        self._logd("Evaluating basis functions and support")
        queries = utils.data.enforce_contiguous(
            queries, dtype="float64", asarray=settings.CHECK_BOUNDS
        )

        if settings.CHECK_BOUNDS:
            self.check.valid_queries(queries)

        return super().basis_and_support(
            queries=queries,
            nthreads=_default_if_none(nthreads, settings.NTHREADS),
        )

    def basis_derivative(self, queries, orders, nthreads=None):
        """
        Returns derivative of basis functions evaluated on the supports
        of given queries.

        Parameters
        ----------
        queries: (n, para_dim) array-like
        orders: (para_dim,) or (n, para_dim) array-like
        nthreads: int

        Returns
        --------
        basis_derivatives: (n, prod(degrees + 1)) np.ndarray
        """
        self._logd("Evaluating basis function derivatives")
        queries = utils.data.enforce_contiguous(
            queries, dtype="float64", asarray=settings.CHECK_BOUNDS
        )

        if settings.CHECK_BOUNDS:
            self.check.valid_queries(queries)

        orders = utils.data.enforce_contiguous(orders, dtype="int32")

        return super().basis_derivative(
            queries=queries,
            orders=orders,
            nthreads=_default_if_none(nthreads, settings.NTHREADS),
        )

    def basis_derivative_and_support(self, queries, orders, nthreads=None):
        """
        Returns derivative of basis functions and their support ids of given
        queries. Same as calling `basis_derivative` and `support` at the same
        time.

        Parameters
        ----------
        queries: (n, para_dim) array-like
        orders: (para_dim,) or (n, para_dim) array-like
        nthreads: int

        Returns
        --------
        basis_derivatives: (n, prod(degrees + 1)) np.ndarray
        supports: (n, prod(degrees + 1)) np.ndarray
        """
        self._logd("Evaluating basis function derivatives")
        queries = utils.data.enforce_contiguous(
            queries, dtype="float64", asarray=settings.CHECK_BOUNDS
        )

        if settings.CHECK_BOUNDS:
            self.check.valid_queries(queries)

        orders = utils.data.enforce_contiguous(orders, dtype="int32")

        return super().basis_derivative_and_support(
            queries=queries,
            orders=orders,
            nthreads=_default_if_none(nthreads, settings.NTHREADS),
        )

    def proximities(
        self,
        queries,
        initial_guess_sample_resolutions=None,
        tolerance=None,
        max_iterations=-1,
        aggressive_search_bounds=False,
        nthreads=None,
        return_verbose=False,
    ):
        """
        Given physical coordinate, finds a parametric coordinate that maps to
        the nearest physical coordinate. Also known as "point inversion".

        Makes initial guess by searching for the nearest points with kd-tree.
        With `initial_guess_sample_resolutions` value,
        physical points will be sampled to plant a kd-tree.

        Probably useful if you have a lot of queries. The planted
        kd-tree is saved in c_spline internally and will re-plant only if
        `initial_guess_sample_resolutions` differ from previous function call.

        Setting any entry in `initial_guess_sample_resolutions` will make sure
        no trees are planted.

        If there's no tree, this will raise runtime error from c++ side.

        Parameters
        -----------
        queries: (n, dim) array-like
        initial_guess_sample_resolutions: (para_dim) array-like
        tolerance: float
          Convergence criteria. Currently for both distance and residual
        max_iterations: int
          Default is (para_dim * 20)
        aggressive_search_bounds: bool
          Default is False.
          Set search bound to aabb of direct neighbors from the sample grid.
        nthreads: int
        return_verbose : bool
          If False, returns only parametric coords

        Returns
        --------
        para_coord: (n, para_dim) np.ndarray
          Parametric coordinates
        phys_coord: (n, dim) np.ndarray
          (only if return_verbose) respective physical coordinates
        phys_diff: (n, dim) np.ndarray
          (only if return_verbose) respective cartesian difference to query
        distance: (n, 1) np.ndarray
          (only if return_verbose) respective 2-norm of difference
        convergence_norm: (n, 1) np.ndarrray
          (only if return_verbose) Newton residual
        first_derivatives: (n, para_dim, dim) np.ndarray
          (only if return_verbose) Convergence information
        second_derivatives: (n, para_dim, para_dim, dim) np.ndarray
          (only if return_verbose) Convergence information
        """
        self._logd("Searching for nearest parametric coord")

        queries = utils.data.enforce_contiguous(queries, dtype="float64")

        verbose_info = super().proximities(
            queries=queries,
            initial_guess_sample_resolutions=_default_if_none(
                initial_guess_sample_resolutions,
                self.control_mesh_resolutions * 2,
            ),
            tolerance=_default_if_none(tolerance, settings.TOLERANCE),
            max_iterations=max_iterations,
            aggressive_search_bounds=aggressive_search_bounds,
            nthreads=_default_if_none(nthreads, settings.NTHREADS),
        )

        if return_verbose:
            return verbose_info
        else:
            if np.any(
                verbose_info[3]
                > _default_if_none(tolerance, settings.TOLERANCE)
            ):
                self._logw(
                    "Remaining distance of some of the points is larger than "
                    "tolerance, please rerun proximity search with "
                    "return_verbose=True to get more information."
                )
            return verbose_info[0]

    def elevate_degrees(self, parametric_dimensions):
        """
        Elevate degree.

        Parameters
        -----------
        parametric_dimensions: list

        Returns
        --------
        None
        """
        super().elevate_degrees(para_dims=parametric_dimensions)
        self._logd(
            f"Elevated {parametric_dimensions}.-dim. " "degree of the spline."
        )
        self._data = dict()

    def reduce_degrees(self, parametric_dimensions, tolerance=None):
        """
        Tries to reduce degree.

        Parameters
        -----------
        parametric_dimensions: list
        tolerance: float

        Returns
        --------
        reduced: list
        """
        reduced = super().reduce_degrees(
            para_dims=parametric_dimensions,
            tolerance=_default_if_none(tolerance, settings.TOLERANCE),
        )

        def meaningful(r):
            """True/False to reduced/failed"""
            return "reduced" if r else "failed"

        self._logd(
            f"Tried to reduce degrees for {parametric_dimensions}.-dims. "
            "Results: ",
            f"{[meaningful(r) for r in reduced]}.",
        )

        if any(reduced):
            self._data = dict()

        return reduced

    def export(self, fname):
        """
        Export spline. Please be aware of the limits of `.iges`.
        vtk exports are excluded here because it is not a true spline export.

        Parameters
        -----------
        fname: str

        Returns
        --------
        None
        """
        # prepare export destination
        fname = str(fname)
        dirname = os.path.dirname(fname)
        if not os.path.isdir(dirname) and dirname != "":
            os.makedirs(dirname)

        # ext tells you what you
        ext = os.path.splitext(fname)[1]

        if ext == ".iges":
            io.iges.export(fname, [self])

        elif ext == ".xml":
            io.xml.export(fname, [self])

        elif ext == ".itd":
            io.irit.export(fname, [self])

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
        Does not check current status.

        Parameters
        -----------
        tolist: bool
          Default is False. Converts `np.ndarray` into lists.

        Returns
        --------
        dict_spline: dict
        """
        self._logd("Preparing dict_spline")
        if core.has_core(self):
            core_properties = self.current_core_properties()
            if tolist:
                new_dict = {}
                for k, v in core_properties.items():
                    if k.startswith("knot_vectors"):
                        new_dict[k] = [vv.tolist() for vv in v]
                    else:
                        new_dict[k] = v.tolist()
                return new_dict
            else:
                return core_properties

        dict_spline = dict()
        # loop and copy entries.
        for p in self.required_properties:
            # no prop? no prob. default is None
            tmp_prop = getattr(self, p, None)
            should_copy = True
            # attr are either list or np.ndarray
            # prepare list if needed.
            if tolist and (tmp_prop is not None):
                if isinstance(tmp_prop, np.ndarray):
                    tmp_prop = tmp_prop.tolist()  # copies
                    should_copy = False
                if p.startswith("knot_vectors"):
                    tmp_prop = [[t for t in tp] for tp in tmp_prop]
                    should_copy = False

            if should_copy:
                tmp_prop = copy.deepcopy(tmp_prop)

            # update
            dict_spline[p] = tmp_prop

        return dict_spline

    def copy(self, saved_data=True):
        """
        Returns deepcopy of stored data and newly initialized self.

        Parameters
        -----------
        saved_data: bool
          Default is True. calls deepcopy on saved data excluding
          coordinate_references

        Returns
        --------
        new_spline: type(self)
        """
        new = type(self)(**self.current_core_properties())
        if saved_data:
            required_properties = self.required_properties
            saved_copy = {}
            for k, v in self._data.items():
                if k in required_properties:
                    continue
                saved_copy[k] = copy.deepcopy(v)

            new._data = saved_copy

        else:
            new._data = dict()

        return new

    def show(self, **kwargs):
        """Equivalent to

        .. code-block:: python

           splinepy.helpme.visualize.show(spline, **kwargs)

        Parameters
        ----------
        kwargs: kwargs
          see splinepy.helpme.visualize.show

        Returns
        -------
        showable_or_plotter: dict or vedo.Plotter
        """
        return visualize.show(self, **kwargs)

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
        return visualize.show(self, return_showable=True, **kwargs)

    # short cuts / alias
    ds = degrees
    kvs = knot_vectors
    cps = control_points
    ws = weights

    # Deprecated - TODO: inform
    knot_vector_bounds = parametric_bounds
