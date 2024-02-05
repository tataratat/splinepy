"""
Abstract Spline
"""

import copy as _copy
import os as _os

import numpy as _np

from splinepy import helpme as _helpme
from splinepy import io as _io
from splinepy import settings as _settings
from splinepy import splinepy_core as _core
from splinepy import utils as _utils
from splinepy._base import SplinepyBase as _SplinepyBase
from splinepy.helpme import visualize as _visualize
from splinepy.helpme.check import Checker as _Checker
from splinepy.helpme.create import Creator as _Creator
from splinepy.helpme.extract import Extractor as _Extractor
from splinepy.helpme.integrate import Integrator as _Integrator
from splinepy.utils.data import SplineData as _SplineData
from splinepy.utils.data import cartesian_product as _cartesian_product


class RequiredProperties(_SplinepyBase):
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
    if not _core.has_core(spl):
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

    cps = cps.view(_utils.data.PhysicalSpaceArray)
    cps._source_ptr = coordinate_pointers[0]
    spl._data["control_points"] = cps

    if spl.is_rational:
        ws = ws.view(_utils.data.PhysicalSpaceArray)
        ws._source_ptr = coordinate_pointers[1]
        spl._data["weights"] = ws.reshape(-1, 1)


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

    # if attr doesn't exist, create one and save to its member defined in
    # slots
    helper_obj = helper_class(spl)
    setattr(spl, attr_name, helper_obj)

    return helper_obj


def _safe_new_core(spl, exclude="?"):
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
    # make sure spl._data has all the properties we need
    if _core.has_core(spl):
        for rp in spl.required_properties:
            if not rp.startswith(exclude):
                _ = getattr(spl, rp, None)

    # try to sync with current status
    spl._new_core(
        keep_properties=True,
        raise_=False,
        **spl._data,
    )


def _safe_array_copy(spl, array, key, coupled_key):
    """
    In splinepy we make safe copies. For control_points and weights,
    to support non-contiguous / homogeneous control points in backend,
    we have a helper class called `PhysicalSpaceArray`. This syncs
    numpy array and backend's control points through `ControlPointPointers`.

    This function checks a non-copy criteria: self assignment.
    This tends to happen during inplace operations.
    For example, spline.cps += 1 becomes
    cp_ref = spline.cps; cp_ref +=1; spline.cps = cp_ref.
    In this case, we don't want to create a new core spline.

    Parameters
    ----------
    spl: Spline
    array: np.ndarray
    key: str
    coupled_key: str

    Returns
    -------
    copied: bool
      Whether copy was made
    """
    # is current value same as given array?
    if spl._data.get(key, None) is array:
        return False

    coupled_array = spl._data.get(coupled_key, None)
    if coupled_array is not None:
        # length check - they need to match
        if len(coupled_array) != len(array):
            raise ValueError(
                f"len({key}) ({len(array)}) "
                f"should match len({coupled_key}) ({len(coupled_array)})."
            )

        # for rational splines, control points and weights are coupled
        # -> in backends, we save weighted control points and they must
        #    belong to a same spline.
        # for example, if we are setting new control points,
        # we need to create a copy of weights so that these new control points
        # are coupled with a fresh weight array.
        #
        # Else, the following assert will raise
        # >>> orig_sample = spl.sample(n)
        # >>> spl.cps = new_cps
        # >>> spl.ws[[1,2,3]] *= .5
        # >>> assert not np.allclose(orig_sample, spl.sample(n))
        if isinstance(coupled_array, _utils.data.PhysicalSpaceArray):
            spl._data[coupled_key] = coupled_array.copy()

    # safety copy
    spl._data[key] = _np.array(array, dtype="float64", copy=True, order="C")
    return True


class Spline(_SplinepyBase, _core.PySpline):
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
        "_show_options",
        "_spline_data",
    )

    __show_option__ = _visualize.SplineShowOption

    def __init__(self, spline=None, **kwargs):
        # return if this is an empty init
        if spline is None and len(kwargs) == 0:
            super().__init__()
            return None

        # current core spline based init if given spline has no local changes
        if spline is not None and isinstance(spline, _core.PySpline):
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
        return _Extractor(self)

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
        checker: Checker
        """
        return _Checker(self)

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
        return _Integrator(self)

    @property
    def create(self):
        """Returns spline creator Can be used to create new splines using
        geometric relations.

        Examples
        --------
        >>> prism = spline.create.extrude(axis=[1, 4, 1])

        Parameters
        ----------
        None

        Returns
        -------
        creator: spline.Creator
        """
        return _Creator(self)

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
        return _get_helper(self, "_spline_data", _SplineData)

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
        _core.annul_core(self)
        self._data = {}

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
        kwargs = _utils.data.enforce_contiguous_values(kwargs)

        # type check for knot vectors
        kvs = kwargs.get("knot_vectors", None)
        if kvs is not None:
            kwargs["knot_vectors"] = [
                kv.numpy() if isinstance(kv, _core.KnotVector) else kv
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
            # BSpline supports viewing-control-points for contiguous arrays
            # such as np.ndarray, so we need to keep it alive.
            if self.name.startswith("BSpline"):
                saved_cps = self._data.get("control_points", None)
                self._data = {}
                if saved_cps is not None:
                    self._logw(
                        "_new_core(keep_properties=False) -",
                        "BSplines need to keep control_points.",
                        "Properties excluding control_points will be cleared.",
                    )
                    self._data["control_points"] = saved_cps
            else:
                self._data = {}

        return _core.has_core(self)

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

        if _core.has_core(self):
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
            _core.annul_core(self)
            self._data = {}
            return None

        # len should match knot_vector's len
        if self.knot_vectors is not None and len(self.knot_vectors) != len(
            degrees
        ):
            raise ValueError(
                f"len(degrees) ({len(degrees)}) should match "
                f"len(self.knot_vectors) ({len(self.knot_vectors)})."
            )

        # set - copies
        self._data["degrees"] = _np.asarray(
            degrees, dtype="int32", order="C"
        ).copy()

        # make sure _new_core call works
        _safe_new_core(self, exclude="degrees")

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

        if _core.has_core(self):
            if not self.has_knot_vectors:
                return None

            # check type and early exit if it is core type
            if kvs is not None and isinstance(kvs[0], _core.KnotVector):
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
            _core.annul_core(self)
            self._data = {}
            return None

        # len should match degrees's len
        if self.degrees is not None and len(knot_vectors) != len(self.degrees):
            raise ValueError(
                f"len(knot_vectors) ({len(knot_vectors)}) should "
                f"match len(self.degrees) ({len(self.degrees)})."
            )

        # set - copies
        new_kvs = []
        for kv in knot_vectors:
            if isinstance(kv, _core.KnotVector):
                new_kvs.append(kv.numpy())
            else:
                new_kvs.append(_utils.data.enforce_contiguous(kv, "float64"))

        self._data["knot_vectors"] = new_kvs

        _safe_new_core(self, exclude="knot_vectors")

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
            self._logd("Retrieving unique knots")
            return [k.unique() for k in self.knot_vectors]

    @property
    def knot_multiplicities(self):
        """
        Returns knot multiplicities.
        Does not store results.

        Parameters
        -----------
        None

        Returns
        --------
        knot_multiplicities: list<np.array>
        """
        if "Bezier" in self.name:
            self._logd(
                "Returning multiplicities of knots if Bezier was BSpline"
            )
            return [_np.array([d + 1, d + 1]) for d in self.degrees]

        else:
            self._logd("Retrieving knot multiplicities")
            return [k.multiplicities() for k in self.knot_vectors]

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

        if _core.has_core(self) and not isinstance(
            cps, _utils.data.PhysicalSpaceArray
        ):
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
            _core.annul_core(self)
            self._data = {}
            return None

        # set - copies if it is not the same value
        if not _safe_array_copy(
            self, control_points, "control_points", "weights"
        ):
            return None

        _safe_new_core(self, exclude="control_points")

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

        return _np.vstack((cps.min(axis=0), cps.max(axis=0)))

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
        return _helpme.mapper.Mapper(self, reference=reference)

    def greville_abscissae(self, duplicate_tolerance=None):
        """
        Returns greville abscissae.

        Parameters
        -----------
        duplicate_tolerance : float
          If negative, duplicates are possible for C^(-1) splines (B-Spline
          family). Otherwise they are filtered out using adjacent abscissae.
          Value represents tolerance between two greville abscissae, to be
          considered as equal

        Returns
        --------
        greville_abscissae: (para_dim) np.ndarray
        """
        return _cartesian_product(
            super()._greville_abscissae_list(
                _default_if_none(duplicate_tolerance, -1.0),
            )
        )

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
        self._data["multi_index"] = _helpme.multi_index.MultiIndex(
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

        if _core.has_core(self) and not isinstance(
            ws, _utils.data.PhysicalSpaceArray
        ):
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
            _core.annul_core(self)
            self._data = {}
            return None

        # set - copies if it is not the same value
        # if it is the same, return
        if not _safe_array_copy(self, weights, "weights", "control_points"):
            return None

        _safe_new_core(self, exclude="weights")

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

        queries = _utils.data.enforce_contiguous(
            queries, dtype="float64", asarray=_settings.CHECK_BOUNDS
        )

        if _settings.CHECK_BOUNDS:
            self.check.valid_queries(queries)

        return super().evaluate(
            queries,
            nthreads=_default_if_none(nthreads, _settings.NTHREADS),
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
        from gustaf.utils import arr as _arr

        resolutions = _arr.enforce_len(resolutions, self.para_dim)

        self._logd(f"Sampling {_np.prod(resolutions)} points from spline.")

        return super().sample(
            resolutions,
            nthreads=_default_if_none(nthreads, _settings.NTHREADS),
        )

    def derivative(self, queries, orders, nthreads=None):
        """
        Evaluates derivatives of spline.

        Parameters
        -----------
        queries: (n, para_dim) array-like
        orders: (para_dim,) or (m, para_dim) array-like
        nthreads: int

        Returns
        --------
        results: (n, m, dim) np.ndarray
          Iff m == 1, it will have (n, dim) shape.
        """
        self._logd("Evaluating derivatives of the spline")

        queries = _utils.data.enforce_contiguous(
            queries, dtype="float64", asarray=_settings.CHECK_BOUNDS
        )

        if _settings.CHECK_BOUNDS:
            self.check.valid_queries(queries)

        orders = _utils.data.enforce_contiguous(orders, dtype="int32")

        return super().derivative(
            queries=queries,
            orders=orders,
            nthreads=_default_if_none(nthreads, _settings.NTHREADS),
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

        queries = _utils.data.enforce_contiguous(
            queries, dtype="float64", asarray=_settings.CHECK_BOUNDS
        )

        if _settings.CHECK_BOUNDS:
            self.check.valid_queries(queries)

        return super().jacobian(
            queries,
            nthreads=_default_if_none(nthreads, _settings.NTHREADS),
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
        queries = _utils.data.enforce_contiguous(
            queries, dtype="float64", asarray=_settings.CHECK_BOUNDS
        )

        if _settings.CHECK_BOUNDS:
            self.check.valid_queries(queries)

        return super().support(
            queries=queries,
            nthreads=_default_if_none(nthreads, _settings.NTHREADS),
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
        queries = _utils.data.enforce_contiguous(
            queries, dtype="float64", asarray=_settings.CHECK_BOUNDS
        )

        if _settings.CHECK_BOUNDS:
            self.check.valid_queries(queries)

        return super().basis(
            queries=queries,
            nthreads=_default_if_none(nthreads, _settings.NTHREADS),
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
        queries = _utils.data.enforce_contiguous(
            queries, dtype="float64", asarray=_settings.CHECK_BOUNDS
        )

        if _settings.CHECK_BOUNDS:
            self.check.valid_queries(queries)

        return super().basis_and_support(
            queries=queries,
            nthreads=_default_if_none(nthreads, _settings.NTHREADS),
        )

    def basis_derivative(self, queries, orders, nthreads=None):
        """
        Returns derivative of basis functions evaluated on the supports
        of given queries.

        Parameters
        ----------
        queries: (n, para_dim) array-like
        orders: (para_dim,) or (m, para_dim) array-like
        nthreads: int

        Returns
        --------
        basis_derivatives: (n, m, prod(degrees + 1)) np.ndarray
          Iff m == 1, it will have (n, prod(degrees + 1)) shape.
        """
        self._logd("Evaluating basis function derivatives")
        queries = _utils.data.enforce_contiguous(
            queries, dtype="float64", asarray=_settings.CHECK_BOUNDS
        )

        if _settings.CHECK_BOUNDS:
            self.check.valid_queries(queries)

        orders = _utils.data.enforce_contiguous(orders, dtype="int32")

        return super().basis_derivative(
            queries=queries,
            orders=orders,
            nthreads=_default_if_none(nthreads, _settings.NTHREADS),
        )

    def basis_derivative_and_support(self, queries, orders, nthreads=None):
        """
        Returns derivative of basis functions and their support ids of given
        queries. Same as calling `basis_derivative` and `support` at the same
        time.

        Parameters
        ----------
        queries: (n, para_dim) array-like
        orders: (para_dim,) or (m, para_dim) array-like
        nthreads: int

        Returns
        --------
        basis_derivatives: (n, m, prod(degrees + 1)) np.ndarray
          Iff m == 1, it will have (n, prod(degrees + 1)) shape.
        supports: (n, prod(degrees + 1)) np.ndarray
        """
        self._logd("Evaluating basis function derivatives")
        queries = _utils.data.enforce_contiguous(
            queries, dtype="float64", asarray=_settings.CHECK_BOUNDS
        )

        if _settings.CHECK_BOUNDS:
            self.check.valid_queries(queries)

        orders = _utils.data.enforce_contiguous(orders, dtype="int32")

        return super().basis_derivative_and_support(
            queries=queries,
            orders=orders,
            nthreads=_default_if_none(nthreads, _settings.NTHREADS),
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

        The planted kd-tree is saved in cpp object internally
        and will re-plant if `initial_guess_sample_resolutions` are positive.
        Default behavior is to re-plant the tree each query.
        If you know that spline doesn't have any inplace changes, and would
        like to avoid it, set negative values.

        If there's no tree, this will raise runtime error.

        Parameters
        -----------
        queries: (n, dim) array-like
        initial_guess_sample_resolutions: (para_dim) array-like
          Default is 2 * control_mesh_resolutions. For complex shapes, we
          recommend setting large values.
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

        queries = _utils.data.enforce_contiguous(queries, dtype="float64")

        # set small tolerance.
        if tolerance is None and _settings.TOLERANCE > 1.0e-18:
            tolerance = 1e-18

        verbose_info = super().proximities(
            queries=queries,
            initial_guess_sample_resolutions=_default_if_none(
                initial_guess_sample_resolutions,
                self.control_mesh_resolutions * 2,
            ),
            tolerance=tolerance,
            max_iterations=max_iterations,
            aggressive_search_bounds=aggressive_search_bounds,
            nthreads=_default_if_none(nthreads, _settings.NTHREADS),
        )

        if return_verbose:
            return verbose_info
        else:
            if _np.any(
                verbose_info[4]
                > _default_if_none(tolerance, _settings.TOLERANCE)
            ):
                self._logw(
                    "Proximity search did not converge within the tolerance "
                    "for some queries. Try to rerun proximity search with "
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
        self._data = {}

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
            tolerance=_default_if_none(tolerance, _settings.TOLERANCE),
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
            self._data = {}

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
        dirname = _os.path.dirname(fname)
        if not _os.path.isdir(dirname) and dirname != "":
            _os.makedirs(dirname)

        # ext tells you what you
        ext = _os.path.splitext(fname)[1]

        if ext == ".iges":
            _io.iges.export(fname, [self])

        elif ext == ".xml":
            _io.xml.export(fname, [self])

        elif ext == ".itd":
            _io.irit.export(fname, [self])

        elif ext == ".npz":
            _io.npz.export(fname, self)

        elif ext == ".json":
            _io.json.export(fname, self)

        elif ext == ".mesh":
            # mfem nurbs mesh.
            _io.mfem.export(fname, self, precision=12)

        else:
            raise ValueError(
                "We can only export "
                "< .iges | .xml | .itd | .npz | .mesh | .json> extensions"
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
        if _core.has_core(self):
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

        dict_spline = {}
        # loop and copy entries.
        for p in self.required_properties:
            # no prop? no prob. default is None
            tmp_prop = getattr(self, p, None)
            should_copy = True
            # attr are either list or np.ndarray
            # prepare list if needed.
            if tolist and (tmp_prop is not None):
                if isinstance(tmp_prop, _np.ndarray):
                    tmp_prop = tmp_prop.tolist()  # copies
                    should_copy = False
                if p.startswith("knot_vectors"):
                    tmp_prop = [list(tp) for tp in tmp_prop]
                    should_copy = False

            if should_copy:
                tmp_prop = _copy.deepcopy(tmp_prop)

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

        # new's _data is populated with current_core_properties(), which they
        # own. So we will add to that
        if saved_data:
            required_properties = self.required_properties
            for k, v in self._data.items():
                if k in required_properties:
                    continue
                new._data[k] = _copy.deepcopy(v)

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
        return _visualize.show(self, **kwargs)

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
        return _visualize.show(self, return_showable=True, **kwargs)

    # short cuts / alias
    ds = degrees
    kvs = knot_vectors
    cps = control_points
    ws = weights
