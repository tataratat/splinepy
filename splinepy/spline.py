"""
Abstract Spline
"""
import copy
import os
from functools import wraps

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


class CoordinateReferences(SplinepyBase):
    """
    Helper class to core.CoordinateReferences. Allows direct
    access to cpp's control points.
    Named `coordinate`, since cpp core splines store weighted control points
    for rational splines.
    Extends functionality to provide numpy's convenient __getitem__.
    """

    __slots__ = ("core", "spline", "id_lookup", "apply_weight")

    def __init__(self, core_obj, spline, apply_weight=False):
        """
        Initialize with core.CoordinateReferences and its owner spline.
        Also initialize id_lookup

        Parameters
        ----------
        core_obj: splinepy_core.CoordinateReferences
        spline: Spline
        apply_weight: bool
          Default is False, Applies weight before assigning
        """
        # hold core
        self.core = core_obj

        # remember coordinate's originating spline
        self.spline = spline

        # for __getitem__ trick
        self.id_lookup = np.arange(len(core_obj), dtype=np.int32).reshape(
            -1, spline.dim
        )

        # what do we do with rational splines?
        self.apply_weight = apply_weight

    def __setitem__(self, key, value):
        """
        Sets items. uses key to get global ids from id_lookup.
        Applies weight if needed.

        Parameters
        ----------
        key: types allowed by numpy.ndarray.__setitem__
        value: array-like or float

        Returns
        -------
        None
        """
        ids = self.id_lookup[key].ravel()

        # in case scalar, and no apply_weight, use direct version
        call_setitem = True
        if isinstance(value, (int, float)):
            if not self.apply_weight:
                self.core.broadcast_scalar(ids, value)
                call_setitem = False
            else:
                value = [value] * ids.size

        # currently we will only support matching size assignment
        value = utils.data.enforce_contiguous(value, "float64")

        # for rational splines, apply weight if wanted
        if self.apply_weight and self.spline.is_rational:
            original_value = value  # keep original here.
            value = value.reshape(-1) * self.spline.weights[
                ids // self.spline.dim
            ].reshape(-1)

        # __setitem__ call checks that array size matches
        if call_setitem:
            self.core[ids] = value

        # if there's local TrackedArray copy, update it too
        saved_cps = self.spline._data.get("properties", dict()).get(
            "control_points", None
        )

        # if None, it doesn't have local copy, so return
        if saved_cps is None:
            return None

        # if spline has local changes, warn!
        if is_modified(self.spline):
            self._logw(
                "trying to update coordinates of a spline",
                "that has local changes.",
                "This may cause unexpected behavior of the spline.",
            )

        # take away the weight if weighted is on.
        if self.apply_weight and self.spline.is_rational:
            value = original_value

        # update and set modified flag to false
        saved_cps.reshape(-1)[ids] = value.reshape(-1)
        saved_cps._modified = False

    def set_with_global_ids(self, global_ids, values):
        """
        Sets values using global indices. len(values) and len(global_ids)
        should match.
        Skips overhead of using __getitem__ to retrieve global ids.
        Does NOT apply any weights

        Parameters
        ----------
        global_ids: (n,) array-like
        values: (n,) array-like

        Returns
        -------
        None
        """
        # checks size match
        self.core[global_ids] = values

    def numpy(self):
        """
        Returns copy of current values as numpy.ndarray

        Returns
        -------
        as_numpy: (n, spline.dim) np.ndarray
        """
        return self.core.numpy().reshape(-1, self.spline.dim)


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
    # take a look at saved values. we don't want this call to create
    # a new copy.
    saved_properties = spl._data.get("properties", dict())
    if len(saved_properties) == 0:
        return False

    for rp in spl.required_properties:
        prop = saved_properties.get(rp, None)
        if prop is not None:
            if utils.data.is_modified(prop):
                return True

    return False


def sync_from_core(spl):
    """
    Clears saved data and syncs given spline with its core spline values.
    Similiar to previously called `_update_p`, but does a bit more.
    Meant for internal use for syncing python exposed splines with cpp splines.
    However, could be also useful for debugging or ensuring purpose.
    Syncing is done inplace.

    Parameters
    ----------
    spl: Spline

    Returns
    -------
    synced: Spline
      Synced input spline.
    """
    # assign default data
    spl._data = _default_data()
    if core.have_core(spl):
        # get current core properies and make them trackable
        _make_core_properties_trackable(spl)
        # set modified flags to false to all (inplace change)
        # trackable properties
        _set_modified_false(spl)
    else:
        utils.log.debug("sync_from_core() - no core to sync. Skipping.")

    return spl


def permute_parametric_axes(spline, permutation_list, inplace=True):
    """
    Permutates the parametric dimensions

    This function can be used, e.g., to  interchange the parametric
    dimensions xi and eta in order to have them in the right orientation
    (applications in boundary condition definition or mfem export)

    Parameters
    ----------
    spline: Spline
    permutation_list : list
        New order of parametric dimensions
    inplace: bool
        Default is True. If True, modifies spline inplace, else, returns
        a modified_spline.

    Returns
    -------
    modified_spline : type(spline)
        spline with reordered parametric dimensions. iff `inplace=True`.
    """
    # Data collector for new spline object
    dict_spline = {}

    # Sanity checks
    if not isinstance(permutation_list, list):
        raise ValueError("Permutation list incomprehensive")
    if not set(range(spline.para_dim)) == set(permutation_list):
        raise ValueError("Permutation list invalid")

    utils.log.debug("Permuting parametric axes...")

    # Update knot_vectors where applicable
    if "knot_vectors" in spline.required_properties:
        dict_spline["knot_vectors"] = [
            spline.knot_vectors[permutation_list[i]]
            for i in range(spline.para_dim)
        ]
    # Update degrees
    dict_spline["degrees"] = [
        spline.degrees[permutation_list[i]] for i in range(spline.para_dim)
    ]

    # Retrieve control mesh resolutions
    ctps_dims = spline.control_mesh_resolutions
    new_ctps_dims = [
        ctps_dims[permutation_list[i]] for i in range(spline.para_dim)
    ]
    n_ctps = spline.control_points.shape[0]

    # Map global to local index
    # i_glob = i + n_i * j  + n_i * n_j * k ...
    local_indices = np.empty([n_ctps, spline.para_dim], dtype=int)
    global_indices = np.arange(n_ctps, dtype=int)
    for i_p in range(spline.para_dim):
        local_indices[:, i_p] = global_indices % ctps_dims[i_p]
        global_indices -= local_indices[:, i_p]
        global_indices = np.floor_divide(global_indices, ctps_dims[i_p])

    # Reorder indices
    local_indices[:] = local_indices[:, permutation_list]

    # Rearange global to local
    global_indices = np.matmul(
        local_indices, np.cumprod([1] + new_ctps_dims)[0:-1]
    )
    # Get inverse mapping
    global_indices = np.argsort(global_indices)
    if "weights" in spline.required_properties:
        dict_spline["weights"] = spline.weights[global_indices]
    dict_spline["control_points"] = spline.control_points[global_indices, :]

    if inplace:
        utils.log.debug("  applying permutation inplace")
        spline.new_core(**dict_spline)

        return None

    else:
        utils.log.debug("  returning permuted spline")
        return type(spline)(**dict_spline)


def _default_data():
    """
    Returns default dict used for splines

    Parameters
    ----------
    None

    Returns
    -------
    data_holder: dict()
    """
    return dict(properties=dict())


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
    # empty dict? skip
    if not spl._data.get("properties", dict()):
        return None

    for rp in spl.required_properties:
        prop = getattr(spl, rp, None)

        if prop is not None:
            if isinstance(prop, list):
                for p in prop:
                    p._modified = False

            else:
                prop._modified = False


def _make_core_properties_trackable(spl):
    """
    Internal use only.
    Get current properties and save them as TrackedArray.
    This marks all properties as modified. You could change this using
    `_set_modified_false(spl)`,
    perhaps a combination with `spl._data = _default_data()`.

    Parameters
    ----------
    spl: Spline

    Returns
    -------
    None
    """
    # copy truely current properties
    current_p = spl.current_core_properties()
    # reset properties
    spl._data["properties"] = dict()

    # loop and make them trackable
    for key, values in current_p.items():
        if key.startswith("knot_vectors"):
            kvs = list()
            for kv in values:
                kvs.append(utils.data.make_tracked_array(kv, copy=False))
            spl._data["properties"][key] = kvs

        else:
            spl._data["properties"][key] = utils.data.make_tracked_array(
                values, copy=False
            )


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


def _new_core_if_modified(func):
    """
    Decorator to sync splines before a function call.

    Parameters
    ----------
    func: callable

    Returns
    -------
    inner: callable
    """

    @wraps(func)
    def inner(*args, **kwargs):
        self = args[0]
        # sync if modified
        # removes saved data
        if is_modified(self):
            self.new_core(
                keep_properties=True, raise_=False, **self._data["properties"]
            )

        # now call
        return func(*args, **kwargs)

    return inner


def _get_property(spl, property_):
    """
    Returns spline property.
    syncs spline from core if needed.

    Parameters
    ----------
    spl: Spline
    property_: str

    Returns
    -------
    spline_property: array-like
    """
    p = spl._data.get("properties", dict()).get(property_, None)

    # property exists, return
    if p is not None:
        return p
    # p is none and core does not exist, return None
    elif not core.have_core(spl) or property_ not in spl.required_properties:
        return None

    # core exists, and just doesn't have local copy
    return sync_from_core(spl)._data["properties"][property_]


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
            self._data = _default_data()
            return None

        # current core spline based init if given spline has no local changes
        if spline is not None and isinstance(spline, core.PySpline):
            # if spline is modified, update this spline first
            if is_modified(spline):
                spline.new_core(
                    keep_properties=True,
                    raise_=False,
                    **self._data["properties"],
                )

            # spline based init - takes SplinepyBase, even the nullptr
            super().__init__(spline)

            # don't share the copy of properties
            self._data = _default_data()
            return None

        else:
            # do they at least contain minimal set of keywards?
            kset = set(kwargs.keys())
            if not RequiredProperties.intersection().issubset(kset):
                raise RuntimeError(
                    f"Given keyword arguments ({kwargs}) don't contain"
                    "minimal set of keywords, "
                    f"{RequiredProperties.intersection()}."
                )
            # we will call new_core to make sure all the array values are
            # contiguous
            super().__init__()  # alloc
            self.new_core(**kwargs, raise_=False)

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
        >>> spline_faces = spline.check.valid_queries()

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
        self._data = _default_data()

    def new_core(self, *, keep_properties=False, raise_=True, **kwargs):
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

        # Spline, you do whatever
        if type(self).__qualname__ == "Spline":
            # maybe minimal set check?
            super().new_core(**kwargs)

        # specified ones needs specific sets of kwargs
        # in case of an incomplete set of kwargs, nothing will happen
        elif set(kwargs.keys()).issuperset(set(self.required_properties)):
            super().new_core(**kwargs)

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
            # old behavior was to remove core here. maybe do so?

        # in case no roundtrip is desired, keep properties alive
        props = {}
        if keep_properties:
            props = self._data["properties"]

        # clear saved data
        self._data = _default_data()
        self._data["properties"] = props
        _set_modified_false(self)

        return core.have_core(self)

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
        degrees = _get_property(self, "degrees")
        if degrees is None:
            return degrees

        # if core spline exists,
        # degrees should not be modifiable
        if core.have_core(self):
            degrees.mutable = False

        return degrees

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
            self._data["properties"]["degrees"] = None
            return None

        # len should match knot_vector's len
        if self.knot_vectors is not None:
            if len(self.knot_vectors) != len(degrees):
                raise ValueError(
                    f"len(degrees) ({len(degrees)}) should match "
                    f"len(self.knot_vectors) ({len(self.knot_vectors)})."
                )

        # set - copies
        self._data["properties"]["degrees"] = utils.data.make_tracked_array(
            degrees, "int32"
        )
        self._logd(f"Degrees set: {self.degrees}")

        # try to sync core with current status
        self.new_core(
            keep_properties=True,
            raise_=False,
            **self._data["properties"],
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
        return _get_property(self, "knot_vectors")

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
            self._data["properties"]["knot_vectors"] = None
            return None

        # len should match degrees's len
        if self.degrees is not None:
            if len(knot_vectors) != len(self.degrees):
                raise ValueError(
                    f"len(knot_vectors) ({len(knot_vectors)}) should "
                    f"match len(self.degrees) ({len(self.degrees)})."
                )

        # set - copies
        self._data["properties"]["knot_vectors"] = [
            utils.data.make_tracked_array(kv, "float64") for kv in knot_vectors
        ]

        self._logd("Knot vectors set:")
        for i, kv in enumerate(self.knot_vectors):
            self._logd(f"  {i}" ". knot vector length: " f"{len(kv)}")

        # try to sync core with current status
        self.new_core(
            keep_properties=True,
            raise_=False,
            **self._data["properties"],
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
            return self.parametric_bounds.T

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
        return _get_property(self, "control_points")

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
            self._data["properties"]["control_points"] = None
            return None

        # len should match weights' len
        if self.weights is not None:
            if len(self.weights) != len(control_points):
                raise ValueError(
                    f"len(control_points) ({len(control_points)}) "
                    "should match len(weights) ({len(weights)})."
                )

        # set - copies
        self._data["properties"][
            "control_points"
        ] = utils.data.make_tracked_array(control_points, "float64")
        self._logd(f"{self.control_points.shape[0]} Control points set.")

        # try to sync core with current status
        self.new_core(
            keep_properties=True,
            raise_=False,
            **self._data["properties"],
        )

    @property
    def control_point_bounds(
        self,
    ):
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
    def coordinate_references(self):
        """
        Returns direct reference of underlying cpp coordinates.

        Parameters
        ----------
        None

        Returns
        -------
        coordinate_references: CoordinateReferences
          See docs of CoordinateReferences.
        """
        # see if it is stored
        cr = self._data.get("coordinate_references", None)

        if cr is not None:
            return cr

        # not stored, create one
        self._data["coordinate_references"] = CoordinateReferences(
            super().coordinate_references(), self
        )

        return self._data["coordinate_references"]

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
        return _get_property(self, "weights")

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
            self._data["properties"]["weights"] = None
            return None

        # len should match control_points' len
        if self.control_points is not None:
            if len(self.control_points) != len(weights):
                raise ValueError(
                    f"len(weights) ({len(weights)}) should match "
                    "len(control_points) ({len(control_points)})."
                )

        # set - copies
        self._data["properties"]["weights"] = utils.data.make_tracked_array(
            weights, "float64"
        ).reshape(-1, 1)

        self._logd(f"{self.weights.shape[0]} Weights set.")

        # try to sync core with current status
        self.new_core(
            keep_properties=True,
            raise_=False,
            **self._data["properties"],
        )

    @_new_core_if_modified
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

    @_new_core_if_modified
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
        self._logd(f"Sampling {np.prod(resolutions)} " "points from spline.")

        return super().sample(
            resolutions, nthreads=_default_if_none(nthreads, settings.NTHREADS)
        )

    @_new_core_if_modified
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

    @_new_core_if_modified
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

    @_new_core_if_modified
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

    @_new_core_if_modified
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

    @_new_core_if_modified
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

    @_new_core_if_modified
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

    @_new_core_if_modified
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

    @_new_core_if_modified
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

    @_new_core_if_modified
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
        self._data = _default_data()

    @_new_core_if_modified
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
            self._data = _default_data()

        return reduced

    @_new_core_if_modified
    def extract_boundaries(self, boundary_ids=None):
        """
        Extracts boundary spline.

        The boundaries deducted from the parametric axis which is normal to the
        boundary (j), if the boundary is at parametric axis position x_j=x_jmin
        the corresponding boundary is 2*j, else at parametric axis position
        x_j=x_jmin the boundary is 2*j+1


        Parameters
        -----------
        boundary_ids: array-like
          Boundary IDs with the enumeration described above

        Returns
        -------
        boundary_spline: type(self)
          boundary spline, which has one less para_dim
        """
        # check if boundaries were saved.
        saved_boundaries = self._data.get("boundaries", None)
        if saved_boundaries is not None:
            if boundary_ids is None:
                return saved_boundaries
            else:
                return [saved_boundaries[bi] for bi in boundary_ids]

        # extract boundaries
        boundaries = [
            type(self)(spline=c)
            for c in core.extract_boundaries(
                self, _default_if_none(boundary_ids, [])
            )
        ]

        # save if this was all-query
        if boundary_ids is None:
            self._data["boundaries"] = boundaries

        return boundaries

    @_new_core_if_modified
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
                    tmp_prop = [t.tolist() for t in tmp_prop]

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
            # shallow copy
            shallow = self._data.copy()

            # don't copy coordinate references
            shallow.pop("coordinate_references", None)

            # call deep copy
            new._data = copy.deepcopy(shallow)

            # copy modified status, as it is created with True flag.
            for k in shallow.get("properties", dict()).keys():
                if k.startswith("knot_vectors"):
                    for kv_new, kv in zip(
                        new._data["properties"]["knot_vectors"],
                        shallow["properties"]["knot_vectors"],
                    ):
                        kv_new._modified = kv._modified
                else:
                    new._data["properties"][k]._modified = shallow[
                        "properties"
                    ][k]._modified

        else:
            new._data = _default_data()

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
