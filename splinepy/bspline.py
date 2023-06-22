import numpy as np

from splinepy import settings, spline, splinepy_core, utils


class BSplineBase(spline.Spline):
    r"""BSpline base. Contains extra operations that are only available for
    bspline families, i.e., :class:`.BSpline` and :class:`.NURBS`.

    Passes all arguments to :code:`super.__init__()`, see :class:`.Spline`.

    All B-Spline-like splines feature knot vectors for each parametric
    dimension. A knot vector can be considered as a set of non-decreasing
    parametric coordinates. If coordinates in the parametric dimension are
    denoted by :math:`u\in\mathbb{R}^{N_{param}}`, then the knot vector
    associated with this dimension contains values

    .. math::
            \Theta_u = \{u_1,u_2,...,u_{l+p+1}\}

    such that

    .. math::
            u_1 \leq u_2 \leq ... \leq u_{l+p+1}

    :math:`u_i` is referred to as the :math:`i`-th *knot*. Intervals of the
    type :math:`[u_i,u_{i+1}]` are known as *knot spans* or *elements*.
    As we can see, the number of knot entries, the number of control points,
    and the polynomial degree of the spline are related to each other by the
    equation :math:`\# knots=l+p+1`. Here, :math:`l` denotes the number of
    control points along a parametric dimension and :math:`p` the polynomial
    order of the spline.
    In the special case of a so-called *open knot vector*, the first and the
    last entry of the knot vector are repeated :math:`p+1` times, that is

    .. math::
          u_1 = ... = u_{p+1} < u_{p+2} < ... < u_l < u_{l+1} = ... = u_{l+p+1}

    Splines with open knot vectors have the property that they exactly
    interpolate the first and the last control point in this parametric
    dimension.

    Based on the knot vectors, the mono-variate basis functions are defined
    recursively via the Cox-de-Boor recursion formula for each parametric
    dimension:

    .. math::
            N^{i;0}(u) &= \begin{cases}
                1 & u\in [u_i,u_{i+1}] \\\\
                0 & \mathrm{else} \end{cases} \\\\
            N^{i;p}(u) &= \frac{u-u_i}{u_{i+p}-u_i} N^{i;p-1}(u)
                + \frac{u_{i+p+1}-u}{u_{i+p+1}-u_{i+1}}N^{i+1;p-1}(u)

    For explanations on how to construct splines using these
    basis functions and usage examples, we refer to the documentation of the
    classes :class:`.BSpline` and :class:`.NURBS`.
    """

    __slots__ = ()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @spline._new_core_if_modified
    def insert_knots(self, parametric_dimension, knots):
        """
        Inserts knots.

        Parameters
        -----------
        parametric_dimension: int
        knots: list or float

        Returns
        --------
        inserted: list
          List of bool. True if the knots are inserted. Otherwise, False
        """
        if parametric_dimension >= self.para_dim:
            raise ValueError("Invalid parametric dimension to remove knots.")

        if isinstance(knots, float):
            knots = [knots]

        if max(knots) > max(self.knot_vectors[parametric_dimension]):
            raise ValueError(
                "One of the query knots not in valid knot range. (Too big)"
            )

        if min(knots) < min(self.knot_vectors[parametric_dimension]):
            raise ValueError(
                "One of the query knots not in valid knot range. "
                "(Too small)"
            )

        inserted = splinepy_core.insert_knots(
            self, parametric_dimension, knots
        )

        self._logd(f"Inserted {len(knots)} knot(s).")

        self._data = spline._default_data()

        return inserted

    def remove_knots(self, parametric_dimension, knots, tolerance=None):
        """
        Tries to removes knots.

        Parameters
        -----------
        parametric_dimension: int
        knots: list or float
        tolerance: float

        Returns
        --------
        removed: list<bool>
          List of bool. True if the knots are Removed. Otherwise, False
        """
        if parametric_dimension >= self.para_dim:
            raise ValueError("Invalid parametric dimension to remove knots.")

        if isinstance(knots, float):
            knots = [knots]

        if max(knots) > max(self.knot_vectors[parametric_dimension]):
            raise ValueError(
                "One of the query knots not in valid knot range. (Too big)"
            )

        if min(knots) < min(self.knot_vectors[parametric_dimension]):
            raise ValueError(
                "One of the query knots not in valid knot range. "
                "(Too small)"
            )

        removed = splinepy_core.remove_knots(
            self,
            parametric_dimension,
            knots,
            tolerance=spline._default_if_none(tolerance, settings.TOLERANCE),
        )

        if any(removed):
            self._data = spline._default_data()

        self._logd(f"Tried to remove {len(knots)} knot(s).")
        self._logd(f"Actually removed {sum(removed)} knot(s).")

        return removed

    def normalize_knot_vectors(self):
        """
        Sets all knot vectors into a range of [0,1], if applicable

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        # check if knot_vectors is already normalized
        # creating new parametric space in backend can be expensive, so
        # let's avoid it if we can.
        para_bounds = self.parametric_bounds
        kv_modified = utils.data.is_modified(
            self._data.get("properties", dict()).get("knot_vectors", [])
        )
        lower_bounds_are_zero = np.allclose(
            para_bounds[0], [0] * self.para_dim
        )
        upper_bounds_are_one = np.allclose(para_bounds[1], [1] * self.para_dim)
        if not kv_modified and lower_bounds_are_zero and upper_bounds_are_one:
            return None

        new_kvs = []
        for i, kv in enumerate(self.knot_vectors):
            offset = kv[0]
            scale = 1 / (kv[-1] - offset)
            new_kvs.append((kv - offset) * scale)

        # use setter to update
        self.knot_vectors = new_kvs

    def extract_bezier_patches(self):
        """
        Extract all knot spans as Bezier patches to perform further operations
        such as compositions and multiplications

        Parameters
        ----------
        None

        Returns
        -------
        extracted Beziers : list
        """
        # Extract bezier patches and create PyRationalBezier objects
        patches = splinepy_core.extract_bezier_patches(self)

        # use core spline based init and name to type conversion to find
        # correct types
        return [settings.NAME_TO_TYPE[p.name](spline=p) for p in patches]


class BSpline(BSplineBase):
    r"""
    BSpline.

    Passes all arguments to :code:`super.__init__()`, see
    :class:`.BSplineBase`.

    With the help of the mono-variate basis functions introduced from the
    B-Spline introduction in :class:`.BSplineBase`, we can now construct
    B-Splines.

    .. note::
        For simplicity, we restrict ourselves to the three most common types
        of splines in the following, namely curves, surfaces and volumes,
        although :code:`splinepy` also supports higher dimensions, see
        the documentation of :class:`.Spline` for more information.

    1. A B-Spline of degree :math:`p` with control points
    :math:`P^i\in\mathbb{R}^{N_{phys}}` and a one-dimensional parameter space
    (:math:`N_{param}=1`) corresponds to a line embedded into the physical
    space:

    .. math::
            C(u) = \sum_{i=1}^{l} N^{i;p}(u) P^i

    2. A B-Spline of degrees :math:`p,q` with control points
    :math:`P^{i,j} \in \mathbb{R}^{N_{phys}}` and a two-dimensional parameter
    space (:math:`N_{param}=2`) corresponds to a surface, embedded into the
    physical space:

    .. math::
            S(u,v) = \sum_{i=1}^{l} \sum_{j=1}^{m} N^{i;p}(u) N^{j;q}(v)
                P^{i,j}

    Due to the tensor-product nature of the B-Spline basis functions, this is
    often rewritten in terms of multi-variate basis functions

    .. math::
            \tilde{N}_{i,j;p,q}(u,v) := N^{i;p}(u) N^{j;q}(v)

    3. A B-Spline of degrees :math:`p,q,r` with control points
    :math:`P^{i,j,k} \in \mathbb{R}^{N_{phys}}` and a three-dimensional
    parameter space (:math:`N_{param}=3`) corresponds to a volume, embedded
    into the physical space:

    .. math::
            V(u,v,w) = \sum_{i=1}^{l} \sum_{j=1}^{m} \sum_{k=1}^{n}
                N^{i;p}(u) N^{j;q}(v) N^{k;r}(w) P^{i,j,k}

    Here, we can introduce the multi-variate basis functions

    .. math::
            \tilde{N}^{i,j,k;p,q,r}(u,v,w) := N^{i;p}(u) N^{j;q}(v) N^{k;r}(w)

    Higher-dimensional instances are constructed accordingly.

    **Usage**:

    .. code-block:: python

        bspline_volume = splinepy.BSpline(
            degrees=[1,1,1],
            knot_vectors=[
                [0.0, 0.0, 1.0, 1.0],
                [0.0, 0.0, 1.0, 1.0],
                [0.0, 0.0, 1.0, 1.0],
            ],
            control_points=[
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, -1.0, 1.0],
                [1.0, 0.0, 1.0],
                [-1.0, 1.0, 2.0],
                [2.0, 2.0, 2.0],
            ]
        )

    Parameters
    ----------
    degrees: (para_dim) array-like
    knot_vectors: (para_dim) list
        list of array-like
    control_points: (n, dim) array-like
    spline: Spline

    Returns
    -------
    None
    """

    __slots__ = "_fitting_queries"

    def __init__(
        self,
        degrees=None,
        knot_vectors=None,
        control_points=None,
        spline=None,
    ):
        super().__init__(
            spline=spline,
            degrees=degrees,
            knot_vectors=knot_vectors,
            control_points=control_points,
        )

    @classmethod
    def interpolate_curve(
        cls,
        query_points,
        degree,
        centripetal=True,
        knot_vector=[],
        save_query=True,
    ):
        """
        Interpolates BSpline Curve through query points.
        Class method.

        Parameters
        -----------
        query_points: (n, 1-10) list-like
        degree: int
        centripetal: bool
          (Optional) Default is True.
        knot_vector: list
          (Optional) Default is an empty list. List of floats.
           If defined, tries to fit curve with given knot vector.
        save_query: bool
          (Optional) Default is True. Saves query points for plotting, or
          whatever.

        Returns
        --------
        fitted: BSpline
        """
        query_points = utils.data.enforce_contiguous(
            query_points, dtype="float64"
        )

        fitted = cls(
            **splinepy_core.interpolate_curve(
                points=query_points,
                degree=degree,
                centripetal=centripetal,
                knot_vector=knot_vector,
            )
        )

        cls._logd("BSpline curve interpolation complete. ")

        if save_query:
            fitted._fitting_queries = query_points

        return fitted

    @classmethod
    def approximate_curve(
        cls,
        query_points,
        degree,
        num_control_points,
        centripetal=True,
        knot_vector=[],
        save_query=True,
    ):
        """
        Approximates BSpline Curve based on query points.

        Parameters
        -----------
        query_points: (n, 1-10) list-like
        degree: int
        num_control_points: int
          Should be smaller than n. If it is same, the result is same as
          interpolate.
        centripetal: bool
          (Optional) Default is True.
        knot_vector: list
          (Optional) Default is an empty list.  List of floats. If
          defined, tries to fit curve with given knot vector.
        save_query: bool
          (Optional) Default is True. Saves query points for plotting, or
          whatever.

        Returns
        --------
        fitted: BSpline
        """
        query_points = utils.data.enforce_contiguous(
            query_points, dtype="float64"
        )

        results = splinepy_core.approximate_curve(
            points=query_points,
            degree=degree,
            n_control_points=num_control_points,
            centripetal=centripetal,
            knot_vector=knot_vector,
        )

        res = results.pop("residual")

        fitted = cls(**results)

        cls._logd("BSpline curve approximation complete. ")
        cls._logd(f"  Approximation residual: {res}")

        if save_query:
            fitted._fitting_queries = query_points

        return fitted

    @classmethod
    def interpolate_surface(
        cls,
        query_points,
        size_u,
        size_v,
        degree_u,
        degree_v,
        centripetal=True,
        reorganize=False,
        save_query=True,
    ):
        """
        Interpolates BSpline Surface through query points.

        Parameters
        -----------
        query_points: (n, 1-10) list-like
        size_u: int
        size_v: int
        degree_u: int
        degree_v: int
        centripetal: bool
          (Optional) Default is True.
        reorganize: bool
          (Optional) Default is False. Reorganize control points, assuming they
          are listed v-direction first, along u-direction.
        save_query: bool
          (Optional) Default is True. Saves query points for plotting, or
          whatever.

        Returns
        --------
        fitted: BSpline
        """
        query_points = utils.data.enforce_contiguous(
            query_points, dtype="float64"
        )

        fitted = cls(
            **splinepy_core.interpolate_surface(
                points=query_points,
                size_u=size_u,
                size_v=size_v,
                degree_u=degree_u,
                degree_v=degree_v,
                centripetal=centripetal,
            )
        )

        cls._logd("BSpline surface interpolation complete. ")

        if save_query:
            fitted._fitting_queries = query_points

        # Reorganize control points.
        if reorganize:
            ri = [v + size_v * u for v in range(size_v) for u in range(size_u)]
            fitted.control_points = fitted.control_points[ri]

        if save_query:
            fitted._fitting_queries = query_points

        return fitted

    @classmethod
    def approximate_surface(
        cls,
        query_points,
        num_points_u,
        num_points_v,
        size_u,
        size_v,
        degree_u,
        degree_v,
        centripetal=True,
        reorganize=False,
        save_query=True,
    ):
        """
        Approximates the BSpline surface in the least-squares sense through
        query points.

        Parameters
        -----------
        query_points: (num_points_u, num_points_v, 1-10) list-like
          The query points must form a rectangular grid along the x and y-axis.
        num_points_u: int
          The number of sampling points along the first parametric direction.
          By default the first parametric direction is along the cartesian
          x-axis, this can be adapted by reorganize.
        num_points_v: int
          The number of sampling points along the first parametric direction.
        size_u: int
          Number of control points along first parametric direction.
        size_v: int
          Number of control points along second parametric direction.
        degree_u: int
        degree_v: int
        centripetal: bool
          (Optional) Default is True.
        reorganize: bool
          (Optional) Default is False. Reorganize control points, assuming they
          are listed v-direction first, along u-direction.
        save_query: bool
          (Optional) Default is True. Saves query points for plotting, or
          whatever.

        Returns
        --------
        fitted: BSpline
        """
        if query_points.shape[0] != int(num_points_u * num_points_v):
            raise ValueError(
                "The query points do not comply the given " "number of points."
            )
        if (size_u > num_points_u) or (size_v > num_points_v):
            raise ValueError(
                """The number of sampling points must be equal
            or larger than the number of control points in each dimension."""
            )

        query_points = utils.data.enforce_contiguous(
            query_points, dtype="float64"
        )

        fitted = cls(
            **splinepy_core.approximate_surface(
                points=query_points,
                num_points_u=num_points_u,
                num_points_v=num_points_v,
                size_u=size_u,
                size_v=size_v,
                degree_u=degree_u,
                degree_v=degree_v,
                centripetal=centripetal,
            )
        )

        cls._logd("BSpline surface approximation complete. ")

        if save_query:
            fitted._fitting_queries = query_points

        # Reorganize control points.
        if reorganize:
            ri = [v + size_v * u for v in range(size_v) for u in range(size_u)]
            fitted.control_points = fitted.control_points[ri]

        if save_query:
            fitted._fitting_queries = query_points

        return fitted

    @property
    def nurbs(self):
        """
        Returns NURBS version of current BSpline by defining all the weights as
        1.

        Parameters
        -----------
        None

        Returns
        --------
        same_nurbs: NURBS
        """
        same_nurbs = settings.NAME_TO_TYPE["NURBS"](
            **self.todict(),
            weights=np.ones(self.control_points.shape[0], dtype=np.float64),
        )

        return same_nurbs
