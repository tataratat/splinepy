import numpy as np

from splinepy import spline
from splinepy import splinepy_core
from splinepy import settings
from splinepy import utils


class BSplineBase(spline.Spline):
    """BSpline base. Contains extra operations that's only available for
    bspline families.
    """
    __slots__ = ()

    def __init__(self, *args, **kwargs):
        """
        BSplineBase. Serves as a base for bspline families.

        Parameters
        -----------
        *args: args
        **kwargs: kwargs

        Returns
        --------
        None
        """
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
        None
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

        splinepy_core.insert_knots(self, parametric_dimension, knots)

        self._logd(f"Inserted {len(knots)} knot(s).")

        spline.sync_from_core(self)

    def remove_knots(self, parametric_dimension, knots, tolerance=None):
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
                tolerance=spline._default_if_none(
                        tolerance, settings.TOLERANCE
                )
        )

        if any(removed):
            spline.sync_from_core(self)

        self._logd(f"Tried to remove {len(knots)} knot(s).")
        self._logd(f"Actually removed {sum(removed)} knot(s).")

    def normalize_knot_vectors(self, ):
        """
        Sets all knot vectors into a range of [0,1], if applicable

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
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
    """
    BSpline.
    """
    __slots__ = ("_fitting_queries")

    def __init__(
            self,
            degrees=None,
            knot_vectors=None,
            control_points=None,
            spline=None,
    ):
        """
        BSpline. Uses Spline.__init__()

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
            save_query=True
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
        query_points = np.ascontiguousarray(query_points, dtype="float64")

        fitted = cls(
                **splinepy_core.interpolate_curve(
                        points=query_points,
                        degree=degree,
                        centripetal=centripetal,
                        knot_vector=knot_vector,
                )
        )

        utils.log.debug("BSpline curve interpolation complete. ")

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
        query_points = np.ascontiguousarray(query_points, dtype="float64")

        results = splinepy_core.approximate_curve(
                points=query_points,
                degree=degree,
                n_control_points=num_control_points,
                centripetal=centripetal,
                knot_vector=knot_vector,
        )
        fitted = cls(**results)
        res = results["residual"]

        utils.log.debug("BSpline curve approximation complete. ")
        utils.log.debug(f"  Approximation residual: {res}")

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
            reorganize=True,
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
        query_points = np.ascontiguousarray(query_points, dtype="float64")

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

        utils.log._logd("BSpline surface interpolation complete. ")

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
    def nurbs(self, ):
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
                weights=np.ones(
                        self.control_points.shape[0], dtype=np.float64
                ),
        )

        return same_nurbs
