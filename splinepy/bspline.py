import logging
import copy

import numpy as np

from splinepy import utils
from splinepy._splinepy import *
from splinepy._spline import Spline
from splinepy.nurbs import NURBS


class BSpline(Spline):

    def __init__(
            self,
            degrees=None,
            knot_vectors=None,
            control_points=None,
    ):
        """
        BSpline.

        Parameters
        -----------
        degrees: (para_dim,) list-like
        knot_vectors: (para_dim, n) list
        control_points: (m, dim) list-like

        Returns
        --------
        None
        """
        super().__init__(
            degrees=degrees,
            knot_vectors=knot_vectors,
            control_points=control_points,
        )

    def _update_c(self,):
        """
        Updates/Init cpp spline, if it is ready to be updated.
        Checks if all the entries are filled before updating.

        Parameters
        -----------
        None

        Returns
        --------
        None
        """
        if (
            self.dim is None
            or self.para_dim is None
            or self.degrees is None
            or self.knot_vectors is None
            or self.control_points is None
        ):
            logging.debug(
                "Spline - Not enough information to update cpp spline. "
                "Skipping update or removing existing backend spline."
            )
            if hasattr(self, "_c_spline"):
                delattr(self, "_c_spline")
            return None

        c_spline_class = f"BSpline{self.para_dim}P{self.dim}D"
        self._c_spline = eval(c_spline_class)(
            degrees=self.degrees,
            knot_vectors=self.knot_vectors,
            control_points=self.control_points,
        )

        logging.debug(f"Spline - Your spline is {self.whatami}.")

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
        # Don't use setters here, since it calls update_c each time
        # and it causes us to lose array references.
        self._degrees = self._c_spline.degrees
        self._knot_vectors = self._c_spline.knot_vectors
        self._control_points = self._c_spline.control_points
        logging.debug(
            "Spline - Updated python spline. CPP spline and python spline are "
            "now identical."
        )

    def interpolate_curve(
            self,
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
          (Optional) Default is an empty list `[]`. List of `float`s.
           If defined, tries to fit curve with given knot vector.
        save_query: bool
          (Optional) Default is True. Saves query points for plotting, or 
          whatever.

        Returns
        --------
        None
        """
        query_points = utils.make_c_contiguous(query_points, dtype="float64")

        dim = query_points.shape[1]
        c_spline_class = f"BSpline1P{dim}D()"
        c_spline = eval(c_spline_class)
        c_spline.knot_vectors = [knot_vector]
        c_spline.interpolate_curve(
            points=query_points,
            degree=degree,
            centripetal=centripetal
        )
        self._c_spline = c_spline

        logging.debug(
            "Spline - BSpline curve interpolation complete. "
            f"Your spline is {self.whatami}."
        )

        if save_query:
            self._fitting_queries = query_points

        self._update_p()

    def approximate_curve(
            self,
            query_points,
            degree,
            num_control_points,
            centripetal=True,
            knot_vector=[],
            save_query=True,
            return_residual=False,
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
          (Optional) Default is an empty list `[]`. List of `float`s. If
          defined, tries to fit curve with given knot vector.
        save_query: bool
          (Optional) Default is True. Saves query points for plotting, or 
          whatever.
        return_residual: bool
          (Optional) Default is False. Returns Approximation residual.

        Returns
        --------
        res: float
          (Optional) Only returned, if `return_residual` is True.
        """
        query_points = utils.make_c_contiguous(query_points, dtype="float64")

        dim = query_points.shape[1]
        c_spline_class = f"BSpline1P{dim}D()"
        c_spline = eval(c_spline_class)
        c_spline.knot_vectors = [knot_vector]
        res = c_spline.approximate_curve(
            points=query_points,
            degree=degree,
            num_control_points=num_control_points,
            centripetal=centripetal
        )
        self._c_spline = c_spline

        logging.debug(
            "Spline - BSpline curve approximation complete. "
            f"Your spline is {self.whatami}."
        )
        logging.debug(f"Spline -   Approximation residual: {res}")

        if save_query:
            self._fitting_queries = query_points

        self._update_p()

        if return_residual:
            return res

    def interpolate_surface(
            self,
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
        None
        """
        query_points = utils.make_c_contiguous(query_points, dtype="float64")

        dim = query_points.shape[1]
        c_spline_class = f"BSpline2P{dim}D()"
        c_spline = eval(c_spline_class)
        c_spline.interpolate_surface(
            points=query_points,
            size_u=size_u,
            size_v=size_v,
            degree_u=degree_u,
            degree_v=degree_v,
            centripetal=centripetal,
        )
        self._c_spline = c_spline

        logging.debug(
            "Spline - BSpline surface interpolation complete. "
            f"Your spline is {self.whatami}."
        )

        if save_query:
            self._fitting_queries = query_points

        self._update_p()

        # Reorganize control points.
        if reorganize:
            ri = [v + size_v * u for v in range(size_v) for u in range(size_u)]
            self.control_points = self._control_points[ri] 

        if save_query:
            self._fitting_queries = query_points

    @property
    def nurbs(self,):
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
        same_nurbs = NURBS(
            **self.todict(),
            weights=np.ones(self.control_points.shape[0], dtype=np.float64),
        )

        return same_nurbs

    def copy(self,):
        """
        Returns freshly initialized BSpline of self.

        Parameters
        -----------
        None

        Returns
        --------
        new_bspline: `BSpline`
        """
        new_bspline = BSpline(**self.todict())

        # add fitting_queries if it has one.
        if hasattr(self, "_fitting_queries"):
            new_bspline._fitting_queries = copy.deepcopy(self._fitting_queries)

        return new_bspline

    def todict(self):
        """
        Returns copy of degrees, knot_vectors, control_points in dict.

        Parameters
        -----------
        None

        Returns
        --------
        dict_spline: dict
          Keys are {degrees, knot_vectors, control_points}.
        """
        return dict(
            degrees=copy.deepcopy(self.degrees),
            knot_vectors=copy.deepcopy(self.knot_vectors),
            control_points=copy.deepcopy(self.control_points),
        )
