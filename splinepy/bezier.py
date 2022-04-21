import logging
import copy

import numpy as np

from splinepy import utils
from splinepy._splinepy import *
from splinepy._spline import Spline

class Bezier(Spline):

    def __init__(self, degrees=None, control_points=None):
        """
        Bezier (Spline).

        Parameters
        -----------
        degrees: (para_dim,) list-like
        control_points: (m, dim) list-like

        Returns
        --------
        None
        """
        super().__init__(
            degrees=degrees,
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
            or self.control_points is None
        ):
            logging.debug(
                "Spline - Not enough information to update cpp spline. "
                "Skipping update or removing existing backend spline."
            )
            if hasattr(self, "_c_spline"):
                delattr(self, "_c_spline")
            return None

        c_spline_class = f"Bezier{self.para_dim}P{self.dim}D"
        self._c_spline = eval(c_spline_class)(
            degrees=self.degrees,
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
        self._control_points = self._c_spline.control_points
        self._dim = self._c_spline.control_points.shape[1]
        self._para_dim = self._c_spline.degrees.shape[0]

        logging.debug(
            "Spline - Updated python spline. CPP spline and python spline are "
            "now identical."
        )

    def recursive_evaluate(self, queries, n_threads=1):
        """
        Evaluates spline, classique way.

        Parameters
        -----------
        queries: (n, para_dim) list-like

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

        logging.debug("Spline - Evaluating spline...")

        return self._c_spline.recursive_evaluate(queries=queries)


    def sample(self, resolution):
        """
        Overwrites basic sample routine by explicitly giving sampling location.

        Parameters
        -----------
        resolution: int or list

        Returns
        --------
        sampled: (prod(resolution), dim) np.ndarray
        """
        if isinstance(resolution, int):
            resolution = [resolution for _ in range(self.para_dim)]

        elif isinstance(resolution, (list or np.ndarray)):
            if len(resolution) != self.para_dim:
                raise ValueError(
                    "Invalid resolution length. Should match para_dim"
                )

        else:
            raise TypeError("resolution only accept list, np.ndarray or int.")

        q = utils.raster_points(
            [
                [0. for _ in range(self.para_dim)],
                [1. for _ in range(self.para_dim)],
            ],
            resolution
        )

        return self.evaluate(q)

    # Operator Overloading

    def __mul__(self, factor):
      """
      Overwrites the * operator to multiply splines with different types of
      data

      The resulting spline fulfils the equation
        factora(t) * factorb(t) = result(t)

      Parameters
      ----------
      factor :  scalar type (float) or spline with same parametric and physical
                dimension or scalar control points
      Returns
      -------
      spline : New spline with updated types
      """
      if isinstance(factor, float):
        return Bezier(
            control_points = self._control_points * factor,
            degrees = self._degrees
          )
      elif isinstance(factor, type(self)):
          if len(factor.degrees) == len(self._degrees):
              if factor.dim == self.dim:
                  result_c_spline = self._c_spline.multiply_with_spline(
                                           factor._c_spline)
              elif factor.dim == 1:
                  result_c_spline = self._c_spline.multiply_with_scalar_spline(
                                           factor._c_spline)
              elif self.dim == 1:
                  result_c_spline = factor._c_spline.\
                                   multiply_with_scalar_spline(self._c_spline)
              else:
                raise NotImplementedError()
              # Copy the c spline to the python object
              result = type(self)()
              result._c_spline = result_c_spline
              result._update_p()
              return result
          else:
              raise TypeError("Multiplication with inequal parametric dimensions")
      else:
          raise TypeError("Multiplication with unknown type requested.")

    def compose(self, inner_function):
      """
      Calculates the spline that formes the composition of the inner function
      spline (function argument), using the caller spline as the outer (or 
      deformation function).

      The resulting spline fulfils the equation
        spline_self(inner_function(t)) = result(t)

      Parameters
      ----------
      inner_function :  spline with parametric dimension <= 3

      Returns
      -------
      spline : New spline that describes composition
      """
      if isinstance(inner_function, type(self)):
          if inner_function.para_dim == 1:
              resulting_c_spline = \
                  self._c_spline.compose_line(inner_function._c_spline)
          elif inner_function.para_dim == 2:
              resulting_c_spline = \
                  self._c_spline.compose_surface(inner_function._c_spline)
          elif inner_function.para_dim == 3:
              resulting_c_spline = \
                  self._c_spline.compose_volume(inner_function._c_spline)
          else :
              raise TypeError("Compositions with high"
                      " parametric dimensions not supported.")          
          # Copy the c spline to the python object
          result = type(self)()
          result._c_spline = resulting_c_spline
          result._update_p()
          return result
      else :
          raise TypeError("Composisiton must be formed with Bezier Splines")
      

    def copy(self):
        """
        """
        pass

    def todict(sef):
        """
        """
        pass
