import numpy as np

from splinepy import utils
from splinepy import settings
from splinepy.spline import Spline
from splinepy import splinepy_core

class BezierBase(Spline):
    """Bezier Base. Contain extra operations that's only
    available for bezier families.
    """
    def __init__(self, *args, **kwargs):
        """
        BezierBase. Serves as a base for bezier families.
        Passes all the args to super.__init__
        """
        super().__init__(*args, **kwargs)

    def __mul__(self, factor):
        """
        Overloads multiplication between splines with different types of
        degrees

        The resulting spline fulfills the equation
          factora(t) * factorb(t) = result(t)

        Parameters
        ----------
        factor :  float or BezierBase
          Spline with compatible dimensionality

        Returns
        -------
        multiplied : Bezier or RationalBezier
          New spline with updated types
        """
        # scalar
        if isinstance(factor, float):
          return type(self)(
              control_points=self.control_points * factor,
              degrees=self.degrees
          )

        # Supports only Bezier families
        if not isinstance(factor, BezierBase):
            raise TypeError(
                    f"Multiplication with {type(factor)}-type is not suppoted."
            )

        # multiply - dimension compatibility is checked in cpp side
        multiplied = splinepy_core.multiply(self, factor)

        # return corresponding type 
        return settings.NAME_TO_TYPE[multiplied.name](multiplied)

    def __add__(self, summand):
        """
        Calculates the spline that formes the sum of the summand and the
        current spline (function argument)

        The resulting spline fulfils the equation
          self(t) + summand(t) = result(t)

        Parameters
        ----------
        summand: type(self) 
          spline with same parametric and physical dimension

        Returns
        -------
        added: type(self)
          New spline that describes sum
        """
        # same type check, same para dim check, same dim check done in cpp
        added = splinepy_core.add(self, summand)

        return type(self)(added)

    def compose(self, inner_function):
        """
        Calculates the spline that formes the composition of the inner function
        spline (function argument), using the caller spline as the outer (or
        deformation function).

        The resulting spline fulfils the equation
          spline_self(inner_function(t)) = result(t)

        Parameters
        ----------
        inner_function : BezierBase

        Returns
        -------
        composed : BezierBase
        """
        # dimension compatibility checked in cpp
        composed = splinepy_core.compose(self, inner_function)

        return settings.NAME_TO_TYPE[composed.name](composed)


class Bezier(BezierBase):

    def __init__(self, degrees=None, control_points=None, spline=None):
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
            spline=spline,
            degrees=degrees,
            control_points=control_points,
        )


