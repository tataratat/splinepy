from splinepy import settings
from splinepy import spline
from splinepy import splinepy_core


class BezierBase(spline.Spline):
    """Bezier Base. Contain extra operations that's only
    available for bezier families.
    """
    __slots__ = ()

    def __init__(self, *args, **kwargs):
        """
        BezierBase. Serves as a base for bezier families.
        Passes all the args to super.__init__
        """
        super().__init__(*args, **kwargs)

    @spline._new_core_if_modified
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
        return settings.NAME_TO_TYPE[multiplied.name](spline=multiplied)

    @spline._new_core_if_modified
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

        return type(self)(spline=added)

    @spline._new_core_if_modified
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

        return settings.NAME_TO_TYPE[composed.name](spline=composed)

    @spline._new_core_if_modified
    def split(self, para_dim, locations):
        """
        Splits spline at given locations along the given para_dim.

        Parameters
        ----------
        para_dim: int
        locations: array-like
          Should be in range of (0, 1)

        Returns
        -------
        splitted: list
          list of splitted splines. Self stays intact.
        """
        if max(locations) > 1 or min(locations) < 0:
            raise ValueError("Invalid split location. Should be in (0, 1).")
        splitted = splinepy_core.split(self, para_dim, locations)

        return [type(self)(spline=s) for s in splitted]

    @spline._new_core_if_modified
    def extract_dim(self, dim):
        """
        Extracts a single physical dimension of a spline.

        Parameters
        ----------
        dim: int

        Returns
        -------
        extracted: type(self)
        """
        if dim >= self.dim:
            raise ValueError(f"Can't extract ({dim}) dim from {self.whatami}.")

        return type(self)(spline=splinepy_core.extract_dim(self, dim))

    @spline._new_core_if_modified
    def composition_derivative(self, inner, inner_derivative):
        """
        Derivative of composition.

        Parameters
        ----------
        inner: BezierBase
        inner_derivative: BezierBase

        Returns
        -------
        composition_der: BezierBase
        """
        # compatibility check is done in cpp
        composition_der = splinepy_core.composition_derivative(
                self, inner, inner_derivative
        )

        return settings.NAME_TO_TYPE[composition_der](spline=composition_der)


class Bezier(BezierBase):

    __slots__ = ()

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
