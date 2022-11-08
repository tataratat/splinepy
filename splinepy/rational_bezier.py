import numpy as np

from splinepy import utils
from splinepy.spline import Spline


class RationalBezier(Spline):

    def __init__(self, degrees=None, control_points=None, weights=None):
        """
        RationalBezier (Spline).

        Parameters
        -----------
        degrees: (para_dim,) list-like
        control_points: (m, dim) list-like
        weights: (m) list-like

        Returns
        --------
        None
        """
        super().__init__(
            degrees=degrees,
            control_points=control_points,
            weights=weights
        )

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

    def __mul__(self, factor):
        """
        Overloads multiplication between splines with different types of
        degrees

        The resulting spline fulfills the equation
          factora(t) * factorb(t) = result(t)

        Parameters
        ----------
        factor :  float or spline with compatible dimensionality
        Returns
        -------
        spline : New spline with updated types
        """
        if isinstance(factor, float):
            return RationalBezier(
                control_points=self._control_points * factor,
                degrees=self._degrees,
                weights=self._weights
            )

        elif isinstance(factor, type(self)):
            if len(factor.degrees) == len(self._degrees):
                if factor.dim == self.dim:
                    res_c_spl = self._c_spline.multiply_with_spline(
                        factor._c_spline
                    )
                elif factor.dim == 1:
                    res_c_spl = self._c_spline.multiply_with_scalar_spline(
                        factor._c_spline
                    )
                elif self.dim == 1:
                    res_c_spl = factor._c_spline.multiply_with_scalar_spline(
                        self._c_spline
                    )
                else:
                    raise NotImplementedError()

                # Copy the c spline to the python object
                result = type(self)()
                result._c_spline = res_c_spl
                result._update_p()
                return result
            else:
                raise TypeError(
                    "Multiplication with inequal parametric dimensions"
                )
        else:
            raise TypeError("Multiplication with unknown type requested.")

    def __add__(self, summand):
        """
        Calculates the spline that formes the sum of the summand and the
        current spline (function argument)

        The resulting spline fulfils the equation
          self(t) + summand(t) = result(t)

        Parameters
        ----------
        summand : RationalBezier
          spline with same parametric and physical dimension

        Returns
        -------
        spline : RationalBezier
          New spline that describes sum
        """
        if isinstance(summand, type(self)):
            if (
                summand.para_dim == self.para_dim
                and summand.dim == self.dim
            ):
                # Calculate Spline sum
                resulting_c_spline = self._c_spline.add_spline(
                    summand._c_spline
                )

                # Copy the c spline to the python object
                result = type(self)()
                result._c_spline = resulting_c_spline
                result._update_p()
                return result

            else:
                raise TypeError("Dimension Mismatch.")

        else:
            raise TypeError(
                "Sum must be formed between RationalBezier Splines")

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
            if not inner_function.dim == self.para_dim:
                raise TypeError(
                    "Dimension Mismatch - Physical Dimension of inner spline "
                    "must match parametric dimension of outer spline"
                )
            if inner_function.para_dim == 1:
                res_c_spl = self._c_spline.compose_line_rr(
                    inner_function._c_spline
                )

            elif inner_function.para_dim == 2:
                res_c_spl = self._c_spline.compose_surface_rr(
                    inner_function._c_spline
                )

            elif inner_function.para_dim == 3:
                res_c_spl = self._c_spline.compose_volume_rr(
                    inner_function._c_spline
                )

            else:
                raise TypeError(
                    "Compositions with high"
                    " parametric dimensions not supported."
                )

            # Copy the c spline to the python object
            result = type(self)()
            result._c_spline = res_c_spl
            result._update_p()
            return result

        if inner_function.whatami.startswith("Bezier"):
            if inner_function.para_dim == 1:
                res_c_spl = self._c_spline.compose_line_rp(
                    inner_function._c_spline
                )

            elif inner_function.para_dim == 2:
                res_c_spl = self._c_spline.compose_surface_rp(
                    inner_function._c_spline
                )

            elif inner_function.para_dim == 3:
                res_c_spl = self._c_spline.compose_volume_rp(
                    inner_function._c_spline
                )

            else:
                raise TypeError(
                    "Compositions with high"
                    " parametric dimensions not supported."
                )

            # Copy the c spline to the python object
            result = type(self)()
            result._c_spline = res_c_spl
            result._update_p()
            return result

        else:
            raise TypeError(
                "Composisiton must be formed with RationalBezier Splines"
            )
