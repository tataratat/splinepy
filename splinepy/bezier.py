import numpy as np

from splinepy import utils
from splinepy.spline import Spline

class Bezier(Spline):

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

    def recursive_evaluate(self, queries, n_threads=1):
        """
        Evaluates spline, classique (de Casteljau) way.

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

        self._logd("Evaluating spline...")

        return self._c_spline.recursive_evaluate(queries=queries)

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
          return Bezier(
              control_points=self._control_points * factor,
              degrees=self._degrees
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
        summand : Bezier
          spline with same parametric and physical dimension

        Returns
        -------
        spline : Bezier
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

          else :
              raise TypeError("Dimension Mismatch.")

        else :
            raise TypeError("Sum must be formed between Bezier Splines")

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
                res_c_spl = self._c_spline.compose_line_pp(
                    inner_function._c_spline
                )

            elif inner_function.para_dim == 2:
                res_c_spl = self._c_spline.compose_surface_pp(
                    inner_function._c_spline
                )

            elif inner_function.para_dim == 3:
                res_c_spl = self._c_spline.compose_volume_pp(
                    inner_function._c_spline
                )

            else :
                raise TypeError(
                    "Compositions with high"
                    " parametric dimensions not supported."
                )

            # Copy the c spline to the python object
            result = type(inner_function)()
            result._c_spline = res_c_spl
            result._update_p()
            return result

        if "RationalBezier" in inner_function.whatami:
            if inner_function.para_dim == 1:
                res_c_spl = self._c_spline.compose_line_pr(
                    inner_function._c_spline
                )

            elif inner_function.para_dim == 2:
                res_c_spl = self._c_spline.compose_surface_pr(
                    inner_function._c_spline
                )

            elif inner_function.para_dim == 3:
                res_c_spl = self._c_spline.compose_volume_pr(
                    inner_function._c_spline
                )

            else :
                raise TypeError(
                    "Compositions with high"
                    " parametric dimensions not supported."
                )

            # Copy the c spline to the python object
            result = type(inner_function)()
            result._c_spline = res_c_spl
            result._update_p()
            return result
            

        else :
            raise TypeError("Composisiton must be formed with Bezier Splines")
