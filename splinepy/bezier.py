from splinepy import settings, spline, splinepy_core


class BezierBase(spline.Spline):
    """Bezier Base. Contain extra operations that's only
    available for bezier families.

    Beziers are a special type of spline where the basis functions have
    a global support. The basis functions are given by the so-called
    *Berstein polynomials*

    .. math::
            N_{i;p}(u) = \\binom{p}{i} u^i (1-u)^{p-i}

    A Bezier spline itself can be seen as a map from the parametric domain
    :math:`[0,1]^{N_{param}}` to the physical domain
    :math:`\\mathbb{R}^{N_{phys}}`. We can distinguish between different types
    of splines depending on the dimension of the parametric space.

    #. A spline of degree :math:`p` with control points
    :math:`P_i\\in\\mathbb{R}^{N_{phys}}` and a one-dimensional parameter space
    corresponds to a line embedded into the physical space:

    .. math::
            C(u) = \\sum_{i=0}^{l} N_{i;p}(u) P_i

    #. A spline of degrees :math:`p,q` with control points
    :math:`P_{i,j}\\in\\mathbb{R}^{N_{phys}}` and a two-dimensional parameter
    space corresponds to a surface, embedded into the physical space:

    .. math::
            C(u,v) = \\sum_{i=0}^{l} \\sum_{j=0}^{m} N_{i;p}(u) N_{j;q}(v)
                P_{i,j}

    Due to the tensor-product nature of the Bezier basis functions, this is
    often rewritten as in terms of multi-variate basis functions

    .. math::
            \\tilde{N}_{i,j;p,q}(u,v) := N_{i;p}(u) N_{j;q}(v)

    #. A spline of degrees :math:`p,q,r` with control points
    :math:`P_{i,j,k}\\in\\mathbb{R}^{N_{phys}}` and a three-dimensional
    parameter space corresponds to a volume, embedded into the physical space:

    .. math::
            C(u,v,w) = \\sum_{i=0}^{l} \\sum_{j=0}^{m} \\sum_{k=0}^{n}
                N_{i;p}(u) N_{j;q}(v) N_{k;r}(w) P_{i,j,k}

    Here, we can introduce the multi-variate basis functions

    .. math::
            \\tilde{N}_{i,j,k;p,q,r}(u,v,w) \\coloneqq N_{i;p}(u) N_{j;q}(v)
                N_{k;r}(w)

    Finally, we like to emphasize that Bezier splines can be seen as a special
    type of B-Splines, where the knot vectors do not feature internal knots
    but only have the first and last entry :math:`p^i+1`-times repeated.
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
                degrees=self.degrees,
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

        return settings.NAME_TO_TYPE[composition_der.name](
            spline=composition_der
        )


class Bezier(BezierBase):
    __slots__ = ()

    def __init__(self, degrees=None, control_points=None, spline=None):
        """
        Bezier (Spline).

        See :class:`.BezierBase` for more information on the theory of Bezier
        splines.

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
