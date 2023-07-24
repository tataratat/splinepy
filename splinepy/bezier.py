import numpy as np

from splinepy import settings, spline, splinepy_core


class BezierBase(spline.Spline):
    r"""Bezier Base. Contains extra operations that are only
    available for bezier families.

    Passes all arguments to :code:`super.__init__()`, see :class:`.Spline`.

    Beziers are a special type of spline where the basis functions have
    a global support. The basis functions are given by the so-called
    *Berstein polynomials*

    .. math::
            B^{i;p}(u) = \binom{p}{i} u^i (1-u)^{p-i}

    The parametric domain of a Bezier spline is always a multi-dimensional
    hypercube, that is :math:`\Omega_{param}=[0,1]^{N_{param}}`.
    We refer to the documentation of the classes :class:`.Bezier` and
    :class:`.RationalBezier` for explanations how to construct splines using
    Berstein polynomials as basis.

    Moreover, we like to make the point that Bezier splines can also be seen
    as a special type of B-Splines with open knot vectors (i.e., the first and
    last entry are repeated :math:`p+1`-times) that do not feature additional
    internal knots.

    For usage examples, we once again refer to the derived classes
    :class:`.Bezier` and :class:`.RationalBezier`.
    """

    __slots__ = ()

    def __init__(self, *args, **kwargs):
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

    def compose(self, inner_function, compute_sensitivities=False):
        r"""
        Calculates the spline that forms the composition of the inner function
        spline (function argument), using the caller spline as the outer (or
        deformation function).

        Given an outer function :math:`\mathcal{T}` and an inner function
        (spline) :math:`\mathcal{K}`, this function returns a new spline
        function :math:`\mathcal{M}` that represents the composition:

        .. math:: \mathcal{M}(\mathbf{u}) =
                  \mathcal{T}\big(\mathcal{K}(\mathbf{u})\big)
                  \quad \forall \mathbf{u}

        If the flag `compute_sensitivities` is set, compose also returns the
        derivative of the composed spline with respect to the outer function's
        control point positions. This corresponds to the basis function
        representation (which is given by a scalar valued spline).
        This functionality differs from :func:`.composition_derivative`, which
        computes the derivatives concerning the inner functions geometric
        parametrization.

        Given an outer function :math:`\mathcal{T}` with control points
        :math:`C_i` and an inner function :math:`\mathcal{K}`. This function
        computes:

        .. math:: \frac{\partial\mathcal{T}(\mathcal{K})}{\partial C_i}
                  (\mathbf{u}) =
                  B_i\big(\mathcal{K}\big)(\mathbf{u}) \quad \forall \mathbf{u}

        where :math:`B_i` represents the basis function associated to control
        point :math:`C_i`.


        Parameters
        ----------
        inner_function : BezierBase
          Bezier-Type Spline that represents inner function
        compute_sensitivities : bool
          Flag for return values. If set to true returns also sensitivities,
          default False

        Returns
        -------
        composed : BezierBase
          Composed function
        sensitivities : list (optional)
          List all splines that represents the derivatives with respect to
          internal control points
        """
        # dimension compatibility checked in cpp
        composed = splinepy_core.compose(self, inner_function)

        if compute_sensitivities:
            composed_sensivities = splinepy_core.compose_sensitivities(
                self, inner_function
            )

            return (
                settings.NAME_TO_TYPE[composed.name](spline=composed),
                [
                    settings.NAME_TO_TYPE[cc.name](spline=cc)
                    for cc in composed_sensivities
                ],
            )
        else:
            return settings.NAME_TO_TYPE[composed.name](spline=composed)

    def composition_derivative(self, inner, inner_derivative):
        r"""
        Derivative of composition when given the differentiated inner function
        with constant outer function. This function differs from
        :func:`.compose_sensitivities`, which computes the derivatives
        concerning the outer function's control point positions.

        Given an outer function :math:`\mathcal{T}` and an inner function
        :math:`\mathcal{K}` with its derivative with respect to some design
        variable :math:`\alpha`, i.e.,
        :math:`\frac{\partial \mathcal{K}}{\partial \alpha}`, this function
        returns the derivative of the composed geometry in the form
        :math:`\frac{\partial \mathcal{T}(\mathcal{K})}{\partial\alpha}` by
        computing.

        .. math:: \frac{\partial\mathcal{T}(\mathcal{K})}{\partial\alpha} =
                  \frac{\partial \mathcal{T}}{\partial u}(\mathcal{K})\cdot
                  \frac{\partial \mathcal{K}}{\partial \alpha}

        Here :math:`u` refers to the parametric coordinates of the deformation
        function.

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


class Bezier(BezierBase):
    r"""
    Bezier (Spline).

    Passes all arguments to :code:`super.__init__()`, see :class:`.BezierBase`.

    We can distinguish between different types of splines depending on the
    dimension of the parametric space.

    .. note::
        For simplicity, we restrict ourselves to the three most common types
        of splines in the following, namely curves, surfaces and volumes,
        although :code:`splinepy` also supports higher dimensions, see
        the documentation of :class:`.Spline` for more information.

    1. A spline of degree :math:`p` with :math:`(p+1)` control points
    :math:`P^i\in\mathbb{R}^{N_{phys}}` and a one-dimensional parameter space
    (i.e., :math:`N_{param}=1`) corresponds to a line embedded into the
    physical space:

    .. math::
            C(u) = \sum_{i=0}^{l} B^{i;p}(u) P^i

    2. A spline of degrees :math:`p,q` with :math:`(p+1)\cdot(q+1)` control
    points :math:`P^{i,j}\in\mathbb{R}^{N_{phys}}` and a two-dimensional
    parameter space (i.e., :math:`N_{param}=2`) corresponds to a surface,
    embedded into the physical space:

    .. math::
            S(u,v) = \sum_{i=0}^{l} \sum_{j=0}^{m} B^{i;p}(u) B^{j;q}(v)
                P^{i,j}

    Due to the tensor-product nature of the Bezier basis functions, this is
    often rewritten as in terms of multi-variate basis functions

    .. math::
            \tilde{B}^{i,j;p,q}(u,v) := B^{i;p}(u) B^{j;q}(v)

    3. A spline of degrees :math:`p,q,r` with :math:`(p+1)\cdot(q+1)\cdot(r+1)`
    control points :math:`P^{i,j,k}\in\mathbb{R}^{N_{phys}}` and a
    three-dimensional parameter space (i.e., :math:`N_{param}=3`) corresponds
    to a volume, embedded into the physical space:

    .. math::
            V^B(u,v,w) = \sum_{i=0}^{l} \sum_{j=0}^{m} \sum_{k=0}^{n}
                B^{i;p}(u) B^{j;q}(v) B^{k;r}(w) P^{i,j,k}

    Here, we can introduce the multi-variate basis functions

    .. math::
            \tilde{B}^{i,j,k;p,q,r}(u,v,w) := B^{i;p}(u) B^{j;q}(v)
                B^{k;r}(w)

    Higher-dimensional instances are constructed accordingly.

    **Usage**:

    .. code-block:: python

        # Polynomial Bezier surface
        polynomial_bezier = splinepy.Bezier(
            degrees=[2,1],
            control_points=[
                [0.0, 0.0],
                [1.0, 0.0],
                [2.0, 1.0],
                [0.0, 2.0],
                [1.0, 1.0],
                [2.0, 2.0]
            ]
        )

    Parameters
    -----------
    degrees: (para_dim,) list-like
    control_points: (m, dim) list-like

    Returns
    --------
    None
    """

    __slots__ = ()

    def __init__(self, degrees=None, control_points=None, spline=None):
        super().__init__(
            spline=spline,
            degrees=degrees,
            control_points=control_points,
        )

    @property
    def bezier(self):
        return self.copy(saved_data=False)

    @property
    def rationalbezier(self):
        return settings.NAME_TO_TYPE["RationalBezier"](
            **self.todict(), weights=np.ones((self.cps.shape[0], 1))
        )

    @property
    def bspline(self):
        return settings.NAME_TO_TYPE["BSpline"](
            spline=splinepy_core.same_spline_with_knot_vectors(self)
        )

    @property
    def nurbs(self):
        return settings.NAME_TO_TYPE["NURBS"](
            spline=splinepy_core.same_spline_with_knot_vectors(
                self.rationalbezier
            )
        )
