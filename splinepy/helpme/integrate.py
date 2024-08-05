from functools import wraps as _wraps

import numpy as _np

from splinepy._base import SplinepyBase as _SplinepyBase
from splinepy.utils.data import cartesian_product as _cartesian_product
from splinepy.utils.data import has_scipy as _has_scipy

if _has_scipy:
    from scipy.sparse import dok_matrix as _dok_matrix
    from scipy.sparse.linalg import spsolve as _spsolve


def _get_integral_measure(spline):
    """
    Determines the appropriate measure to be used in integration

    If the spline dimension matches its parametric dimension it will return a
    Callable in the form

    .. math::
        \\mathcal{J}_S = det(\\mathbf(J))

    If the physical dimension is bigger then the paramtric dimension it will
    return

    .. math::
        \\mathcal{J}_S = \big( det(\\mathbf(J)^T \\mathbf(J)\big)^{0.5}

    Parameters
    ----------
    spline : Spline / Multipatch
      For parametric and physical dimension

    Returns
    -------
    measure : Callable
      single patch only
    """
    # Check dimensionality
    if spline.dim == spline.para_dim:

        def measure(spline_patch, positions):
            return _np.linalg.det(spline_patch.jacobian(positions))

        return measure

    elif spline.dim > spline.para_dim:

        def measure(spline_patch, positions):
            jacs = spline_patch.jacobian(positions)
            return _np.sqrt(
                _np.linalg.det(_np.einsum("oji,ojk->oik", jacs, jacs))
            )

        return measure

    else:
        raise ValueError("`Volume` not supported if para_dim > dim")


def _default_quadrature_orders(spline):
    expected_degree = spline.degrees * spline.para_dim + 1
    if spline.is_rational:
        expected_degree += 2
        spline._logd("Integration on rational spline is only approximation")

    # Gauss-Legendre is exact for polynomials 2*n-1
    return _np.ceil((expected_degree + 1) * 0.5).astype("int")


def _get_quadrature_information(spline, orders=None):
    """
    Select appropriate integration order (gauss-legendre)

    Determinante of a polynomial spline with para_dim==dim has degree

    .. math::
        p_i^{det} = n_{dim} \\cdot p_i - 1

    cf. [Mantzaflaris et al., 2017,
    DOI:http://dx.doi.org/10.1016/j.cma.2016.11.013]

    The same order approximation will also be used for the metric

    .. math::
        \\mathcal{J}_S = \big( det(\\mathbf(J)^T \\mathbf(J)\big)^{0.5}

    Parameters
    ----------
    spline : Spline
      Spline for integration
    orders : array-like (optional)
      Orders along every parametric dimension

    Returns
    -------
    positions : np.ndarray
      quadrature position in unit-square
    weights : np.ndarray
      quadrature weights
    """

    # Determine integration points
    positions = []
    weights = []

    # Determine integration orders
    if orders is None:
        quad_orders = _default_quadrature_orders(spline)

    else:
        quad_orders = _np.ascontiguousarray(orders, dtype=int).flatten()
        if quad_orders.size != spline.para_dim:
            raise ValueError(
                "Integration order must be array of size para_dim"
            )

    for order in quad_orders:
        # Get legendre quadratuer points
        pos, ws = _np.polynomial.legendre.leggauss(order)
        # from (-1,1) to (0,1)
        positions.append(pos * 0.5 + 0.5)
        weights.append(ws * 0.5)

    # summarize integration points
    return (
        _cartesian_product(positions),
        _cartesian_product(weights).prod(axis=1),
    )


def volume(spline, orders=None):
    r"""Compute volume of a given spline

    Parameters
    ----------
    spline : Spline
        (self if called via integrator)
      splinepy - spline type
    orders : array-like (optional)
      order for gauss quadrature

    Returns
    -------
    volume : float
      Integral of dim-dimensional object
    """
    from splinepy.spline import Spline as _Spline

    # Check i_nput type
    if not isinstance(spline, _Spline):
        raise NotImplementedError("integration only works for splines")

    # Retrieve aux info
    meas = _get_integral_measure(spline)
    positions, weights = _get_quadrature_information(spline, orders)

    # Calculate Volume
    if spline.has_knot_vectors:
        volume = _np.sum(
            [
                _np.sum(meas(b, positions) * weights)
                for b in spline.extract.beziers()
            ]
        )
    else:
        volume = _np.sum(meas(spline, positions) * weights)
    return volume


def parametric_function(
    spline,
    function,
    orders=None,
):
    """Integrate a function defined within the parametric domain

    Parameters
    ----------
    spline : Spline
        (self if called via integrator)
    function : Callable
    orders : optional

    Returns
    -------
    integral : np.ndarray
    """
    from splinepy.spline import Spline as _Spline

    # Check i_nput type
    if not isinstance(spline, _Spline):
        raise NotImplementedError("integration only works for splines")

    # Retrieve aux info
    meas = _get_integral_measure(spline)
    positions, weights = _get_quadrature_information(spline, orders)

    # Calculate Volume
    if spline.has_knot_vectors:
        # positions must be mapped into each knot-element
        para_view = spline.create.parametric_view(axes=False)

        # get initial shape
        initial = function([positions[0]])
        result = _np.zeros(initial.shape[1])
        for bezier_element in para_view.extract.beziers():
            quad_positions = bezier_element.evaluate(positions)
            result += _np.einsum(
                "i...,i,i->...",
                function(quad_positions),
                meas(spline, quad_positions),
                weights,
                optimize=True,
            )

    else:
        result = _np.sum(
            function(positions) * meas(spline, positions) * weights, axis=1
        )
    return result


def physical_function(
    function,  # noqa ARG001
    orders,  # noqa ARG001
):
    """Integrate a function defined within the physical domain"""
    raise NotImplementedError(
        "Function not implemented yet. Please feel free to write an issue, if "
        "you need it: github.com/tatarata/splinepy/issues"
    )


class Integrator:
    """Helper class to integrate some values on a given spline geometry

    Examples
    --------

    .. code-block:: python

        splinepy.helpme.integrate.volume()
        # Equivalent to
        spline.integrate.volume()

    Parameters
    ----------
    spline : Spline
      Spline parent
    """

    __slots__ = ("_helpee",)

    def __init__(self, spl):
        self._helpee = spl

    @_wraps(volume)
    def volume(self, *args, **kwargs):
        return volume(self._helpee, *args, **kwargs)

    @_wraps(parametric_function)
    def parametric_function(self, *args, **kwargs):
        return parametric_function(self._helpee, *args, **kwargs)


class Transformation:
    __slots__ = (
        "_spline",
        "_solution_field",
        "_mapper",
        "_para_dim",
        "_ukv",
        "_n_elems",
        "_quad_positions",
        "_quad_weights",
        "_grid_ids",
        "_all_supports",
        "_all_element_quad_points",
        "_all_jacobians",
        "_all_jacobian_inverses",
        "_all_jacobian_determinants",
    )

    def __init__(self, spline, solution_field=None, orders=None):
        self._spline = spline
        self._solution_field = solution_field
        if solution_field is not None:
            self._mapper = self._solution_field.mapper(reference=self._spline)

        self._para_dim = spline.para_dim
        if self._para_dim == 3:
            raise NotImplementedError("Not yet tested for 3D")

        if solution_field is None:
            self._ukv = spline.unique_knots
        else:
            self._ukv = self._solution_field.unique_knots
        n_elems_per_dim = [len(kv) - 1 for kv in self._ukv]
        self._n_elems = _np.prod(n_elems_per_dim)

        # Gauss-Legendre quadrature points and weights
        spline_for_quad = spline if solution_field is None else solution_field

        if orders is None:
            quad_positions = []
            quad_weights = []
            for dim_quadrature_order in _default_quadrature_orders(
                spline_for_quad
            ):
                quad_position, quad_weight = _np.polynomial.legendre.leggauss(
                    deg=dim_quadrature_order
                )
                # Make quadrature points go from [0,1] instead of [-1,1]
                quad_positions.append((quad_position + 1) / 2)
                # Adjust weights accordingly
                quad_weights.append(quad_weight / 2)

            self._quad_positions = _cartesian_product(quad_positions)
            self._quad_weights = _np.prod(
                _cartesian_product(quad_weights), axis=1
            )
        else:
            self._quad_positions, self._quad_weights = (
                _get_quadrature_information(spline_for_quad, orders)
            )

        # Precompute grid IDs
        self._grid_ids = _cartesian_product(
            [_np.arange(n_elems) for n_elems in n_elems_per_dim],
            reverse=True,
        )

        self._all_supports = None
        self._all_element_quad_points = None
        self._all_jacobians = None
        self._all_jacobian_inverses = None
        self._all_jacobian_determinants = None

    def check_element_id_validity(self, element_id):
        """Check if given element ID is valid

        Parameters
        -----------
        element_id: int
            ID of element in spline's element. ID-array is 1D
        """
        assert element_id >= 0
        assert element_id < self._n_elems

    @property
    def all_supports(self):
        """Supports of all quadrature points.
        List of <n_elements> entries of support"""
        return self._all_supports

    @property
    def all_quad_points(self):
        """Quadrature points of all elements.
        Dimensions [<n_elements>, <n_quad_pts>, 2]"""
        return self._all_element_quad_points

    @property
    def all_jacobians(self):
        """Jacobians of all elements.
        Dimensions [<n_elements> <n_quad_pts>, <para_dim>, <para_dim>]"""
        return self._all_jacobians

    @property
    def all_jacobian_inverses(self):
        """Inverses of Jacobians of all elements.
        Dimensions [<n_elements> <n_quad_pts>, <para_dim>, <para_dim>]"""
        return self._all_jacobian_inverses

    @property
    def all_jacobian_determinants(self):
        """Determinants of Jacobians of all elements.
        Dimensions [<n_elements> <n_quad_pts>]"""
        return self._all_jacobian_determinants

    @property
    def quadrature_weights(self):
        return self._quad_weights

    def get_element_grid_id(self, element_id):
        """Compute element ID in grid

        Parameters
        ----------
        element_id: int
            ID of element of spline

        Returns
        ---------
        element_grid_id: list<int>
            ID of element in grid
        """
        if self._para_dim == 3:
            raise NotImplementedError(
                "Element grid ID not yet implemented for 3D"
            )

        return self._grid_ids[element_id, :]

    def get_element_quad_points(self, element_id):
        """For given element computes quad points

        Parameters
        -----------
        element_id: int
            ID of element in spline's element

        Returns
        -----------
        element_quad_points: np.ndarray
            Quadrature points for element
        """
        self.check_element_id_validity(element_id)

        if self._all_element_quad_points is not None:
            return self._all_element_quad_points[element_id]

        element_grid_id = self.get_element_grid_id(element_id)

        element_corner_points = _np.vstack(
            [
                ukv_dim[e_dim_id : (e_dim_id + 2)]
                for ukv_dim, e_dim_id in zip(self._ukv, element_grid_id)
            ]
        )
        element_lengths = _np.diff(element_corner_points, axis=1).ravel()
        element_midpoints = _np.mean(element_corner_points, axis=1)

        # Bring center to origin and scale
        element_quad_points = (self._quad_positions - 0.5) * element_lengths
        # Apply offset
        element_quad_points += element_midpoints

        return element_quad_points

    def get_element_support(self, element_id):
        """Get support for quadrature points in element

        Parameters
        ------------
        element_id: int
            ID of element

        Returns
        ---------
        support: np.ndarray
            Support for element. All quadrature points have same support
        """
        element_quad_points = self.get_element_quad_points(element_id)

        # All quad points in element have same support, therefore take arbitrary
        # one
        relevant_quad_point = element_quad_points[0, :]

        if self._solution_field is None:
            return self._spline.support(relevant_quad_point)
        else:
            return self._solution_field.support(relevant_quad_point)

    def jacobian(self, element_id):
        """Return Jacobian of single element at quadrature points

        Parameters
        ------------
        element_id: list<int>
            ID of element in grid

        Returns
        ---------
        element_jacobian: np.ndarray
            Jacobian of element evaluated at quadrature points
        """
        if self._all_jacobians is not None:
            return self._all_jacobians[element_id]

        element_quad_points = self.get_element_quad_points(element_id)

        return self._spline.jacobian(element_quad_points)

    def jacobian_inverse(self, element_id):
        """Return inverse of Jacobian of single element, evaluated at quadrature points

        Parameters
        ------------
        element_id: list<int>
            ID of element in grid

        Returns
        ---------
        element_inverse_jacobian: np.ndarray
            Inverse of Jacobian of element evaluated at quadrature points
        """
        if self._all_jacobian_inverses is not None:
            return self._all_jacobian_inverses[element_id]

        element_jacobians = self.jacobian(element_id)
        element_jacobian_inverse = _np.stack(
            [
                _np.linalg.inv(element_jacobian)
                for element_jacobian in element_jacobians
            ]
        )

        return element_jacobian_inverse

    def jacobian_determinant(self, element_id):
        """Return determinant of Jacobian of single element, evaluated at
        quadrature points

        Parameters
        ------------
        element_id: list<int>
            ID of element in grid

        Returns
        ---------
        element_jacobian_determinant: np.ndarray
            Determinant of Jacobian of element evaluated at quadrature points
        """
        if self._all_jacobian_determinants is not None:
            return self._all_jacobian_determinants[element_id]

        element_jacobians = self.jacobian(element_id)
        return _np.array(
            [
                _np.linalg.det(element_jacobian)
                for element_jacobian in element_jacobians
            ]
        )

    def compute_all_element_quad_points(self, recompute=False):
        """Compute the quadrature points of all elements

        Parameters
        ----------
        recompute: bool
            Recompute quadrature points
        """
        if self._all_element_quad_points is not None and not recompute:
            return

        # compute element lengths and center points
        element_lengths = _cartesian_product(
            [_np.diff(dim_ukv) for dim_ukv in self._ukv]
        )
        element_midpoints = (
            _cartesian_product([dim_ukv[:-1] for dim_ukv in self._ukv])
            + element_lengths / 2
        )
        # Scale quad points for each element
        quad_points_centered = _np.einsum(
            "ij,hj->hij", self._quad_positions, element_lengths
        )
        # apply offset to quad points
        n_elements, n_quad_points, _ = quad_points_centered.shape
        offsets = element_midpoints - element_lengths / 2
        self._all_element_quad_points = quad_points_centered + _np.repeat(
            (offsets).reshape(n_elements, 1, -1), n_quad_points, 1
        )

    def compute_all_supports(self, recompute=False):
        """Compute the support for all quadrature points

        Parameters
        --------------
        recompute: bool
            Recompute quadrature points
        """
        if self._all_supports is not None and not recompute:
            return

        self.compute_all_element_quad_points(recompute=recompute)
        relevant_spline = (
            self._spline
            if self._solution_field is None
            else self._solution_field
        )
        self._all_supports = [
            relevant_spline.support(quad_points[:1, :]).ravel()
            for quad_points in self._all_element_quad_points
        ]

    def compute_all_element_jacobians(self, recompute=False):
        """Compute Jacobians of each element at each quadrature point

        Parameters
        ----------
        recompute: bool
            Recompute Jacobians
        """
        if self._all_jacobians is not None and not recompute:
            return

        self.compute_all_element_quad_points(recompute=recompute)
        self._all_jacobians = _np.stack(
            [
                self._spline.jacobian(quad_points)
                for quad_points in self._all_element_quad_points
            ]
        )

    def compute_all_element_jacobian_inverses(self, recompute=False):
        """Compute Jacobians' inverses of each element at each quadrature point

        Parameters
        ----------
        recompute: bool
            Recompute Jacobians' inverses
        """
        if self._all_jacobian_inverses is not None and not recompute:
            return

        self.compute_all_element_jacobians(recompute=recompute)

        self._all_jacobian_inverses = _np.stack(
            [
                _np.stack(
                    [
                        _np.linalg.inv(element_jacobian)
                        for element_jacobian in element_jacobians
                    ]
                )
                for element_jacobians in self._all_jacobians
            ]
        )

    def compute_all_element_jacobian_determinants(self, recompute=False):
        """Compute Jacobians' determinants of each element at each quadrature point

        Parameters
        ----------
        recompute: bool
            Recompute Jacobians' determinants
        """
        if self._all_jacobian_determinants is not None and not recompute:
            return

        self.compute_all_element_jacobians(recompute=recompute)

        self._all_jacobian_determinants = _np.stack(
            [
                _np.stack(
                    [
                        _np.linalg.det(element_jacobian)
                        for element_jacobian in element_jacobians
                    ]
                )
                for element_jacobians in self._all_jacobians
            ]
        )


class FieldIntegrator(_SplinepyBase):
    __slots__ = (
        "_helpee",
        "_solution_field",
        "_mapper",
        "_global_rhs",
        "_global_system_matrix",
        "_trafo",
        "_supports",
        "_ndofs",
    )

    def __init__(self, geometry, solution_field=None, orders=None):
        """ """
        self._helpee = geometry
        if solution_field is None:
            self._solution_field = geometry.copy()
            self._solution_field.control_points = _np.ones(
                (geometry.cps.shape[0], 1)
            )
        else:
            self._solution_field = solution_field
        self._ndofs = int(
            _np.prod(
                [
                    len(kv) - 1 - deg
                    for deg, kv in zip(
                        self._solution_field.degrees,
                        self._solution_field.knot_vectors,
                    )
                ]
            )
        )
        self._mapper = self._solution_field.mapper(reference=self._helpee)

        self.reset(orders)

    def reset(self, orders=None):
        """ """
        self._trafo = Transformation(
            self._helpee, self._solution_field, orders
        )
        self.precompute_transformation()

        self._supports = None
        self._global_rhs = None
        self._global_system_matrix = None

    def precompute_transformation(self):
        """Computes the quadrature points, jacobians and their determinants
        of all elements in spline"""
        self._trafo.compute_all_supports()
        self._trafo.compute_all_element_jacobian_determinants()

    @property
    def supports(self):
        """ """
        if self._supports is None:
            self._supports = self._helpee.supports(self._trafo.all_quad_points)

        return self._supports

    def assemble_matrix(self, function, matrixout=None):
        """Assemble the system matrix for a given function. If system matrix is
        already assembled, it will add values on top of existing matrix.

        Parameters
        ------------
        function: callable
            Function which defines how to assemble the system matrix
        matrixout: np.ndarray
            Assembled matrix will be stored there. Default is global system matrix
        """
        # Initialize system matrix if not already
        if self._global_system_matrix is None and matrixout is None:
            global_size = (self._ndofs, self._ndofs)
            if _has_scipy:
                self._global_system_matrix = _dok_matrix(global_size)
            else:
                self._global_system_matrix = _np.zeros(global_size)

        # If other matrix is used, check if it compatible
        if matrixout is not None:
            if _has_scipy:
                assert isinstance(
                    matrixout, _dok_matrix
                ), "Matrixout must be scipy sparse dok matrix"
            else:
                assert isinstance(matrixout, _np.ndarray)
            assert matrixout.shape == (self._ndofs, self._ndofs)

        system_matrix = (
            self._global_system_matrix if matrixout is None else matrixout
        )

        quad_weights = self._trafo.quadrature_weights

        # Element loop
        for element_jacobian_det, element_support, element_quad_points in zip(
            self._trafo.all_jacobian_determinants,
            self._trafo.all_supports,
            self._trafo.all_quad_points,
        ):
            element_matrix = function(
                mapper=self._mapper,
                quad_points=element_quad_points,
                quad_weights=quad_weights,
                jacobian_det=element_jacobian_det,
            )
            matrix_element_support = _cartesian_product(
                [element_support, element_support]
            )
            system_matrix[
                matrix_element_support[:, 0], matrix_element_support[:, 1]
            ] += element_matrix

    def assemble_vector(self, function, current_sol=None, vectorout=None):
        """Assemble the rhs for a given function. If rhs is already assembled,
        it will add values on top of existing rhs.

        Parameters
        ------------
        function: callable
            Function which defines how to assemble the rhs vector
        current_sol: np.ndarray
            Current solution vector. Needed for nonlinear forms
        """
        # Initialize rhs vector
        if self._global_rhs is None:
            self._global_rhs = _np.zeros(self._ndofs)

        # Ensure that current solution has right dimensions
        if current_sol is not None:
            assert len(current_sol) == self._ndofs

        if vectorout is not None:
            assert isinstance(vectorout, _np.ndarray)
            assert len(vectorout) == self._ndofs

        rhs_vector = self._global_rhs if vectorout is None else vectorout

        # Prepare function arguments, which stay the same for all elements
        function_args = {
            "mapper": self._mapper,
            "quad_weights": self._trafo.quadrature_weights,
        }

        # Element loop
        for element_jacobian_det, element_support, element_quad_points in zip(
            self._trafo.all_jacobian_determinants,
            self._trafo.all_supports,
            self._trafo.all_quad_points,
        ):
            # Assemble element vector
            function_args["quad_points"] = element_quad_points
            function_args["jacobian_det"] = element_jacobian_det
            if current_sol is not None:
                function_args["current_sol"] = current_sol[element_support]
            element_vector = function(**function_args)

            # Add element vector to global rhs vector
            rhs_vector[element_support] += element_vector

    def check_if_assembled(self):
        """
        Check if system matrix and rhs are already assembled
        """
        if self._global_system_matrix is None or self._global_rhs is None:
            raise ValueError("System is not yet fully assembled")

    def get_boundary_dofs(self):
        """
        Get indices of boundary dofs.

        Returns
        -------------
        indices: np.ndarray
            Indices of relevant boundary dofs
        """
        relevant_spline = (
            self._helpee
            if self._solution_field is None
            else self._solution_field
        )

        multi_index = relevant_spline.multi_index

        indices = _np.unique(
            _np.hstack(
                (
                    multi_index[0, :],
                    multi_index[-1, :],
                    multi_index[:, 0],
                    multi_index[:, -1],
                )
            )
        )

        return indices

    def assign_homogeneous_dirichlet_boundary_conditions(self):
        """
        Assembles homogeneous Dirichlet boundary conditions
        """
        self.check_if_assembled()

        indices = self.get_boundary_dofs()

        self._global_system_matrix[indices, :] = 0
        self._global_system_matrix[indices, indices] = 1
        self._global_rhs[indices] = 0

    def apply_dirichlet_boundary_conditions(self, function):
        """
        Applies Dirichlet boundary conditions via :math:`L^2`-projection.

        For a given function g, it solves the following equation

        .. math:: \\sum\\limits_{j=1}^n \alpha_j (N_i, N_j) = (g, N_i) \\quad
        \forall i = 1, \\dots, n.

        Then, the :math:`Pg`, the :math:`L^2`-projection of the function, is
        given by

        .. math:: Pg(x) = \\sum\\limits_{j=1}^n \alpha_j \\phi_j(x)

        For the Dirichlet boundary conditions, only the DoFs corresponding to
        the boundaries are taken from the :math:`L^2`-projection.

        Parameters
        -------------
        function: callable
            Function to apply. Input are points, output is scalar
        """
        self.check_if_assembled()

        # Assemble mass matrix
        global_size = (self._ndofs, self._ndofs)
        mass_matrix = (
            _dok_matrix(global_size) if _has_scipy else _np.zeros(global_size)
        )

        def mass_lhs(mapper, quad_points, quad_weights, jacobian_det):
            bf_values = mapper._field_reference.basis(quad_points)
            element_matrix = _np.einsum(
                "qi,qj,q,q->ij",
                bf_values,
                bf_values,
                quad_weights,
                jacobian_det,
                optimize=True,
            )
            return element_matrix.ravel()

        self.assemble_matrix(mass_lhs, matrixout=mass_matrix)

        # Assemble rhs: f(x) * N
        rhs_vector = _np.zeros(self._ndofs)

        def rhs_function(mapper, quad_points, quad_weights, jacobian_det):
            bf = mapper._field_reference.basis(quad_points)
            quad_points_forward = mapper._geometry_reference.evaluate(
                quad_points
            )
            function_values = function(quad_points_forward)
            element_vector = _np.einsum(
                "qj,q,q,q->j",
                bf,
                function_values,
                quad_weights,
                jacobian_det,
                optimize=True,
            )
            return element_vector

        self.assemble_vector(rhs_function, vectorout=rhs_vector)

        # Solve system to get all dofs
        dof_vector = _np.empty(self._ndofs)
        if _has_scipy:
            dof_vector = _spsolve(mass_matrix.tocsr(), rhs_vector)
        else:
            dof_vector = _np.linalg.solve(mass_matrix, rhs_vector)

        # Get relevant dofs
        indices = self.get_boundary_dofs()

        # Apply Dirichlet to relevant boundary dofs
        self._global_system_matrix[indices, :] = 0
        self._global_system_matrix[indices, indices] = 1
        self._global_rhs[indices] = dof_vector[indices]

    def solve_linear_system(self):
        """
        Solve linear system for system matrix and rhs
        """
        self.check_if_assembled()

        if _has_scipy:
            solution_vector = _spsolve(
                self._global_system_matrix.tocsr(), self._global_rhs
            )
        else:
            solution_vector = _np.linalg.solve(
                self._global_system_matrix, self._global_rhs
            )
        self._solution_field.control_points = solution_vector.reshape(-1, 1)
