"""
Helps you map derivatives of fields and basis function into the physical
(image) domain
"""
import numpy as np

from splinepy import utils
from splinepy._base import SplinepyBase


class Mapper(SplinepyBase):
    r"""Map expressions and derivatives into the physical domain using a
    geometric mapper

    The mapper has access to two fields: the solution field and the underlying
    geometry representation, both are assumed to be spline-functions. Hence, we
    assume a geometry representation in the form:

    .. math:: x_i (\mathbf{u}) = \sum_b N^b (\mathbf{u}) x^b_i

    as well as a field representation in the form :

    .. math:: f_i (\mathbf{u}) = \sum_a N^a (\mathbf{u}) f^a_i

    With :math:`x_i`, :math:`f_i` being (vector valued) geometry and solution
    field functions, with their respective shape/basis functions :math:`N^b`,
    :math:`N^a` and their associated coefficients/control points, :math:`x^b_i`
    and :math:`f^a_i`. Further, using the definition of the Jacobian and it's
    inverse we can write

    .. math:: \mathbf{J}^{-1} = \bar{\mathbf{J}}
    .. math:: J_{ij} = \frac{\partial x_i}{\partial u_j}
    .. math:: \bar{J}_{ij} = \frac{\partial u_i}{\partial x_j}

    Basis function derivatives of the first order can be mapped into the
    physical domain, using the following equation:

    .. math:: \frac{\partial N^a}{\partial x_i} =
      \frac{ \partial N^a}{\partial u_j} \frac{ \partial u_j}{\partial x_i} =
      \frac{\partial N^a}{ \partial u_j} \bar{J}_{ji}

    Using the definition of the solution field, field gradients (and
    divergences) can be computed in a similar fashion

    .. math:: \frac{\partial f_i}{\partial x_j} =
        \sum_a \frac{\partial N^a}{\partial x_j} f^a_i

    Divergences are computed by summing up the respective components.

    Second order derivatives of the basis function values with respect to the
    reference coordinates are written as (similar for :math:`H^b`):

    .. math:: H^a_{ij} = \frac{\partial^2 N^a}{\partial u_i\partial u_j}

    The hessians can be mapped into physical space using the following
    equation:

    .. math:: \frac{\partial^2 N^a}{\partial x_i\partial x_j} =
      \bar{J}_{li} \left(
        H^a_{kl} - \frac{\partial N^a}{\partial u_n}\bar{J}_{nm} x^b_m H^b_{lk}
      \right) \bar{J}_{kj}

    Here, the expression :math:`x^b_mH^b_{lk}` represents the Hessian of the
    geometric reference spline :math:`\tilde{H}_{mlk}`, which can sometimes be
    calculated more efficiently using spline methods implemented in SplineLib
    and Bezman.

    Here we used the following identity:

    .. math:: \frac{\partial \bar{J}_{li} }{\partial u_k} =
      - \bar{J}_{lm} \frac{\partial J_{mn}}{\partial u_k} \bar{J}_{ni}\\


    Parameters
    ----------
    field : spline
      Field that is to be mapped
    reference : spline
      Reference for mapping (geometry)
    """

    def __init__(self, field, reference):
        self._field_reference = field
        self._geometry_reference = reference

        # Check the parametric dimensions
        if not reference.para_dim == field.para_dim:
            raise ValueError("Parametric dimension mismatch")
        if not reference.para_dim == reference.dim:
            raise ValueError(
                "Mismatch between physical and parametric dimension for "
                "geometry representation"
            )
        # Easy access
        self._para_dim = reference.para_dim

    def basis_function_derivatives(
        self,
        queries,
        gradient=False,
        hessian=False,
        laplacian=False,
        nthreads=None,
    ):
        """Function to retrieve more than one basis function derivative

        More efficient implementation if more than one derivative is required,
        as many values can be precalculated. See class documentation for more
        details.

        Parameters
        ----------
        queries: (n, para_dim) array-like
          Evaluation points in the parametric domain
        gradient : bool
          Evaluate Basis Function Gradient mapped into the physical domain
        hessian : bool
          Evaluate Basis Function Hessian mapped into the physical domain
        laplacian : bool
          Evaluate Basis Function Laplacian mapped into the physical domain
        nthreads: int
          Number of threads available for the computation

        Returns
        -------
        results : dict
          Dictionary with required values stored with same name as function
          arguments
        """
        # Naming convention for indices (more indices follow from
        # implementation details in class documentation):
        # q : query-point
        # s : support-index
        # a : solution field shape_function
        # b : geometry field shape_function
        self._logd("Evaluating basis function gradients in physical space")
        queries = utils.data.enforce_contiguous(queries, dtype="float64")

        # Precompute required values
        invjacs = np.linalg.inv(
            self._geometry_reference.jacobian(queries, nthreads)
        )
        results = {}

        bf_gradients = np.empty(
            (
                queries.shape[0],
                np.prod(self._field_reference.degrees + 1),
                self._para_dim,
            )
        )
        (
            bf_gradients[:, :, 0],
            support,
        ) = self._field_reference.basis_derivative_and_support(
            queries=queries,
            # Array size M with [0, ..., 0, 1, 0, ..., 0] with 1 at position k
            orders=np.eye(1, M=self._para_dim, k=0),
            nthreads=nthreads,
        )
        results["support"] = support
        for i in range(1, self._para_dim):
            bf_gradients[:, :, i] = self._field_reference.basis_derivative(
                queries=queries,
                orders=np.eye(1, M=self._para_dim, k=i),
                nthreads=nthreads,
            )
        if gradient:
            results["gradient"] = np.einsum(
                "qsi,qij->qsj", bf_gradients, invjacs, optimize=True
            )
            results["support"] = support

        if hessian or laplacian:
            # Retrieve basis function hessians from both the geometry as from
            # the field
            bf_hessians = np.empty(
                (
                    queries.shape[0],
                    np.prod(self._field_reference.degrees + 1),
                    self._para_dim,
                    self._para_dim,
                )
            )
            for i in range(self._para_dim):
                for j in range(i, self._para_dim):
                    bf_hessians[
                        :, :, i, j
                    ] = self._field_reference.basis_derivative(
                        queries=queries,
                        orders=np.eye(1, M=self._para_dim, k=i)
                        + np.eye(1, M=self._para_dim, k=j),
                        nthreads=nthreads,
                    )
                    if i != j:
                        bf_hessians[:, :, j, i] = bf_hessians[:, :, i, j]
            # This is unnecessary if isoparametric (but if only field is high
            # order, this is more efficient)
            geo_hessians = np.empty(
                (
                    queries.shape[0],
                    self._para_dim,  # dim==para_dim
                    self._para_dim,
                    self._para_dim,
                )
            )
            for i in range(self._para_dim):
                for j in range(i, self._geometry_reference.para_dim):
                    geo_hessians[
                        :, :, i, j
                    ] = self._geometry_reference.derivative(
                        queries=queries,
                        orders=np.eye(1, M=self._para_dim, k=i)
                        + np.eye(1, M=self._para_dim, k=j),
                        nthreads=nthreads,
                    )
                    if i != j:
                        geo_hessians[:, :, j, i] = geo_hessians[:, :, i, j]

            # Overwrite bf_hessians (see documentation for indices)
            bf_hessians -= np.einsum(
                "qan,qnm,qmlk->qalk",
                bf_gradients,
                invjacs,
                geo_hessians,
                optimize=True,
            )
        if hessian:
            results["hessian"] = np.einsum(
                "eli,ealk,ekj->eaij",
                invjacs,
                bf_hessians,
                invjacs,
                optimize=True,
            )
            results["support"] = support
        if laplacian:
            if hessian:
                results["laplacian"] = np.einsum(
                    "eaii->ea", results["hessian"], optimize=True
                )
            else:
                results["laplacian"] = np.einsum(
                    "eli,ealk,eki->ea",
                    invjacs,
                    bf_hessians,
                    invjacs,
                    optimize=True,
                )
                results["support"] = support

        return results

    def field_derivatives(
        self,
        queries,
        gradient=False,
        divergence=False,
        hessian=False,
        laplacian=False,
        basis_function_values=False,
        nthreads=None,
    ):
        """Function to retrieve more than one field derivative

        More efficient implementation if more than one derivative is required.
        Can also return basis function values if both are required, e.g., for
        some assembly. See class documentation for more details.

        Parameters
        ----------
        queries: (n, para_dim) array-like
          Evaluation points in the parametric domain
        gradient : bool
          Evaluate Gradient mapped into the physical domain
        divergence : bool
          Evaluate Divergenec of a vector field
        hessian : bool
          Evaluate Hessian mapped into the physical domain
        laplacian : bool
          Evaluate Laplacian mapped into the physical domain
        basis_function_values : bool
          Return basis function derivatives in dictionary
        nthreads: int
          Number of threads available for the computation

        Returns
        -------
        results : dict
          Dictionary with required values stored with same name as function
          arguments (basis function derivatives as dictionary in dictionary)
        """
        # Remark : In terms of computations, it might actually be more
        # efficient, to calculate the derivatives of the spline values
        # directly, rather than to precompute al basis function derivatives, as
        # more efficient methods might be implemented in both SplineLib and
        # Bezman, however, for sake of simplicity, we will stick to the basis
        # function based implemntation here, as basis function derivatives are
        # mostly required in IGA applications.

        # Compute required basis function values
        as_dictionary = self.basis_function_derivatives(
            queries=queries,
            gradient=(gradient or divergence),
            hessian=hessian,
            laplacian=laplacian,
            nthreads=nthreads,
        )

        results = {}
        supports = as_dictionary["support"]
        if basis_function_values:
            results["basis_function_values"] = as_dictionary
        # Start computation
        if gradient:
            results["gradient"] = np.einsum(
                "qsd,qsv->qvd",
                as_dictionary["gradient"],
                self._field_reference.control_points[supports, :],
                optimize=True,
            )
        if divergence:
            # Can only be performed on vector fields with para_dim=dim
            if self._field_reference.para_dim == self._field_reference.dim:
                if gradient:
                    results["divergence"] = np.einsum(
                        "qii->q", results["gradient"], optimize=True
                    )
                else:
                    results["divergence"] = np.einsum(
                        "qsd,qsd->q",
                        as_dictionary["gradient"],
                        self._field_reference.control_points[supports, :],
                        optimize=True,
                    )
            else:
                raise ValueError(
                    "Divergence can only be performed on vector fields with "
                    "para_dim = dim"
                )

        if hessian:
            results["hessian"] = np.einsum(
                "qsij,qsv->qvij",
                as_dictionary["hessian"],
                self._field_reference.control_points[supports, :],
                optimize=True,
            )
        if laplacian:
            if hessian:
                results["laplacian"] = np.einsum(
                    "qdii->qd", results["hessian"], optimize=True
                )
            else:
                results["laplacian"] = np.einsum(
                    "qs,qsd->qd",
                    as_dictionary["laplacian"],
                    self._field_reference.control_points[supports, :],
                    optimize=True,
                )

        return results

    def basis_gradient_and_support(self, queries, nthreads=None):
        """Map gradient of basis functions into the physical domain

        Calls self.basis_function_derivatives internally, see class
        documentation for more details.

        Parameters
        ----------
        queries: (n, para_dim) array-like
        nthreads: int

        Returns
        --------
        gradient: (n, prod(degrees + 1), para_dim) np.ndarray
        support: (n, prod(degrees + 1)) np.ndarray
        """
        as_dict = self.basis_function_derivatives(
            queries,
            gradient=True,
            hessian=False,
            laplacian=False,
            nthreads=nthreads,
        )
        return (as_dict["gradient"], as_dict["support"])

    def gradient(self, queries, nthreads=None):
        r"""Map gradient field into the physical domain

        Gradient is in form

        .. math:: J^k_i = \frac{\partial f^k}{\partial x_i}

        See class documentation for more details.

        Parameters
        ----------
        queries: (n, para_dim) array-like
          Points for evaluation
        nthreads: int
          Threads used for calculation

        Returns
        --------
        gradient: (n, dim, para_dim) np.ndarray
        """
        return self.field_derivatives(
            queries,
            gradient=True,
            divergence=False,
            hessian=False,
            laplacian=False,
            basis_function_values=False,
            nthreads=nthreads,
        )["gradient"]

    def divergence(self, queries, nthreads=None):
        """Map field divergence (where applicable) into the physical domain

        Returns scalar values for each query. See class documentation for more
        details.

        Parameters
        ----------
        queries: (n, para_dim) array-like
          Points for evaluation
        nthreads: int
          Threads used for calculation

        Returns
        --------
        gradient: (n) np.ndarray
        """
        if self._field_reference.para_dim != self._field_reference.dim:
            raise ValueError(
                "Divergence can only be performed on vector fields with "
                "para_dim = dim"
            )
        return self.field_derivatives(
            queries,
            gradient=False,
            divergence=True,
            hessian=False,
            laplacian=False,
            basis_function_values=False,
            nthreads=nthreads,
        )["divergence"]

    def basis_hessian_and_support(self, queries, nthreads=None):
        """Map hessian of basis functions into the physical domain

        Calls self.basis_function_derivatives internally, see class
        documentation for more details.

        Parameters
        ----------
        queries: (n, para_dim) array-like
        nthreads: int

        Returns
        --------
        hessians: (n, prod(degrees + 1), para_dim, para_dim) np.ndarray
        support: (n, prod(degrees + 1)) np.ndarray
        """
        as_dict = self.basis_function_derivatives(
            queries,
            gradient=False,
            hessian=True,
            laplacian=False,
            nthreads=nthreads,
        )
        return (as_dict["hessian"], as_dict["support"])

    def hessian(self, queries, nthreads=None):
        r"""Map hessian field into the physical domain

        Hessian is in form

        .. math:: H^k_{ij} = \frac{\partial^2 f_k}{\partial x_i \partial x_j}

        Parameters
        ----------
        queries: (n, para_dim) array-like
          Points for evaluation
        nthreads: int
          Threads used for calculation

        Returns
        --------
        gradient: (n, prod(degrees + 1), paradim) np.ndarray
        """
        return self.field_derivatives(
            queries,
            gradient=False,
            divergence=False,
            hessian=True,
            laplacian=False,
            basis_function_values=False,
            nthreads=nthreads,
        )["hessian"]

    def basis_laplacian_and_support(self, queries, nthreads=None):
        """Map laplacian of basis functions into the physical domain

        Calls self.basis_function_derivatives internally, see class
        documentation for more details.

        Parameters
        ----------
        queries: (n, para_dim) array-like
        nthreads: int

        Returns
        --------
        laplacian: (n, prod(degrees + 1)) np.ndarray
        support: (n, prod(degrees + 1)) np.ndarray
        """
        as_dict = self.basis_function_derivatives(
            queries,
            gradient=False,
            hessian=False,
            laplacian=True,
            nthreads=nthreads,
        )
        return (as_dict["laplacian"], as_dict["support"])

    def laplacian(self, queries, nthreads=None):
        """Map laplacian field into the physical domain

        Parameters
        ----------
        queries: (n, para_dim) array-like
          Points for evaluation
        nthreads: int
          Threads used for calculation

        Returns
        --------
        laplacian: (n) np.ndarray
        """
        return self.field_derivatives(
            queries,
            gradient=False,
            divergence=False,
            hessian=False,
            laplacian=True,
            basis_function_values=False,
            nthreads=nthreads,
        )["laplacian"]
