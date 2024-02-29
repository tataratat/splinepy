from _typeshed import Incomplete
from splinepy._base import SplinepyBase as _SplinepyBase

class Mapper(_SplinepyBase):
    def __init__(self, field, reference) -> None: ...
    def basis_function_derivatives(self, queries, gradient: bool = False, hessian: bool = False, laplacian: bool = False, nthreads: Incomplete | None = None): ...
    def field_derivatives(self, queries, gradient: bool = False, divergence: bool = False, hessian: bool = False, laplacian: bool = False, basis_function_values: bool = False, nthreads: Incomplete | None = None): ...
    def basis_gradient_and_support(self, queries, nthreads: Incomplete | None = None): ...
    def gradient(self, queries, nthreads: Incomplete | None = None): ...
    def divergence(self, queries, nthreads: Incomplete | None = None): ...
    def basis_hessian_and_support(self, queries, nthreads: Incomplete | None = None): ...
    def hessian(self, queries, nthreads: Incomplete | None = None): ...
    def basis_laplacian_and_support(self, queries, nthreads: Incomplete | None = None): ...
    def laplacian(self, queries, nthreads: Incomplete | None = None): ...
