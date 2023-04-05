# splinepy
[![workflow](https://github.com/tataratat/splinepy/actions/workflows/main.yml/badge.svg)](https://github.com/tataratat/splinepy/actions)
[![PyPI version](https://badge.fury.io/py/splinepy.svg)](https://badge.fury.io/py/splinepy)

splinepy is a python library for splines of arbitrary dimensions and degrees.
The library supports Bezier, Rational Bezier, BSpline and NURBS with fast and easy-to-use APIs.

## Install guide
splinepy wheels are available for python3.6+ for MacOS, Linux, and Windows:
```bash
pip install --upgrade pip
pip install splinepy
```

It is also possible to install current development version using `pip`. It requires a compiler that supports C++17 or higher (C++20 for debug mode - tested with gcc-10.3 and clang-12). Here are two variants:
1) Fast build - minimal and debug mode
```bash
SPLINEPY_MINIMAL_DEBUG_BUILD=True pip install git+https://github.com/tataratat/splinepy.git@main -vvv
```

2) Same build as in PyPI - full set of splines and optimized build
```bash
pip install git+https://github.com/tataratat/splinepy.git@main -vvv
```
`-vvv` is not necessary, but we suggest using it, since you can see the build progress. Full build (the second option) may take a while.

Of course, you can install directly from the source.
In addition to aforementioned compilers, this requires a cmake3.16+. If you don't have cmake, easiest way to install it would be: `pip install cmake`.
```bash
git clone git@github.com:tataratat/splinepy.git
cd splinepy
git submodule update --init --recursive
python3 setup.py install
```
For visualization and extended functionalities, please take a look at [gustaf](https://github.com/tataratat/gustaf)!

## Quick start
```python
import splinepy
import numpy as np

# Initialize bspline with any array-like input
bspline = splinepy.BSpline(
    degrees=[2, 1],
    knot_vectors=[
        [0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        [0.0, 0.0, 1.0, 1.0],
    ],
    control_points=[
        [0. , 0. ],
        [0.5, 0. ],
        [1. , 0. ],
        [0. , 1. ],
        [0.5, 1. ],
        [1. , 1. ],
    ],
)

# We always store control points in 2D arrays with
# (total_number_of_control_points, physical_dimension) shape.
# They fill control grid by iterating lower-indexed dimensions first.
# But if you prefer the grid-like structure, this should hold
grid_cps = np.empty(2, 3, 2)  # (dim, n_cps_u, n_cps_v)
gird_cps[:, 0, 0] = [0. , 0. ]
gird_cps[:, 0, 1] = [0.5, 0. ]
gird_cps[:, 0, 2] = [1. , 0. ]
gird_cps[:, 1, 0] = [0. , 1. ]
gird_cps[:, 1, 1] = [0.5, 1. ]
gird_cps[:, 1, 2] = [1. , 1. ]

assert np.allclose(
    bspline.control_points,
    grid_cps.reshape(-1, 2)  # (-1, dim)
)

# Evaluate spline mapping.
# First, let's form parametric coordinate queries
queries = [
    [0.1, 0.2], # first query
    [0.4, 0.5], # second query
    [0.1156, 0.9091],  # third query
]
physical_coords = bspline.evaluate(queries)

# we can also execute this in parallel using multithread executions in c++ side
# (probably makes more sense if you have giant queries)
physical_coords_parallel = bspline.evaluate(queries, nthreads=2)

# this holds
assert np.allclose(physical_coords, physical_coords_parallel)

```

## Feature Summary
For details, please take a look at the [documentation](https://tataratat.github.io/splinepy).
Most of the functions are vectorized and capable of multithread executions.

#### Common features
| Method                         | Description                                                                                                                                           |
| ------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| evaluate()                     | Given parametric coordinates, returns physical (i.e., mapped / evaluated) coordinate                                                                  |
| derivative()                   | Given parametric cooridnates and order of partial derivatives, returns physical derivatives                                                           |
| jacobian()                     | Given a number of query points, returns the jacobian of the geometry or field                                                                         |
| sample()                       | Given number of sampling points per parametric dimension,  returns evaluated physical coordinates, which are equally distributed in parametric space. |
| basis_and_support()            | Given parametric coordinates, returns basis function values and their support. Satisfies partition of unity.                                          |
| basis_derivative_and_support() | Given parametric coordinates and order of partial derivatives, returns basis function derivative values.                                              |
| proximities()                  | Given physical coordinates, returns parametric coordinates that maps to the nearest physical coordinate. Often referred as "point inversion".         |
| elevate_degrees()              | Elevates Spline degrees along specified parametric dimensions                                                                                         |
| extract_boundaries()           | Given boundary ids, returns extracted boundary splines.                                                                                               |

#### BSpline, NURBS
| Method                   | Description                                                                                                       |
| ------------------------ | ----------------------------------------------------------------------------------------------------------------- |
| reduce_degrees()         | Reduces spline degrees along specified parametric dimensionans, as long as the spline stays under given tolerance |
| insert_knots()           | Insert knots at given locations                                                                                   |
| remove_knots()           | Removes knots at given locations, as long as the spline stays under given tolerance                               |
| extract_bezier_patches() | Extracts each knot spans as a Bezier spline                                                                       |

#### Bezier, Rational Bezier
| Method                   | Description                                                                                                                                                                                                                        |
| ------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| multiply (*)             | Given two Beziers `a` and `b`, returns a Bezier `c` that satisfies: c.evaluate(query) = a.evalute(query) * b.evaluate(query)                                                                                                       |
| add (+)                  | Given two Beziers `a` and `b`, returns a Bezier `c` that satisfies:  c.evaluate(query) = a.evaluate(query) + b.evaluate(query)                                                                                                     |
| derivative_spline()      | Given order or partial derivatives, returns a bezier `c` that satisfies: c.derivative(query, orders) = a.derivative_spline(orders).evaluate(query)                                                                                 |
| split()                  | Splits Bezier into multiple patches at defined locations                                                                                                                                                                           |
| extract_dim()            | Extract a single physical dimension                                                                                                                                                                                                |
| compose()                | Given an inner function spline `a` and an outer function spline `b`, it returns their functional composition that satisfies c.evaluate(query) = b.evaluate(a.evaluate(query))                                                      |
| composition_derivative() | Given an outer function `a`, an inner function `b` and its derivative `b'` with respect to some variable `s`, it returns the derivative with respect to the same variable `s` of the composition `a(b)` by applying the chain-rule |


### Multipatch
Splinepy offers a common interface for multipatch geometries, i.e., geometries consisting of multiple, individual splines of arbitrary types. This concept is both used for complex geometries as for Isogeometric Analysis. Multipatch objects have the following functionalities:
 - determine patch-interfaces automatically
 - identification of boundary faces
 - boundary assignement using different techniques, relying either on the boundary position or on the continuity inbetween patches
 - Boundary extraction

### IO
Available in `splinepy.io`.

| Formats | Description                                                                                                                                          |
| ------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- |
| iges    | Loads/Exports splines from an [IGES](https://en.wikipedia.org/wiki/IGES) file                                                                        |
| irit    | [IRIT](https://www.cs.technion.ac.il/~irit/) compatible format                                                                                       |
| json    | (Custom) easy-to-read format, supports base64 encoding                                                                                               |
| mfem    | [MFEM](https://mfem.org) compatible `.mesh` format. Supports structured multi-patch splines in `controlpoints_cartesian` and 2D single-patch splines |
| gismo   | [GISMO](https://gismo.github.io) compatible `.xml` format                                                                                            |
| npz     | Based on np.savez()                                                                                                                                  |
| xml     | [RWTH CATS](https://www.cats.rwth-aachen.de/) spline format                                                                                          |


### BSpline fitting
Implements fitting routines from the [The NURBS Book](https://link.springer.com/book/10.1007/978-3-642-97385-7). Available as BSpline's classmethod.
- interpolate_curve()
- approximate_curve()
- interpolate_surface()
- approximate_surface()


## Dependencies
Followings are direct dependencies for splinepy. If something is missing, please feel free to extend!
#### c++
- [pybind11](https://github.com/pybind/pybind11)
- [SplineLib](https://github.com/tataratat/SplineLib)
- [bezman](https://github.com/tataratat/bezman)
- [napf](https://github.com/tataratat/napf) - wraps / based on [nanoflann](https://github.com/jlblancoc/nanoflann)
#### python
- [numpy](https://numpy.org)
#### build
- [cmake](https://cmake.org)
- [setuptools](https://setuptools.pypa.io/en/latest/)
- [wheel](https://wheel.readthedocs.io/en/stable/)
