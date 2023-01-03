# splinepy
[![workflow](https://github.com/tataratat/splinepy/actions/workflows/main.yml/badge.svg)](https://github.com/tataratat/splinepy/actions)
[![PyPI version](https://badge.fury.io/py/splinepy.svg)](https://badge.fury.io/py/splinepy)
  
splinepy is a python library for splines of arbitrary dimensions and degrees.
The library supports Bezier, Rational Bezier, BSpline and NURBS with fast and easy-to-use APIs.

## Install guide
For python3.6+ splinepy is available through `pip`:
```
pip install --upgrade pip
pip install splinepy
```

It is also possible to install current development version using `pip`. It requires a compiler that supports C++17 or higher (C++20 for debug mode - tested with gcc-10.3 and clang-12). Here are two variants:
```
# Fast build - minimal and debug mode
SPLINEPY_MINIMAL_DEBUG_BUILD=True pip install git+https://github.com/tataratat/splinepy.git@main -vvv

# Same build as in PyPI - full set of splines and optimized build
pip install git+https://github.com/tataratat/splinepy.git@main -vvv
```
The last part `-vvv` is not necessary, but we suggest using it, since you can see the build progress. Full build (the second option) may take a while.

Of course, you can install directly from the source.
In addition to aforementioned compilers, this requires a cmake3.16+. If you don't have cmake, easiest way to install it would be: `pip install cmake`.
```
git clone git@github.com:tataratat/splinepy.git
cd splinepy
git submodule update --init --recursive
python3 setup.py install
```

## Feature Summary
For details, please take a look at the [documentation](https://tataratat.github.io/splinepy).
Most of the functions are vectorized and capable of multithread executions.

### Common features:
| Method | Description |
| ------ | ----------- |
| evaluate() | Given parametric coordinates, returns physical (i.e., mapped / evaluated) coordinate |
| derivative() | Given parametric cooridnates and order of partial derivatives, returns physical derivatives |
| sample() | Given number of sampling points per parametric dimension,  returns evaluated physical coordinates, which are equally distributed in parametric space. |
| basis_and_support() | Given parametric coordinates, returns basis function values and their support. Satisfies partition of unity. |
| basis_derivative_and_support() | Given parametric coordinates and order of partial derivatives, returns basis function derivative values. |
| proximities() | Given physical coordinates, returns parametric coordinates that maps to the nearest physical coordinate. Often referred as "point inversion". |
| elevate_degrees() | Elevates Spline degrees along specified parametric dimensions |

### BSpline, NURBS
| Method | Description |
| ------ | ----------- |
| reduce_degrees() | Reduces spline degrees along specified parametric dimensionans, as long as the spline stays under given tolerance |
| insert_knots() | Insert knots at given locations |
| remove_knots() | Removes knots at given locations, as long as the spline stays under given tolerance |

### Bezier, Rational Bezier
| Method | Description |
| ------ | ----------- |
| multiply (*) | Given two Beziers, returns a Bezier, which evalutation equals the product of two inputs |
| add (+) | Given two Beziers, returns a Bezier, which evalutation equals the product of two inputs |



### Multipatch:

### IO:


## Quick start
```
Coming Soon!
```
Test version of documentations are available [here](https://tataratat.github.io/splinepy)


### Dependencies
