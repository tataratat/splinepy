# splinepy
[![workflow](https://github.com/tataratat/splinepy/actions/workflows/main.yml/badge.svg)](https://github.com/tataratat/splinepy/actions)
[![PyPI version](https://badge.fury.io/py/splinepy.svg)](https://badge.fury.io/py/splinepy)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/tataratat/try-splinepy/main)

splinepy is a python library for splines of arbitrary dimensions and degrees.
The library supports Bezier, Rational Bezier, BSpline and NURBS with fast and easy-to-use APIs.

## Install guide
splinepy wheels are available for python3.8+ for MacOS, Linux, and Windows:
```bash
# including all optional dependencies
pip install "splinepy[all]"  # quotation marks required for some shells
# or
pip install splinepy
```


Of course, you can install it directly from the source.
In addition to the aforementioned compilers, this requires a cmake3.16+. If you don't have `cmake`, the easiest way to install it would be: `pip install cmake`.
```bash
git clone git@github.com:tataratat/splinepy.git
cd splinepy
git submodule update --init --recursive
pip install -e .
```

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
        [0.0, 0.0],  # [0, 0] (control grid index)
        [0.5, 0.0],  # [1, 0]
        [1.0, 0.0],  # [2, 0]
        [0.0, 1.0],  # [0, 1]
        [0.5, 1.0],  # [1, 1]
        [1.0, 1.0],  # [2, 1]
    ],
)

# We always store control points in 2D arrays with shape
# (total_number_of_control_points, physical_dimension).
# The indexing of the control grid is defined by iterating
# lower-indexed dimensions first. But if you prefer a
# grid-like structure, try
multi_index = bspline.multi_index
grid_cps = np.empty(bspline.control_points.shape)
grid_cps[multi_index[0, 0]] = [0.0, 0.0]
grid_cps[multi_index[1, 0]] = [0.5, 0.0]
grid_cps[multi_index[2, 0], 0] = 1.0
# which also supports ranges
grid_cps[multi_index[:, 0], 1] = 0.0
grid_cps[multi_index[:, 1], 1] = 1.0
grid_cps[multi_index[:, 1], 0] = [0.0, 0.5, 1.0]

assert np.allclose(bspline.control_points, grid_cps)

# Evaluate spline mapping.
# First, let's form parametric coordinate queries
queries = [
    [0.1, 0.2],  # first query
    [0.4, 0.5],  # second query
    [0.1156, 0.9091],  # third query
]
physical_coords = bspline.evaluate(queries)

# we can also execute this in parallel using multithread
# executions on c++ side (for heavy multi-queries scenarios)
physical_coords_parallel = bspline.evaluate(queries, nthreads=2)

# this holds
assert np.allclose(physical_coords, physical_coords_parallel)
```
You can also try `splinepy` online by clicking the [Binder](https://mybinder.org/v2/gh/tataratat/try-splinepy/main) badge above!

## Feature Summary
For details, please take a look at the [documentation](https://tataratat.github.io/splinepy).
Most of the functions are vectorized and capable of multithread executions.

### Splines
__Any type of spline is capable of:__
- computing spline mappings, derivatives, partial derivatives, jacobian, basis functions, basis function derivatives, basis function partial derivatives, and proximity (point inversion, nearest mapping search),
- degree elevation,
- extracting boundary splines, and
- visualization (see [visualizing with splinepy](docs/markdown/spline_plotting.md)).

In addition to the common features, __Bezier and Rational Bezier__ can:
- add/multiply two splines,
- split itself into multiple patches,
- create derivative splines, and
- compose an inner spline into an outer spline and compute its composition derivative

and __BSpline and NURBS__ can:
- reduce degrees,
- insert and remove knots, and
- extract bezier patches.

Some __BSpline fitting__ routines from [The NURBS Book](https://link.springer.com/book/10.1007/978-3-642-97385-7):
- curve interpolation/approximation
- surface interpolation/approximation

### Multipatch
Splinepy offers a common interface for multipatch geometries, i.e., geometries consisting of multiple, individual splines of arbitrary types. This concept is used for complex geometries and for Isogeometric Analysis. __Multipatch__ objects have the following functionalities:
 - determine patch-interfaces automatically
 - identification of boundary faces
 - boundary assignment using different techniques, relying either on the boundary position or on the continuity in between patches
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


## Dependencies
The following are direct dependencies for splinepy. Please feel free to check out the repositories linked.

| Package | Description                                             | Python | C++ |
| ------- | ------------------------------------------------------- | :----: | :---: |
| [pybind11](https://github.com/pybind/pybind11) | Binds c++ and python | X | X |
| [SplineLib](https://github.com/tataratat/SplineLib) | Main functionalities for BSplines and NURBS |    | X |
| [bezman](https://github.com/tataratat/bezman)       | Main functionalities for Beziers and rational Beziers |    | X |
| [napf](https://github.com/tataratat/napf)           | Creates k-d trees that provide initial guess for proximity search. Wraps [nanoflann](https://github.com/jlblancoc/nanoflann) |   | X |
| [numpy](https://numpy.org) | Fast array data storage and manipulation | X |   |
| [gustaf](https://github.com/tataratat/gustaf) | Conversion to mesh representation, visualization, and helpers | X |  |
| [scipy](https://scipy.org) | (Optional) Creates sparse matrices, where applicable | X |   |
| [cmake](https://cmake.org) | Platform independent build system for c++ implementations |   | X |
| [setuptools](https://setuptools.pypa.io/en/latest/) | Build python package  | X |  |
| [wheel](https://wheel.readthedocs.io/en/stable/)    | Implementation of python binary packaging standard | X | X |
| [cibuildwheel](https://cibuildwheel.readthedocs.io/en/stable/) | Builds and tests wheels on CI server | X | X |
