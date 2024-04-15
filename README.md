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


You can install it directly from the source:
```bash
git clone git@github.com:tataratat/splinepy.git
cd splinepy
git submodule update --init --recursive
pip install -e .
```

## Quick start
### 1. Create a spline
Here, we will create a [NURBS](https://tataratat.github.io/splinepy/_generated/splinepy.nurbs.NURBS.html#splinepy.nurbs.NURBS) for the following example. Alternatively, we can also create [Bezier](https://tataratat.github.io/splinepy/_generated/splinepy.bezier.Bezier.html#splinepy.bezier.Bezier), [RationalBezier](https://tataratat.github.io/splinepy/_generated/splinepy.rational_bezier.RationalBezier.html#splinepy.rational_bezier.RationalBezier), and [BSpline](https://tataratat.github.io/splinepy/_generated/splinepy.bspline.BSpline.html#splinepy.bspline.BSpline).

```python
import splinepy

# Initialize nurbs with any array-like input
nurbs = splinepy.NURBS(
    degrees=[2, 1],
    knot_vectors=[
        [0, 0, 0, 1, 1, 1],
        [0, 0, 1, 1],
    ],
    control_points=[
        [-1.0, 0.0],
        [-1.0, 1.0],
        [0.0, 1.0],
        [-2.0, 0.0],
        [-2.0, 2.0],
        [0.0, 2.0],
    ],
    weights=[
        [1.0],
        [2**-0.5],
        [1.0],
        [1.0],
        [2**-0.5],
        [1.0],
    ],
)

# vizusalize
nurbs.show()
```

### 2. Modifications
All the splines can be modified. For example, by:
1. directly accessing properties,
2. [elevating degrees](https://tataratat.github.io/splinepy/_generated/splinepy.spline.Spline.elevate_degrees.html#splinepy.spline.Spline.elevate_degrees),
3. [inserting knots](https://tataratat.github.io/splinepy/_generated/splinepy.bspline.BSplineBase.insert_knots.html#splinepy.bspline.BSplineBase.insert_knots),
4. [reducing degrees](https://tataratat.github.io/splinepy/_generated/splinepy.bspline.BSplineBase.reduce_degrees.html) and [removing knots](https://tataratat.github.io/splinepy/_generated/splinepy.bspline.BSplineBase.remove_knots.html) with a specified tolerance

*Note: currently {2, 3, 4} are limited to BSpline families.*
```python
# start with a copy of the original spline
modified = nurbs.copy()

# manipulate control points
# 1. all at once
modified.control_points /= 2.0
# 2. indexwise (flat indexing)
modified.control_points[[3, 4, 5]] *= [1.3, 2.]
# 3. with grid-like indexing using multi_index helper
multi_index = modified.multi_index
modified.control_points[multi_index[0, 1]] = [-.1, -.1]
modified.control_points[multi_index[2, :]] += [2., .1]

modified.show()  # visualize Nr. 1

# elevate degrees and insert knots
modified.elevate_degrees([0, 1])
modified.show()  # visualize Nr. 2

modified.insert_knots(1, [.5])
modified.show()  # visualize Nr. 3
```

### 3. Evaluate
You can evaluate spline's basis functions, mapping, and their derivatives by giving parametric coordinate queries.
They should be 2D array-like objects and functions return 2D np.ndarray.
```python
# first, create parametric coordinate queries
queries = [
    [0.1, 0.2],  # first query
    [0.4, 0.5],  # second query
    [0.1156, 0.9091],  # third query
]

# evaluate basis, spline and derivatives.
# for derivatives, specify order per parametric dimension.
basis = nurbs.basis(queries)
basis_derivative = nurbs.basis_derivative(queries, [1, 1])
physical_coordinates = nurbs.evaluate(queries)
physical_derivatives = nurbs.derivative(queries, [2, 0])
```
Many of `splinepy`'s multi-query functions can be executed in parallal usiung multithread executions on c++ side. For that, set either global flag or pass `nthreads` argument.
```python
p_basis0 = nurbs.basis(queries, nthreads=2)
# or
splinepy.settings.NTHREADS = 3
p_basis1 = nurbs.basis(queries)
```

### 4. Helper Modules
There's a list of helper modules under the namespace `splinepy.helpme` to boost prototyping efficiencies. Please checkout the full list [here](https://tataratat.github.io/splinepy/_generated/splinepy.helpme.html)!

#### 4.1 Create
[splinepy.helpme.create](https://tataratat.github.io/splinepy/_generated/splinepy.helpme.create.html#module-splinepy.helpme.create) module can help you create several primitive shapes and another spline based on existing spline. For the latter, you can directly access such functions through [spline.create](https://tataratat.github.io/splinepy/_generated/splinepy.spline.Spline.create.html#splinepy.spline.Spline.create)
```python
# basic shapes
splinepy.show(
    ["arc", splinepy.helpme.create.arc(radius=3, angle=70)],
    ["box", splinepy.helpme.create.box(1, 2, 3)],  # length per dim
    ["circle", splinepy.helpme.create.circle(radius=2)],
    ["sphere", splinepy.helpme.create.sphere(outer_radius=2)],
    [
        "disk",
        splinepy.helpme.create.disk(outer_radius=3, inner_radius=2, angle=256),
    ],
    [
        "torus",
        splinepy.helpme.create.torus(torus_radius=3, section_outer_radius=1.5),
    ],
)
```

```python
# derived shapes
splinepy.show(
    ["extruded", nurbs.create.extruded(extrusion_vector=[1, 2, 3])],
    [
        "revolved",
        nurbs.create.revolved(
            axis=[1, 0, 0],
            center=[-1, -1, 0],
            angle=50,
        ),
    ],
)
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
| json    | (Custom) easy-to-read format, supports base64 encoding.                                                                                              |
| mfem    | [MFEM](https://mfem.org) compatible `.mesh` format. Supports structured multi-patch splines in `controlpoints_cartesian` and 2D single-patch splines |
| gismo   | [GISMO](https://gismo.github.io) compatible `.xml` format                                                                                            |
| npz     | Based on np.savez()                                                                                                                                  |
| svg     | Exports spline to svg representation. Limited to 2D geometric splines. Adheres to many of the show_options.                                          |
| xml     | [RWTH CATS](https://www.cats.rwth-aachen.de/) spline format                                                                                          |

## Dependencies
The following are direct dependencies for `splinepy`. Please feel free to check out the repositories linked.

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
