# splinepy (c++)

Library for prototyping spline geometries of arbitrary dimensions and degrees, and IGA.
`splinepy` implements core types and functions in c++ for better runtime performance.
Using [pybind11](https://github.com/pybind/pybind11), these implementations are available to use in python.
You can, for example, either use this part of the project to create pure c++ executables or to interoperate with other `c++`/`c`/`fortran` libraries.


## Install guide
There are two main ways to install the c++:

### 1) using `pip`
```bash
# sample format
# -> pip install -e. -v --config-settings=cmake.args="<cmake-build-args>" --config-settings=cmake.build-type="<build-type>"

# for example, minimal (up to spline dim 3), optimized build
pip install -e. -v --config-settings=cmake.args=-DSPLINEPY_MORE=OFF --config-settings=cmake.build-type="Release"
```
This will install compiled library to `site-packges` path.
You can find this path with:
```bash
python -c "import site; print(site.getsitepackages())"
```

### 2) using `cmake` only
```
mkdir build
cd build
cmake .. # along with other cmake options
```
