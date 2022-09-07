# splinepy
[![workflow](https://github.com/tataratat/splinepy/actions/workflows/main.yml/badge.svg)](https://github.com/tataratat/splinepy/actions)

[![PyPI version](https://badge.fury.io/py/splinepy.svg)](https://badge.fury.io/py/splinepy)  

Python Bezier, BSpline and NURBS library with C++ backend, [bezman](https://github.com/tataratat/bezman) and [SplineLib](https://github.com/tataratat/SplineLib).


## Install guide
Option 1: `pip`
```
pip install --upgrade pip
pip install splinepy
```
It is also possible to install current development version. It requires that your system has a compiler that supports C++17 or higher (It is tested with gcc-10.3 and clang-12).
```
pip install git+https://github.com/tataratat/splinepy.git@main -vvv
```
The last part `-vvv` is not necessary, but we suggest using it, since you can see the build progress. In other words, source builds takes quite long and we are working to shorten this!


Option 2: `setup.py`  
You can locally compile and install using `setup.py`.
This requires C++20 (C++17 for release mode) compatible C++ compiler
and cmake version 3.16 or higher.
```
git clone git@github.com:tataratat/splinepy.git
cd splinepy
git submodule update --init --recursive
python3 setup.py install
```

## Quick start
```
Coming Soon!
```
Test version of documentations are available [here](https://tataratat.github.io/splinepy)
